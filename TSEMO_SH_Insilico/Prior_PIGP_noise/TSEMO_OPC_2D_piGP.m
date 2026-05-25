clear all
close all
clc

%% =========================================================================
%  TSEMO_OPC_2D_piGP.m  —  Physics-informed TSEMO with concentration denoising
%
%  KEY ADDITION vs previous version:
%  After each MIR measurement, raw concentrations are passed through
%  denoise_concentration.m before objectives are computed.
%  The denoiser blends physics prediction with measurement using a weight w
%  learned from historical residuals — no explicit noise level needed.
%
%  Files needed in MATLAB path:
%    TSEMO_V4_1a_OPC_piGP.m    main algorithm (self-contained)
%    piGP_prior_mean.m          physics prior mean
%    NLikelihood_piGP.m         extended MAP likelihood
%    TrainingOfGP_piGP.m        joint hyperparameter optimisation
%    posterior_sample_piGP.m    spectral sampling
%    denoise_concentration.m    concentration denoiser  <-- NEW
%    TSEMO_options.m            unchanged
%    f_PFR_2nd_order.m          simulation dummy only
%    f_PFR_kinetics.m           simulation dummy only
%    + TSEMO package files (invChol, paretofront, nsga2, Direct, etc.)
% =========================================================================

%% Reactor specification
V_R  = 2;     c1_0 = 0.2;   c2_0 = 0.2;
tau_min = 0.25;  tau_max = 2;
T_min   = 30;    T_max   = 50;

%% OPC-UA (uncomment for real experiment)
%uaClient = opcua('localhost',4841); connect(uaClient)
%LabvisionNode_T_IN  = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.T_IN');
%T_IN  = findNodeByName(LabvisionNode_T_IN,'Value');
%LabvisionNode_Q_IN  = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q_IN');
%Q_IN  = findNodeByName(LabvisionNode_Q_IN,'Value');
%LabvisionNode_Q2_IN = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q2_IN');
%Q2_IN = findNodeByName(LabvisionNode_Q2_IN,'Value');
%LabvisionNode_T_OUT  = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.T_OUT');
%T_OUT  = findNodeByName(LabvisionNode_T_OUT,'Value');
%LabvisionNode_Q_OUT  = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q_OUT');
%Q_OUT  = findNodeByName(LabvisionNode_Q_OUT,'Value');
%LabvisionNode_Q2_OUT = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q2_OUT');
%Q2_OUT = findNodeByName(LabvisionNode_Q2_OUT,'Value');
%LabvisionNode_C_MIR  = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.C_MIR');
%C_MIR  = findNodeByName(LabvisionNode_C_MIR,'Value');
%writeValue(uaClient, C_MIR, [1, 3, 4, 5])

%% MIR / SQLite (uncomment for real experiment)
%load betaPLS_fullDoE_withplusEtOH_10b.mat
%betaPLS = betaPLS(:,[1 3 4 5]);
t_meas = 3;
%{
t = timer;
t.UserData = struct('N_Meas',0,'timestamp',[],'newMeas',1:128,'AllMeas',[]);
t.TimerFcn = @(~,~) readSqlDatabase2(t,betaPLS,uaClient,C_MIR,T_IN,Q_IN,Q2_IN,V_R);
t.Period = 5; t.ExecutionMode = 'fixedRate'; start(t);
%}

%% Step 1: Problem specification
no_outputs = 2;  no_inputs = 2;
lb = [tau_min, T_min];  ub = [tau_max, T_max];
M_1=159.09; M_2=87.12; M_F=18.99; M_H=1.00; M_3=M_1+M_2-M_F-M_H;

%% Denoising state (persists across all experiments)
w_blend           = 0.3;    % starting blend weight — no history yet
X_hist            = [];     % [tau, T] of all experiments so far (raw)
C_hist            = [];     % raw measured concentrations of all experiments
theta_STY_current = [-2.5; 35000; 1.0; 0.0];   % prior params, updated each iter
theta_SEL_current = [-2.5; 35000; 1.0; 0.0];

%% Step 2: Initial dataset
dataset_size = no_inputs * 2;
X = lhsdesign(dataset_size, no_inputs);
X = sortrows(X, 2);
Y = zeros(dataset_size, no_outputs);

%% Step 3: Initial experiments
iter = 0;
for k = 1:size(X,1)
    X(k,:)  = X(k,:).*(ub-lb) + lb;
    tau     = X(k,1);   tau_sec = tau*60;
    T       = round(X(k,2),2);
    Q_tot   = V_R/tau;  Q_i = Q_tot/2;
    iter    = iter+1;

    disp('---------------------------------------------')
    disp(['Exp.No: ' num2str(iter) '  tau=' num2str(tau) 'min  T=' num2str(T) 'C'])

    %{
    %% Hardware block (uncomment for real experiment)
    writeValue(uaClient,T_OUT,T);
    writeValue(uaClient,Q_OUT,0.02); writeValue(uaClient,Q2_OUT,0.02);
    tic; start_time=tic; temperature_data=[];
    while abs(readValue(T_IN)-T)>0.1
        temperature_data=[temperature_data; toc(start_time), readValue(T_IN)];
    end
    temperature_data=[temperature_data; toc(start_time), readValue(T_IN)];
    time_thermostat(k)=toc/60;
    writeValue(uaClient,Q_OUT,Q_i); writeValue(uaClient,Q2_OUT,Q_i); tic
    while abs(Q_i-round(readValue(Q_IN),2))>0.02||abs(Q_i-round(readValue(Q2_IN),2))>0.02; end
    time_pumps(k)=toc/60;
    time_SS=0.1*tau; tic; pause(time_SS*60); timestamp_SS=datetime('now');
    time_steadyState(k)=toc/60;
    tic; NoMeas=3;
    while t.UserData.timestamp-seconds(t_meas)<timestamp_SS; end
    N_1st=t.UserData.N_Meas;
    for j=1:NoMeas
        while t.UserData.N_Meas<N_1st+j
            if t.UserData.N_Meas==N_1st+NoMeas-1; break; end
        end
    end
    c_out_raw=mean(t.UserData.C(end-2:end,:));
    time_Measurement(k)=toc/60;
    %}

    %%%% Simulation dummy %%%%
    c_out_raw = f_PFR_2nd_order(X(k,:), c1_0, c2_0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %% [KEY STEP] Denoise concentrations
    [c_out, w_blend] = denoise_concentration(c_out_raw, tau, T, c1_0, c2_0, ...
                                              theta_STY_current, theta_SEL_current, ...
                                              w_blend, X_hist, C_hist);

    % Store raw measurement in history (denoiser learns from raw data)
    X_hist = [X_hist; tau, T];
    C_hist = [C_hist; c_out_raw];

    % Objectives computed from DENOISED concentrations
    Y(k,1) = -log(M_3 .* c_out(3) ./ tau_sec);
    Y(k,2) = -log(c_out(3) ./ (c_out(3)+c_out(4)+c_out(5)));

    disp(['  Raw c_out:      ' num2str(c_out_raw,'%.4f ')])
    disp(['  Denoised c_out: ' num2str(c_out,'%.4f ') ' (w=' num2str(w_blend,'%.2f') ')'])
    disp(['  STY=' num2str(exp(-Y(k,1)),'%.4f') '  SEL=' num2str(exp(-Y(k,2)),'%.4f')])
end

%% Step 4: Run physics-informed TSEMO
%
%  The algorithm proposes new [tau,T] and calls the objective function.
%  We pass a wrapper that applies denoising at each new evaluation,
%  using the state (w_blend, X_hist, C_hist, theta priors) built up
%  during initial experiments and updated as optimisation proceeds.
%
%  After each TSEMO iteration, theta_history carries the updated prior
%  params so the denoiser can be kept in sync.

opt = TSEMO_options;
opt.maxeval          = 20;
opt.NoOfBachSequential = 1;

% Pack denoiser state into a struct for the objective wrapper
ds = struct('w',         w_blend, ...
            'X_hist',    X_hist, ...
            'C_hist',    C_hist, ...
            'theta_STY', theta_STY_current, ...
            'theta_SEL', theta_SEL_current, ...
            'c1_0',      c1_0, ...
            'c2_0',      c2_0, ...
            'M_3',       M_3);

% Objective wrapper — denoises each proposed point before computing Y
f_denoised = @(x) obj_with_denoising(x, ds);

[Xpareto, Ypareto, X, Y, XparetoGP, YparetoGP, YparetoGPstd, hypf, theta_history] = ...
    TSEMO_V4_1a_OPC_piGP(f_denoised, X, Y, lb, ub, opt, c1_0, c2_0, V_R, t_meas);

%% Step 5: Visualise
Y1=exp(-Y(:,1)); Y2=exp(-Y(:,2));
Ypareto1=exp(-Ypareto(:,1)); Ypareto2=exp(-Ypareto(:,2));
YparetoGP1=exp(-YparetoGP(:,1)); YparetoGP2=exp(-YparetoGP(:,2));
YparetoGPstd1=exp(-YparetoGP(:,1))-exp(-(YparetoGP(:,1)+YparetoGPstd(:,1)));
YparetoGPstd2=exp(-YparetoGP(:,2))-exp(-(YparetoGP(:,2)+YparetoGPstd(:,2)));

figure; hold on
plot(Y1(1:dataset_size),     Y2(1:dataset_size),     '.','MarkerSize',14,'color',[0.8500 0.3250 0.0980])
plot(Y1(dataset_size+1:end), Y2(dataset_size+1:end), 'x','MarkerSize', 8,'LineWidth',2,'color',[0.8500 0.3250 0.0980])
plot(Ypareto1,  Ypareto2,  'O','MarkerSize',8,'LineWidth',2,'color',"#EDB120")
errorbar(YparetoGP1,YparetoGP2,YparetoGPstd2,YparetoGPstd2,...
         YparetoGPstd1,YparetoGPstd1,'.','MarkerSize',8,'color',"#7E2F8E",'LineWidth',0.5)
legend('Initial LHC','Algorithm','Pareto front','GP Pareto front','Location','Northeast')
grid on
xlabel('STY (g \cdot L^{-1} \cdot s^{-1})')
ylabel('Selectivity (-)')
title('Physics-informed TSEMO + concentration denoising — SNAr')

T_unc=table((1:length(YparetoGPstd1))',YparetoGPstd1,YparetoGPstd2,...
    'VariableNames',{'Point','STY_GP_STD','SEL_GP_STD'});
disp('--- GP Pareto standard deviations ---'); disp(T_unc)

n_iter_done=size(theta_history,1);
param_labels={'log k_{ref}','E_a (kJ/mol)','\alpha','\beta'};
colors={[0.2 0.5 0.8],[0.8 0.3 0.1]};
figure('Name','Prior parameter evolution')
for p=1:4
    subplot(2,2,p); hold on; grid on
    for obj=1:2
        vals=theta_history(:,p,obj);
        if p==2, vals=vals/1000; end
        plot(1:n_iter_done,vals,'o-','Color',colors{obj},...
             'MarkerFaceColor',colors{obj},'LineWidth',1.5,'MarkerSize',5)
    end
    xlabel('Iteration'); ylabel(param_labels{p}); title(param_labels{p})
    if p==1, legend({'STY','SEL'},'Location','best'); end
end
sgtitle('Physics prior parameters learned from data')


%% =========================================================================
%  Objective wrapper — applies denoising at each TSEMO-proposed point
% =========================================================================
function Y_out = obj_with_denoising(x, ds)
% Called by TSEMO_V4_1a_OPC_piGP at each proposed experiment.
% x is [1 x 2]: [tau(min), T(degC)]

tau     = x(1);
T       = round(x(2), 2);
tau_sec = tau * 60;

%%%% Simulation dummy — replace with OPC-UA + MIR hardware call %%%%
c_out_raw = f_PFR_2nd_order(x, ds.c1_0, ds.c2_0);
%%%%

% Denoise using current state
[c_out, ~] = denoise_concentration(c_out_raw, tau, T, ds.c1_0, ds.c2_0, ...
                                    ds.theta_STY, ds.theta_SEL, ...
                                    ds.w, ds.X_hist, ds.C_hist);

Y_out = zeros(1,2);
Y_out(1,1) = -log(ds.M_3 .* c_out(3) ./ tau_sec);
Y_out(1,2) = -log(c_out(3) ./ (c_out(3)+c_out(4)+c_out(5)));
end
