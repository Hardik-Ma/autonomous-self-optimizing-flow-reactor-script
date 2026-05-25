clear all
close all
clc

%% =========================================================================
%  TSEMO_OPC_2D_piGP.m  —  Physics-informed TSEMO for SNAr optimisation
%
%  Run this script.  It calls TSEMO_V4_1a_OPC_piGP which uses a
%  Damkohler-based physics prior mean jointly fitted with the GP.
%
%  Files needed in MATLAB path:
%    TSEMO_V4_1a_OPC_piGP.m    (main algorithm, self-contained)
%    piGP_prior_mean.m          (physics prior mean function)
%    NLikelihood_piGP.m         (extended MAP likelihood)
%    TrainingOfGP_piGP.m        (joint hyperparameter optimisation)
%    posterior_sample_piGP.m    (spectral sampling + mean_sample_piGP)
%    TSEMO_options.m            (unchanged from original TSEMO package)
%    f_PFR_2nd_order.m          (simulation dummy only)
%    f_PFR_kinetics.m           (simulation dummy only)
%    + all original TSEMO package files (invChol, paretofront, nsga2,
%      Direct, hypervolume2D, hypervolume3D mex files etc.)
% =========================================================================

%% Reactor specification
V_R  = 2;     % reactor volume [mL]
c1_0 = 0.2;   % initial DFNB concentration [mol/L]
c2_0 = 0.2;   % initial pyrrolidine concentration [mol/L]

tau_min = 0.25;  tau_max = 2;    % residence time bounds [min]
T_min   = 30;    T_max   = 50;   % temperature bounds [degC]

%% OPC-UA connection (uncomment for real experiment)
%uaClient = opcua('localhost',4841);
%connect(uaClient)
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
t_meas = 3;   % MIR measurement time [sec]
%{
t = timer;
t.UserData   = struct('N_Meas',0,'timestamp',[],'newMeas',1:128,'AllMeas',[]);
t.TimerFcn   = @(~,~) readSqlDatabase2(t,betaPLS,uaClient,C_MIR,T_IN,Q_IN,Q2_IN,V_R);
t.Period     = 5;
t.ExecutionMode = 'fixedRate';
start(t);
%}

%% Step 1: Problem specification
no_outputs = 2;
no_inputs  = 2;
lb = [tau_min, T_min];
ub = [tau_max, T_max];

% Black-box handle for simulation. Not called in real experiments.
f = @(x) f_PFR_2nd_order(x, c1_0, c2_0);

%% Step 2: Initial dataset (Latin hypercube)
dataset_size = no_inputs * 2;
X = lhsdesign(dataset_size, no_inputs);
X = sortrows(X, 2);           % sort by T, low to high
Y = zeros(dataset_size, no_outputs);

%% Step 3: Run initial experiments
M_1=159.09; M_2=87.12; M_F=18.99; M_H=1.00;
M_3=M_1+M_2-M_F-M_H;  M_4=M_3;  M_5=M_3+M_2-M_F-M_H;

iter = 0;
for k = 1:size(X,1)
    X(k,:)  = X(k,:).*(ub-lb) + lb;
    tau     = X(k,1);
    tau_sec = X(k,1)*60;
    T       = round(X(k,2), 2);
    Q_tot   = V_R/tau;
    Q_i     = Q_tot/2;
    iter    = iter + 1;

    disp('---------------------------------------------')
    disp(['Exp.No: ' num2str(iter) '  tau=' num2str(tau) 'min  T=' num2str(T) 'C'])
    disp(['  Q1=Q2=' num2str(round(Q_i,3)) ' mL/min'])

    %{
    %% Thermostat (uncomment for real experiment)
    writeValue(uaClient, T_OUT, T);
    writeValue(uaClient, Q_OUT, 0.02); writeValue(uaClient, Q2_OUT, 0.02);
    tic; start_time=tic; temperature_data=[];
    while abs(readValue(T_IN)-T)>0.1
        temperature_data=[temperature_data; toc(start_time), readValue(T_IN)];
    end
    temperature_data=[temperature_data; toc(start_time), readValue(T_IN)];
    time_thermostat(k)=toc/60;
    filename=sprintf('Temperature_Profile_Exp%d_T%.1fC.csv',k,T);
    writetable(table(temperature_data(:,1),temperature_data(:,2),...
        repmat(T,size(temperature_data,1),1),repmat(k,size(temperature_data,1),1),...
        'VariableNames',{'Time_sec','Temperature_C','SetTemperature_C','ExperimentNumber'}),filename);

    %% Pumps (uncomment for real experiment)
    writeValue(uaClient,Q_OUT,Q_i); writeValue(uaClient,Q2_OUT,Q_i); tic
    while abs(Q_i-round(readValue(Q_IN),2))>0.02||abs(Q_i-round(readValue(Q2_IN),2))>0.02; end
    time_pumps(k)=toc/60;

    %% Steady state (uncomment for real experiment)
    time_SS=0.1*tau; tic; pause(time_SS*60); timestamp_SS=datetime('now');
    time_steadyState(k)=toc/60;

    %% MIR measurement (uncomment for real experiment)
    tic; NoMeas=3;
    while t.UserData.timestamp-seconds(t_meas)<timestamp_SS; end
    N_1st=t.UserData.N_Meas;
    for j=1:NoMeas
        while t.UserData.N_Meas<N_1st+j
            if t.UserData.N_Meas==N_1st+NoMeas-1; break; end
        end
    end
    c_out=mean(t.UserData.C(end-2:end,:));
    time_Measurement(k)=toc/60;
    %}

    %%%% Simulation dummy — replace with hardware block above %%%%
    c_out = f_PFR_2nd_order(X(k,:), c1_0, c2_0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Y(k,1) = -log(M_3.*c_out(3)./tau_sec);
    Y(k,2) = -log(c_out(3)./(c_out(3)+c_out(4)+c_out(5)));
    disp(['  STY=' num2str(exp(-Y(k,1))) '  SEL=' num2str(exp(-Y(k,2)))])
    disp(['  c_out=' num2str(c_out)])
end

%% Step 4: Run physics-informed TSEMO
opt = TSEMO_options;
opt.maxeval          = 5;
opt.NoOfBachSequential = 1;

[Xpareto, Ypareto, X, Y, XparetoGP, YparetoGP, YparetoGPstd, hypf, theta_history] = ...
    TSEMO_V4_1a_OPC_piGP(f, X, Y, lb, ub, opt, c1_0, c2_0, V_R, t_meas);

% OUTPUTS:
%   Xpareto, Ypareto      Pareto set/front from experimental data
%   X, Y                  Full dataset (all experiments)
%   XparetoGP, YparetoGP  GP-predicted Pareto (use this as final result)
%   YparetoGPstd          GP uncertainty at each Pareto point
%   hypf                  Final kernel hyperparameters
%   theta_history         [n_iter x 4 x 2] prior params per iteration
%                         (:,:,1)=STY, (:,:,2)=SEL
%                         [log_kref, Ea(J/mol), alpha, beta]

%% Step 5: Visualise results

% Back-transform from log-space
Y1          = exp(-Y(:,1));
Y2          = exp(-Y(:,2));
Ypareto1    = exp(-Ypareto(:,1));
Ypareto2    = exp(-Ypareto(:,2));
YparetoGP1  = exp(-YparetoGP(:,1));
YparetoGP2  = exp(-YparetoGP(:,2));
YparetoGPstd1 = exp(-YparetoGP(:,1)) - exp(-(YparetoGP(:,1)+YparetoGPstd(:,1)));
YparetoGPstd2 = exp(-YparetoGP(:,2)) - exp(-(YparetoGP(:,2)+YparetoGPstd(:,2)));

% Pareto front plot
figure; hold on
plot(Y1(1:dataset_size),     Y2(1:dataset_size),     '.','MarkerSize',14,'color',[0.8500 0.3250 0.0980])
plot(Y1(dataset_size+1:end), Y2(dataset_size+1:end), 'x','MarkerSize', 8,'LineWidth',2,'color',[0.8500 0.3250 0.0980])
plot(Ypareto1,  Ypareto2,  'O','MarkerSize',8,'LineWidth',2,'color',"#EDB120")
errorbar(YparetoGP1, YparetoGP2, YparetoGPstd2, YparetoGPstd2, ...
         YparetoGPstd1, YparetoGPstd1, '.','MarkerSize',8,'color',"#7E2F8E",'LineWidth',0.5)
legend('Initial LHC','Algorithm','Pareto front','GP Pareto front (piGP)','Location','Northeast')
grid on
xlabel('STY (g \cdot L^{-1} \cdot s^{-1})')
ylabel('Selectivity (-)')
title('Physics-informed TSEMO — SNAr reaction')

% GP uncertainty table
T_unc = table((1:length(YparetoGPstd1))', YparetoGPstd1, YparetoGPstd2, ...
    'VariableNames',{'Point','STY_GP_STD','SEL_GP_STD'});
disp('--- GP Pareto standard deviations ---')
disp(T_unc)

% Prior parameter evolution (diagnostic)
n_iter_done  = size(theta_history,1);
param_labels = {'log k_{ref}','E_a (kJ/mol)','\alpha','\beta'};
colors       = {[0.2 0.5 0.8],[0.8 0.3 0.1]};
figure('Name','Prior parameter evolution')
for p = 1:4
    subplot(2,2,p); hold on; grid on
    for obj = 1:2
        vals = theta_history(:,p,obj);
        if p==2, vals=vals/1000; end
        plot(1:n_iter_done, vals, 'o-','Color',colors{obj},...
             'MarkerFaceColor',colors{obj},'LineWidth',1.5,'MarkerSize',5)
    end
    xlabel('Iteration'); ylabel(param_labels{p}); title(param_labels{p})
    if p==1, legend({'STY','SEL'},'Location','best'); end
end
sgtitle('Physics prior parameters learned from data')
