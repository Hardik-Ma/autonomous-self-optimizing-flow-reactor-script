clear all
close all
clc

%% Reactor+design spezification
V_R = 2; %2ml
c1_0 = 0.2; %mol/L
c2_0 = 0.2; %mol/L

% operating limits
tau_min = 0.25; %min
tau_max = 2; %min (I have changed it to 2 min)
T_min = 30; %°C
T_max = 50; %°C

%% OPC UA connection - LabManager
%uaClient = opcua('localhost',4841); % number for LabManager (You nedd to check if opc_toolbox is available: license('inuse')
%connect(uaClient)

%% find nodes of LabManager (cycle send!!)
% Actual Values (Inputs)
% Thermostat
%LabvisionNode_T_IN = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.T_IN'); %'Projektname:OPC_UA_SERVER(device).Var'
%T_IN = findNodeByName(LabvisionNode_T_IN,'Value');
% SyrDos1
%LabvisionNode_Q_IN = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q_IN'); 
%Q_IN = findNodeByName(LabvisionNode_Q_IN,'Value');
% SyrDos2
%LabvisionNode_Q2_IN = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q2_IN'); 
%Q2_IN = findNodeByName(LabvisionNode_Q2_IN,'Value');

% Setpoints (Outputs)
% Thermostat
%LabvisionNode_T_OUT = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.T_OUT'); 
%T_OUT = findNodeByName(LabvisionNode_T_OUT,'Value');
% SyrDos1
%LabvisionNode_Q_OUT = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q_OUT'); 
%Q_OUT = findNodeByName(LabvisionNode_Q_OUT,'Value');
% SyrDos2
%LabvisionNode_Q2_OUT = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.Q2_OUT'); 
%Q2_OUT = findNodeByName(LabvisionNode_Q2_OUT,'Value');
% MIR
%LabvisionNode_C_MIR = findNodeByName(uaClient.Namespace,'FLOWREACTOR:OPC_UA_SERVER.C_MIR'); 
%C_MIR = findNodeByName(LabvisionNode_C_MIR,'Value');
% writeValue(uaClient ,C_MIR, [1, 3, 4, 5]) %[1, 3, 4, 5]

%% SqliteDatabase connection - MIR
% load PLS model
load betaPLS_fullDoE_withplusEtOH_10b.mat     %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
betaPLS = betaPLS(:,[1 3 4 5]); %!!!!!!!!!!!!!!!
%measurment time
t_meas = 3; %sec !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Initiale Anzahl der Datensätze
lastRowCount = 0;  % Setze initialen Wert, wird später aktualisiert

% Timer konfigurieren
%{
t = timer;
% t.UserData = struct('count', 0, 'newMeasurement', 1:128); % UserData definieren
t.UserData = struct('N_Meas', 0,'timestamp',[], 'newMeas', 1:128,'AllMeas',[]); % sql2
% t.UserData = struct('count', 0, 'newMeasurement', 1:128,'measurements',[],'runningAverageMeas',1:128); % sql3
t.TimerFcn = @(~,~) readSqlDatabase2(t, betaPLS, uaClient, C_MIR, T_IN, Q_IN, Q2_IN, V_R); % Callback-Funktion definieren
t.Period = 5;  % every 5 seconds
t.ExecutionMode = 'fixedRate';

% Starte den Timer
start(t);
%}
%% Step 1: Specify problem
no_outputs = 2;              % number of objectives
no_inputs  = 2;              % number of decision variables
lb = [tau_min,T_min];        % define lower bound on decision variables, [lb1,lb2,...] [c2 equivalent,tau,T]
ub =  [tau_max,T_max];       % define upper bound on decision variables, [ub1,ub2,...]

f = @(x)f_PFR_2nd_order(x,c1_0,c2_0);   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Step 2: Generate initial dataset
dataset_size = no_inputs*2;             % initial dataset size 2*no_inputs
X = lhsdesign(dataset_size,no_inputs);  % Latin hypercube design
X = sortrows(X,2);                      % sort by temperature from smaller to higher
Y = zeros(dataset_size,no_outputs);     % corresponding matrix of response data

%% Step 3: Run initial Experiments
iter = 0;
for k = 1:size(X,1)
    X(k,:) = X(k,:).*(ub-lb)+lb;        % adjustment of bounds to real units

    tau = X(k,1); %min
    tau_sec = X(k,1).*60; %tau in sec
    T = X(k,2); %°C
    T = round(T,2); %round to 2nd digit

    % Transfer tau to flow
    Q_tot = V_R./tau; %mL/min
    Q_i = Q_tot./2;  % this is because we have two main syrdos pumps

    iter = iter+1;
    disp('---------------------------------------------')
    disp(['Exp.No: ' num2str(iter)])
T
%{
    %% 3.1a) Setpoint Temperature
    writeValue(uaClient ,T_OUT, T);
    % while waiting for T-Setpoint, set pumps to minimum flow
    writeValue(uaClient ,Q_OUT, 0.02);
    writeValue(uaClient ,Q2_OUT, 0.02);
    tic
    temperature_data = [];  % Initialize array to store [time, temperature]
    start_time = tic;
    disp('1b) Setpoint of Thermostat')
    disp(['     T = ' num2str(T) ' °C']) 
    % Waiting Loop for Thermostat
    while abs(readValue(T_IN)-T)>0.1
        current_time = toc(start_time);
        current_temp = readValue(T_IN);
        temperature_data = [temperature_data; current_time, current_temp];
    end
    % Record final point
    current_time = toc(start_time);
    current_temp = readValue(T_IN);
    temperature_data = [temperature_data; current_time, current_temp];

    time_thermostat(k) = toc./60;  % Save time of Thermostat
    disp(['Thermostat time: ' num2str(time_thermostat(k)) ' min']);
    %% Save temperature profile data to CSV
    filename = sprintf('Temperature_Profile_Exp%d_T%.1fC.csv', k, T);
    % Create table with time and temperature data
    data_table = table();
    data_table.Time_sec = temperature_data(:,1);
    data_table.Temperature_C = temperature_data(:,2);
    data_table.SetTemperature_C = repmat(T, size(temperature_data,1), 1);
    data_table.ExperimentNumber = repmat(k, size(temperature_data,1), 1);
    % Save to CSV
    writetable(data_table, filename);
    
    disp(['Temperature profile saved to: ' filename]);
   %% 3.1b) Setpoint Pumps
    writeValue(uaClient ,Q_OUT, Q_i);
    writeValue(uaClient ,Q2_OUT, Q_i);
    tic
    disp('1a) Setpoint of Pumps')
    disp(['     Q1 = ' num2str(round(Q_i,2)) ' ml/min'])
    disp(['     Q2 = ' num2str(round(Q_i,2)) ' ml/min'])
    % Waiting Loop for Pumps
    % Important for PID-controlled Pumps (Syringe pumps: Q_IST = Q_SOLL in seconds)
    % while readValue(Q_IN) ~= Q_i || readValue(Q2_IN) ~= Q_i
    while abs(Q_i-round(readValue(Q_IN),2))>0.02 || abs(Q_i-round(readValue(Q2_IN),2))>0.02
    end
    time_pumps(k) = toc./60;  % Save time of pumps
    disp(['Pumps time: ' num2str(time_pumps(k)) ' min']);

    %% 3.2) SteadyState-Time
    % Waiting Loop for SS
    time_SS = 0.1.*tau; %min    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAKE SIMILAR in the main TSEMO Code
    tic
    disp('2) SteadyStateTime')
    disp(['     tau = ' num2str(tau) ' min'])
    pause(time_SS*60) %sec.
    timestamp_SS = datetime('now');
    time_steadyState(k) = toc./60;  % Time for Steady-State
    disp(['Steady-state time: ' num2str(time_steadyState(k)) ' min']);
%}
    %% 3.3) Online MIR measurment (include PLS Regression, etc.)
    %{
    tic
    NoMeas = 3;
    disp(['d) Take ' num2str(NoMeas) ' Measurements'])
    % Waiting Loop for Measurement in SteadyState
    while t.UserData.timestamp - seconds(t_meas) < timestamp_SS  % ???Attention???
    end
    N_1st = t.UserData.N_Meas;
    % Waiting Loop for x Measurements
    NoMeas = 3;
    for j=1:NoMeas
        while t.UserData.N_Meas < N_1st+j 
            if t.UserData.N_Meas == N_1st + NoMeas-1
                break
            end
        end
    end
    % calculation of response data from mean of last 3 measurments 
    c_out = mean(t.UserData.C(end-2:end,:)); %  Test!! save as vector %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %}
    %%%% Dummy Experiment %%%%%
    c_out = f_PFR_2nd_order(X(k,:),c1_0,c2_0); %(k,:)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %%%%

    %time_Measurement(k) = toc./60;  % Time for MIR measurement
    %disp(['Measurment time: ' num2str(time_Measurement(k)) ' min']);

    % Molecular Mass [g/mol]
    M_1 = 159.09;
    M_2 = 87.12;
    M_F = 18.99;
    M_H = 1.00;
    M_3 = M_1 + M_2 - M_F - M_H;
    M_4 = M_3;
    M_5 = M_3 + M_2 - M_F - M_H;
 % c_out(2) is component 3
    Y(k,1) = -log(M_3.*c_out(3)./(tau_sec));                 %STY (g.L-3.sec-1)
    %Y(k,2) = -log(c_out(2)./(c_out(2)+c_out(3)+c_out(4)));      %neu selectivity
    Y(k,2) = -log(c_out(3)./(c_out(3)+c_out(4)+c_out(5)));      %neu selectivity (just based on product concentration)
    %Y(k,2) = -log(c_out(2)./(c1_0-c_out(1)));      %selectivity--> Not appropriate (change to Purity or E-factor for 3 DOF with 3rd FeedPump)
    %Y(k,2) = log((M_4.*c_out(4)+ M_5.*c_out(5))./(M_1.*c_out(1) + M_3.*c_out(3) + M_4.*c_out(4)+ M_5.*c_out(5))); % impurities 
    %Y(k,1) = -log((c1_0-c_out(1))./c1_0);      %Converion 

    disp(['Meas.No:' num2str(k) ])
    disp(['    Y = ' num2str(exp(-Y(k,1)))])
    disp(['    SEL = ' num2str(exp(-Y(k,2)))])
    disp(['    c_out = ' num2str(c_out)])
    disp('---------------------------------------------')
    disp('---------------------------------------------')

    %% Visulization
    %{
    figure(2)
    %colororder("sail")
    subplot(2,1,1)
    bar(k,[time_thermostat(k) time_pumps(k) time_steadyState(k) time_Measurement(k) 0],'stacked');
    hold on
    xlabel('Exp.No')
    ylabel('time [min]')
    legend('Thermostat', 'Pumpen', 'Steady-State','Measurement','AutoUpdate','off');
    title('Initial Exp.')
%}
end


%% Step 4: Start algorithm to find Pareto front
opt = TSEMO_options;             % call options for solver, see TSEMO_options file to adjust
opt.maxeval = 100;                % number of function evaluations before termination
opt.NoOfBachSequential = 1;      % number of function evaluations per iteration
% Total number of iterations = opt.maxeval/opt.NoOfBachSequential

[Xpareto,Ypareto,X,Y,XparetoGP,YparetoGP,YparetoGPstd,hypf] = TSEMO_V4_1a_OPC(f,X,Y,lb,ub,opt,c1_0,c2_0,V_R,t_meas); %

% INPUTS
% f denotes the function to be optimized
% X and Y   are the initial datasets to create a surrogate model
% lb and ub are the lower and upper bound of the decision variables
% opt is the option structure of the algorithm

% OUTPUTS
%   Xpareto and Ypareto correspond to the current best Pareto set and Pareto
%   front respectively. 
%   X and Y are the complete dataset of the decision variables and 
%   the objectives respectively. 
%   XparetoGP and YparetoGP represent the Pareto set and Pareto front of the 
%   final Gaussian process model within the algorithm. It is recommended to
%   use these as final result for problems with measurement noise. 
%   YparetoGPstd denotes the standard deviations of the predictions of
%   the GP pareto front YparetoGP  
%   hypf represents the final hyperparameters found for analysis

% For each iteration the current iteration number is displayed, the
% predicted hypervolume improvement and the time taken.

% TS-EMO creates a log file named "TSEMO_log.txt" that contains all relevant information 
% over the entire algorithm run. 

%% Step 5: Visualise results
% adapt results from logarithmic form to normal.
% ADAPT always for obective functions!!!!!!!!!!!

Y1 = exp(-Y(:,1));
Y2 = exp(-Y(:,2));
Ypareto1 = exp(-Ypareto(:,1));
Ypareto2 = exp(-Ypareto(:,2));
YparetoGP1 = exp(-YparetoGP(:,1));
YparetoGP2 = exp(-YparetoGP(:,2));
YparetoGPstd1 = exp(-YparetoGP(:,1))-exp(-(YparetoGP(:,1)+YparetoGPstd(:,1))); 
YparetoGPstd2 = exp(-YparetoGP(:,2))-exp(-(YparetoGP(:,2)+YparetoGPstd(:,2))); 
%{

Y1 = exp(-Y(:,1));
Y2 = exp(Y(:,2));
Ypareto1 = exp(-Ypareto(:,1));
Ypareto2 = exp(Ypareto(:,2));
YparetoGP1 = exp(-YparetoGP(:,1));
YparetoGP2 = exp(YparetoGP(:,2));
YparetoGPstd1 = exp(-YparetoGP(:,1))-exp(-(YparetoGP(:,1)+YparetoGPstd(:,1))); 
YparetoGPstd2 = -exp(YparetoGP(:,2))-exp((YparetoGP(:,2)+YparetoGPstd(:,2))); 
%}
figure
hold on
plot(Y1(1:dataset_size),Y2(1:dataset_size),'.','MarkerSize',14,'color',[0.8500 0.3250 0.0980])
plot(Y1(dataset_size+1:end),Y2(dataset_size+1:end),'x','MarkerSize',8,'LineWidth',2,'color',[0.8500 0.3250 0.0980])
plot(Ypareto1,Ypareto2,'O','MarkerSize',8,'LineWidth',2,'color',"#EDB120")
errorbar(YparetoGP1,YparetoGP2,YparetoGPstd2,YparetoGPstd2,YparetoGPstd1,YparetoGPstd1,'.','MarkerSize',8,'color',"#7E2F8E",'LineWidth',0.5) %,'color',"#7E2F8E"
% plot(YparetoGP1,YparetoGP2,'.','MarkerSize',8,'color',"#7E2F8E",'LineWidth',0.5) %,'color',"#7E2F8E"
legend('Initial LHC','Algorithm','Pareto front','GP Pareto front','Location','Northeast')
grid on
xlabel('STY (3) [(g.L-3.sec-1)]')
ylabel('SEL')
% set(gca, 'XScale', 'log', 'YScale', 'log')

% Display GP prediction uncertainties in a table
T = table((1:length(YparetoGPstd1))', ...
          YparetoGPstd1, ...
          YparetoGPstd2, ...
          'VariableNames', {'Point','STY_GP_STD','sel_GP_STD'});

disp('--- GP Prediction Standard Deviations (original units) ---')
disp(T)


% figure
% hold on
% plot(Y(1:dataset_size,1),Y(1:dataset_size,2),'.','MarkerSize',14,'color',[0.8500 0.3250 0.0980])
% plot(Y(dataset_size+1:end,1),Y(dataset_size+1:end,2),'x','MarkerSize',8,'LineWidth',2,'color',[0.8500 0.3250 0.0980])
% plot(Ypareto(:,1),Ypareto(:,2),'O','MarkerSize',8,'LineWidth',2,'color',"#EDB120")
% % plot(YparetoGP(:,1),YparetoGP(:,2),'.','MarkerSize',8,'color',"#7E2F8E",'LineWidth',0.5) %,'color',"#7E2F8E"
% errorbar(YparetoGP(:,1),YparetoGP(:,2),YparetoGPstd(:,2),YparetoGPstd(:,2),YparetoGPstd(:,1),YparetoGPstd(:,1),'.','MarkerSize',8,'color',"#7E2F8E",'LineWidth',0.5) %,'color',"#7E2F8E"

