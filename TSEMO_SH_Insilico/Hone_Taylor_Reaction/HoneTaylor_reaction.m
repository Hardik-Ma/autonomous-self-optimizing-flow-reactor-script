clear all
close all
clc
%%%%% SImulation of nucleophilic aromatic substitution reaction of 
% DFNB + Pyrrolidine

% decision variables
tau = 100; % residence time [s]
T = 50+273.15; % temperature [K]K
c0 = [1000; 1000; 0; 0; 0]; % inlet concentration [mol/m3!!] 

% Definition components
LK = 1; % limiting component (Unterschuss-Komponente Edukt)
ZK = 3; % key component

%% ODE-function for PFR
tspan = [0 tau]; %s
[t x] = ode45(@(t,x) f_PFR_kinetics(t,x,T), tspan, c0);

%% objective parameters
%conversion of component 1
X = 1-x(:,1)/c0(1);

%yield of component 3
Y = (x(:,3)-c0(3))./c0(1);

%selectivity of component 3
S = (x(:,3)-c0(3))./(c0(1)-x(:,1));

%output concentrations
c_out = x(end,:);

%% Diagrams
%concentration
figure
plot(t,x./1000)
xlabel('Residence time [s]')
ylabel('Concentration [mol/L]')
legend('A','B','C','D','E')
grid on

%conversion of component 1
figure
plot(t,X)
xlabel('Residence time [s]')
ylabel('Conversion of A [mol/mol]')
grid on
ylim([0 1])

%yield of component 3
figure
plot(t,Y)
xlabel('Residence time [s]')
ylabel('Yield of 3 [mol/mol]')
grid on
ylim([0 1])

%selectivity of component 3
figure
plot(t,S)
xlabel('Residence time [s]')
ylabel('Selectivity of 3 [mol/mol]')
grid on
ylim([0 1])

