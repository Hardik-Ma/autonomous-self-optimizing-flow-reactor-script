function f = fun_PFR_SNAr_4D(x,rd)
% Calculation of objective functions for SNAr-Reaction in Plug Flow Reactor
% experimental variables: 
% [concentration 1; molar equivalent 2; residence time; temperature]

% exp.variables (to be manipulated for optimization)
x(:,1) = x(:,1).*1000;     %mol/m3
x(:,3) = x(:,3).*60;       %s
x(:,4) = x(:,4) + 273.15;  %Kelvin
tspan = [0 x(:,3)]; 

% Solve of differential equation of PFR (see help for ode45)
% Integration from 0 to tau
% Start concentration: [c10; c10*eq; 0; 0; 0]
[t,c] = ode45(@(t,c) f_PFR_kinetics(t,c,x(:,4)), tspan, [x(:,1); x(:,1).*x(:,2); 0; 0; 0]);

% Outlet concentration
c_end = c(end,:) + rd*c(end,:).*randn(1,5); 

% Molecular Mass [g/mol]
M_1 = 159.09;
M_2 = 87.12;
M_F = 18.99;
M_H = 1.00;
M_3 = M_1 + M_2 - M_F - M_H;
M_4 = M_3;
M_5 = M_3 + M_2 - M_F - M_H;

% objectives
% f(:,1) = -c(end,3)./x(:,1);                 %yield
% f(:,1) = -(M_3.*c(end,3)./x(:,3)).*3.6;     %space-time yield
% f(:,2) = -c(end,3)./(x(:,1)-c(end,1));      %selectivity
% f(:,2) = ((M_1.*c(end,1) + M_2.*c(end,2) + M_4.*c(end,4 )+ M_5.*c(end,5))./(M_3.*c(end,3)));    %e-factor

% log-objectives
% f(:,1) = -log(c(end,3)./x(:,1));                 %yield
f(:,1) = -log((M_3.*c_end(3)./x(:,3)));     %space-time yield
% f(:,2) = -log(c(end,3)./(x(:,1)-c(end,1)));      %selectivity
f(:,2) = log(((M_1.*c_end(1) + M_2.*c_end(2) + M_4.*c_end(4)+ M_5.*c_end(5))./(M_3.*c_end(3))));    %e-factor