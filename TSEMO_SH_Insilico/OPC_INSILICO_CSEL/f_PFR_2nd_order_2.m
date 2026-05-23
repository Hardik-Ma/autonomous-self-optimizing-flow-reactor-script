% exp.variables (to be manipulated for optimization)
x(:,1) = x(:,1) + 273.15;  %Kelvin
tau = x(:,2);       %min 
tau_sec = x(:,2).*60; %Residence Time in S 
tspan = [0 tau_sec]; 

% Solve of differential equation of PFR (see help for ode45)
% Integration from 0 to tau
% Start concentration: [c10; c10*eq; 0; 0; 0]
[t,c] = ode45(@(t,c) f_PFR_kinetics(t,c,x(:,1)), tspan, [c1_0; c2_0; 0; 0; 0]);

% Outlet concentration
c_end = c(end,:) + rd*c(end,:).*randn(1,5); 
c_end = c_end./1000; %converting it to mol/L
% Normalization + min/max

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

% log-objectives (scaling for TSEMO)
%f(:,1) = -log(c(end,3)./x(:,2));                 %yield
f(:,1) = -log((M_3.*c_end(3)./(tau_sec)));     %max. (min-(-STY)) space-time yield 
c1s = (c1_0./1000);
%f(:,2) = -log(c_end(3)./(c1s-c_end(1)));      %alt selectivity
%f(:,2) = -log(c_end(3)./(c_end(3)+c_end(4)+c_end(5)));      %neu selectivity
%f(:,2) = -log((M_3.*c_end(3))./(M_3.*c_end(3)+M_4.*c_end(4)+M_5.*c_end(5))); % Mass Based selectivity
%f(:,2) = log(((M_1.*c_end(1) + M_2.*c_end(2) + M_4.*c_end(4)+ M_5.*c_end(5))./(M_3.*c_end(3))));    % min. e-factor 
%f(:,2) = log(((M_1.*c_end(1) + M_4.*c_end(4)+ M_5.*c_end(5))./(M_3.*c_end(3))));    % min. e-factor without (comp 2)
f(:,2) = log((M_4.*c_end(4)+ M_5.*c_end(5))./(M_1.*c_end(1) + M_3.*c_end(3) + M_4.*c_end(4)+ M_5.*c_end(5))); % impurities