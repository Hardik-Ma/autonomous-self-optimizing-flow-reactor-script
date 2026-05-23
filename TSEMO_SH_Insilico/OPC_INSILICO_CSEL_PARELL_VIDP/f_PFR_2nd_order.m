function c_end = f_PFR_2nd_order(X,c1_0,c2_0)

% Constant DFNB concentration
c1_0 = c1_0*1000; %mol/m3 c1_0
c2_0 = c2_0*1000; %c2_0

% exp.variables (to be manipulated for optimization)
tau = X(:,1).*60;       %s
T = X(:,2) + 273.15;  %Kelvin
tspan = [0 (X(:,1).*60)]; 

% Solve of differential equation of PFR (see help for ode45)
% Integration from 0 to tau
% Start concentration: [c10; c10*eq; 0; 0; 0]
[t,c] = ode45(@(t,c) f_PFR_kinetics(t,c,T), tspan, [c1_0; c2_0; 0; 0; 0]);

% random error
rd = 0.00; % 0% random error
c_end = c(end,:) + rd*c(end,:).*randn(1,5); %% min/max
%c_end = c(end,:)
%c_end = c_end./sum(c_end([1 3 4 5])).*c1_0; %mol/L
c_end = c_end./1000; %mol/L
%c_end = c_end([1 3 4 5]); % in mol/L
end

