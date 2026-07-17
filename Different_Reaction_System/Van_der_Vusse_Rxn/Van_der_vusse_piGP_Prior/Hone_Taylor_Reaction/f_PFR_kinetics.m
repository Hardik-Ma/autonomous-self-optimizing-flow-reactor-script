function dcdz = f_PFR_kinetics(z,x,T)
% Molbalance of Plug Flow reactor with reaction kinetics for 2nd order
% Similar to MBDoE function
% DOI: 10.1039/c6re00109b

R = 8.314; % ideal gas constants [JK-1mole-1]
Tm = 90 + 273.15; % reference temperture [K]

% stochiometric matrix
 N = [-1, -1, 0, 0;
     -1, -1, -1, -1;
      1, 0, -1, 0;
      0, 1, 0, -1;
      0, 0, 1, 1];    
    
% components
c = x(1:5); 

% activation energies
E_A = [33.3e3; 35.3e3; 38.9e3; 44.8e3];
% pre-exp factors
k_0 = [57.9e-2; 2.7e-2; 0.865e-2; 1.63e-2]*1e-3;
% reaction constant
k = k_0.*exp(-E_A/R*(1/T-1/Tm));
% concentration
c_A = c(1);
c_B = c(2);
c_C = c(3);
c_D = c(4);
c_E = c(5);

% reaction rates
r_1 = k(1)*c_A*c_B;
r_2 = k(2)*c_A*c_B;
r_3 = k(3)*c_C*c_B;
r_4 = k(4)*c_D*c_B;
    
r = [r_1; r_2; r_3; r_4];

% molbalance PFR
dcdz = N*r; 

end
