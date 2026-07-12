function dcdt = f_van_de_vusse_kinetics(t, c, T, k_ref, Ea, Tref)
% F_VAN_DE_VUSSE_KINETICS  Mole balance RHS for an isothermal PFR running
% the Van de Vusse reaction network:
%
%   A --k1--> B --k2--> C      (both first order)
%   2A --k3--> D                (second order in c_A)
%
% Mirrors the structural pattern of f_PFR_kinetics.m (stoichiometric
% matrix + Arrhenius rates), adapted for a 3-reaction / 4-species network.
%
% UNIT CONVENTION (deliberate deviation from f_PFR_kinetics.m):
%   Concentrations stay in mol/L throughout — no mol/m3 internal
%   conversion. Van de Vusse literature rate constants (Bequette,
%   Chen/Kremling/Allgower) are natively reported on an L/mol basis, so
%   working in mol/L avoids a second, easy-to-miss unit conversion on
%   top of the h->s one already needed for k3.
%
% CONVENTION WARNING (read before changing k3 or the stoichiometric row
% for D): the control literature (Chen et al. 1995, Bequette, Doyle et
% al. 1995) universally writes the A-balance as
%   dcA/dt  = ... - k1*cA - k3*cA^2
% i.e. k3 as reported ALREADY is the coefficient of the A-loss term, not
% a "reaction extent" rate that gets multiplied by a stoichiometric
% factor of 2. To stay consistent with literature k3 values while still
% conserving atoms, D is produced at HALF the cA^2 rate (since 2 mol A
% is consumed per mol D formed). If you instead build a "clean" integer
% stoichiometric matrix with -2 on the A row for reaction 3, you will
% silently double the side-reaction rate relative to what k3 was fit to.
%
% INPUTS
%   t      - integration variable (residence time, s) [unused by the
%            RHS itself, required by ode45's calling convention]
%   c      - [cA; cB; cC; cD], concentrations [mol/L]
%   T      - reactor temperature [K] (isothermal, passed through)
%   k_ref  - [k1_ref; k2_ref; k3_ref], rate constants AT Tref
%            k1_ref, k2_ref [1/s] ; k3_ref [L/(mol s)]
%   Ea     - [Ea1; Ea2; Ea3], activation energies [J/mol]
%   Tref   - reference temperature [K] (centre of operating range,
%            same role as Tm=363.15 K in f_PFR_kinetics.m)
%
% OUTPUT
%   dcdt   - [dcA; dcB; dcC; dcD] / dt [mol/L/s]

R = 8.314; % ideal gas constant [J/K/mol]

% stoichiometric matrix (rows: A,B,C,D | cols: r1, r2, r3)
% r3 here is defined as k3*cA^2, i.e. already the A-consumption rate
% for the dimerization (see CONVENTION WARNING above) -> D row is 0.5,
% not the "clean" -2/+1 pattern you'd write from stoichiometry alone.
N = [-1,  0, -1;
      1, -1,  0;
      0,  1,  0;
      0,  0,  0.5];

cA = c(1);
cB = c(2);

% Arrhenius, referenced to Tref (same centred-reference convention as
% f_PFR_kinetics.m, so this drops directly into a piGP prior-mean
% structure later: k(T) = k_ref * exp(-Ea/R * (1/T - 1/Tref)) )
k = k_ref .* exp(-Ea / R .* (1 / T - 1 / Tref));

r1 = k(1) * cA;
r2 = k(2) * cB;
r3 = k(3) * cA^2;

r = [r1; r2; r3];

dcdt = N * r;

end
