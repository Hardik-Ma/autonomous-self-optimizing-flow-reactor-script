function Y = f_van_de_vusse_objective(X, cA_0, k_ref, Ea, Tref, rd)
% F_VAN_DE_VUSSE_OBJECTIVE  Full TSEMO objective function for the Van de
% Vusse network -- this is what you hand to TSEMO_V4_2_generic.m as the
% "f" argument, NOT f_van_de_vusse_pfr directly (which only returns raw
% concentrations, not the final objective vector).
%
% Objectives (locked in piGP_Validation_Project_Context.md, Section 5):
%   yield of B         = cB / cA_0                          in [0,1]
%   selectivity         = cB / (cB + cC + cD)                in [0,1]
% Internal minimisation form, matching the SNAr TSEMO convention so
% exp(-Y) recovers the physical quantity for plotting:
%   Y(:,1) = -log(yield of B)
%   Y(:,2) = -log(selectivity)
%
% INPUTS  (k_ref/Ea/Tref/rd: same convention as f_van_de_vusse_pfr.m)
%   X      - [n x 2], col1 = tau [min], col2 = T [degC]
%   cA_0   - inlet concentration of A [mol/L]
%   k_ref  - [k1_ref; k2_ref; k3_ref] at Tref [1/s; 1/s; L/(mol s)]
%   Ea     - [Ea1; Ea2; Ea3] [J/mol]
%   Tref   - reference temperature [K]
%   rd     - relative noise level (default 0.01)
%
% OUTPUT
%   Y      - [n x 2]: Y(:,1) = -log(yield of B), Y(:,2) = -log(selectivity)

if nargin < 6
    rd = 0.01;
end

c = f_van_de_vusse_pfr(X, cA_0, k_ref, Ea, Tref, rd);  % [n x 4]: cA cB cC cD

yield_B     = c(:,2) ./ cA_0;
selectivity = c(:,2) ./ (c(:,2) + c(:,3) + c(:,4));

Y = [-log(max(yield_B, 1e-12)), -log(max(selectivity, 1e-12))];

end
