function c_end = f_van_de_vusse_pfr(X, cA_0, k_ref, Ea, Tref, rd)
% F_VAN_DE_VUSSE_PFR  Isothermal PFR simulator for the Van de Vusse
% reaction network, mirroring the structural pattern of
% f_PFR_2nd_order.m (ode45 call + multiplicative noise at the outlet).
%
%   A --k1--> B --k2--> C      (both first order)
%   2A --k3--> D                (second order in c_A)
%
% Decision variables: X = [tau (min), T (degC)] -- same convention as
% the thesis TSEMO scripts (tau in minutes, T in Celsius on the
% TSEMO-facing side; converted to SI internally).
%
% INPUTS
%   X      - [n x 2], col1 = tau [min], col2 = T [degC]
%   cA_0   - inlet concentration of A [mol/L], cB_0 = cC_0 = cD_0 = 0
%            (fresh A-only feed, standard van de Vusse setup)
%   k_ref  - [k1_ref; k2_ref; k3_ref], rate constants AT Tref
%            k1_ref, k2_ref [1/s] ; k3_ref [L/(mol s)]
%            NOTE: literature k10/k20/k30 are typically reported per
%            HOUR (Chen et al. 1995) or per MINUTE (Bequette). Convert
%            to per-SECOND before calling this function -- see the
%            example instantiation in f_van_de_vusse_pfr_validate.m.
%   Ea     - [Ea1; Ea2; Ea3], activation energies [J/mol]
%   Tref   - reference temperature [K]
%   rd     - relative noise level (default 0.01, matches thesis
%            in-silico convention; set rd = 0 for the epsilon-constraint
%            / verification runs, exactly as in f_PFR_2nd_order.m)
%
% OUTPUT
%   c_end  - [n x 4] outlet concentrations [cA cB cC cD], mol/L

if nargin < 6
    rd = 0.01; % 1% Noise injected into the system
end

n = size(X, 1);
c_end = zeros(n, 4);

for i = 1:n
    tau = X(i, 1) * 60;        % min -> s
    T   = X(i, 2) + 273.15;    % degC -> K
    tspan = [0, tau];

    c0 = [cA_0; 0; 0; 0];      % fresh feed: only A present at inlet

    [~, c] = ode45(@(t, c) f_van_de_vusse_kinetics(t, c, T, k_ref, Ea, Tref), ...
                   tspan, c0);

    c_raw = max(c(end, :), 0);            % clip solver undershoot before noise
    c_meas = c_raw .* (1 + rd * randn(1, 4));
    c_end(i, :) = max(c_meas, 0);         % guard against noise pushing <0
end

end
