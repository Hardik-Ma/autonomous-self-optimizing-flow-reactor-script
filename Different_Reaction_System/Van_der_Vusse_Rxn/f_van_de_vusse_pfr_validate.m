%% f_van_de_vusse_pfr_validate.m
%
% Verification suite for f_van_de_vusse_pfr.m / f_van_de_vusse_kinetics.m,
% to run BEFORE any epsilon-constraint or TSEMO sweep. Checks four
% analytical properties that the numerical solver must reproduce.
%
% Parameter set used below: Chen, Kremling & Allgower (1995), the
% standard non-isothermal van de Vusse CSTR benchmark. Tref and the
% (tau, T) test ranges are PROVISIONAL -- revisit once operating bounds
% are chosen from the epsilon-constraint front.

clear; close all; clc;

%% --- Literature parameters (Chen, Kremling & Allgower, 1995) ---
R = 8.314; % J/mol/K

k0_h = [1.287e12; 1.287e12; 9.043e9];   % [1/h; 1/h; L/(mol h)]
Ea   = [9758.3; 9758.3; 8560] * R;      % J/mol (converted from the
                                          % k=k0*exp(Ea_K/T) convention:
                                          % Ea_J = |Ea_K| * R)

Tref = 398.15;                          % K, PROVISIONAL reference temp
                                          % (centre of an assumed 100-150 degC
                                          % operating range -- revisit)

k_ref_h = k0_h .* exp(-Ea / (R * Tref)); % rate constants at Tref, per-hour basis
k_ref   = k_ref_h / 3600;                % -> per-second basis (uniform for
                                          % k1, k2 [1/s] and k3 [L/(mol s)]
                                          % since only the time unit changes)

cA_0 = 5.1; % mol/L (literature feed concentration)

fprintf('k_ref @ Tref=%.2f K: k1=%.5g /s  k2=%.5g /s  k3=%.5g L/(mol s)\n', ...
        Tref, k_ref(1), k_ref(2), k_ref(3));

%% --- Check 1: atom balance (side reaction ON) ---
% cA + cB + cC + 2*cD = cA_0 must hold exactly for ALL tau, T, with or
% without the side reaction, since it's a mass/atom identity, not a
% consequence of any particular parameter choice. rd=0 (noise-free) so
% this isolates solver error from measurement noise.

fprintf('\n--- Check 1: atom balance (cA+cB+cC+2cD = cA_0) ---\n');
tau_test = [0.05, 0.25, 1.0, 2.0, 5.0];   % min
T_test   = [80, 100, 120, 140];           % degC
max_err  = 0;
for tau = tau_test
    for T = T_test
        c = f_van_de_vusse_pfr([tau, T], cA_0, k_ref, Ea, Tref, 0);
        balance = c(1) + c(2) + c(3) + 2*c(4);
        err = abs(balance - cA_0);
        max_err = max(max_err, err);
    end
end
fprintf('  max |atom balance error| over grid = %.3e mol/L  (tol 1e-6)\n', max_err);
assert(max_err < 1e-6, 'Atom balance violated -- check stoichiometric matrix N.');
fprintf('  PASS\n');

%% --- Check 2: closed-form series solution, side reaction OFF (k3=0) ---
% With the dimerization switched off, A->B->C reduces to the classic
% first-order series reaction. Two closed forms depending on whether
% k1 == k2:
%   general (k1 != k2):
%     cB(tau) = cA0 * k1/(k2-k1) * [exp(-k1*tau) - exp(-k2*tau)]
%   degenerate (k1 == k2) -- this is NOT a hypothetical: under the
%   Chen/Kremling/Allgower set, k10=k20 and Ea1=Ea2 EXACTLY, so
%   k1(T) == k2(T) at every T, confirmed numerically, not assumed:
%     cB(tau) = cA0 * k1 * tau * exp(-k1*tau)
% The script detects which case it's in, so it stays correct if you
% later switch to a parameter set with genuinely distinct Ea1, Ea2.

fprintf('\n--- Check 2: closed-form series solution (k3 = 0) ---\n');
T_check = 110; % degC, fixed
k_at_T  = k_ref .* exp(-Ea/R .* (1/(T_check+273.15) - 1/Tref));
k1 = k_at_T(1); k2 = k_at_T(2);
degenerate = abs(k1 - k2) < 1e-9 * max(k1, k2);
fprintf('  k1 = %.6g /s, k2 = %.6g /s  (k1==k2: %d)\n', k1, k2, degenerate);

Ea_series = Ea;
k_ref_series = k_ref; k_ref_series(3) = 0; % side reaction off

tau_grid_s = linspace(1, 600, 25); % s (up to 10 min)
max_rel_err_cB = 0;
for tau_s = tau_grid_s
    tau_min = tau_s / 60;
    c_num = f_van_de_vusse_pfr([tau_min, T_check], cA_0, k_ref_series, Ea_series, Tref, 0);

    if degenerate
        cB_ana = cA_0 * k1 * tau_s * exp(-k1*tau_s);
    else
        cB_ana = cA_0 * k1/(k2-k1) * (exp(-k1*tau_s) - exp(-k2*tau_s));
    end

    rel_err_cB = abs(c_num(2) - cB_ana) / max(cB_ana, 1e-9);
    max_rel_err_cB = max(max_rel_err_cB, rel_err_cB);
end
fprintf('  max relative error in cB vs. closed form = %.3e  (tol 1e-4)\n', max_rel_err_cB);
assert(max_rel_err_cB < 1e-4, 'PFR solver does not match series closed form -- check kinetics RHS.');
fprintf('  PASS\n');

%% --- Check 3: peak-cB location and value (side reaction OFF) ---
% Analytical optimum of the series reaction:
%   general (k1 != k2):    tau_opt = ln(k2/k1) / (k2 - k1)
%   degenerate (k1 == k2): tau_opt = 1/k1,  cB_max = cA0/e

fprintf('\n--- Check 3: peak cB location/value (k3 = 0) ---\n');
if degenerate
    tau_opt_s  = 1 / k1;
    cB_max_ana = cA_0 / exp(1);
else
    tau_opt_s  = log(k2/k1) / (k2 - k1);
    cB_max_ana = cA_0 * k1/(k2-k1) * (exp(-k1*tau_opt_s) - exp(-k2*tau_opt_s));
end

tau_opt_min = tau_opt_s / 60;
c_at_opt = f_van_de_vusse_pfr([tau_opt_min, T_check], cA_0, k_ref_series, Ea_series, Tref, 0);

fprintf('  tau_opt = %.4f min (%.2f s)\n', tau_opt_min, tau_opt_s);
fprintf('  cB_max analytical = %.5f mol/L, numerical = %.5f mol/L\n', ...
        cB_max_ana, c_at_opt(2));
rel_err_peak = abs(c_at_opt(2) - cB_max_ana) / cB_max_ana;
fprintf('  relative error = %.3e (tol 1e-4)\n', rel_err_peak);
assert(rel_err_peak < 1e-4, 'Peak cB mismatch -- check tau discretisation/ODE tolerances.');
fprintf('  PASS\n');

%% --- Check 4: short-tau differential limit ---
% As tau -> 0, cB starts at 0 and rises with initial slope k1*cA_0
% (the -k2*cB loss term is second-order small since cB is itself small).
% Full network (side reaction ON) -- this checks the RHS directly, not
% just the reduced series case.

fprintf('\n--- Check 4: short-tau differential limit (full network) ---\n');
dtau_s = 0.01; % s, small
c_small = f_van_de_vusse_pfr([dtau_s/60, T_check], cA_0, k_ref, Ea, Tref, 0);
slope_num = c_small(2) / dtau_s;
k1_full = k_ref(1) * exp(-Ea(1)/R * (1/(T_check+273.15) - 1/Tref));
slope_ana = k1_full * cA_0;
rel_err_slope = abs(slope_num - slope_ana) / slope_ana;
fprintf('  numerical dcB/dtau|0 = %.5f, analytical k1*cA_0 = %.5f\n', slope_num, slope_ana);
fprintf('  relative error = %.3e (tol 1e-3, finite-dtau discretisation)\n', rel_err_slope);
assert(rel_err_slope < 1e-3, 'Short-tau slope mismatch -- check inlet condition / RHS.');
fprintf('  PASS\n');

fprintf('\nAll checks passed. Solver is ready for epsilon-constraint / TSEMO use.\n');
