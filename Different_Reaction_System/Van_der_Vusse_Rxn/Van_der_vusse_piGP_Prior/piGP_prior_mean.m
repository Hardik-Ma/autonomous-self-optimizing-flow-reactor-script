function [m_yield, m_sel] = piGP_prior_mean(X_scaled, lb, ub, cA_0, theta_yield, theta_sel)
% PIGP_PRIOR_MEAN  Physics-motivated prior mean for Van de Vusse TSEMO.
%
% SAME FILENAME as the SNAr version by necessity: NLikelihood_piGP.m and
% posterior_sample_piGP.m call this function BY NAME (not via a handle),
% so this file must be named piGP_prior_mean.m to be picked up. Keep
% this in a SEPARATE working directory from your SNAr piGP files -- do
% NOT overwrite the thesis-validated SNAr copy with this one.
%
% Encodes ONLY the dominant A->B series pathway (the 2A->D side reaction
% is deliberately omitted here -- same simplification philosophy as the
% SNAr prior omitting 3 of its 4 reaction steps; the gap is absorbed by
% the learned alpha/beta correction, not by adding more mechanism):
%
%   First-order consecutive step A->B, closed form EXACT under k1==k2,
%   which is exactly true for the locked Chen/Kremling/Allgower kinetic
%   set (k10=k20, Ea1=Ea2 -- confirmed numerically, see
%   f_van_de_vusse_pfr_validate.m Check 2/3):
%     yield_B(tau)     = k1*tau*exp(-k1*tau)                    [cB/cA0]
%     selectivity(tau) = k1*tau*exp(-k1*tau) / (1 - exp(-k1*tau))
%   Arrhenius:  k1(T) = k_ref * exp(-Ea/R * (1/T_K - 1/T_ref))
%   Da = k1(T)*tau_sec   -- dimensionless. NOTE: first-order reaction,
%   so Da has NO concentration term, unlike the SNAr prior's 2nd-order
%   Da = k*c0*tau. This is why cA_0 is accepted below but unused.
%
%   yield ~ alpha*[Da*exp(-Da)]                    + beta
%   sel   ~ alpha*[Da*exp(-Da)/(1-exp(-Da))]        + beta
%
% Parameters {log_kref, Ea, alpha, beta} are LEARNED from data, exactly
% as in the SNAr prior -- the CLOSED FORM is physics (derived from the
% true ODE structure, exact in the side-reaction-off / k1=k2 limit), the
% NUMBERS are not hard-coded from the literature Chen/Kremling/Allgower
% values used in the simulator itself.
%
% Architecture note: theta_yield and theta_sel are INDEPENDENT 4-vectors
% (mirroring the SNAr prior's independent theta_STY/theta_SEL), even
% though physically both quantities are governed by the same true k1.
% This was a deliberate choice to keep the piGP methodology identical
% across reaction systems for a fair validation comparison, at the cost
% of not exploiting the shared-k1 structure. Revisit if you want a
% more physically faithful (but less directly comparable) design.
%
% INPUTS
%   X_scaled    [n x 2]  inputs scaled [0,1]: col1=tau, col2=T
%   lb          [1 x 2]  [tau_min(min), T_min(degC)]
%   ub          [1 x 2]  [tau_max(min), T_max(degC)]
%   cA_0        scalar   initial A concentration [mol/L] -- accepted for
%                        signature parity with the SNAr version, NOT
%                        used in the Da definition (see note above)
%   theta_yield [4 x 1]  [log_kref, Ea(J/mol), alpha, beta] for yield of B
%   theta_sel   [4 x 1]  [log_kref, Ea(J/mol), alpha, beta] for selectivity
%
% OUTPUTS
%   m_yield     [n x 1]  prior mean values for -log(yield of B)
%   m_sel       [n x 1]  prior mean values for -log(selectivity)

R     = 8.314;
T_ref = 125 + 273.15;   % K -- PROVISIONAL, matches Tref in
                          % f_van_de_vusse_pfr_validate.m and
                          % TSEMO_VdV_Multiple_Runs.m. Keep all three in
                          % sync if you recentre after locking real bounds.

n       = size(X_scaled, 1);
m_yield = zeros(n, 1);
m_sel   = zeros(n, 1);

k_ref_yield = exp(theta_yield(1));   Ea_yield   = theta_yield(2);
alpha_yield = theta_yield(3);        beta_yield = theta_yield(4);

k_ref_sel   = exp(theta_sel(1));     Ea_sel     = theta_sel(2);
alpha_sel   = theta_sel(3);          beta_sel   = theta_sel(4);

for i = 1:n
    tau   = X_scaled(i,1)*(ub(1)-lb(1)) + lb(1);   % min
    T_C   = X_scaled(i,2)*(ub(2)-lb(2)) + lb(2);   % degC
    tau_s = tau * 60;                                % sec
    T_K   = T_C + 273.15;                           % K

    arr_yield = exp(-Ea_yield/R * (1/T_K - 1/T_ref));
    arr_sel   = exp(-Ea_sel/R   * (1/T_K - 1/T_ref));

    Da_yield = k_ref_yield * arr_yield * tau_s;   % dimensionless
    Da_sel   = k_ref_sel   * arr_sel   * tau_s;

    yield_phys = alpha_yield * (Da_yield * exp(-Da_yield)) + beta_yield;
    sel_phys   = alpha_sel   * series_ratio(Da_sel)         + beta_sel;

    m_yield(i) = -log(max(yield_phys, 1e-12));
    m_sel(i)   = -log(max(sel_phys,   1e-12));
end
end

function r = series_ratio(Da)
% Da*exp(-Da)/(1-exp(-Da)) with a Taylor branch near Da=0 to avoid 0/0.
% The ratio -> 1 as Da -> 0 (selectivity -> 1 at tau -> 0: no time for
% either B or C to have formed yet, so whatever tiny amount of A has
% reacted is still essentially all B). Verified numerically against the
% closed-form direct evaluation across Da in [1e-8, 50] before shipping.
if Da < 1e-4
    r = 1 - Da/2 + Da^2/12;
else
    r = Da*exp(-Da) / (1 - exp(-Da));
end
end
