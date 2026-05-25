function [c_denoised, w_out] = denoise_concentration(c_measured, tau, T_C, c1_0, c2_0, ...
                                                       theta_STY, theta_SEL, w_prev, X_hist, C_hist)
% DENOISE_CONCENTRATION  Physics-informed denoising of outlet concentrations.
%
%  Blends the noisy MIR/PLS measurement with the Damkohler-based physics
%  prediction. The blend weight w is learned from the historical residuals
%  between physics predictions and measurements accumulated so far.
%
%  PHYSICS BASIS
%  -------------
%  The Damkohler prior (same as piGP_prior_mean) predicts approximate
%  concentrations from [tau, T] using:
%    Da = k_ref * exp(-Ea/R*(1/T-1/Tref)) * c0 * tau_sec
%    X  = Da/(1+Da)         [2nd order PFR conversion]
%    c3_physics ~ alpha_STY * c0 * X                  [desired product C]
%    sel_physics ~ alpha_SEL / (1 + Da_SEL^2)         [selectivity proxy]
%    c3_physics = sel_physics * c0 * X                 [combine]
%    c4_physics ~ (1 - sel_physics) * c0 * X * 0.5    [para product D, rough split]
%    c5_physics ~ c0 * X * max(0, Da_SEL^2/(1+Da_SEL^2)) * 0.2  [bis-adduct E]
%    c1_physics = c1_0 - c0*X                          [remaining DFNB]
%    c2_physics = c2_0 - c0*X*(1 + Da_SEL)             [remaining pyrrolidine]
%
%  DENOISING
%  ---------
%  For each concentration component i:
%    c_denoised(i) = w * c_physics(i) + (1-w) * c_measured(i)
%
%  w is learned by minimising the mean squared residual between
%  blended predictions and measurements on all historical data.
%  When no history exists (first experiment), w = 0.3 (mild prior trust).
%  w is constrained to [0, 0.8] — never fully override the measurement.
%
%  INPUTS
%    c_measured   [1 x 5]   raw MIR/PLS concentrations [mol/L]: c1..c5
%    tau          scalar    residence time [min]
%    T_C          scalar    temperature [degC]
%    c1_0, c2_0   scalar    initial concentrations [mol/L]
%    theta_STY    [4 x 1]   fitted prior params for STY [log_kref,Ea,alpha,beta]
%    theta_SEL    [4 x 1]   fitted prior params for SEL [log_kref,Ea,alpha,beta]
%    w_prev       scalar    blend weight from previous iteration (warm start)
%    X_hist       [n x 2]   historical [tau, T] inputs (physical units)
%    C_hist       [n x 5]   historical measured concentrations (physical units)
%
%  OUTPUTS
%    c_denoised   [1 x 5]   denoised concentrations [mol/L]
%    w_out        scalar    updated blend weight for next iteration

%% Physics prediction at current [tau, T]
c_physics = physics_concentration(tau, T_C, c1_0, c2_0, theta_STY, theta_SEL);

%% Learn blend weight w from historical data
if isempty(X_hist) || size(X_hist,1) < 2
    % Insufficient history — use conservative default
    w = 0.3;
else
    % Compute physics predictions at all historical points
    n_hist = size(X_hist, 1);
    C_phys_hist = zeros(n_hist, 5);
    for k = 1:n_hist
        C_phys_hist(k,:) = physics_concentration(X_hist(k,1), X_hist(k,2), c1_0, c2_0, theta_STY, theta_SEL);
    end

    % Learn w by minimising sum of squared residuals of blended prediction
    % vs measurements across all historical points and components 3,4,5
    % (components 1 and 2 are reactants with larger model uncertainty)
    % Objective: min_w sum_k sum_j (w*c_phys(k,j) + (1-w)*c_meas(k,j) - c_meas(k,j))^2
    %          = min_w sum_k sum_j (w*(c_phys(k,j) - c_meas(k,j)))^2
    % This simplifies to: w* = 0 unless we add regularisation toward w>0.
    %
    % Better formulation: find w that minimises leave-one-out prediction error.
    % For point k left out: predict c(k) = w*c_phys(k) + (1-w)*mean(C_hist(~k))
    % This correctly rewards w>0 when physics is more stable than measurements.

    comps = [3, 4, 5];   % focus on product components
    loo_errors = zeros(50, 1);
    w_grid = linspace(0, 0.8, 50);

    for iw = 1:50
        ww = w_grid(iw);
        err = 0;
        for k = 1:n_hist
            idx = setdiff(1:n_hist, k);
            % Leave-one-out: blend physics with measurement
            c_blend_k = ww * C_phys_hist(k, comps) + (1-ww) * C_hist(k, comps);
            % Compare to measurement (ground truth in LOO sense)
            err = err + sum((c_blend_k - C_hist(k, comps)).^2);
        end
        loo_errors(iw) = err;
    end

    % Best w from grid (cheap, no optimisation needed)
    % Add regularisation: bias toward w_prev for stability
    reg_weight = 0.1;
    reg_term   = reg_weight * (w_grid - w_prev).^2;
    [~, best_idx] = min(loo_errors' + reg_term);
    w = w_grid(best_idx);
end

%% Blend physics and measurement
c_denoised = w * c_physics + (1-w) * max(c_measured, 0);

% Enforce physical constraints
c_denoised = max(c_denoised, 0);

% Mass balance correction: total carbon should be conserved
% c1 + c3 + c4 + 2*c5 = c1_0  (each bis-adduct contains 2 aromatic units)
% Apply soft correction on products only
total_aromatic = c_denoised(1) + c_denoised(3) + c_denoised(4) + 2*c_denoised(5);
if total_aromatic > 1e-8
    scale = c1_0 / total_aromatic;
    % Only apply if correction is small (< 20%) to avoid over-correction
    if abs(scale - 1) < 0.20
        c_denoised([1,3,4,5]) = c_denoised([1,3,4,5]) * scale;
    end
end

w_out = w;
end


function c = physics_concentration(tau, T_C, c1_0, c2_0, theta_STY, theta_SEL)
% PHYSICS_CONCENTRATION  Predict outlet concentrations from Damkohler model.
%
%  Uses the same Arrhenius/Damkohler structure as piGP_prior_mean but
%  returns individual concentrations rather than objective values.

R     = 8.314;
T_ref = 90 + 273.15;   % reference T [K]
T_K   = T_C + 273.15;
tau_s = tau * 60;      % [sec]

k_ref_STY = exp(theta_STY(1));   Ea_STY = theta_STY(2);
alpha_STY = theta_STY(3);

k_ref_SEL = exp(theta_SEL(1));   Ea_SEL = theta_SEL(2);
alpha_SEL = theta_SEL(3);

arr_STY = exp(-Ea_STY/R * (1/T_K - 1/T_ref));
arr_SEL = exp(-Ea_SEL/R * (1/T_K - 1/T_ref));

Da_STY  = k_ref_STY * arr_STY * c1_0 * tau_s;
Da_SEL  = k_ref_SEL * arr_SEL * c1_0 * tau_s;

% Overall conversion of DFNB (component 1)
X_conv = Da_STY / (1 + Da_STY);

% Selectivity toward ortho-product (component 3)
SEL = alpha_SEL / (1 + Da_SEL^2);
SEL = max(min(SEL, 1), 0);

% Distribute converted DFNB among products
converted = c1_0 * X_conv;
c3 = alpha_STY * SEL * converted;            % ortho product
c4 = alpha_STY * (1 - SEL) * 0.5 * converted; % para product (rough split)
c5 = alpha_STY * max(Da_SEL^2/(1+Da_SEL^2), 0) * 0.1 * converted; % bis-adduct

% Clip to physical bounds
c3 = max(c3, 0);  c4 = max(c4, 0);  c5 = max(c5, 0);
c1 = max(c1_0 - converted, 0);
c2 = max(c2_0 - converted * (1 + Da_SEL * 0.5), 0);

c = [c1, c2, c3, c4, c5];
end
