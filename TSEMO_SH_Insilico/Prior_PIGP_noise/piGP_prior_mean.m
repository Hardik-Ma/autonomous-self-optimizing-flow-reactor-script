function [m_STY, m_SEL] = piGP_prior_mean(X_scaled, lb, ub, c1_0, theta_STY, theta_SEL)
% PIGP_PRIOR_MEAN  Physics-motivated prior mean for SNAr TSEMO.
%
%  Encodes ONLY what is known without fitted kinetics:
%    - 2nd order bimolecular reaction (from reaction scheme)
%    - Approximate PFR: X = Da/(1+Da),  Da = k(T)*c0*tau_sec
%    - Arrhenius: k(T) = k_ref * exp(-Ea/R * (1/T_K - 1/T_ref))
%    - STY  ~ alpha*(c0*X*M_C/tau_sec) + beta
%    - SEL  ~ alpha/(1 + Da^2) + beta
%
%  Parameters {log_kref, Ea, alpha, beta} are LEARNED from data.
%  No kinetic values from literature are hard-coded.
%
%  INPUTS
%    X_scaled   [n x 2]  inputs scaled [0,1]: col1=tau, col2=T
%    lb         [1 x 2]  [tau_min(min), T_min(degC)]
%    ub         [1 x 2]  [tau_max(min), T_max(degC)]
%    c1_0       scalar   initial DFNB concentration [mol/L]
%    theta_STY  [4 x 1]  [log_kref, Ea(J/mol), alpha, beta] for STY
%    theta_SEL  [4 x 1]  [log_kref, Ea(J/mol), alpha, beta] for SEL
%
%  OUTPUTS
%    m_STY      [n x 1]  prior mean values for -log(STY) [physical units]
%    m_SEL      [n x 1]  prior mean values for -log(SEL) [physical units]

R     = 8.314;
T_ref = 90 + 273.15;      % reference T [K] - centre of operating range
M_C   = 159.09 + 87.12 - 18.99 - 1.00;  % MW of ortho-product [g/mol]

n     = size(X_scaled, 1);
m_STY = zeros(n, 1);
m_SEL = zeros(n, 1);

k_ref_STY = exp(theta_STY(1));   Ea_STY = theta_STY(2);
alpha_STY = theta_STY(3);        beta_STY = theta_STY(4);

k_ref_SEL = exp(theta_SEL(1));   Ea_SEL = theta_SEL(2);
alpha_SEL = theta_SEL(3);        beta_SEL = theta_SEL(4);

for i = 1:n
    tau   = X_scaled(i,1)*(ub(1)-lb(1)) + lb(1);   % min
    T_C   = X_scaled(i,2)*(ub(2)-lb(2)) + lb(2);   % degC
    tau_s = tau * 60;                                % sec
    T_K   = T_C + 273.15;                           % K

    arr_STY = exp(-Ea_STY/R * (1/T_K - 1/T_ref));
    arr_SEL = exp(-Ea_SEL/R * (1/T_K - 1/T_ref));

    Da_STY = k_ref_STY * arr_STY * c1_0 * tau_s;
    Da_SEL = k_ref_SEL * arr_SEL * c1_0 * tau_s;

    X_conv   = Da_STY / (1 + Da_STY);   % 2nd-order PFR conversion
    STY_phys = alpha_STY * (c1_0 * X_conv * M_C / tau_s) + beta_STY;
    SEL_phys = alpha_SEL / (1 + Da_SEL^2) + beta_SEL;

    m_STY(i) = -log(max(STY_phys, 1e-12));
    m_SEL(i) = -log(max(SEL_phys, 1e-12));
end
end
