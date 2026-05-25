function f = posterior_sample_piGP(Xnew, Yresidual, OptGP, theta_prior, lb, ub, c1_0, obj_idx)
% POSTERIOR_SAMPLE_PIGP  Draw one Thompson spectral sample from the
%   physics-informed GP posterior.
%
%   GP is trained on residuals r = Ynew - m_scaled(X|theta_prior).
%   The returned function handle evaluates the sample at test points
%   and adds back the prior mean m(x*|theta_prior).

nSp = OptGP.nSpectralpoints;
[n, D] = size(Xnew);
ell = exp(OptGP.hyp.cov(1:D));
sf2 = exp(2*OptGP.hyp.cov(D+1));
sn2 = exp(2*OptGP.hyp.lik);

sW1 = lhsdesign(nSp, D, 'criterion', 'none');
sW2 = lhsdesign(nSp, D, 'criterion', 'none');
if OptGP.cov ~= inf
    W = repmat(1./ell', nSp,1) .* norminv(sW1) .* sqrt(OptGP.cov./chi2inv(sW2, OptGP.cov));
else
    W = randn(nSp, D) .* repmat(1./ell', nSp, 1);
end
b   = 2*pi*lhsdesign(nSp, 1, 'criterion', 'none');
phi = sqrt(2*sf2/nSp) * cos(W*Xnew' + repmat(b,1,n));

A        = phi*phi' + sn2*eye(nSp);
invA     = invChol(A);
mu_th    = invA * phi * Yresidual;
cov_th   = sn2*invA; cov_th = (cov_th+cov_th')/2;
theta_sp = mvnrnd(mu_th, cov_th)';

% Prior mean standardisation constants (fixed from training set)
dummy = zeros(4,1);
if obj_idx == 1
    m_tr = piGP_prior_mean(Xnew, lb, ub, c1_0, theta_prior, dummy);
else
    m_tr = piGP_prior_mean(Xnew, lb, ub, c1_0, dummy, theta_prior);
end
m_mu = mean(m_tr); m_sig = max(std(m_tr), 1e-8);

f = @(x) sample_eval(x, theta_sp, W, b, sf2, nSp, ...
                      theta_prior, lb, ub, c1_0, obj_idx, m_mu, m_sig);
end

function vals = sample_eval(x, theta_sp, W, b, sf2, nSp, ...
                             theta_prior, lb, ub, c1_0, obj_idx, m_mu, m_sig)
phi_x  = sqrt(2*sf2/nSp) * cos(W*x' + repmat(b,1,size(x,1)));
resid  = (theta_sp' * phi_x)';
dummy  = zeros(4,1);
if obj_idx == 1
    m_t = piGP_prior_mean(x, lb, ub, c1_0, theta_prior, dummy);
else
    m_t = piGP_prior_mean(x, lb, ub, c1_0, dummy, theta_prior);
end
vals = resid + (m_t - m_mu)/m_sig;
end

function [f, varf] = mean_sample_piGP(Xnew, Yresidual, OptGP, theta_prior, lb, ub, c1_0, obj_idx)
% MEAN_SAMPLE_PIGP  Deterministic GP posterior mean + variance.
%   Used for final Pareto front estimation after optimisation completes.

nSp = OptGP.nSpectralpoints;
[n, D] = size(Xnew);
ell = exp(OptGP.hyp.cov(1:D));
sf2 = exp(2*OptGP.hyp.cov(D+1));
sn2 = exp(2*OptGP.hyp.lik);

sW1 = lhsdesign(nSp, D, 'criterion', 'none');
sW2 = lhsdesign(nSp, D, 'criterion', 'none');
if OptGP.cov ~= inf
    W = repmat(1./ell', nSp,1) .* norminv(sW1) .* sqrt(OptGP.cov./chi2inv(sW2, OptGP.cov));
else
    W = randn(nSp, D) .* repmat(1./ell', nSp, 1);
end
b   = 2*pi*lhsdesign(nSp, 1, 'criterion', 'none');
phi = sqrt(2*sf2/nSp) * cos(W*Xnew' + repmat(b,1,n));

A      = phi*phi' + sn2*eye(nSp);
invA   = invChol(A);
mu_th  = invA * phi * Yresidual;

dummy = zeros(4,1);
if obj_idx == 1
    m_tr = piGP_prior_mean(Xnew, lb, ub, c1_0, theta_prior, dummy);
else
    m_tr = piGP_prior_mean(Xnew, lb, ub, c1_0, dummy, theta_prior);
end
m_mu = mean(m_tr); m_sig = max(std(m_tr), 1e-8);

phi_fn = @(x) sqrt(2*sf2/nSp)*cos(W*x' + repmat(b,1,size(x,1)));

f = @(x) mean_eval(x, mu_th, phi_fn, theta_prior, lb, ub, c1_0, obj_idx, m_mu, m_sig);
varf = @(x) sn2 + sn2*phi_fn(x)'*invA*phi_fn(x);
end

function vals = mean_eval(x, mu_th, phi_fn, theta_prior, lb, ub, c1_0, obj_idx, m_mu, m_sig)
dummy = zeros(4,1);
if obj_idx == 1
    m_t = piGP_prior_mean(x, lb, ub, c1_0, theta_prior, dummy);
else
    m_t = piGP_prior_mean(x, lb, ub, c1_0, dummy, theta_prior);
end
vals = (mu_th' * phi_fn(x))' + (m_t - m_mu)/m_sig;
end
