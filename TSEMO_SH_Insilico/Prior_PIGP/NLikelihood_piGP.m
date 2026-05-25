function [NLL, dNLL] = NLikelihood_piGP(hypVar, Xnew, Ynew_raw, K_M, OptGP, lb, ub, c1_0, obj_idx)
% NLIKELIHOOD_PIGP  Joint MAP negative log-likelihood for physics-informed GP.
%
%  hypVar layout: [log_lambda1, log_lambda2, log_sf, log_sn,
%                  log_kref, Ea(J/mol), alpha, beta]
%  The GP is trained on residuals r = Ynew_raw - m_scaled(X|theta_prior).
%  MAP priors are placed on all parameters.

[n, D]  = size(Xnew);
h1      = OptGP.h1;   % D+1 cov hyperparams
h2      = OptGP.h2;   % 1  lik hyperparam
h_gp    = h1 + h2;
h_prior = 4;          % [log_kref, Ea, alpha, beta]

hyp_cov      = hypVar(1:h1);
hyp_lik      = hypVar(h1+1:h_gp);
theta_prior  = hypVar(h_gp+1:h_gp+h_prior);

%% Physics prior mean at training points
dummy = zeros(4,1);
if obj_idx == 1
    [m_phys, ~] = piGP_prior_mean(Xnew, lb, ub, c1_0, theta_prior, dummy);
else
    [~, m_phys] = piGP_prior_mean(Xnew, lb, ub, c1_0, dummy, theta_prior);
end

m_mu  = mean(m_phys);
m_sig = std(m_phys);
if m_sig < 1e-8, m_sig = 1; end
m_scaled = (m_phys - m_mu) / m_sig;
r = Ynew_raw - m_scaled;   % residuals the GP models

%% Covariance matrix (Matern, same as original TSEMO)
if OptGP.cov ~= inf, d = OptGP.cov; else, d = 1; end
ell = exp(hyp_cov(1:D));
sf2 = exp(2*hypVar(D+1));
K   = zeros(n,n);
for i = 1:D
    K = K_M(:,(i-1)*n+1:i*n) * d/ell(i)^2 + K;
end
if OptGP.cov ~= inf
    sqrtK = sqrt(K); expnK = exp(-sqrtK);
else
    expnK = exp(-1/2*K); sqrtK = [];
end
if     OptGP.cov == 3,   t=sqrtK; Km=(1+t).*expnK;
elseif OptGP.cov == 1,            Km=expnK;
elseif OptGP.cov == 5,   t=sqrtK; Km=(1+t.*(1+t/3)).*expnK;
elseif OptGP.cov == inf,          Km=expnK;
end
K_full = sf2*Km + eye(n)*exp(hyp_lik*2);
K_full = (K_full+K_full')/2;
try
    CH = chol(K_full); invK = CH\(CH'\eye(n));
catch
    CH = chol(K_full+eye(n)*1e-4); invK = CH\(CH'\eye(n));
end
logDetK = 2*sum(log(abs(diag(CH))));

%% MAP priors on GP hyperparameters
logprior = 0;
dlogpriorcov = zeros(1,h1);
for i = 1:h1
    [A, dlogpriorcov(i)] = priorGauss(OptGP.priorcov(1), OptGP.priorcov(2), hyp_cov(i));
    logprior = logprior + A;
end
dlogpriorlik = zeros(1,h2);
for i = 1:h2
    [A, dlogpriorlik(i)] = priorGauss(OptGP.priorlik(1), OptGP.priorlik(2), hyp_lik(i));
    logprior = logprior + A;
end

%% MAP priors on physics prior parameters
% Centred on physically reasonable SNAr values, wide enough to learn freely
[lp_logk,  dlp_logk]  = priorGauss(-2.5,   9.0,      theta_prior(1));  % log(k_ref)
[lp_Ea,    dlp_Ea]    = priorGauss(40000,   4e8,      theta_prior(2));  % Ea [J/mol]
[lp_alpha, dlp_alpha] = priorGauss(1.0,     1.0,      theta_prior(3));  % alpha
[lp_beta,  dlp_beta]  = priorGauss(0.0,     0.25,     theta_prior(4));  % beta
logprior = logprior + lp_logk + lp_Ea + lp_alpha + lp_beta;

%% Negative log-posterior
NLL = n/2*log(2*pi) + 1/2*logDetK + 1/2*(r'*invK*r) - logprior;

%% Gradients
if nargout == 2
    c = invK * r;
    dsq_M = zeros(n, n*D);
    for i = 1:D
        dsq_M(:,(i-1)*n+1:i*n) = K_M(:,(i-1)*n+1:i*n)*d/ell(i)^2;
    end

    % Analytical gradients for GP hyperparameters (on residuals r)
    dNLL_cov = zeros(h1,1);
    for i = 1:h1
        dK = covMaternanisotropic(OptGP.cov, hyp_cov, sqrtK, expnK, dsq_M, Xnew, [], i);
        b  = invK*dK;
        dNLL_cov(i) = 1/2*trace(b) - 1/2*(r'*b*c);
    end
    dNLL_lik = zeros(h2,1);
    for i = 1:h2
        dK = 2*exp(hyp_lik(i))*eye(n)*exp(hyp_lik(i));
        b  = invK*dK;
        dNLL_lik(i) = 1/2*trace(b) - 1/2*(r'*b*c);
    end

    % Numerical gradients for prior mean parameters (central differences)
    eps_fd = 1e-5;
    dNLL_prior = zeros(h_prior,1);
    for ip = 1:h_prior
        tp_f = theta_prior; tp_f(ip) = tp_f(ip)+eps_fd;
        tp_b = theta_prior; tp_b(ip) = tp_b(ip)-eps_fd;
        dummy = zeros(4,1);
        if obj_idx == 1
            [mf,~] = piGP_prior_mean(Xnew,lb,ub,c1_0,tp_f,dummy);
            [mb,~] = piGP_prior_mean(Xnew,lb,ub,c1_0,tp_b,dummy);
        else
            [~,mf] = piGP_prior_mean(Xnew,lb,ub,c1_0,dummy,tp_f);
            [~,mb] = piGP_prior_mean(Xnew,lb,ub,c1_0,dummy,tp_b);
        end
        mf_s = (mf-mean(mf))/max(std(mf),1e-8);
        mb_s = (mb-mean(mb))/max(std(mb),1e-8);
        dm = (mf_s - mb_s)/(2*eps_fd);
        dNLL_prior(ip) = c' * dm;   % dr/dtheta = -dm/dtheta, double negative cancels
    end

    dNLL = [dNLL_cov - dlogpriorcov'; ...
            dNLL_lik - dlogpriorlik'; ...
            dNLL_prior - [dlp_logk; dlp_Ea; dlp_alpha; dlp_beta]];
end
end

function [lp,dlp] = priorGauss(mu,s2,x)
% Copyright (c) 2005-2017 Carl Edward Rasmussen & Hannes Nickisch. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%    1. Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%    2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY CARL EDWARD RASMUSSEN & HANNES NICKISCH ``AS IS''
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL CARL EDWARD RASMUSSEN & HANNES NICKISCH OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The views and conclusions contained in the software and documentation
% are those of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of Carl Edward Rasmussen & Hannes Nickisch.
%
% The code and associated documentation is available from http://gaussianprocess.org/gpml/code.</pre>

lp  = -(x-mu).^2/(2*s2) - log(2*pi*s2)/2;
dlp = -(x-mu)/s2;

end

function K = covMaternanisotropic(d, hyp, sqrtK,expnK, dsq_M, x, z, i)
% Copyright (c) 2005-2017 Carl Edward Rasmussen & Hannes Nickisch. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%    1. Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%    2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY CARL EDWARD RASMUSSEN & HANNES NICKISCH ``AS IS''
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL CARL EDWARD RASMUSSEN & HANNES NICKISCH OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The views and conclusions contained in the software and documentation
% are those of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of Carl Edward Rasmussen & Hannes Nickisch.
%
% The code and associated documentation is available from http://gaussianprocess.org/gpml/code.</pre>

[n,D] = size(x);
sf2 = exp(2*hyp(D+1));

if nargin<7                                                        % covariances
    if      d == 3, t = sqrtK ; m =  (1 + t).*expnK;
    elseif  d == 1,             m =  expnK;
    elseif  d == 5, t = sqrtK ; m =  (1 + t.*(1+t/3)).*expnK;
    elseif  d == inf, m = expnK;
    end
    K = sf2*m;
else                                                               % derivatives
    if i<=D                                               % length scale parameter
        Ki = dsq_M(:,(i-1)*n+1:i*n) ;
        if     d == 3,             dm = expnK;
        elseif d == 1, t = sqrtK ; dm = (1./t).*expnK;
        elseif d == 5, t = sqrtK ; dm = ((1+t)/3).*expnK;
        elseif d == inf; dm = -1/2*expnK;
        end
        
        K = sf2*dm.*Ki;
        K(Ki<1e-12) = 0;                                    % fix limit case for d=1
    elseif i==D+1                                            % magnitude parameter
        if      d == 3, t = sqrtK ; m =  (1 + t).*expnK;
        elseif  d == 1,             m =  expnK;
        elseif  d == 5, t = sqrtK ; m =  (1 + t.*(1+t/3)).*expnK;
        elseif  d == inf,           m = expnK;
        end
        K = 2*sf2*m;
    end
end

end
