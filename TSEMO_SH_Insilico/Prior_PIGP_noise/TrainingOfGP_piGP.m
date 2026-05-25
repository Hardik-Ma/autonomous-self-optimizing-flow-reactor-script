function [OptGPhyp, theta_prior_out] = TrainingOfGP_piGP(Xnew, Ynew, OptGP, lb, ub, c1_0, obj_idx, theta_prior_init)
% TRAININGOFGP_PIGP  Joint MAP training of GP hyperparameters and physics
%   prior mean parameters via DIRECT global search + fmincon local refine.
%
%  Optimisation variable: hypVar = [log_lam1, log_lam2, log_sf, log_sn,
%                                    log_kref, Ea(J/mol), alpha, beta]
%
%  INPUTS
%    Xnew             [n x D]   scaled inputs [0,1]
%    Ynew             [n x 1]   standardised outputs
%    OptGP            struct    GP options (cov type, priors, h1, h2, etc.)
%    lb, ub           [1 x 2]   physical bounds [tau_min/max, T_min/max]
%    c1_0             scalar    initial concentration [mol/L]
%    obj_idx          1 or 2    which objective (1=STY, 2=SEL)
%    theta_prior_init [4 x 1]   warm-start for [log_kref,Ea,alpha,beta]
%
%  OUTPUTS
%    OptGPhyp         struct    fitted GP hyperparameters (.cov, .lik)
%    theta_prior_out  [4 x 1]   fitted prior mean parameters

[n, D] = size(Xnew);
h1     = OptGP.h1;
h2     = OptGP.h2;
h_gp   = h1 + h2;
h_prior = 4;
h_total = h_gp + h_prior;

%% Precompute squared-distance matrix
a   = Xnew';
K_M = zeros(n, n*D);
for i = 1:D
    K_M(:,(i-1)*n+1:i*n) = sqdist(a(i,:), a(i,:));
end

%% Objective
obj_fun.f = @(hypVar) NLikelihood_piGP(hypVar, Xnew, Ynew, K_M, OptGP, lb, ub, c1_0, obj_idx);

%% Bounds
lb_gp        = ones(h_gp,1) * log(sqrt(1e-3));
ub_gp        = ones(h_gp,1) * log(sqrt(1e3));
lb_gp(h_gp)  = -6;
ub_gp(h_gp)  = OptGP.noiselimit;
lb_prior = [-14;  1000; 0.01; -2.0];
ub_prior = [  5; 1e5;  10.0;  2.0];
lb_all   = [lb_gp; lb_prior];
ub_all   = [ub_gp; ub_prior];
bounds   = [lb_all, ub_all];

%% Initial point
x0_gp = [OptGP.hyp.cov(:); OptGP.hyp.lik(:)];
if nargin < 8 || isempty(theta_prior_init)
    x0_prior = [-2.5; 35000; 1.0; 0.0];
else
    x0_prior = theta_prior_init(:);
end
x0_prior = max(min(x0_prior, ub_prior), lb_prior);
x0_all   = [x0_gp; x0_prior];

%% Stage 1: DIRECT global search
opts_d.maxevals = OptGP.fun_eval * h_total;
opts_d.maxits   = 100000 * h_total;
opts_d.maxdeep  = 100000 * h_total;
opts_d.showits  = 0;
[~, x0_ref] = Direct(obj_fun, bounds, opts_d);

%% Stage 2: fmincon local refinement
LSopt = optimoptions('fmincon', ...
    'Algorithm',               'interior-point', ...
    'Display',                 'off', ...
    'SpecifyObjectiveGradient', true, ...
    'TolFun',                  1e-10, ...
    'TolX',                    1e-12, ...
    'MaxFunctionEvaluations',  5000*h_total);
hypResult = fmincon(obj_fun.f, x0_ref, [], [], [], [], lb_all, ub_all, [], LSopt);

%% Extract results
OptGPhyp.cov    = hypResult(1:h1);
OptGPhyp.lik    = hypResult(h1+1:h_gp);
theta_prior_out = hypResult(h_gp+1:h_gp+h_prior);
end

function D = sqdist(X1, X2)
% Pairwise squared Euclidean distance. Copyright (c) 2016 Mo Chen.
D = bsxfun(@plus, dot(X2,X2,1), dot(X1,X1,1)') - 2*(X1'*X2);
D(D<0) = 0;
end
