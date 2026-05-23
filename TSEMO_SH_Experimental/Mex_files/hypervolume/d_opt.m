% Define parameters
nfactors = 3;     % Number of factors (e.g., N, P, K)
nruns = 10;       % Number of runs (design points)

% Generate D-optimal design for a quadratic model
[dCE, X] = cordexch(nfactors, nruns, 'quadratic', 'tries', 10);

% Display results
disp('D-Optimal Design Points (Factor Levels):');
disp(dCE);
disp('Design Matrix (Model Terms):');
disp(X);
