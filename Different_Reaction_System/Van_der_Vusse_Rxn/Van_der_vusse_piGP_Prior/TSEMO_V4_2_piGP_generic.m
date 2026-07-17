function [Xpareto,Ypareto,X,Y,XParetoGP,YParetoGP,YParetoGPstd,hypf,theta_history] = ...
    TSEMO_V4_2_piGP_generic(f, X, Y, lb, ub, opt, cA_0, varargin)
% TSEMO_V4_2_PIGP_GENERIC  Physics-informed TSEMO, reaction-agnostic.
%
%  Forked from TSEMO_V4_1a_OPC_piGP.m (left untouched on disk for thesis
%  reproducibility) -- same fix as the plain-TSEMO fork
%  (TSEMO_V4_2_generic.m): f is now used for the true re-evaluation of
%  every proposed point, instead of a hardcoded f_PFR_2nd_order call +
%  SNAr STY/selectivity formulas that ignored f entirely. f must be a
%  full objective function (X -> Y, already log-transformed), e.g.
%  f_van_de_vusse_objective.m.
%
%  Also fixed: the same parfor-incompatible bare tic/toc pattern found
%  in the plain-TSEMO core file (see iter_tic below).
%
%  Also recentred: cold-start theta_yield/theta_sel below, and reaction-
%  specific trailing args (c2_0, V_R, t_meas) dropped -- absorbed into
%  varargin for backward compatibility with old-style calls, unused.
%
%  Drop-in replacement for TSEMO_V4_1a_OPC_piGP.m.
%  Injects a learnable Damkohler-based prior mean m(x|theta) into each GP.
%  theta = {log_kref, Ea, alpha, beta} is optimised jointly with GP
%  hyperparameters via MAP at every iteration (TrainingOfGP_piGP).
%
%  Extra output: theta_history [n_iter x 4 x 2]
%    (:,:,1) = yield-of-B prior params per iteration
%    (:,:,2) = selectivity prior params per iteration
%
%  cA_0 is still a required argument (not folded into the f closure)
%  because piGP_prior_mean, TrainingOfGP_piGP, and posterior_sample_piGP
%  all need it directly for the physics prior -- it isn't just a
%  simulator input, so it can't be hidden inside f the way k_ref/Ea/rd
%  can for f_van_de_vusse_objective.

it     = 1;
Opt    = set_option_structure(opt, X, Y);
n_iter = ceil(Opt.maxeval / Opt.NoOfBachSequential);

% Cold-start prior parameters [log_kref, Ea(J/mol), alpha, beta]
% Recentred for Van de Vusse (was [-2.5; 35000; 1.0; 0.0], SNAr-scale) --
% keep in sync with NLikelihood_piGP.m's MAP prior centers and
% TrainingOfGP_piGP.m's default warm-start if you retune any of them.
theta_yield   = [-4.8; 75000; 1.0; 0.0];
theta_sel     = [-4.8; 75000; 1.0; 0.0];
theta_history = zeros(n_iter, 4, 2);

create_log_file(X, Y, Opt, f, lb, ub)

for i = 1:n_iter
    iter_tic = tic; % ID-based timer -- parfor-safe, see TSEMO_V4_2_generic.m
                     % for the full explanation of why bare tic/toc breaks
                     % under parfor.

    %% Scale variables
    [Xnew, Ynew, Y_mean, Y_std] = ScaleVariables(X, Y, lb, ub, Opt);

    %% Joint MAP: GP hyperparameters + prior mean parameters
    fprintf('  [piGP] Iter %d: fitting GP + prior params...\n', it)
    [Opt.GP(1).hyp, theta_yield] = TrainingOfGP_piGP(Xnew, Ynew(:,1), Opt.GP(1), lb, ub, cA_0, 1, theta_yield);
    [Opt.GP(2).hyp, theta_sel]   = TrainingOfGP_piGP(Xnew, Ynew(:,2), Opt.GP(2), lb, ub, cA_0, 2, theta_sel);
    theta_history(i,:,1) = theta_yield';
    theta_history(i,:,2) = theta_sel';
    fprintf('  Yield: k_ref=%.3g, Ea=%.1fkJ/mol, a=%.3f, b=%.3f\n', exp(theta_yield(1)),theta_yield(2)/1000,theta_yield(3),theta_yield(4))
    fprintf('  Sel:   k_ref=%.3g, Ea=%.1fkJ/mol, a=%.3f, b=%.3f\n', exp(theta_sel(1)),theta_sel(2)/1000,theta_sel(3),theta_sel(4))

    %% Compute residuals
    dummy = zeros(4,1);
    m1 = piGP_prior_mean(Xnew, lb, ub, cA_0, theta_yield, dummy);
    m2 = piGP_prior_mean(Xnew, lb, ub, cA_0, dummy, theta_sel);
    m1s = (m1-mean(m1))/max(std(m1),1e-8);
    m2s = (m2-mean(m2))/max(std(m2),1e-8);
    Yres = [Ynew(:,1)-m1s, Ynew(:,2)-m2s];

    %% Draw spectral posterior samples (physics-informed)
    Opt.Sample(1).f = posterior_sample_piGP(Xnew, Yres(:,1), Opt.GP(1), theta_yield, lb, ub, cA_0, 1);
    Opt.Sample(2).f = posterior_sample_piGP(Xnew, Yres(:,2), Opt.GP(2), theta_sel, lb, ub, cA_0, 2);

    %% Standard TSEMO: Pareto of samples -> HV improvement -> xNew
    [Sample_pareto, Sample_xpareto, Sample_nadir] = Find_sample_pareto(Opt, i);
    Opt.warmstart_pareto = Sample_xpareto;
    [index, hv_imp] = hypervolume_improvement_index(Ynew, Sample_nadir, Sample_pareto, Opt);
    xNew  = Sample_xpareto(index,:);
    Xnew  = [Xnew; xNew];
    for j = 1:Opt.Gen.NoOfInputDim
        xnewtrue(:,j) = xNew(:,j)*(ub(j)-lb(j)) + lb(j);
    end

    %% Evaluate true system via the SAME f handle used for the initial
    %% dataset -- this is the fix. The original file hardcoded
    %% f_PFR_2nd_order + SNAr STY/selectivity math here regardless of f.
    for l = 1:size(xnewtrue,1)
        ytrue(l,:) = f(xnewtrue(l,:));
    end

    X = [X; xnewtrue];
    Y = [Y; ytrue];
    disp('---------------------------------------------')
    disp(['Opt.Exp.No: ' num2str(it) '  tau=' num2str(xnewtrue(:,1)') 'min  T=' num2str(xnewtrue(:,2)') 'C'])
    disp(['  Y(1)=' num2str(exp(-Y(end,1))) '  Y(2)=' num2str(exp(-Y(end,2)))])
    disp('---------------------------------------------')

    front   = paretofront(Y);
    Xpareto = X(front,:);
    Ypareto = Y(front,:);
    update_log_file(it, hv_imp, toc(iter_tic), xnewtrue, ytrue, Opt, Y, ub, lb)
    if it==1, fprintf('%10s %10s %10s\n','Iteration','HypImp','Time(s)'); end
    fprintf('%10d %10.4g %10.3g\n', it, hv_imp, toc(iter_tic));
    it = it+1;

    %% Final iteration: GP Pareto
    if i == n_iter
        [Xnew, Ynew, Y_mean, Y_std] = ScaleVariables(X, Y, lb, ub, Opt);
        [Opt.GP(1).hyp, theta_yield] = TrainingOfGP_piGP(Xnew,Ynew(:,1),Opt.GP(1),lb,ub,cA_0,1,theta_yield);
        [Opt.GP(2).hyp, theta_sel]   = TrainingOfGP_piGP(Xnew,Ynew(:,2),Opt.GP(2),lb,ub,cA_0,2,theta_sel);

        dummy = zeros(4,1);
        m1f=(piGP_prior_mean(Xnew,lb,ub,cA_0,theta_yield,dummy));
        m2f=(piGP_prior_mean(Xnew,lb,ub,cA_0,dummy,theta_sel));
        m1fs=(m1f-mean(m1f))/max(std(m1f),1e-8);
        m2fs=(m2f-mean(m2f))/max(std(m2f),1e-8);
        Yres_f=[Ynew(:,1)-m1fs, Ynew(:,2)-m2fs];

        hypf = zeros(Opt.Gen.NoOfInputDim+2, Opt.Gen.NoOfGPs);
        for j = 1:Opt.Gen.NoOfGPs
            ch=exp(Opt.GP(j).hyp.cov);
            hypf(:,j)=[ch(1:Opt.Gen.NoOfInputDim).*(ub-lb)'; ch(end)*std(Y(:,j)); exp(Opt.GP(j).hyp.lik)*std(Y(:,j))];
        end

        [Opt.Mean(1).f, Opt.Mean(1).varf] = mean_sample_piGP(Xnew,Yres_f(:,1),Opt.GP(1),theta_yield,lb,ub,cA_0,1);
        [Opt.Mean(2).f, Opt.Mean(2).varf] = mean_sample_piGP(Xnew,Yres_f(:,2),Opt.GP(2),theta_sel,lb,ub,cA_0,2);
        [Mean_pareto, Mean_xpareto] = Find_mean_pareto(Opt);

        for j = 1:Opt.Gen.NoOfInputDim
            XParetoGP(:,j) = Mean_xpareto(:,j)*(ub(j)-lb(j)) + lb(j);
        end
        for j = 1:Opt.Gen.NoOfGPs
            YParetoGP(:,j) = Mean_pareto(:,j)*Y_std(j) + Y_mean(j);
        end
        for j = 1:Opt.Gen.NoOfGPs
            for k = 1:Opt.pop
                YParetoGPstd(k,j) = sqrt(Opt.Mean(j).varf(Mean_xpareto(k,:)))*Y_std(j);
            end
        end
        final_log_update(Xpareto,Ypareto,X,Y,XParetoGP,YParetoGP,hypf,Opt)
        fprintf('\n--- Final prior parameters ---\n')
        fprintf('Yield: k_ref=%.4g, Ea=%.1fkJ/mol, alpha=%.4f, beta=%.4f\n', exp(theta_yield(1)),theta_yield(2)/1000,theta_yield(3),theta_yield(4))
        fprintf('Sel:   k_ref=%.4g, Ea=%.1fkJ/mol, alpha=%.4f, beta=%.4f\n', exp(theta_sel(1)),theta_sel(2)/1000,theta_sel(3),theta_sel(4))
    end
end
end

function [Xnew,Ynew,Y_mean,Y_std] = ScaleVariables(X,Y,lb,ub,Opt)
Xnew=zeros(size(X)); Ynew=zeros(size(Y));
Y_mean=zeros(Opt.Gen.NoOfGPs,1); Y_std=zeros(Opt.Gen.NoOfGPs,1);
for i=1:size(X,2), Xnew(:,i)=(X(:,i)-lb(i))/(ub(i)-lb(i)); end
for i=1:size(Y,2)
    Y_mean(i)=mean(Y(:,i)); Y_std(i)=std(Y(:,i));
    if Y_std(i)<1e-10, Y_std(i)=1; end
    Ynew(:,i)=(Y(:,i)-Y_mean(i))/Y_std(i);
end
end

function Opt=set_option_structure(Opt,X,Y)
Opt.Gen.NoOfGPs=size(Y,2); Opt.Gen.NoOfInputDim=size(X,2);
for i=1:Opt.Gen.NoOfGPs
    Opt.GP(i).cov=Opt.GP(i).matern; Opt.GP(i).noiselimit=0; Opt.GP(i).var=10;
    Opt.GP(i).h1=Opt.Gen.NoOfInputDim+1; Opt.GP(i).h2=1;
    Opt.GP(i).priorlik=[-6,Opt.GP(i).var]; Opt.GP(i).priorcov=[0,Opt.GP(i).var];
    Opt.GP(i).hyp.cov=zeros(1,Opt.GP(i).h1); Opt.GP(i).hyp.lik=log(1e-2);
end
end


function create_log_file(X,Y,Opt,f,lb,ub)
try
    function_name = func2str(f);
    string1 = '';
    string2 = {};
    string3 = '';
    string4 = '';
    string5 = {};
    string6 = '';
    string7 = '';
    for i = 1:size(X,2)
        string1 = strcat(string1,'%8.4f ');
        string2 = {string2{:},strcat('x',num2str(i))};
        string3 = strcat(string3,'%+8s ');
    end
    for i = 1:size(Y,2)
        string4 = strcat(string4,'%8.4f ');
        string5 = {string5{:},strcat('f',num2str(i))};
        string6 = strcat(string6,'%+8s ');
        string7 = strcat(string7,'%8d ');
    end
    
    TSEMO_log = fopen( 'TSEMO_log.txt', 'w');
    fprintf(TSEMO_log,'\n %s %s \n', 'TSEMO log file created on',date);
    fprintf(TSEMO_log,'\n %s','This file shows the initial specifications of TSEMO and logs the output.');
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'License information');
    fprintf(TSEMO_log,'\n %s \n', 'BSD 2-Clause License');
    fprintf(TSEMO_log,'\n %s', 'Copyright (c) 2017, Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin');
    fprintf(TSEMO_log,'\n %s \n', 'All rights reserved.');
    fprintf(TSEMO_log,'\n %s', 'Redistribution and use in source and binary forms, with or without');
    fprintf(TSEMO_log,'\n %s \n', 'modification, are permitted provided that the following conditions are met:');
    fprintf(TSEMO_log,'\n %s   ', '*Redistributions of source code must retain the above copyright notice, this');
    fprintf(TSEMO_log,'\n %s \n', ' list of conditions and the following disclaimer.');
    fprintf(TSEMO_log,'\n %s   ', '*Redistributions in binary form must reproduce the above copyright notice,');
    fprintf(TSEMO_log,'\n %s   ', ' this list of conditions and the following disclaimer in the documentation');
    fprintf(TSEMO_log,'\n %s \n', ' and/or other materials provided with the distribution.');
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'Problem specifications');
    fprintf(TSEMO_log,'\n %s %s \n', 'Function used:  ',function_name);
    fprintf(TSEMO_log,'\n %s %d', 'Number of inputs:  ',size(X,2));
    fprintf(TSEMO_log,'\n %s %d \n', 'Number of outputs: ',size(Y,2));
    fprintf(TSEMO_log,'\n %s', 'Lower bounds of decision variables:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,string1,lb);
    fprintf(TSEMO_log,'\n \n %s', 'Upper bounds of decision variables:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,string1,ub);
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'Algorithm options');
    fprintf(TSEMO_log,'\n %s %d', 'Maximum number of function evaluations: ',Opt.maxeval);
    fprintf(TSEMO_log,'\n %s %d', 'Sample batch size:                      ',Opt.NoOfBachSequential);
    fprintf(TSEMO_log,'\n %s %d \n', 'Number of algorithm iterations:         ',ceil(Opt.maxeval/Opt.NoOfBachSequential));
    fprintf(TSEMO_log,'\n %s %d', 'Genetic algorithm population size:       ',Opt.pop);
    fprintf(TSEMO_log,'\n %s %d \n', 'Genetic algorithm number of generations: ',Opt.Generation);
    fprintf(TSEMO_log,strcat('\n','%s',string6,'\n'),'                                         ',string5{:});
    fprintf(TSEMO_log,strcat('%s',string7,'\n'), ' Number of spectral sampling points:     ',Opt.GP(1:size(Y,2)).nSpectralpoints);
    fprintf(TSEMO_log,strcat('%s',string7), ' Type of matern function:                ',Opt.GP(1:size(Y,2)).matern);
    fprintf(TSEMO_log,strcat('\n','%s',string7), ' Direct evaluations per input dimension: ',Opt.GP(1:size(Y,2)).fun_eval);
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'Initial data set');
    fprintf(TSEMO_log,'\n %s %d \n', 'Number of initial data points: ',size(X,1));
    fprintf(TSEMO_log,'\n %s', 'Initial input data matrix:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),X');
    fprintf(TSEMO_log,'\n %s', 'Initial output data matrix:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),Y');
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    fclose(TSEMO_log) ;
catch
    warning('There was an error when writing the log file. Maybe the file was currently open in another program. The algorithm continues to run but some parts of the log file are maybe not correct.') ;
    fclose('all') ;
end


end

function update_log_file(it,hv_imp,elapsed_time,xnewtrue,ytrue,Opt,Y,ub,lb)
try
    string1 = '';
    string2 = {};
    string3 = '';
    string4 = '';
    string5 = {};
    string6 = '';
    string7 = '';
    string8 = {};
    for i = 1:Opt.Gen.NoOfInputDim
        string1 = strcat(string1,'%8.4f ');
        string2 = {string2{:},strcat('x',num2str(i))};
        string3 = strcat(string3,'%+8s ');
        string8 = {string8{:},strcat('lambda',num2str(i))};
    end
    string8 = {string8{:},'sigmaf'};
    string8 = {string8{:},'sigman'};
    for i = 1:Opt.Gen.NoOfGPs
        string4     = strcat(string4,'%8.4f ');
        string5     = {string5{:},strcat('f',num2str(i))};
        string6     = strcat(string6,'%+8s ');
        string7     = strcat(string7,'%8d ');
        hypcov      = exp(Opt.GP(i).hyp.cov);
        hypmat(:,i) = [hypcov(1:end-1).*(1./(ub-lb))';hypcov(end)*std(Y(:,i));exp(Opt.GP(i).hyp.lik)*std(Y(:,i))];
    end
    
    TSEMO_log = fopen( 'TSEMO_log.txt', 'a');
    fprintf(TSEMO_log,'\n %s %d \n', 'Algorithm iteration',it);
    fprintf(TSEMO_log,'\n %s %8.4f', 'Predicted hypervolume improvement: ',hv_imp);
    fprintf(TSEMO_log,'\n %s %8.4f \n', 'Time taken: ',elapsed_time);
    fprintf(TSEMO_log,'\n %s', 'Proposed evaluation point(s): ');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),xnewtrue');
    fprintf(TSEMO_log,'\n %s', 'Corresponding observation(s): ');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),ytrue);
    fprintf(TSEMO_log,'\n %s', 'Current hyperparameter values: ');
    fprintf(TSEMO_log,strcat('\n','%+16s',string6,'\n'),'Hyperparameter',string5{:});
    for i = 1:Opt.Gen.NoOfInputDim+2
        fprintf(TSEMO_log,strcat('%+16s',string4,'\n'),string8{i},hypmat(i,:));
    end
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    fclose(TSEMO_log) ;
catch
    warning('There was an error when writing the log file. Maybe the file was currently open in another program. The algorithm continues to run but some parts of the log file are maybe not correct.') ;
    fclose('all') ;
end


end

function final_log_update(Xpareto,Ypareto,X,Y,XParetoGP,YParetoGP,hypf,Opt)
try
    string1 = '';
    string2 = {};
    string3 = '';
    string4 = '';
    string5 = {};
    string6 = '';
    string7 = '';
    string8 = {};
    for i = 1:Opt.Gen.NoOfInputDim
        string1 = strcat(string1,'%8.4f ');
        string2 = {string2{:},strcat('x',num2str(i))};
        string3 = strcat(string3,'%+8s ');
        string8 = {string8{:},strcat('lambda',num2str(i))};
    end
    string8 = {string8{:},'sigmaf'};
    string8 = {string8{:},'sigman'};
    for i = 1:Opt.Gen.NoOfGPs
        string4     = strcat(string4,'%8.4f ');
        string5     = {string5{:},strcat('f',num2str(i))};
        string6     = strcat(string6,'%+8s ');
        string7     = strcat(string7,'%8d ');
    end
    
    TSEMO_log = fopen( 'TSEMO_log.txt', 'a');
    fprintf(TSEMO_log,'\n %s \n', 'Final algorithm output');
    fprintf(TSEMO_log,'\n %s', 'Final input data matrix:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),X');
    fprintf(TSEMO_log,'\n %s', 'Final output data matrix:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),Y');
    fprintf(TSEMO_log,'\n %s', 'Input data matrix of corresponding Pareto front:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),Xpareto');
    fprintf(TSEMO_log,'\n %s', 'Output data matrix of corresponding Pareto front:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),Ypareto');
    fprintf(TSEMO_log,'\n %s', 'Input data matrix of Pareto front of final Gaussian process model:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),XParetoGP');
    fprintf(TSEMO_log,'\n %s', 'Output data matrix of Pareto front of final Gaussian process model:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),YParetoGP');
    fprintf(TSEMO_log,'\n %s', 'Final hyperparameter values: ');
    fprintf(TSEMO_log,strcat('\n','%+16s',string6,'\n'),'Hyperparameter',string5{:});
    for i = 1:Opt.Gen.NoOfInputDim+2
        fprintf(TSEMO_log,strcat('%+16s',string4,'\n'),string8{i},hypf(i,:));
    end
    fclose(TSEMO_log) ;
catch
    warning('There was an error when writing the log file. Maybe the file was currently open in another program. The algorithm continues to run but some parts of the log file are maybe not correct.') ;
    fclose('all') ;
end


end

function v = hypervolumemonte(P,r,N)
% Copyright (c) 2009, Yi Cao
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% HYPERVOUME    Hypervolume indicator as a measure of Pareto front estimate.
%   V = HYPERVOLUME(P,R,N) returns an estimation of the hypervoulme (in
%   percentage) dominated by the approximated Pareto front set P (n by d)
%   and bounded by the reference point R (1 by d). The estimation is doen
%   through N (default is 1000) uniformly distributed random points within
%   the bounded hyper-cuboid.
%
%   V = HYPERVOLUMN(P,R,C) uses the test points specified in C (N by d).
%
% See also: paretofront, paretoGroup

% Version 1.0 by Yi Cao at Cranfield University on 20 April 2008

% Example
%{
% an random exmaple
F=(randn(100,3)+5).^2;
% upper bound of the data set
r=max(F);
% Approximation of Pareto set
P=paretofront(F);
% Hypervolume
v=hypervolume(F(P,:),r,100000);
%}
% https://se.mathworks.com/matlabcentral/fileexchange/19651-hypervolume-indicator

% Check input and output
error(nargchk(2,3,nargin));
error(nargoutchk(0,1,nargout));

P=P*diag(1./r);
[n,d]=size(P);
if nargin<3
    N=1000;
end
if ~isscalar(N)
    C=N;
    N=size(C,1);
else
    C=rand(N,d);
end

fDominated=false(N,1);
lB=min(P);
fcheck=all(bsxfun(@gt, C, lB),2);

for k=1:n
    if any(fcheck)
        f=all(bsxfun(@gt, C(fcheck,:), P(k,:)),2);
        fDominated(fcheck)=f;
        fcheck(fcheck)=~f;
    end
end

v=sum(fDominated)/N;


end

function [f,const] = pareto_objective(x,OptSample)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

Opt.Sample = OptSample;
f = zeros(size(x,1),size(Opt.Sample,2));
for i = 1:size(Opt.Sample,2)
    f(:,i) = Opt.Sample(i).f(x);
end
const = [];


end

function [Sample_pareto,Sample_xpareto,Sample_nadir] = Find_sample_pareto(Opt,it)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

D = Opt.Gen.NoOfInputDim;               % Number of input dimensions
options = nsgaopt();                    % create default options structure
options.popsize = Opt.pop;              % populaion size
options.maxGen  = Opt.Generation;       % max generation
options.numObj = Opt.Gen.NoOfGPs;       % number of objectives
options.numVar = D;                     % number of design variables
options.numCons = 0;                    % number of constraints
options.outputfuns = [];                % saving pop
options.lb = zeros(1,D);                % lower bound of x
options.ub = ones(1,D);                 % upper bound of x
options.objfun = @(x) pareto_objective(x,Opt.Sample);        % objective function handle
options.useParallel = 'no';             % parallel computation is non-essential here
if it > 1
options.initfun = {@(opt,pop) initpop_new(opt,pop,Opt)};    
end
[~,result] = evalc('nsga2(options);');  % begin the optimization!

Sample_xpareto = zeros(Opt.pop,D);
Sample_pareto = zeros(Opt.pop,Opt.Gen.NoOfGPs);
result = result.pops(Opt.Generation,:);

for k = 1:Opt.pop
    Sample_xpareto(k,:) = result(k).var;
    Sample_pareto(k,:) = result(k).obj;
end

for k = 1:Opt.Gen.NoOfGPs
    Sample_nadir(k) = max(Sample_pareto(:,k));
end


end

function pop = initpop_new(opt,pop,Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2020-22-05.

for p = 1:opt.popsize
    pop(p).var = Opt.warmstart_pareto(p,:);
end



end

function [Mean_pareto,Mean_xpareto] = Find_mean_pareto(Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

D = Opt.Gen.NoOfInputDim;               % Number of input dimensions
options = nsgaopt();                    % create default options structure
options.popsize = Opt.pop;              % populaion size
options.maxGen  = Opt.Generation;       % max generation
options.numObj = Opt.Gen.NoOfGPs;       % number of objectives
options.numVar = D;                     % number of design variables
options.numCons = 0;                    % number of constraints
options.outputfuns = [];                % saving pop
options.lb = zeros(1,D);                % lower bound of x
options.ub = ones(1,D);                 % upper bound of x
options.objfun = @(x) pareto_objective(x,Opt.Mean); % objective function handle
options.useParallel = 'no';             % parallel computation is non-essential here
[~,result] = evalc('nsga2(options);');  % begin the optimization!

Mean_xpareto = zeros(Opt.pop,D);
Mean_pareto  = zeros(Opt.pop,Opt.Gen.NoOfGPs);
result       = result.pops(Opt.Generation,:);

for k = 1:Opt.pop
    Mean_xpareto(k,:) = result(k).var;
    Mean_pareto(k,:) = result(k).obj;
end



end

function hv = hypervolume_2D(Yfront,r)
% The mex-file used follows the description of:
% 'M. Emmerich, K. Yang, A. Deutz, H. Wang and C. M. Fonseca. A
% Multicriteria Generalization of Bayesian Global Optimization.'
% Group website: http://liacs.leidenuniv.nl/~csmoda/index.php?page=code

AYfront = remove_points_above_reference(Yfront,r);
if isempty(AYfront)
    hv = 0;
else
    normvec = min(AYfront,[],1);
    A =  AYfront-repmat(normvec,size(AYfront,1),1);
    A =  A * diag(1./(r-normvec));
    A = -A + ones(size(A));
    A = sortrows(A,2);
    hyp_percentage = hypervolume2D(A,[0,0]);
    hv = prod(r-normvec)*hyp_percentage;
end


end

function hv = hypervolume_3D(Yfront,r)
% The mex-file used follows the description of:
% 'K. Yang, M. Emmerich, A. Deutz and C. M. Fonseca. A
% Computing 3-D Expected Hypervolume Improvement and Related Integrals in
% Asymptotically Optimal Time.'
% Group website: http://liacs.leidenuniv.nl/~csmoda/index.php?page=code

AYfront = remove_points_above_reference(Yfront,r);
if isempty(AYfront)
    hv = 0;
else
    normvec = min(AYfront,[],1);
    A = AYfront-repmat(normvec,size(AYfront,1),1) ;
    A = A *diag(1./(r-normvec));
    A = -A + ones(size(A)) ;
    A = sortrows(A,3);
    hyp_percentage = hypervolume3D(A,[0,0,0],[1,1,1]);
    hv = prod(r-normvec)*hyp_percentage;
end


end

function [index,hv_imp] = hypervolume_improvement_index(Ynew,Sample_nadir,Sample_pareto,Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

r = Sample_nadir + 0.01*(max(Sample_pareto)-min(Sample_pareto));
index = [];

for i = 1 : Opt.NoOfBachSequential
    
    Yfront = Ynew(paretofront(Ynew),:);
    
    if size(Ynew,2) == 2
        hvY = hypervolume_2D(Yfront,r);
        
        for k = 1:size(Sample_pareto,1)
            A = [Ynew;Sample_pareto(k,:)];
            Afront = A(paretofront(A),:);
            hv = hypervolume_2D(Afront,r);
            hv_improvement(k) = hv-hvY;
        end
        
    elseif size(Ynew,2) == 3
        hvY = hypervolume_3D(Yfront,r);
        
        for k = 1:size(Sample_pareto,1)
            A = [Ynew;Sample_pareto(k,:)];
            Afront = A(paretofront(A),:);
            hv = hypervolume_3D(Afront,r);
            hv_improvement(k) = hv-hvY;
        end
        
    else
        AYfront = remove_points_above_reference(Yfront,r);
        normvec = min(AYfront,[],1);
        hyp_percentage = hypervolumemonte(AYfront-repmat(normvec,size(AYfront,1),1),r-normvec,3000);
        hvY = prod(r-normvec)*hyp_percentage;
        
        hv_improvement = zeros(size(Sample_pareto,1),1);
        for k = 1:size(Sample_pareto,1)
            B = [Ynew;Sample_pareto(k,:)];
            Bfront = B(paretofront(B),:);
            ABfront = remove_points_above_reference(Bfront,r);
            if isempty(ABfront)
                hv_improvement(k) = 0;
            else
                normvec = min(ABfront,[],1);
                hyp_percentage = hypervolumemonte(ABfront-repmat(normvec,size(ABfront,1),1),r-normvec,10000);
                hv = prod(r-normvec)*hyp_percentage;
                hv_improvement(k) = hv-hvY;
            end
        end
    end
    
    if i == 1
        hvY0 = hvY;
    end
    
    [~,Currentindex] = max(hv_improvement);
    Ynew = [Ynew;Sample_pareto(Currentindex,:)];
    index = [index;Currentindex];
end
hv_imp = hv_improvement(index(end))+hvY-hvY0;


end

function A = remove_points_above_reference(Afront,r)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

[A,~] = sortrows(Afront);
for p = 1:size(Afront,2)
    A = A(A(:,p)<=r(p),:);
end


end

function [f, varf] = mean_sample_piGP(Xnew, Yresidual, OptGP, theta_prior, lb, ub, cA_0, obj_idx)
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
    m_tr = piGP_prior_mean(Xnew, lb, ub, cA_0, theta_prior, dummy);
else
    m_tr = piGP_prior_mean(Xnew, lb, ub, cA_0, dummy, theta_prior);
end
m_mu = mean(m_tr); m_sig = max(std(m_tr), 1e-8);

phi_fn = @(x) sqrt(2*sf2/nSp)*cos(W*x' + repmat(b,1,size(x,1)));

f = @(x) mean_eval(x, mu_th, phi_fn, theta_prior, lb, ub, cA_0, obj_idx, m_mu, m_sig);
varf = @(x) sn2 + sn2*phi_fn(x)'*invA*phi_fn(x);
end

function vals = mean_eval(x, mu_th, phi_fn, theta_prior, lb, ub, cA_0, obj_idx, m_mu, m_sig)
dummy = zeros(4,1);
if obj_idx == 1
    m_t = piGP_prior_mean(x, lb, ub, cA_0, theta_prior, dummy);
else
    m_t = piGP_prior_mean(x, lb, ub, cA_0, dummy, theta_prior);
end
vals = (mu_th' * phi_fn(x))' + (m_t - m_mu)/m_sig;
end
