% TSEMO_piGP_Multiple_Runs_Different_Meval_and_DatasetSize.m
%
% Physics-informed TSEMO (piGP) — multi-batch parallel runner.
%
% Mirrors TSEMO_Multiple_Runs_Different_Meval_and_DatasetSize.m exactly,
% but calls TSEMO_V4_1a_OPC_piGP instead of TSEMO_V4_1a_OPC and saves the
% extra output theta_history (learned prior parameters per iteration).
%
% Objectives: STY and Selectivity (SNAr / Hone-Taylor reaction system)
%
% Required folders on MATLAB path (add via addpath before running):
%   Direct                  - DIRECT global optimiser
%   Hone_Taylor_Reaction    - f_PFR_2nd_order, f_PFR_kinetics
%   Mex_files               - hypervolume2D/3D MEX binaries, invChol, etc.
%   NGPM_v1.4               - nsga2, nsgaopt, paretofront
%
% Required .m files (must be findable on path):
%   TSEMO_V4_1a_OPC_piGP.m
%   piGP_prior_mean.m
%   NLikelihood_piGP.m
%   TrainingOfGP_piGP.m
%   posterior_sample_piGP.m
%   TSEMO_options.m

clear all
close all
clc

%% ========================================================================
%% CONFIGURATION - Define your batches here
%% ========================================================================

% Each row: [maxeval, number_of_runs, initial_dataset_size]
batches = [
    % maxeval, runs, initial_points
    10,  50,  4;   % Batch 1: 2 runs, maxeval=10, 4 initial LHS points
    20,  50,  4;   % Batch 1: 2 runs, maxeval=10, 4 initial LHS points
    50,  50,  4;   % Batch 1: 2 runs, maxeval=10, 4 initial LHS points
    80,  50,  4;   % Batch 1: 2 runs, maxeval=10, 4 initial LHS points
    100,  50,  4;   % Batch 1: 2 runs, maxeval=10, 4 initial LHS points
];

% Number of parallel workers.  [] = use all available cores.
num_workers = [];

%% ========================================================================
%% CREATE MASTER RESULTS FOLDER
%% ========================================================================

timestamp           = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
master_folder_name  = ['TSEMO_piGP_MultiMeval_MultiDataset_', timestamp];
current_dir         = pwd;
master_folder       = fullfile(current_dir, master_folder_name);
mkdir(master_folder);

fprintf('========================================================================\n');
fprintf('TSEMO piGP OPTIMISATION - VARYING MAXEVAL AND INITIAL DATASET SIZE\n');
fprintf('========================================================================\n');
fprintf('Number of batches: %d\n', size(batches, 1));
fprintf('Current directory: %s\n', current_dir);
fprintf('Master results folder: %s\n', master_folder);
fprintf('\nBatch configuration:\n');
fprintf('%6s %10s %8s %15s\n', 'Batch', 'MaxEval', 'Runs', 'Initial Points');
fprintf('------------------------------------------------------\n');
for i = 1:size(batches, 1)
    fprintf('%6d %10d %8d %15d\n', i, batches(i,1), batches(i,2), batches(i,3));
end
fprintf('========================================================================\n\n');

%% ========================================================================
%% CHECK AND START PARALLEL POOL
%% ========================================================================

use_parallel = true;

if license('test', 'Distrib_Computing_Toolbox')
    fprintf('Parallel Computing Toolbox license detected - Setting up parallel pool...\n');

    try
        poolobj = gcp('nocreate');

        if isempty(poolobj)
            if isempty(num_workers)
                poolobj = parpool();
            else
                poolobj = parpool(num_workers);
            end
            fprintf('Parallel pool started with %d workers\n', poolobj.NumWorkers);
        else
            fprintf('Using existing parallel pool with %d workers\n', poolobj.NumWorkers);
        end

        % Attach all .m files in the current directory to the workers
        fprintf('Adding required files to parallel workers...\n');
        try
            current_files = dir('*.m');
            files_to_add  = {};
            for i = 1:length(current_files)
                files_to_add{end+1} = fullfile(pwd, current_files(i).name); %#ok<SAGROW>
            end
            if ~isempty(files_to_add)
                addAttachedFiles(poolobj, files_to_add);
                fprintf('Successfully added %d files to parallel workers\n', length(files_to_add));
            end
        catch ME_attach
            fprintf('Warning: Could not attach all files. Error: %s\n', ME_attach.message);
        end

        fprintf('========================================================================\n\n');

    catch ME
        fprintf('WARNING: Could not start parallel pool.\n');
        fprintf('Error message: %s\n', ME.message);
        fprintf('\nFalling back to SEQUENTIAL execution...\n');
        fprintf('========================================================================\n\n');
        use_parallel = false;
    end

else
    fprintf('WARNING: Parallel Computing Toolbox not available.\n');
    fprintf('Running in SEQUENTIAL mode instead.\n');
    fprintf('========================================================================\n\n');
    use_parallel = false;
end

%% ========================================================================
%% RUN ALL BATCHES
%% ========================================================================

all_batch_results = cell(size(batches, 1), 1);
total_start_time  = tic;

for batch_idx = 1:size(batches, 1)

    current_maxeval      = batches(batch_idx, 1);
    current_num_runs     = batches(batch_idx, 2);
    current_dataset_size = batches(batch_idx, 3);

    fprintf('\n\n');
    fprintf('************************************************************************\n');
    fprintf('                    STARTING BATCH %d of %d\n', batch_idx, size(batches, 1));
    fprintf('          maxeval=%d, runs=%d, initial points=%d\n', ...
        current_maxeval, current_num_runs, current_dataset_size);
    fprintf('************************************************************************\n\n');

    % Create batch folder (name encodes all three parameters)
    batch_folder = fullfile(master_folder, ...
        sprintf('Batch_%d_Meval_%d_InitPts_%d', batch_idx, current_maxeval, current_dataset_size));
    mkdir(batch_folder);

    batch_start_time = tic;

    if use_parallel
        fprintf('Using PARALLEL execution with %d workers\n\n', poolobj.NumWorkers);
        parfor run_num = 1:current_num_runs
            run_single_piGP_optimization(run_num, batch_folder, current_maxeval, current_dataset_size);
        end
    else
        fprintf('Using SEQUENTIAL execution (one run at a time)\n\n');
        for run_num = 1:current_num_runs
            run_single_piGP_optimization(run_num, batch_folder, current_maxeval, current_dataset_size);
        end
    end

    batch_elapsed_time = toc(batch_start_time);

    % Store batch metadata
    all_batch_results{batch_idx}.maxeval      = current_maxeval;
    all_batch_results{batch_idx}.num_runs     = current_num_runs;
    all_batch_results{batch_idx}.dataset_size = current_dataset_size;
    all_batch_results{batch_idx}.elapsed_time = batch_elapsed_time;
    all_batch_results{batch_idx}.folder       = batch_folder;

    fprintf('\n');
    fprintf('************************************************************************\n');
    fprintf('                    BATCH %d COMPLETED\n', batch_idx);
    fprintf('          Time: %.2f seconds (%.2f minutes)\n', batch_elapsed_time, batch_elapsed_time/60);
    fprintf('************************************************************************\n');

end

total_elapsed_time = toc(total_start_time);

%% ========================================================================
%% FINAL SUMMARY
%% ========================================================================

fprintf('\n\n');
fprintf('========================================================================\n');
fprintf('              ALL BATCHES COMPLETED SUCCESSFULLY\n');
fprintf('========================================================================\n\n');

summary_table = table();
for i = 1:length(all_batch_results)
    summary_table.Batch(i)               = i;
    summary_table.MaxEval(i)             = all_batch_results{i}.maxeval;
    summary_table.InitialPoints(i)       = all_batch_results{i}.dataset_size;
    summary_table.NumRuns(i)             = all_batch_results{i}.num_runs;
    summary_table.Time_min(i)            = all_batch_results{i}.elapsed_time / 60;
    summary_table.AvgTime_per_run_sec(i) = all_batch_results{i}.elapsed_time / all_batch_results{i}.num_runs;
end

disp(summary_table);

fprintf('\nTotal time for all batches: %.2f seconds (%.2f minutes, %.2f hours)\n', ...
    total_elapsed_time, total_elapsed_time/60, total_elapsed_time/3600);
fprintf('Results location: %s\n', master_folder);
fprintf('========================================================================\n');

% Save master summary
writetable(summary_table, fullfile(master_folder, 'Batch_Summary.csv'));
save(fullfile(master_folder, 'All_Batch_Results.mat'), ...
    'all_batch_results', 'batches', 'summary_table');

fprintf('\nSummary saved to: %s\n', fullfile(master_folder, 'Batch_Summary.csv'));

fprintf('\n');
fprintf('========================================================================\n');
fprintf('                        FOLDER STRUCTURE\n');
fprintf('========================================================================\n');
fprintf('%s/\n', master_folder_name);
for i = 1:size(batches, 1)
    fprintf('├── Batch_%d_Meval_%d_InitPts_%d/\n', i, batches(i,1), batches(i,3));
    fprintf('│   ├── Run_1/\n');
    fprintf('│   ├── Run_2/\n');
    fprintf('│   └── ...\n');
end
fprintf('├── Batch_Summary.csv\n');
fprintf('└── All_Batch_Results.mat\n');
fprintf('========================================================================\n');


%% ========================================================================
%% HELPER FUNCTION — run one physics-informed TSEMO optimisation
%% ========================================================================

function run_single_piGP_optimization(run_num, batch_folder, maxeval, dataset_size)
% Runs a single piGP-TSEMO trial.  Designed to execute on a parallel worker.

    fprintf('[Worker %d] Starting Run %d (maxeval=%d, %d initial pts)\n', ...
        getCurrentWorker_piGP(), run_num, maxeval, dataset_size);

    run_start_time = tic;

    %% ------------------------------------------------------------------
    %% Create run folder (with retry logic for parallel file-system races)
    %% ------------------------------------------------------------------
    run_folder  = fullfile(batch_folder, sprintf('Run_%d', run_num));
    max_retries = 5;
    retry_count = 0;
    folder_created = false;

    while ~folder_created && retry_count < max_retries
        try
            if ~exist(run_folder, 'dir')
                mkdir(run_folder);
            end
            if exist(run_folder, 'dir')
                folder_created = true;
            else
                retry_count = retry_count + 1;
                pause(0.5);
            end
        catch ME_mkdir
            retry_count = retry_count + 1;
            fprintf('[Worker %d] Warning: Attempt %d to create folder failed: %s\n', ...
                getCurrentWorker_piGP(), retry_count, ME_mkdir.message);
            pause(1);
        end
    end

    if ~folder_created
        error('Failed to create folder %s after %d attempts', run_folder, max_retries);
    end

    original_dir = pwd;
    try
        cd(run_folder);
    catch ME_cd
        error('Cannot change to folder %s: %s', run_folder, ME_cd.message);
    end

    try
        %% ------------------------------------------------------------------
        %% Reactor + design specification
        %% ------------------------------------------------------------------
        V_R  = 2;    % reactor volume [mL]
        c1_0 = 0.2;  % initial DFNB concentration [mol/L]
        c2_0 = 0.2;  % initial pyrrolidine concentration [mol/L]

        tau_min = 0.25;  tau_max = 2;   % residence time [min]
        T_min   = 30;    T_max   = 50;  % temperature [degC]

        %% ------------------------------------------------------------------
        %% Step 1: Specify problem
        %% ------------------------------------------------------------------
        no_outputs = 2;
        no_inputs  = 2;
        lb = [tau_min, T_min];
        ub = [tau_max, T_max];

        f = @(x) f_PFR_2nd_order(x, c1_0, c2_0);

        %% ------------------------------------------------------------------
        %% Step 2: Generate LHS initial dataset
        %% ------------------------------------------------------------------
        X = lhsdesign(dataset_size, no_inputs);
        X = sortrows(X, 2);               % sort by temperature column
        Y = zeros(dataset_size, no_outputs);

        %% ------------------------------------------------------------------
        %% Step 3: Evaluate initial experiments
        %% ------------------------------------------------------------------
        % Molecular masses [g/mol]
        M_1 = 159.09;  M_2 = 87.12;  M_F = 18.99;  M_H = 1.00;
        M_3 = M_1 + M_2 - M_F - M_H;   % ortho-product MW

        for k = 1:size(X, 1)
            X(k,:)   = X(k,:) .* (ub - lb) + lb;
            tau      = X(k,1);
            tau_sec  = tau * 60;
            T        = round(X(k,2), 2);

            Q_tot = V_R / tau;
            Q_i   = Q_tot / 2;          %#ok<NASGU>  (kept for OPC-UA use)

            c_out = f_PFR_2nd_order(X(k,:), c1_0, c2_0);

            Y(k,1) = -log(M_3 .* c_out(3) ./ tau_sec);              % STY
            Y(k,2) = -log(c_out(3) ./ (c_out(3)+c_out(4)+c_out(5)));% Selectivity
        end

        %% ------------------------------------------------------------------
        %% Step 4: Run physics-informed TSEMO
        %% ------------------------------------------------------------------
        opt = TSEMO_options;
        opt.maxeval           = maxeval;
        opt.NoOfBachSequential = 1;

        t_meas = 3;  % MIR measurement time [sec] — unused in simulation mode

        % TSEMO_V4_1a_OPC_piGP returns the extra output theta_history
        % compared to the standard TSEMO_V4_1a_OPC:
        %   theta_history [n_iter x 4 x 2]
        %     (:,:,1) = STY prior params  [log_kref, Ea(J/mol), alpha, beta]
        %     (:,:,2) = SEL prior params  [log_kref, Ea(J/mol), alpha, beta]
        [Xpareto, Ypareto, X, Y, XparetoGP, YparetoGP, YparetoGPstd, hypf, theta_history] = ...
            TSEMO_V4_1a_OPC_piGP(f, X, Y, lb, ub, opt, c1_0, c2_0, V_R, t_meas);

        %% ------------------------------------------------------------------
        %% Step 5: Rename output files with descriptive suffix
        %% ------------------------------------------------------------------
        run_elapsed_time = toc(run_start_time);
        current_folder   = pwd;
        file_suffix      = sprintf('Run_%d_Meval_%d_InitPts_%d', run_num, maxeval, dataset_size);

        % TSEMO log
        log_src = fullfile(current_folder, 'TSEMO_log.txt');
        if exist(log_src, 'file') == 2
            try
                movefile(log_src, fullfile(current_folder, ['TSEMO_log_' file_suffix '.txt']));
            catch ME_log
                fprintf('[Worker %d] Warning: Could not rename TSEMO_log.txt: %s\n', ...
                    getCurrentWorker_piGP(), ME_log.message);
            end
        end

        % GP Pareto front — inputs CSV
        csv_in_src = fullfile(current_folder, 'GP_Pareto_Front_Inputs.csv');
        if exist(csv_in_src, 'file') == 2
            try
                movefile(csv_in_src, fullfile(current_folder, ['GP_Pareto_Front_Inputs_' file_suffix '.csv']));
            catch; end
        end

        % GP Pareto front — outputs CSV
        csv_out_src = fullfile(current_folder, 'GP_Pareto_Front_Outputs.csv');
        if exist(csv_out_src, 'file') == 2
            try
                movefile(csv_out_src, fullfile(current_folder, ['GP_Pareto_Front_Outputs_' file_suffix '.csv']));
            catch; end
        end

        %% ------------------------------------------------------------------
        %% Step 6: Save workspace
        %%   theta_history is the key extra variable vs the standard runner.
        %%   It records how the physics prior parameters evolved each iteration.
        %% ------------------------------------------------------------------
        save(['Workspace_' file_suffix '.mat'], ...
            'Xpareto', 'Ypareto', 'X', 'Y', ...
            'XparetoGP', 'YparetoGP', 'YparetoGPstd', 'hypf', ...
            'theta_history', ...              % <-- piGP-specific
            'lb', 'ub', 'c1_0', 'c2_0', 'V_R', ...
            'run_elapsed_time', 'maxeval', 'dataset_size');

        %% ------------------------------------------------------------------
        %% Step 7: Save theta_history as CSV (one table per objective)
        %%   Rows = iterations, columns = [log_kref, Ea_kJmol, alpha, beta]
        %% ------------------------------------------------------------------
        n_iter_done = size(theta_history, 1);
        iter_col    = (1:n_iter_done)';
        col_names   = {'Iteration', 'log_kref', 'Ea_kJmol', 'alpha', 'beta'};

        % STY prior parameter history
        T_STY = array2table([iter_col, theta_history(:,1,1), ...
                             theta_history(:,2,1)/1000, ...   % Ea: J/mol -> kJ/mol
                             theta_history(:,3,1), ...
                             theta_history(:,4,1)], ...
                            'VariableNames', col_names);
        writetable(T_STY, ['PriorParams_STY_' file_suffix '.csv']);

        % SEL prior parameter history
        T_SEL = array2table([iter_col, theta_history(:,1,2), ...
                             theta_history(:,2,2)/1000, ...
                             theta_history(:,3,2), ...
                             theta_history(:,4,2)], ...
                            'VariableNames', col_names);
        writetable(T_SEL, ['PriorParams_SEL_' file_suffix '.csv']);

        %% ------------------------------------------------------------------
        %% Step 8: Create Pareto front visualisation
        %% ------------------------------------------------------------------

        % Back-transform from log-space to physical units
        Y1            = exp(-Y(:,1));
        Y2            = exp(-Y(:,2));
        Ypareto1      = exp(-Ypareto(:,1));
        Ypareto2      = exp(-Ypareto(:,2));
        YparetoGP1    = exp(-YparetoGP(:,1));
        YparetoGP2    = exp(-YparetoGP(:,2));
        YparetoGPstd1 = exp(-YparetoGP(:,1)) - exp(-(YparetoGP(:,1) + YparetoGPstd(:,1)));
        YparetoGPstd2 = exp(-YparetoGP(:,2)) - exp(-(YparetoGP(:,2) + YparetoGPstd(:,2)));

        fig = figure('Visible', 'off');
        hold on
        plot(Y1(1:dataset_size),      Y2(1:dataset_size), ...
            '.', 'MarkerSize', 14, 'color', [0.8500 0.3250 0.0980])
        plot(Y1(dataset_size+1:end),  Y2(dataset_size+1:end), ...
            'x', 'MarkerSize', 8, 'LineWidth', 2, 'color', [0.8500 0.3250 0.0980])
        plot(Ypareto1, Ypareto2, ...
            'O', 'MarkerSize', 8, 'LineWidth', 2, 'color', "#EDB120")
        errorbar(YparetoGP1, YparetoGP2, ...
            YparetoGPstd2, YparetoGPstd2, YparetoGPstd1, YparetoGPstd1, ...
            '.', 'MarkerSize', 8, 'color', "#7E2F8E", 'LineWidth', 0.5)
        legend('Initial LHS', 'Algorithm', 'Pareto front', ...
            'GP Pareto front (piGP)', 'Location', 'Northeast')
        grid on
        xlabel('STY (g \cdot L^{-1} \cdot s^{-1})')
        ylabel('Selectivity (-)')
        title(sprintf('Run %d  —  maxeval=%d  —  %d Init Pts  (piGP)', ...
            run_num, maxeval, dataset_size))

        saveas(fig, ['Pareto_Plot_' file_suffix '.png']);
        saveas(fig, ['Pareto_Plot_' file_suffix '.fig']);
        close(fig);

        %% ------------------------------------------------------------------
        %% Step 9: Prior parameter evolution diagnostic plot
        %% ------------------------------------------------------------------
        param_labels = {'log k_{ref}', 'E_a (kJ/mol)', '\alpha', '\beta'};
        colors_p     = {[0.2 0.5 0.8], [0.8 0.3 0.1]};

        fig2 = figure('Visible', 'off', 'Name', 'Prior param evolution');
        for p = 1:4
            subplot(2, 2, p); hold on; grid on
            for obj = 1:2
                vals = theta_history(:, p, obj);
                if p == 2, vals = vals / 1000; end   % J/mol -> kJ/mol for display
                plot(1:n_iter_done, vals, 'o-', ...
                    'Color', colors_p{obj}, ...
                    'MarkerFaceColor', colors_p{obj}, ...
                    'LineWidth', 1.5, 'MarkerSize', 5)
            end
            xlabel('Iteration');
            ylabel(param_labels{p});
            title(param_labels{p})
            if p == 1
                legend({'STY', 'SEL'}, 'Location', 'best');
            end
        end
        sgtitle(sprintf('Prior params — Run %d, maxeval=%d, %d init pts', ...
            run_num, maxeval, dataset_size))

        saveas(fig2, ['PriorParam_Evolution_' file_suffix '.png']);
        saveas(fig2, ['PriorParam_Evolution_' file_suffix '.fig']);
        close(fig2);

        fprintf('[Worker %d] Run %d (maxeval=%d, %d pts) done in %.2f min\n', ...
            getCurrentWorker_piGP(), run_num, maxeval, dataset_size, run_elapsed_time/60);

    catch ME
        fprintf('[Worker %d] ERROR in Run %d: %s\n', ...
            getCurrentWorker_piGP(), run_num, ME.message);
        try; cd(original_dir); catch; end
        rethrow(ME);
    end

    % Always return to original directory
    try; cd(original_dir); catch; end

end  % run_single_piGP_optimization


%% ========================================================================
%% HELPER FUNCTION — get current parallel worker ID (0 = main thread)
%% ========================================================================

function worker_id = getCurrentWorker_piGP()
    try
        t         = getCurrentTask();
        worker_id = t.ID;
    catch
        worker_id = 0;
    end
end
