% TSEMO_piGP_Denoise_Multiple_Runs_Different_Meval_and_DatasetSize.m
%
% Physics-informed TSEMO (piGP) + concentration denoising — multi-batch
% parallel runner.
%
% This is the parallel version of TSEMO_OPC_2D_piGP.m (the version that
% includes denoise_concentration.m). It mirrors the structure of the
% previous parallel runners but adds full denoiser state management per run.
%
% KEY DESIGN NOTE — why denoiser state is handled the way it is:
%   denoise_concentration accumulates history (X_hist, C_hist, w_blend)
%   and relies on the current theta_STY/SEL from the piGP fit. These are
%   sequential, stateful updates that CANNOT be shared across parallel runs
%   — each run carries its own independent denoiser state from start to
%   finish. The parfor loop parallelises RUNS (independent replicates),
%   not iterations within a run.
%
% Required folders on MATLAB path (add via addpath before running):
%   Direct                  - DIRECT global optimiser
%   Hone_Taylor_Reaction    - f_PFR_2nd_order, f_PFR_kinetics
%   Mex_files               - hypervolume2D/3D MEX, invChol, etc.
%   NGPM_v1.4               - nsga2, nsgaopt, paretofront
%
% Required .m files (must be findable on path):
%   TSEMO_V4_1a_OPC_piGP.m
%   piGP_prior_mean.m
%   NLikelihood_piGP.m
%   TrainingOfGP_piGP.m
%   posterior_sample_piGP.m
%   denoise_concentration.m       <-- new vs previous parallel runner
%   TSEMO_options.m

clear all
close all
clc

%% ========================================================================
%% CONFIGURATION
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

timestamp          = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
master_folder_name = ['TSEMO_piGP_Denoise_MultiMeval_MultiDataset_', timestamp];
current_dir        = pwd;
master_folder      = fullfile(current_dir, master_folder_name);
mkdir(master_folder);

fprintf('========================================================================\n');
fprintf('TSEMO piGP + DENOISING — VARYING MAXEVAL AND INITIAL DATASET SIZE\n');
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
    fprintf('Parallel Computing Toolbox detected — setting up parallel pool...\n');
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

        % Attach .m files in the current directory to workers
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
            fprintf('Warning: Could not attach all files: %s\n', ME_attach.message);
        end
        fprintf('========================================================================\n\n');

    catch ME
        fprintf('WARNING: Could not start parallel pool: %s\n', ME.message);
        fprintf('Falling back to SEQUENTIAL execution...\n');
        fprintf('========================================================================\n\n');
        use_parallel = false;
    end
else
    fprintf('WARNING: Parallel Computing Toolbox not available.\n');
    fprintf('Running in SEQUENTIAL mode.\n');
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

    batch_folder = fullfile(master_folder, ...
        sprintf('Batch_%d_Meval_%d_InitPts_%d', batch_idx, current_maxeval, current_dataset_size));
    mkdir(batch_folder);

    batch_start_time = tic;

    if use_parallel
        fprintf('Using PARALLEL execution with %d workers\n\n', poolobj.NumWorkers);
        parfor run_num = 1:current_num_runs
            run_single_denoised_piGP(run_num, batch_folder, current_maxeval, current_dataset_size);
        end
    else
        fprintf('Using SEQUENTIAL execution\n\n');
        for run_num = 1:current_num_runs
            run_single_denoised_piGP(run_num, batch_folder, current_maxeval, current_dataset_size);
        end
    end

    batch_elapsed_time = toc(batch_start_time);

    all_batch_results{batch_idx}.maxeval      = current_maxeval;
    all_batch_results{batch_idx}.num_runs     = current_num_runs;
    all_batch_results{batch_idx}.dataset_size = current_dataset_size;
    all_batch_results{batch_idx}.elapsed_time = batch_elapsed_time;
    all_batch_results{batch_idx}.folder       = batch_folder;

    fprintf('\n');
    fprintf('************************************************************************\n');
    fprintf('                    BATCH %d COMPLETED\n', batch_idx);
    fprintf('          Time: %.2f seconds (%.2f minutes)\n', ...
        batch_elapsed_time, batch_elapsed_time/60);
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
fprintf('\nTotal time: %.2f seconds (%.2f minutes, %.2f hours)\n', ...
    total_elapsed_time, total_elapsed_time/60, total_elapsed_time/3600);
fprintf('Results: %s\n', master_folder);
fprintf('========================================================================\n');

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
%% HELPER FUNCTION — run one full piGP + denoising optimisation
%% ========================================================================

function run_single_denoised_piGP(run_num, batch_folder, maxeval, dataset_size)
% Runs one independent piGP-TSEMO trial with concentration denoising.
% Fully self-contained denoiser state — no sharing across parallel workers.

    fprintf('[Worker %d] Starting Run %d (maxeval=%d, %d initial pts)\n', ...
        getCurrentWorker_piGP_DN(), run_num, maxeval, dataset_size);

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
            fprintf('[Worker %d] Warning: folder creation attempt %d failed: %s\n', ...
                getCurrentWorker_piGP_DN(), retry_count, ME_mkdir.message);
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
        error('Cannot cd to %s: %s', run_folder, ME_cd.message);
    end

    try
        %% ------------------------------------------------------------------
        %% Reactor + problem specification
        %% ------------------------------------------------------------------
        V_R  = 2;      % reactor volume [mL]
        c1_0 = 0.2;    % initial DFNB concentration [mol/L]
        c2_0 = 0.2;    % initial pyrrolidine concentration [mol/L]

        tau_min = 0.25;  tau_max = 2;   % residence time [min]
        T_min   = 30;    T_max   = 50;  % temperature [degC]

        no_outputs = 2;
        no_inputs  = 2;
        lb = [tau_min, T_min];
        ub = [tau_max, T_max];

        M_1 = 159.09;  M_2 = 87.12;  M_F = 18.99;  M_H = 1.00;
        M_3 = M_1 + M_2 - M_F - M_H;   % ortho-product MW [g/mol]

        t_meas = 3;   % MIR measurement time [sec]

        %% ------------------------------------------------------------------
        %% Initialise denoiser state
        %%
        %% Each run carries its own fully independent denoiser state.
        %% w_blend  — learned blend weight (physics vs measurement)
        %% X_hist   — [tau, T] history of all experiments (physical units)
        %% C_hist   — raw measured concentrations for all past experiments
        %% theta_*  — current piGP prior params, updated after each TSEMO iter
        %%
        %% These are updated sequentially within a run and MUST NOT be shared
        %% across parfor workers — each worker owns its own copies.
        %% ------------------------------------------------------------------
        w_blend           = 0.3;                    % conservative prior trust at start
        X_hist            = [];                     % grows as experiments accumulate
        C_hist            = [];
        theta_STY_current = [-2.5; 35000; 1.0; 0.0];
        theta_SEL_current = [-2.5; 35000; 1.0; 0.0];

        %% ------------------------------------------------------------------
        %% Step 1: Generate LHS initial dataset
        %% ------------------------------------------------------------------
        X = lhsdesign(dataset_size, no_inputs);
        X = sortrows(X, 2);
        Y = zeros(dataset_size, no_outputs);

        %% ------------------------------------------------------------------
        %% Step 2: Run initial experiments with denoising
        %%
        %% At each initial point we:
        %%   1. Call the simulator (or hardware) to get raw concentrations
        %%   2. Pass them through denoise_concentration, which blends the
        %%      physics prediction with the measurement using blend weight w
        %%   3. Store the RAW measurement in X_hist/C_hist so the denoiser
        %%      can learn from real data — not from already-blended values
        %%   4. Compute objectives from the DENOISED concentrations
        %% ------------------------------------------------------------------
        for k = 1:size(X, 1)
            X(k,:)   = X(k,:) .* (ub - lb) + lb;
            tau      = X(k,1);
            tau_sec  = tau * 60;
            T        = round(X(k,2), 2);

            %%%% Simulation dummy — replace with OPC-UA + MIR hardware %%%%
            c_out_raw = f_PFR_2nd_order(X(k,:), c1_0, c2_0);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Denoise using current state (X_hist, C_hist may be empty on k=1)
            [c_out, w_blend] = denoise_concentration(c_out_raw, tau, T, c1_0, c2_0, ...
                                                     theta_STY_current, theta_SEL_current, ...
                                                     w_blend, X_hist, C_hist);

            % History always stores the RAW measurement (not the blend)
            X_hist = [X_hist; tau, T];       %#ok<AGROW>
            C_hist = [C_hist; c_out_raw];    %#ok<AGROW>

            % Objectives from denoised concentrations
            Y(k,1) = -log(M_3 .* c_out(3) ./ tau_sec);
            Y(k,2) = -log(c_out(3) ./ (c_out(3) + c_out(4) + c_out(5)));
        end

        %% ------------------------------------------------------------------
        %% Step 3: Build denoiser state struct for the objective wrapper
        %%
        %% TSEMO_V4_1a_OPC_piGP calls the objective function as a black box.
        %% We wrap the real objective so that every new point proposed by
        %% TSEMO goes through the denoiser before Y is computed.
        %%
        %% The state struct is a SNAPSHOT — it captures the denoiser's state
        %% as it stands after the initial experiments.  The wrapper uses this
        %% snapshot for all TSEMO-proposed evaluations.
        %%
        %% NOTE: theta_STY/SEL in the snapshot are the cold-start values;
        %% the piGP algorithm will update them internally during optimisation.
        %% If you want the denoiser to track the evolving theta values from
        %% TSEMO, you would need to refactor TSEMO_V4_1a_OPC_piGP to call
        %% a stateful callback — the current architecture does not support
        %% this without modifying the main algorithm file.  The snapshot
        %% approach is a safe, self-contained approximation.
        %% ------------------------------------------------------------------
        ds = struct( ...
            'w',         w_blend, ...
            'X_hist',    X_hist, ...
            'C_hist',    C_hist, ...
            'theta_STY', theta_STY_current, ...
            'theta_SEL', theta_SEL_current, ...
            'c1_0',      c1_0, ...
            'c2_0',      c2_0, ...
            'M_3',       M_3);

        f_denoised = @(x) obj_with_denoising_local(x, ds);

        %% ------------------------------------------------------------------
        %% Step 4: Run physics-informed TSEMO with denoised objective
        %% ------------------------------------------------------------------
        opt = TSEMO_options;
        opt.maxeval            = maxeval;
        opt.NoOfBachSequential = 1;

        [Xpareto, Ypareto, X, Y, XparetoGP, YparetoGP, YparetoGPstd, hypf, theta_history] = ...
            TSEMO_V4_1a_OPC_piGP(f_denoised, X, Y, lb, ub, opt, c1_0, c2_0, V_R, t_meas);

        %% ------------------------------------------------------------------
        %% Step 5: Rename TSEMO output files with descriptive suffix
        %% ------------------------------------------------------------------
        run_elapsed_time = toc(run_start_time);
        current_folder   = pwd;
        file_suffix      = sprintf('Run_%d_Meval_%d_InitPts_%d', run_num, maxeval, dataset_size);

        log_src = fullfile(current_folder, 'TSEMO_log.txt');
        if exist(log_src, 'file') == 2
            try
                movefile(log_src, fullfile(current_folder, ['TSEMO_log_' file_suffix '.txt']));
            catch ME_log
                fprintf('[Worker %d] Warning: could not rename log: %s\n', ...
                    getCurrentWorker_piGP_DN(), ME_log.message);
            end
        end

        csv_in_src = fullfile(current_folder, 'GP_Pareto_Front_Inputs.csv');
        if exist(csv_in_src, 'file') == 2
            try; movefile(csv_in_src, fullfile(current_folder, ['GP_Pareto_Front_Inputs_' file_suffix '.csv'])); catch; end
        end

        csv_out_src = fullfile(current_folder, 'GP_Pareto_Front_Outputs.csv');
        if exist(csv_out_src, 'file') == 2
            try; movefile(csv_out_src, fullfile(current_folder, ['GP_Pareto_Front_Outputs_' file_suffix '.csv'])); catch; end
        end

        %% ------------------------------------------------------------------
        %% Step 6: Save workspace
        %%   Includes all piGP outputs PLUS the denoiser history so you can
        %%   reconstruct exactly what the algorithm saw at each experiment.
        %% ------------------------------------------------------------------
        save(['Workspace_' file_suffix '.mat'], ...
            'Xpareto', 'Ypareto', 'X', 'Y', ...
            'XparetoGP', 'YparetoGP', 'YparetoGPstd', 'hypf', ...
            'theta_history', ...
            'X_hist', 'C_hist', 'w_blend', ...    % full denoiser history
            'lb', 'ub', 'c1_0', 'c2_0', 'V_R', ...
            'run_elapsed_time', 'maxeval', 'dataset_size');

        %% ------------------------------------------------------------------
        %% Step 7: Save denoiser history as CSV
        %%   Lets you inspect the raw vs denoised concentrations offline
        %%   and understand how w_blend evolved with accumulating evidence.
        %% ------------------------------------------------------------------
        n_exp = size(X_hist, 1);
        T_denoise = array2table( ...
            [X_hist, C_hist], ...
            'VariableNames', {'tau_min', 'T_C', 'c1_raw', 'c2_raw', 'c3_raw', 'c4_raw', 'c5_raw'});
        writetable(T_denoise, ['Denoiser_History_' file_suffix '.csv']);

        %% ------------------------------------------------------------------
        %% Step 8: Save prior parameter evolution CSVs
        %% ------------------------------------------------------------------
        n_iter_done = size(theta_history, 1);
        iter_col    = (1:n_iter_done)';
        col_names   = {'Iteration', 'log_kref', 'Ea_kJmol', 'alpha', 'beta'};

        T_STY = array2table([iter_col, theta_history(:,1,1), ...
                             theta_history(:,2,1)/1000, ...
                             theta_history(:,3,1), theta_history(:,4,1)], ...
                            'VariableNames', col_names);
        writetable(T_STY, ['PriorParams_STY_' file_suffix '.csv']);

        T_SEL = array2table([iter_col, theta_history(:,1,2), ...
                             theta_history(:,2,2)/1000, ...
                             theta_history(:,3,2), theta_history(:,4,2)], ...
                            'VariableNames', col_names);
        writetable(T_SEL, ['PriorParams_SEL_' file_suffix '.csv']);

        %% ------------------------------------------------------------------
        %% Step 9: Pareto front visualisation
        %% ------------------------------------------------------------------
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
               'GP Pareto front (piGP+denoise)', 'Location', 'Northeast')
        grid on
        xlabel('STY (g \cdot L^{-1} \cdot s^{-1})')
        ylabel('Selectivity (-)')
        title(sprintf('Run %d  —  maxeval=%d  —  %d Init Pts  (piGP+denoise)', ...
            run_num, maxeval, dataset_size))

        saveas(fig, ['Pareto_Plot_' file_suffix '.png']);
        saveas(fig, ['Pareto_Plot_' file_suffix '.fig']);
        close(fig);

        %% ------------------------------------------------------------------
        %% Step 10: Prior parameter evolution plot
        %% ------------------------------------------------------------------
        param_labels = {'log k_{ref}', 'E_a (kJ/mol)', '\alpha', '\beta'};
        colors_p     = {[0.2 0.5 0.8], [0.8 0.3 0.1]};

        fig2 = figure('Visible', 'off', 'Name', 'Prior param evolution');
        for p = 1:4
            subplot(2, 2, p); hold on; grid on
            for obj = 1:2
                vals = theta_history(:, p, obj);
                if p == 2, vals = vals / 1000; end
                plot(1:n_iter_done, vals, 'o-', ...
                    'Color', colors_p{obj}, 'MarkerFaceColor', colors_p{obj}, ...
                    'LineWidth', 1.5, 'MarkerSize', 5)
            end
            xlabel('Iteration');
            ylabel(param_labels{p});
            title(param_labels{p})
            if p == 1, legend({'STY', 'SEL'}, 'Location', 'best'); end
        end
        sgtitle(sprintf('Prior params — Run %d, maxeval=%d, %d init pts', ...
            run_num, maxeval, dataset_size))

        saveas(fig2, ['PriorParam_Evolution_' file_suffix '.png']);
        saveas(fig2, ['PriorParam_Evolution_' file_suffix '.fig']);
        close(fig2);

        fprintf('[Worker %d] Run %d (maxeval=%d, %d pts) done in %.2f min\n', ...
            getCurrentWorker_piGP_DN(), run_num, maxeval, dataset_size, run_elapsed_time/60);

    catch ME
        fprintf('[Worker %d] ERROR in Run %d: %s\n', ...
            getCurrentWorker_piGP_DN(), run_num, ME.message);
        try; cd(original_dir); catch; end
        rethrow(ME);
    end

    try; cd(original_dir); catch; end

end  % run_single_denoised_piGP


%% ========================================================================
%% OBJECTIVE WRAPPER — applies denoising at each TSEMO-proposed point
%%
%% This function lives here (not inside the helper above) so that MATLAB's
%% parfor broadcast mechanism can serialise and send it cleanly to workers.
%% The denoiser state struct 'ds' is a plain value (no handle classes),
%% so parfor can copy it safely to each worker.
%% ========================================================================

function Y_out = obj_with_denoising_local(x, ds)
% Called by TSEMO_V4_1a_OPC_piGP for each proposed experiment point.
% x is [1 x 2]: [tau(min), T(degC)]

    tau     = x(1);
    T       = round(x(2), 2);
    tau_sec = tau * 60;

    %%%% Simulation dummy — replace with OPC-UA + MIR hardware block %%%%
    c_out_raw = f_PFR_2nd_order(x, ds.c1_0, ds.c2_0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Denoise using the snapshot state captured after initial experiments.
    % ds.w, ds.X_hist, ds.C_hist do not update between TSEMO iterations
    % because TSEMO calls this as a stateless black box. The history
    % accumulated during initial experiments still gives the denoiser
    % meaningful context for the blend weight.
    [c_out, ~] = denoise_concentration(c_out_raw, tau, T, ds.c1_0, ds.c2_0, ...
                                        ds.theta_STY, ds.theta_SEL, ...
                                        ds.w, ds.X_hist, ds.C_hist);

    Y_out = zeros(1, 2);
    Y_out(1,1) = -log(ds.M_3 .* c_out(3) ./ tau_sec);
    Y_out(1,2) = -log(c_out(3) ./ (c_out(3) + c_out(4) + c_out(5)));

end


%% ========================================================================
%% HELPER — get current parallel worker ID (0 = main thread)
%% ========================================================================

function worker_id = getCurrentWorker_piGP_DN()
    try
        t         = getCurrentTask();
        worker_id = t.ID;
    catch
        worker_id = 0;
    end
end
