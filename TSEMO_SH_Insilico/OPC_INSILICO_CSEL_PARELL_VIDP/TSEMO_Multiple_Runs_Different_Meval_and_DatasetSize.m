
clear all
close all
clc


%% CONFIGURATION - Define your batches here

% Define batches: [maxeval, number_of_runs, initial_dataset_size]
% Each row is one batch with THREE parameters
batches = [
    % maxeval, runs, initial_points
    10,  2,  4;   % Batch 1: 2 runs with maxeval=10, 4 initial points
];

% Number of workers (CPU cores) to use
num_workers = [];  % Use all available cores

% Create master results folder with timestamp
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
master_folder_name = ['TSEMO_MultiMeval_MultiDataset_', timestamp];

% Get absolute path for results folder
current_dir = pwd;
master_folder = fullfile(current_dir, master_folder_name);
mkdir(master_folder);

fprintf('========================================================================\n');
fprintf('TSEMO OPTIMIZATION - VARYING MAXEVAL AND INITIAL DATASET SIZE\n');
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


%% CHECK AND START PARALLEL POOL

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
        
        % Add all required files to parallel workers
        fprintf('Adding required files to parallel workers...\n');
        try
            current_files = dir('*.m');
            files_to_add = {};
            for i = 1:length(current_files)
                files_to_add{end+1} = fullfile(pwd, current_files(i).name);
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


%% RUN ALL BATCHES

% Storage for all results
all_batch_results = cell(size(batches, 1), 1);
total_start_time = tic;

for batch_idx = 1:size(batches, 1)
    
    current_maxeval = batches(batch_idx, 1);
    current_num_runs = batches(batch_idx, 2);
    current_dataset_size = batches(batch_idx, 3);
    
    fprintf('\n\n');
    fprintf('************************************************************************\n');
    fprintf('                    STARTING BATCH %d of %d\n', batch_idx, size(batches, 1));
    fprintf('          maxeval=%d, runs=%d, initial points=%d\n', ...
        current_maxeval, current_num_runs, current_dataset_size);
    fprintf('************************************************************************\n');
    fprintf('\n');
    
    % Create folder for this batch (with all three parameters in name)
    batch_folder = fullfile(master_folder, ...
        sprintf('Batch_%d_Meval_%d_InitPts_%d', batch_idx, current_maxeval, current_dataset_size));
    mkdir(batch_folder);
    
    % Start timer for this batch
    batch_start_time = tic;
    
    if use_parallel
        fprintf('Using PARALLEL execution with %d workers\n\n', poolobj.NumWorkers);
        parfor run_num = 1:current_num_runs
            run_single_optimization_combined(run_num, batch_folder, current_maxeval, current_dataset_size);
        end
    else
        fprintf('Using SEQUENTIAL execution (one run at a time)\n\n');
        for run_num = 1:current_num_runs
            run_single_optimization_combined(run_num, batch_folder, current_maxeval, current_dataset_size);
        end
    end
    
    batch_elapsed_time = toc(batch_start_time);
    
    % Store batch results
    all_batch_results{batch_idx}.maxeval = current_maxeval;
    all_batch_results{batch_idx}.num_runs = current_num_runs;
    all_batch_results{batch_idx}.dataset_size = current_dataset_size;
    all_batch_results{batch_idx}.elapsed_time = batch_elapsed_time;
    all_batch_results{batch_idx}.folder = batch_folder;
    
    fprintf('\n');
    fprintf('************************************************************************\n');
    fprintf('                    BATCH %d COMPLETED\n', batch_idx);
    fprintf('          Time: %.2f seconds (%.2f minutes)\n', batch_elapsed_time, batch_elapsed_time/60);
    fprintf('************************************************************************\n');
    
end

total_elapsed_time = toc(total_start_time);

%% FINAL SUMMARY

fprintf('\n\n');
fprintf('========================================================================\n');
fprintf('              ALL BATCHES COMPLETED SUCCESSFULLY\n');
fprintf('========================================================================\n\n');

% Create summary table
summary_table = table();
for i = 1:length(all_batch_results)
    summary_table.Batch(i) = i;
    summary_table.MaxEval(i) = all_batch_results{i}.maxeval;
    summary_table.InitialPoints(i) = all_batch_results{i}.dataset_size;
    summary_table.NumRuns(i) = all_batch_results{i}.num_runs;
    summary_table.Time_min(i) = all_batch_results{i}.elapsed_time / 60;
    summary_table.AvgTime_per_run_sec(i) = all_batch_results{i}.elapsed_time / all_batch_results{i}.num_runs;
end

disp(summary_table);

fprintf('\nTotal time for all batches: %.2f seconds (%.2f minutes, %.2f hours)\n', ...
    total_elapsed_time, total_elapsed_time/60, total_elapsed_time/3600);
fprintf('Results location: %s\n', master_folder);
fprintf('========================================================================\n');

% Save summary
writetable(summary_table, fullfile(master_folder, 'Batch_Summary.csv'));
save(fullfile(master_folder, 'All_Batch_Results.mat'), 'all_batch_results', 'batches', 'summary_table');

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


%% HELPER FUNCTION - Runs a single optimization

function run_single_optimization_combined(run_num, batch_folder, maxeval, dataset_size)
    % This function runs on each parallel worker
    
    fprintf('[Worker %d] Starting Run %d (maxeval=%d, %d initial pts)\n', ...
        getCurrentWorker(), run_num, maxeval, dataset_size);
    
    run_start_time = tic;
    
    % Create folder for this run with retry logic
    run_folder = fullfile(batch_folder, sprintf('Run_%d', run_num));
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
                getCurrentWorker(), retry_count, ME_mkdir.message);
            pause(1);
        end
    end
    
    if ~folder_created
        error('Failed to create folder %s after %d attempts', run_folder, max_retries);
    end
    
    % Change to run folder with error handling
    original_dir = pwd;
    try
        cd(run_folder);
    catch ME_cd
        error('Cannot change to folder %s: %s', run_folder, ME_cd.message);
    end
    
    try
       
        %% Reactor + design specification
        
        V_R = 2; %2ml
        c1_0 = 0.2; %mol/L
        c2_0 = 0.2; %mol/L

        % operating limits
        tau_min = 0.25; %min
        tau_max = 2; %min
        T_min = 30; %°C
        T_max = 50; %°C
        
       
        %% Step 1: Specify problem
       
        no_outputs = 2;
        no_inputs  = 2;
        lb = [tau_min, T_min];
        ub = [tau_max, T_max];

        f = @(x)f_PFR_2nd_order(x,c1_0,c2_0);
        
        
        %% Step 2: Generate LHS with specified dataset size
        X = lhsdesign(dataset_size, no_inputs);
        X = sortrows(X, 2);  % Sort by temperature
        Y = zeros(dataset_size, no_outputs);
        
        %% Step 3: Run initial Experiments
        for k = 1:size(X,1)
            X(k,:) = X(k,:).*(ub-lb)+lb;
            tau = X(k,1);
            tau_sec = X(k,1) * 60;
            T = X(k,2);
            T = round(T, 2);

            Q_tot = V_R / tau;
            Q_i = Q_tot / 2;

            c_out = f_PFR_2nd_order(X(k,:), c1_0, c2_0);

            % Molecular Mass [g/mol]
            M_1 = 159.09;
            M_2 = 87.12;
            M_F = 18.99;
            M_H = 1.00;
            M_3 = M_1 + M_2 - M_F - M_H;
            M_4 = M_3;
            M_5 = M_3 + M_2 - M_F - M_H;

            % OBJECTIVES: STY and Selectivity
            Y(k,1) = -log(M_3.*c_out(3)./(tau_sec));  % STY
            Y(k,2) = -log(c_out(3)./(c_out(3)+c_out(4)+c_out(5)));  % Selectivity
        end
        

        %% Step 4: Start TSEMO algorithm

        opt = TSEMO_options;
        opt.maxeval = maxeval;
        opt.NoOfBachSequential = 1;
        
        t_meas = 3;
        
        [Xpareto,Ypareto,X,Y,XparetoGP,YparetoGP,YparetoGPstd,hypf] = ...
            TSEMO_V4_1a_OPC(f,X,Y,lb,ub,opt,c1_0,c2_0,V_R,t_meas);
        

        %% Step 5: Rename files with run number

        
        run_elapsed_time = toc(run_start_time);
        
        % Get current directory
        current_folder = pwd;
        
        % Create filename suffix
        file_suffix = sprintf('Run_%d_Meval_%d_InitPts_%d', run_num, maxeval, dataset_size);
        
        % Rename output files with full error handling
        log_file = fullfile(current_folder, 'TSEMO_log.txt');
        if exist(log_file, 'file') == 2
            try
                movefile(log_file, fullfile(current_folder, ['TSEMO_log_' file_suffix '.txt']));
            catch ME_log
                fprintf('[Worker %d] Warning: Could not rename TSEMO_log.txt\n', getCurrentWorker());
            end
        end
        
        csv_input_file = fullfile(current_folder, 'GP_Pareto_Front_Inputs.csv');
        if exist(csv_input_file, 'file') == 2
            try
                movefile(csv_input_file, fullfile(current_folder, ['GP_Pareto_Front_Inputs_' file_suffix '.csv']));
            catch
            end
        end
        
        csv_output_file = fullfile(current_folder, 'GP_Pareto_Front_Outputs.csv');
        if exist(csv_output_file, 'file') == 2
            try
                movefile(csv_output_file, fullfile(current_folder, ['GP_Pareto_Front_Outputs_' file_suffix '.csv']));
            catch
            end
        end
        
        % Save workspace
        save(['Workspace_' file_suffix '.mat'], ...
            'Xpareto', 'Ypareto', 'X', 'Y', 'XparetoGP', 'YparetoGP', 'YparetoGPstd', 'hypf', ...
            'lb', 'ub', 'c1_0', 'c2_0', 'V_R', 'run_elapsed_time', 'maxeval', 'dataset_size');
        

        %% Step 6: Create visualization

        
        % Transform from logarithmic form to normal
        Y1 = exp(-Y(:,1));
        Y2 = exp(-Y(:,2));
        Ypareto1 = exp(-Ypareto(:,1));
        Ypareto2 = exp(-Ypareto(:,2));
        YparetoGP1 = exp(-YparetoGP(:,1));
        YparetoGP2 = exp(-YparetoGP(:,2));
        YparetoGPstd1 = exp(-YparetoGP(:,1))-exp(-(YparetoGP(:,1)+YparetoGPstd(:,1))); 
        YparetoGPstd2 = exp(-YparetoGP(:,2))-exp(-(YparetoGP(:,2)+YparetoGPstd(:,2))); 
        
        fig = figure('Visible', 'off');
        hold on
        plot(Y1(1:dataset_size),Y2(1:dataset_size),'.','MarkerSize',14,'color',[0.8500 0.3250 0.0980])
        plot(Y1(dataset_size+1:end),Y2(dataset_size+1:end),'x','MarkerSize',8,'LineWidth',2,'color',[0.8500 0.3250 0.0980])
        plot(Ypareto1,Ypareto2,'O','MarkerSize',8,'LineWidth',2,'color',"#EDB120")
        errorbar(YparetoGP1,YparetoGP2,YparetoGPstd2,YparetoGPstd2,YparetoGPstd1,YparetoGPstd1,'.','MarkerSize',8,'color',"#7E2F8E",'LineWidth',0.5)
        legend('Initial LHS','Algorithm','Pareto front','GP Pareto front','Location','Northeast')
        grid on
        xlabel('STY (3) [(g.L^{-3}.sec^{-1})]')
        ylabel('SEL (Selectivity)')
        title(sprintf('Run %d - maxeval=%d - %d Init Pts', run_num, maxeval, dataset_size))
        
        saveas(fig, ['Pareto_Plot_' file_suffix '.png']);
        saveas(fig, ['Pareto_Plot_' file_suffix '.fig']);
        close(fig);
        
        fprintf('[Worker %d] Run %d (maxeval=%d, %d pts) done: %.2f min\n', ...
            getCurrentWorker(), run_num, maxeval, dataset_size, run_elapsed_time/60);
        
    catch ME
        fprintf('[Worker %d] ERROR in Run %d: %s\n', getCurrentWorker(), run_num, ME.message);
        try
            cd(original_dir);
        catch
        end
        rethrow(ME);
    end
    
    % Return to original directory
    try
        cd(original_dir);
    catch
    end
    
end

%% HELPER FUNCTION - Get current worker ID


function worker_id = getCurrentWorker()
    try
        t = getCurrentTask();
        worker_id = t.ID;
    catch
        worker_id = 0;
    end
end
