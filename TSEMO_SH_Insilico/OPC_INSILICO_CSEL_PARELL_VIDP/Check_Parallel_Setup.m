% Check_Parallel_Setup.m
% Diagnostic script to check if Parallel Computing Toolbox is properly set up

fprintf('\n========================================================================\n');
fprintf('           PARALLEL COMPUTING TOOLBOX DIAGNOSTIC\n');
fprintf('========================================================================\n\n');

%% Check 1: License
fprintf('1. Checking license...\n');
if license('test', 'Distrib_Computing_Toolbox')
    fprintf('   ✓ Parallel Computing Toolbox license found\n\n');
else
    fprintf('   ✗ Parallel Computing Toolbox license NOT found\n');
    fprintf('   → You need to install/activate the toolbox\n\n');
    return;
end

%% Check 2: Installed toolboxes
fprintf('2. Checking installed toolboxes...\n');
v = ver;
pct_installed = false;
for i = 1:length(v)
    if contains(v(i).Name, 'Parallel', 'IgnoreCase', true)
        fprintf('   ✓ Found: %s (Version %s)\n', v(i).Name, v(i).Version);
        pct_installed = true;
    end
end
if ~pct_installed
    fprintf('   ✗ Parallel Computing Toolbox not in installed toolboxes\n');
    fprintf('   → Install it from MATLAB Add-Ons\n\n');
    return;
end
fprintf('\n');

%% Check 3: Number of cores
fprintf('3. Checking CPU cores...\n');
try
    num_cores = feature('numcores');
    fprintf('   ✓ Number of physical cores: %d\n\n', num_cores);
catch
    fprintf('   ? Could not determine number of cores\n\n');
end

%% Check 4: Try to get current pool
fprintf('4. Checking existing parallel pool...\n');
try
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        fprintf('   ○ No parallel pool currently running\n\n');
    else
        fprintf('   ✓ Parallel pool exists with %d workers\n\n', poolobj.NumWorkers);
    end
catch ME
    fprintf('   ✗ ERROR calling gcp: %s\n', ME.message);
    fprintf('   → This is the problem! See solutions below.\n\n');
    gcp_error = true;
end

%% Check 5: Try to create a pool
if ~exist('gcp_error', 'var') || ~gcp_error
    fprintf('5. Attempting to create parallel pool...\n');
    try
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            fprintf('   Creating new pool with 2 workers (test)...\n');
            poolobj = parpool(2);
            fprintf('   ✓ SUCCESS! Parallel pool created with %d workers\n', poolobj.NumWorkers);
            fprintf('   Closing test pool...\n');
            delete(poolobj);
            fprintf('   ✓ Test complete\n\n');
        else
            fprintf('   Pool already exists, skipping creation test\n\n');
        end
    catch ME
        fprintf('   ✗ ERROR creating pool: %s\n', ME.message);
        fprintf('   → See solutions below\n\n');
        pool_error = true;
    end
end

%% Summary and Solutions
fprintf('========================================================================\n');
fprintf('                              SUMMARY\n');
fprintf('========================================================================\n\n');

if ~exist('gcp_error', 'var') && ~exist('pool_error', 'var')
    fprintf('✓ ALL CHECKS PASSED!\n');
    fprintf('  Your parallel computing setup is working correctly.\n');
    fprintf('  You can use TSEMO_Multiple_Runs_Parallel.m\n\n');
else
    fprintf('✗ ISSUES DETECTED\n\n');
    
    fprintf('POSSIBLE SOLUTIONS:\n\n');
    
    fprintf('Solution 1: Restart MATLAB\n');
    fprintf('  Sometimes MATLAB needs a restart for toolbox changes to take effect.\n\n');
    
    fprintf('Solution 2: Check MATLAB version\n');
    fprintf('  Run: version\n');
    fprintf('  Parallel pool (parpool) requires MATLAB R2013b or later.\n');
    fprintf('  If you have an older version, you need to use ''matlabpool'' instead.\n\n');
    
    fprintf('Solution 3: Reinstall/Update Parallel Computing Toolbox\n');
    fprintf('  1. Go to Home → Add-Ons → Get Add-Ons\n');
    fprintf('  2. Search for "Parallel Computing Toolbox"\n');
    fprintf('  3. Install or update it\n\n');
    
    fprintf('Solution 4: Check cluster/profile settings\n');
    fprintf('  Run: parallel.defaultClusterProfile\n');
    fprintf('  Should return ''local'' or ''Processes''\n\n');
    
    fprintf('Solution 5: Delete existing preferences\n');
    fprintf('  Run: parallel.internal.settings.clear\n');
    fprintf('  Then restart MATLAB\n\n');
    
    fprintf('Solution 6: Use Sequential Version Instead\n');
    fprintf('  If parallel doesn''t work, use TSEMO_Multiple_Runs_Simple.m\n');
    fprintf('  It will run sequentially but still organize your results.\n\n');
end

fprintf('========================================================================\n');
fprintf('For more help, run: doc parpool\n');
fprintf('========================================================================\n');

%% Additional diagnostic info
fprintf('\nADDITIONAL INFORMATION:\n\n');
fprintf('MATLAB Version: %s\n', version);
fprintf('MATLAB Root: %s\n', matlabroot);

try
    profile = parallel.defaultClusterProfile;
    fprintf('Default Cluster Profile: %s\n', profile);
catch
    fprintf('Default Cluster Profile: Could not determine\n');
end

fprintf('\n');
