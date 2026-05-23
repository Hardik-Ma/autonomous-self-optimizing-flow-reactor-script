# README: TSEMO_Multiple_Runs_Different_Meval_and_DatasetSize.m

## Overview

This MATLAB script automates multi-batch Bayesian optimization experiments using the **TSEMO (Thompson Sampling Efficient Multi-objective Optimization)** algorithm applied to a **2nd-order Plug Flow Reactor (PFR)** problem.

The script systematically varies two key experimental parameters across batches:
- **`maxeval`** — the maximum number of objective function evaluations TSEMO is allowed to perform
- **Initial dataset size** — the number of Latin Hypercube Sampling (LHS) points used to initialize the surrogate model

Each batch can run multiple independent optimization trials (runs), optionally in parallel. Results for every run are saved with structured file naming and organized into a clearly defined folder hierarchy.

The two objectives being optimized are:
- **STY** — Space-Time Yield of product 3
- **Selectivity** — Selectivity towards product 3

---

## Requirements

### MATLAB Version
- MATLAB R2019b or later (recommended for full `table`, `parfor`, and `saveas` compatibility)

### Required MATLAB Toolboxes
| Toolbox | Required? | Purpose |
|---|---|---|
| **Parallel Computing Toolbox** | Optional (but recommended) | Enables `parfor` parallel execution of runs within each batch |
| **Statistics and Machine Learning Toolbox** | Required | Provides `lhsdesign` for Latin Hypercube Sampling |

> **Note:** If the Parallel Computing Toolbox is not available or licensed, the script automatically falls back to sequential execution — no manual changes are needed.

### Required External Dependencies (Folders to Add to MATLAB Path)

Before running this script, the following folders **must be added to the MATLAB path**. These contain the TSEMO algorithm, the reactor model, compiled MEX files, and supporting utilities.

| Folder | Contents |
|---|---|
| `Direct` | The DIRECT global optimization algorithm used internally by TSEMO |
| **`Hone_Taylor_Reaction`** | The reactor model and objective function `f_PFR_2nd_order` for the Hone–Taylor reaction system |
| `Mex_files` | Compiled MEX binaries that accelerate Gaussian Process operations |
| `NGPM_v1.4` | Non-dominated Genetic Programming Module, used by the Pareto front computation |

#### How to Add Folders to Path

**Option A — Temporarily (in the MATLAB Command Window before running):**
```matlab
addpath(genpath('C:\path\to\Direct'));
addpath(genpath('C:\path\to\Hone_Taylor_Reaction'));
addpath(genpath('C:\path\to\Mex_files'));
addpath(genpath('C:\path\to\NGPM_v1.4'));
```

**Option B — Permanently (via MATLAB GUI):**  
Go to **Home → Set Path → Add with Subfolders**, add each folder, then click **Save**.

**Option C — Via a startup script (`startup.m`):**  
Place a `startup.m` file in your MATLAB working directory or userpath containing the `addpath` calls above. MATLAB will execute it automatically on startup.

> **Important for Parallel Workers:** The script automatically attempts to attach all `.m` files in the current directory to the parallel pool workers. However, the four dependency folders above must also be on the path of each worker. If you encounter errors in parallel mode related to missing functions, add the paths inside the `run_single_optimization_combined` function directly, or configure your parallel pool's `PathDependencies`.

---

## Configuration

Open the script and edit the **CONFIGURATION** section at the top:

```matlab
batches = [
    % maxeval, runs, initial_points
    10,  2,  4;   % Batch 1: 2 runs with maxeval=10, 4 initial points
    50,  3,  8;   % Batch 2: 3 runs with maxeval=50, 8 initial points
    100, 5, 16;   % Batch 3: 5 runs with maxeval=100, 16 initial points
];
```

Each row defines one batch with three parameters:

| Column | Parameter | Description |
|---|---|---|
| 1 | `maxeval` | Maximum number of TSEMO evaluations (algorithm iterations) |
| 2 | `runs` | Number of independent optimization trials in this batch |
| 3 | `initial_points` | Number of LHS points for initializing the surrogate model |

You may also configure:
```matlab
num_workers = [];   % [] = use all available CPU cores
                    % Set to a number, e.g. 4, to limit workers
```

---

## How to Run

1. Add the required dependency folders to your MATLAB path (see above).
2. Set your working directory to the folder containing this script.
3. Edit the `batches` matrix to define your experimental plan.
4. Run the script:
   ```matlab
   run('TSEMO_Multiple_Runs_Different_Meval_and_DatasetSize.m')
   ```
   or open the file in the MATLAB Editor and press **Run (F5)**.

---

## What the Script Does (Step by Step)

For each batch and each run within it, the script:

1. **Defines the reactor problem** — sets reactor volume (`V_R = 2 mL`), feed concentrations (`c1_0 = c2_0 = 0.2 mol/L`), and operating bounds (residence time τ: 0.25–2 min; temperature T: 30–50 °C).
2. **Generates an LHS initial dataset** — creates `initial_points` space-filling samples using Latin Hypercube Sampling over the input domain.
3. **Evaluates the initial experiments** — calls `f_PFR_2nd_order` for each LHS point to compute STY and Selectivity, stored as log-transformed negative values (minimization convention).
4. **Runs TSEMO** — calls `TSEMO_V4_1a_OPC` with the configured `maxeval` budget to iteratively find the Pareto-optimal front.
5. **Renames output files** — appends a descriptive suffix (`Run_N_Meval_M_InitPts_P`) to all output files generated by TSEMO for traceability.
6. **Saves the workspace** — saves all key variables to a `.mat` file.
7. **Generates and saves plots** — produces a Pareto front visualization saved as both `.png` and `.fig`.

---

## Output Folder Structure

When the script runs, it creates a **master results folder** named with a timestamp in the current working directory:

```
TSEMO_MultiMeval_MultiDataset_YYYY-MM-DD_HH-MM-SS/
│
├── Batch_1_Meval_10_InitPts_4/
│   ├── Run_1/
│   │   ├── TSEMO_log_Run_1_Meval_10_InitPts_4.txt
│   │   ├── GP_Pareto_Front_Inputs_Run_1_Meval_10_InitPts_4.csv
│   │   ├── GP_Pareto_Front_Outputs_Run_1_Meval_10_InitPts_4.csv
│   │   ├── Workspace_Run_1_Meval_10_InitPts_4.mat
│   │   ├── Pareto_Plot_Run_1_Meval_10_InitPts_4.png
│   │   └── Pareto_Plot_Run_1_Meval_10_InitPts_4.fig
│   │
│   └── Run_2/
│       ├── TSEMO_log_Run_2_Meval_10_InitPts_4.txt
│       ├── GP_Pareto_Front_Inputs_Run_2_Meval_10_InitPts_4.csv
│       ├── GP_Pareto_Front_Outputs_Run_2_Meval_10_InitPts_4.csv
│       ├── Workspace_Run_2_Meval_10_InitPts_4.mat
│       ├── Pareto_Plot_Run_2_Meval_10_InitPts_4.png
│       └── Pareto_Plot_Run_2_Meval_10_InitPts_4.fig
│
├── Batch_2_Meval_50_InitPts_8/
│   ├── Run_1/
│   │   └── ... (same file structure as above)
│   ├── Run_2/
│   └── Run_3/
│
├── Batch_Summary.csv
└── All_Batch_Results.mat
```

### Folder and File Naming Conventions

#### Master Folder
```
TSEMO_MultiMeval_MultiDataset_YYYY-MM-DD_HH-MM-SS
```
- Timestamp format ensures each experiment run creates a unique, non-overwriting folder.
- Example: `TSEMO_MultiMeval_MultiDataset_2025-06-15_09-30-00`

#### Batch Folders
```
Batch_{batch_index}_Meval_{maxeval}_InitPts_{initial_points}
```
- `batch_index` — sequential batch number (1, 2, 3, …)
- `maxeval` — the maximum number of TSEMO evaluations for this batch
- `initial_points` — number of LHS initialization points for this batch
- Example: `Batch_1_Meval_50_InitPts_8`

#### Run Folders
```
Run_{run_number}
```
- Each run gets its own subfolder inside its batch folder.
- Example: `Run_3`

#### Files Inside Each Run Folder

All files share a common suffix: `Run_{N}_Meval_{M}_InitPts_{P}`, where N = run number, M = maxeval, P = initial points.

| File | Description |
|---|---|
| `TSEMO_log_Run_N_Meval_M_InitPts_P.txt` | Text log generated by TSEMO during the optimization run. Contains iteration-by-iteration progress, hyperparameter updates, and convergence information. |
| `GP_Pareto_Front_Inputs_Run_N_Meval_M_InitPts_P.csv` | Input variable values (τ, T) corresponding to the GP-predicted Pareto-optimal points. Each row is one Pareto point in the input space. |
| `GP_Pareto_Front_Outputs_Run_N_Meval_M_InitPts_P.csv` | Objective values (STY, Selectivity) for the GP-predicted Pareto-optimal points. Each row corresponds to the same Pareto point as in the Inputs CSV. |
| `Workspace_Run_N_Meval_M_InitPts_P.mat` | Full MATLAB workspace saved at the end of the run. Contains: `Xpareto`, `Ypareto` (true evaluated Pareto front), `X`, `Y` (all evaluated points), `XparetoGP`, `YparetoGP`, `YparetoGPstd` (GP-predicted Pareto front with uncertainty), `hypf` (GP hyperparameters), `lb`, `ub`, `c1_0`, `c2_0`, `V_R`, `run_elapsed_time`, `maxeval`, `dataset_size`. |
| `Pareto_Plot_Run_N_Meval_M_InitPts_P.png` | Static image of the Pareto front plot (suitable for reports and quick inspection). |
| `Pareto_Plot_Run_N_Meval_M_InitPts_P.fig` | MATLAB figure file of the Pareto front plot (can be reopened and edited in MATLAB). |

#### Master-Level Summary Files

| File | Description |
|---|---|
| `Batch_Summary.csv` | A summary table with one row per batch. Columns: `Batch`, `MaxEval`, `InitialPoints`, `NumRuns`, `Time_min` (total batch time in minutes), `AvgTime_per_run_sec` (average time per individual run in seconds). Useful for comparing computational cost across configurations. |
| `All_Batch_Results.mat` | MATLAB `.mat` file containing `all_batch_results` (cell array with metadata for each batch), `batches` (the original configuration matrix), and `summary_table`. Useful for programmatic post-processing across all batches. |

---

## Pareto Plot Legend

Each `Pareto_Plot_*.png/.fig` uses the following colour coding:

| Symbol | Colour | Meaning |
|---|---|---|
| `.` (dots) | Orange | Initial LHS points (evaluated before TSEMO started) |
| `x` (crosses) | Orange | Points evaluated by the TSEMO algorithm |
| `O` (circles) | Yellow | True Pareto front (from evaluated data) |
| `. ` with error bars | Purple | GP-predicted Pareto front (with uncertainty bounds) |

---

## Troubleshooting

| Issue | Likely Cause | Solution |
|---|---|---|
| `Undefined function 'f_PFR_2nd_order'` | `Hone_Taylor_Reaction` folder not on path | Add it via `addpath` |
| `Undefined function 'TSEMO_V4_1a_OPC'` | TSEMO source folder not on path | Add the TSEMO folder via `addpath` |
| `Undefined function 'lhsdesign'` | Statistics and Machine Learning Toolbox missing | Install or license the toolbox |
| `parfor` running sequentially | Parallel Computing Toolbox not licensed | Expected — script falls back automatically |
| Workers fail with missing function errors | Dependencies not accessible to workers | Pass paths explicitly inside `run_single_optimization_combined`, or use `addAttachedFiles` with full paths |
| MEX file errors on a new machine | MEX files compiled for different OS/architecture | Recompile MEX sources on the target machine |

---

## Key Parameters (Reactor Model)

These are hardcoded in `run_single_optimization_combined` and reflect the physical reactor setup:

| Parameter | Value | Description |
|---|---|---|
| `V_R` | 2 mL | Reactor volume |
| `c1_0` | 0.2 mol/L | Initial concentration of reactant 1 |
| `c2_0` | 0.2 mol/L | Initial concentration of reactant 2 |
| `tau_min` / `tau_max` | 0.25 / 2 min | Residence time bounds |
| `T_min` / `T_max` | 30 / 50 °C | Temperature bounds |
| `t_meas` | 3 | Measurement time parameter passed to TSEMO |
| `NoOfBachSequential` | 1 | Number of sequential batch candidates per TSEMO iteration |
