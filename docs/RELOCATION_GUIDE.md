# Relocation Guide

## Purpose
Use this guide when moving the MSM analysis code to a new machine or path while keeping behavior reproducible for both humans and AI agents.

## 1. Repository Contract
This repository is code-and-doc only.

It should not store:
- raw source data
- generated condition packs
- generated figures/tables
- binary model/result artifacts

## 2. Required External Inputs
See `docs/DATA_DEPENDENCIES.md` for details.
At minimum, most workflows depend on:
- `dataTabFULL.xlsx` (or equivalent canonical trial table)
- `MSM_neuronSortMASTER_recoveryFile*.xlsx` for some extraction/reconciliation steps
- `.nex` spike files for neural workflows

## 3. Path Strategy (Recommended)
Adopt two local roots on each machine:
- `MSM_CODE_ROOT`: path to this repo
- `MSM_DATA_ROOT`: path to local raw/input datasets (not versioned)
- `MSM_OUTPUT_ROOT`: path for generated outputs (not versioned)

Example (Windows):
- `MSM_CODE_ROOT = C:\Users\<user>\MSM`
- `MSM_DATA_ROOT = D:\MSM_data`
- `MSM_OUTPUT_ROOT = D:\MSM_outputs`

## 4. Before Running Scripts
1. Open target script.
2. Locate hardcoded absolute paths.
3. Replace with current local paths.
4. Confirm input file existence.
5. Confirm output folder exists and is non-versioned.

## 5. Expected Relocation Failure Modes
- stale hardcoded path from previous machine
- script expecting condition packs that are not yet generated
- mismatch between expected and actual sheet variable names
- role/condition labels drifting across source sheets

## 6. Recommended First Validation Sequence
1. `plot_condition_by_date_251204.m`
2. `plot_choice_distributions_251205.m`
3. `GLM_behav_only_251210.m`
4. `RL_basline_and_saline_separated_251217.m`
5. active strategy suite (listed in `README.md`)

## 7. AI/Handoff Notes
When handing off to a new collaborator or AI:
- provide current `MSM_DATA_ROOT` and `MSM_OUTPUT_ROOT`
- state which condition-pack generation script was last run
- state which cohort/session exclusions were applied
- state which script/date suffix is considered canonical for each analysis family
