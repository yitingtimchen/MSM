# MSM Preflight Report

Generated: 2026-02-16
Project root: `C:\Users\plattlab\MSM`

## Summary
- Path refactor status: complete for `.m` scripts (no remaining `Stock_market_experiment` path references).
- Canonical behavior extractor: `d20251205_extract_behavior_condition_packs_with_split.m`.
- Neural extractor retained: `d20251118_extract_spike_event_condition_packs.m`.

## Core Input Checks
- OK `data/tooling/MATLAB_Copy/dataTabFULL.xlsx`
- OK `data/metadata/MSM_neuronSortMASTER_recoveryFile.xlsx`
- OK `data/metadata/MSM_neuronSortMASTER_recoveryFile_modified.xlsx`
- OK `data/metadata/MSM_neuronSortMASTER_recoveryFile_modified_filtered.xlsx`
- OK `data/neural/raw_external/MSM_Sorted`

## Output Targets
- Behavior packs: `outputs_local/condition_packs_behavior_only_from_dataTabFULL`
- Neural packs: `outputs_local/condition_packs_251118`

## Runnable Now
- `d20251205_extract_behavior_condition_packs_with_split.m`
- `d20251118_extract_spike_event_condition_packs.m`

## Blocked Until Packs Exist
Behavior and RL/strategy scripts that load `*_condition_pack.mat` are blocked until the extractor(s) are run.

## Recommended Unblock Sequence
1. Run `d20251205_extract_behavior_condition_packs_with_split.m`
2. Confirm `outputs_local/condition_packs_behavior_only_from_dataTabFULL/*.mat` exists
3. Run behavior analyses
4. If needed, run `d20251118_extract_spike_event_condition_packs.m` and then neural analyses
