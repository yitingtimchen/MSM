# MSM Workflow Map

Project root: `C:\Users\plattlab\MSM`

## Active Script Locations

- Current scripts: `scripts/active/`
- Legacy provenance scripts: `scripts/legacy/am_orig_codes_250602/`
- Non-canonical experimental scripts: `scripts/legacy/bad_codes/`

## Active Branches

1. Behavior-only branch (canonical)
- Extractor: `scripts/active/extract_and_save_behav_only_by_cond_and_split_251205.m`
- Input table: `C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx`
- Output folder: `C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL`

2. Neural branch (kept)
- Extractor: `scripts/active/extract_and_save_spike_and_event_by_cond_251118.m`
- Output folder: `C:\Users\plattlab\MSM\outputs_local\condition_packs_251118`

## Recommended Run Order

1. Run `scripts/active/extract_and_save_behav_only_by_cond_and_split_251205.m`
2. Run behavior analyses in `scripts/active/`
3. If needed, run `scripts/active/extract_and_save_spike_and_event_by_cond_251118.m`
4. Run neural analyses in `scripts/active/`
