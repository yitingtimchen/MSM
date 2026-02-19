# MSM Workflow Map

Project root: `C:\Users\plattlab\MSM`

## Active Script Locations

- Current scripts: `scripts/active/`
- Legacy provenance scripts: `scripts/legacy/am_orig_codes_250602/`
- Non-canonical experimental scripts: `scripts/legacy/bad_codes/`

## Active Branches

1. Behavior-only branch (canonical)
- Extractor: `scripts/active/d20251205_extract_behavior_condition_packs_with_split.m`
- Input table: `C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx`
- Output folder: `C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL`

2. Neural branch (kept)
- Extractor: `scripts/active/d20251118_extract_spike_event_condition_packs.m`
- Output folder: `C:\Users\plattlab\MSM\outputs_local\condition_packs_251118`

## Recommended Run Order

1. Run `scripts/active/d20251205_extract_behavior_condition_packs_with_split.m`
2. Run behavior analyses in `scripts/active/`
3. If needed, run `scripts/active/d20251118_extract_spike_event_condition_packs.m`
4. Run neural analyses in `scripts/active/`
