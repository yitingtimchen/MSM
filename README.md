# Monkey Stock Market (MSM)

Last updated: 2026-02-16

This project is now organized in a repository-style layout with clear separation of code, data, docs, and outputs.

## Repository Layout

- `scripts/active/`: canonical and current MATLAB scripts
- `scripts/legacy/am_orig_codes_250602/`: preserved legacy scripts from original author
- `scripts/legacy/bad_codes/`: non-canonical experimental scripts
- `data/`: copied local data organized by modality
- `outputs_local/`: generated local outputs and condition packs
- `docs/`: runbooks, catalogs, and preflight reports

## Canonical Branches

1. Behavior branch
- Extractor: `scripts/active/extract_and_save_behav_only_by_cond_and_split_251205.m`
- Input: `C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx`
- Output: `C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL`

2. Neural branch
- Extractor: `scripts/active/extract_and_save_spike_and_event_by_cond_251118.m`
- Output: `C:\Users\plattlab\MSM\outputs_local\condition_packs_251118`

## Quick Start

1. Open MATLAB in `C:\Users\plattlab\MSM`.
2. Run `scripts/SETUP_PATHS.m` once per session.
3. Follow `docs/PIPELINE_RUNBOOK.md`.

## Core Docs

- `docs/PIPELINE_RUNBOOK.md`
- `docs/SCRIPT_CATALOG.md`
- `docs/STRICT_PREFLIGHT_2026-02-16.md`
- `docs/BRANCH_CLASSIFICATION_2026-02-16.md`
- `WORKFLOW_MAP.md`
