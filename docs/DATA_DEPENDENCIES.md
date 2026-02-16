# Data Dependencies

## Scope
This file documents external data expected by scripts in this repository.
Data is intentionally excluded from git and must be provisioned locally.

## Core Inputs

### 1) Behavioral Trial Table
Typical filename:
- `dataTabFULL.xlsx`

Used for:
- behavior-only condition pack extraction
- training trends and choice analyses
- market/session indexing

Critical fields referenced in docs/scripts:
- `option` (1=Buy, 2=Hold, 3=Sell)
- `SessionOfDay`
- `session`
- monkey identity fields
- condition/drug/context fields
- timing and reward-related fields (`fb1`, `fb2` often used as timing anchors)

### 2) Recovery Workbook(s)
Typical filenames:
- `MSM_neuronSortMASTER_recoveryFile.xlsx`
- modified variants used in historical scripts

Used for:
- event extraction workflows
- neuron/session reconciliation

Caution:
- treat as a clue trail for reconciliation, not guaranteed source-of-truth for reward magnitude modeling

### 3) Neural Files
Typical format:
- `.nex`
- `.pl2`

Used for:
- spike/event extraction
- neural GLM and PSTH workflows

Current local copy in this project:
- `data/neural/raw_external/MSM_Sorted` (copied from `E:\MSM_Sorted` on 2026-02-16)

## Generated-but-Required Intermediates
Many downstream scripts expect condition packs generated upstream.
Typical generated files:
- `*_condition_pack.mat`

These should be generated locally and stored outside git.

## Data Quality Risks To Track
- inconsistent condition labeling across historical files
- mixed opponent policy behavior in AI condition
- condition rollout confounded with training date
- monkey-level behavioral heterogeneity

## Minimal Data Validation Checklist
1. Validate required input files exist at configured local paths.
2. Confirm key columns/variables exist with expected names.
3. Confirm condition labels and OT/saline encodings match script assumptions.
4. Confirm session indexing logic remains stable after relocation.
