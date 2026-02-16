# MSM Project Root Policy

Last updated: 2026-02-16

## Canonical Locations
- Canonical code root: `C:\Users\plattlab\MSM`
- Legacy archive root (read-only): `C:\Users\plattlab\Tim\Stock_market_experiment`

## Rule of Operation
1. All new code changes must happen only in `C:\Users\plattlab\MSM`.
2. Do not edit scripts in legacy locations.
3. Legacy locations are preserved as provenance snapshots.
4. Raw data are copied into `MSM/data` by modality for one-stop execution.
5. Generated outputs should go to non-versioned output folders.

## Why
This prevents split-brain development (different code versions in multiple roots) and keeps future human/AI handoffs deterministic.

## Session Start Checklist
1. Open repo at `C:\Users\plattlab\MSM`.
2. Confirm working directory is `MSM`.
3. Confirm data roots in `docs/PATHS.local.example.md`.
4. Run analyses only from this root.

## Session End Checklist
1. Save script edits in `MSM`.
2. Record what was run and where outputs were written.
3. Do not back-port edits to legacy paths.
