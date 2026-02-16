# Local Paths (Example Template)

Copy this to `docs/PATHS.local.md` and edit for your machine.

## Canonical Project Paths
- `MSM_CODE_ROOT`: `C:\Users\plattlab\MSM`
- `MSM_DATA_ROOT`: `C:\Users\plattlab\MSM\data`
- `MSM_OUTPUT_ROOT`: `C:\Users\plattlab\MSM\outputs_local`

## Data By Modality (Current MSM Layout)
- Behavioral/raw logs: `MSM_DATA_ROOT\\behavioral`
- Neural/raw or neural-adjacent sources: `MSM_DATA_ROOT\\neural`
- Recovery workbooks and metadata: `MSM_DATA_ROOT\\metadata`
- Legacy MATLAB resources/toolchains: `MSM_DATA_ROOT\\tooling`
- Human-study branch data: `MSM_DATA_ROOT\\human`
- Reference docs/poster/provenance materials: `MSM_DATA_ROOT\\reference`

## Notes
- `MSM_DATA_ROOT` is intentionally ignored by git.
- Keep large outputs out of git as well.
