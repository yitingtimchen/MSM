# Data Layout in MSM

Last updated: 2026-02-16

## Canonical local data root
- `C:\Users\plattlab\MSM\data`

## Modality folders
- `data/behavioral/AWH2`
- `data/human/Human data`
- `data/metadata`
- `data/neural/raw`
- `data/neural/raw_external/MSM_Sorted`
  - `data/neural/raw_external/MSM_Sorted/nex`
  - `data/neural/raw_external/MSM_Sorted/pl2`
  - `data/neural/raw_external/MSM_Sorted/avi`
  - `data/neural/raw_external/MSM_Sorted/mat`
- `data/tooling/MATLAB_Copy`
- `data/reference`
- `data/provenance/MSM_092017`

## Sizing snapshot
- `data/behavioral`: ~290 MB
- `data/human`: ~1.2 GB
- `data/metadata`: ~140 KB
- `data/neural`: ~212 GB
- `data/tooling`: ~673 MB
- `data/reference`: ~277 MB
- `data/provenance`: ~728 KB
- `data` total: ~214 GB

## Raw neural files detected in legacy root
Detected explicit neural/ephys formats (`.nex`, `.plx`, `.nev`, `.ns*`, `.edf`, `.pl2`, `.smr`, `.h5`, `.hdf5`):
- 3 files total
- mirrored under `data/neural/raw/...`

If additional raw neural data exist elsewhere (other drives/paths or non-standard extensions), add them under:
- `data/neural/raw_external/<source_name>/`

and document in `docs/DATA_DEPENDENCIES.md`.
