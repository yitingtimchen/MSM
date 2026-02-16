# Data Copy Manifest (2026-02-16)

## Operation
Copied data from legacy source:
- `C:\Users\plattlab\Tim\Stock_market_experiment`

Destination:
- `C:\Users\plattlab\MSM\data`

Mode:
- copy only (no move/delete in legacy source)

## Source -> Destination mapping
- `AWH2` -> `data/behavioral/AWH2`
- `Human data` -> `data/human/Human data`
- `MATLAB_Copy` -> `data/tooling/MATLAB_Copy`
- `MSM notes` -> `data/reference/MSM notes`
- `MSM_092017` -> `data/provenance/MSM_092017`
- `MSM_neuronSortMASTER_recoveryFile.xlsx` -> `data/metadata/`
- `MSM_neuronSortMASTER_recoveryFile_modified.xlsx` -> `data/metadata/`
- `MSM_neuronSortMASTER_recoveryFile_modified_filtered.xlsx` -> `data/metadata/`
- `SfN_Poster2019_AWH.pdf` -> `data/reference/`
- `E:\MSM_Sorted` -> `data/neural/raw_external/MSM_Sorted`

## Neural raw extension scan
Extensions searched:
- `.nex`, `.nev`, `.ns*`, `.edf`, `.plx`, `.pl2`, `.smr`, `.ced`, `.h5`, `.hdf5`

Detected files (3):
- `MATLAB_Copy/chronux_2_11/dataio/HowToReadNexFilesInMatlab/test.nex`
- `MATLAB_Copy/sigTOOL/sigTOOL Neuroscience Toolkit/File/menu_Import/group_NeuroScience File Formats/plx/explore1.plx`
- `MATLAB_Copy/sigTOOL/sigTOOL Neuroscience Toolkit/File/menu_Import/group_NeuroScience File Formats/plx/postfwd5.plx`

Mirrored to:
- `data/neural/raw/...` with source-relative paths preserved.

## External neural archive copy (2026-02-16)
Source:
- `E:\MSM_Sorted`

Destination:
- `data/neural/raw_external/MSM_Sorted`

Copy method:
- `rsync -a`

Result:
- 316 files copied
- approx. 212 GB
- extension breakdown: 195 `*.nex`, 118 `*.pl2`

Notes:
- Source drive content was copied only; no source deletion/modification was performed.

## Integrity spot-check
- Copy completed without reported rsync errors.
- Legacy source location was not modified by this operation.
