# Neural Source Comparison

Date: 2026-02-16  
Compared:
- `D:\Annamarie\StockMarketTask1_Monkeys\Plexon`
- `C:\Users\plattlab\MSM\data\neural\raw_external\MSM_Sorted`

## Executive Summary
These two locations are **not identical**. They have substantial overlap, but each contains files the other does not. The Annamarie source is larger and includes many video files.

## Counts
- Files in Annamarie source: **615**
- Files in MSM source: **316**
- Filename overlap: **149**
- Exact overlap (same filename + same size): **133**
- Same filename but different size: **16**
- Annamarie-only filenames: **466**
- MSM-only filenames: **167**

## File Type Profile
Annamarie source (top extensions):
- `.AVI`: 304
- `.pl2`: 244
- `.nex`: 67

MSM source (top extensions):
- `.nex`: 195
- `.pl2`: 118
- `.mat`: 2
- `.AVI`: 1

## Integrity Spot-Check
A hash check on overlapping files was performed:
- 20 sampled exact filename+size matches
- **20/20 hashes matched** (byte-identical in sample)

## Important Difference Class
There are 16 files with the same filename but different sizes (all `.nex` in the sampled list), indicating probable alternate sorting/export versions. Example names:
- `051218L1-sorted-01.nex`
- `052818B1-sorted-01.nex`
- `060818L1-sorted-01.nex`
- `080318B2-sorted.nex`

## Practical Interpretation
- Treat `D:\Annamarie\StockMarketTask1_Monkeys\Plexon` as an additional source, not a duplicate.
- If publication completeness is the goal, preserve both sources and maintain provenance.
- Resolve same-name/different-size `.nex` files by explicit version policy before merging.
