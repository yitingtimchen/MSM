# Stock Market Experiment Workflow Map

This file is a practical map for `C:\Users\plattlab\Tim\Stock_market_experiment\Tim`.

## 1) What This Folder Contains

You have three major generations of analysis:

1. Event/spike extraction and neuron-linked analyses (older core):
- `extract_and_save_spike_and_event_by_cond_251118.m`
- `GLM_with_neur_251123.m`
- `RL_and_GLM_251124.m`
- Uses `condition_packs_251118`

2. Behavior-only pack building and GLM/RL:
- `extract_and_save_event_only_by_cond_from_recovFile_250923.m`
- `extract_and_save_behav_only_by_cond_from_dataTabFULL_251201.m`
- `extract_and_save_behav_only_by_cond_and_split_251205.m`
- `GLM_behav_only_251210.m`
- `RL_subject_251126.m`
- `RL_basline_and_saline_separated_251217.m`

3. New MSM strategy analyses (current/latest scripts):
- `find_best_strategy_per_market_260209.m`
- `pBuy_vs_priceBuy_x_training_live_split_roles_260212.m`
- `pBuy_vs_portfolio_size_260212.m`
- `pBuy_control_trialNum_marketOrig_260212.m`
- `optimal_choice_analysis_260212.m`
- `theoretical_tot_rew1_vs_rew2_vs_training_260212.m`
- `more_buy_or_less_sell_260213.m`
- `repeat_opp_analysis_251214.m`
- Mostly use `condition_packs_behavior_only_from_dataTabFULL`

## 2) Canonical Directories (Use These First)

Primary behavior-only analysis directory:
- `condition_packs_behavior_only_from_dataTabFULL`

Older snapshots/alternatives (keep for provenance, not default):
- `condition_packs_behavior_only_v2`
- `condition_packs_behavior_only_from_dataTabFULL_manualSessions`
- `condition_packs_251118`
- `condition_packs_250922`
- `condition_packs`

## 3) Recommended Run Order (Default in 2026)

If packs already exist and you only want analyses:

1. `plot_condition_by_date_251204.m` (sanity check session spread).
2. `plot_choice_distributions_251205.m` and `GLM_behav_only_251210.m` (behavior baseline).
3. `RL_basline_and_saline_separated_251217.m` (RL parameter/model comparisons).
4. MSM strategy suite:
- `find_best_strategy_per_market_260209.m`
- `optimal_choice_analysis_260212.m`
- `pBuy_vs_priceBuy_x_training_live_split_roles_260212.m`
- `pBuy_vs_portfolio_size_260212.m`
- `pBuy_control_trialNum_marketOrig_260212.m`
- `theoretical_tot_rew1_vs_rew2_vs_training_260212.m`
- `more_buy_or_less_sell_260213.m`
- `repeat_opp_analysis_251214.m`

## 4) Scripts To Treat As Historical / Not Default

Likely legacy snapshots:
- `AM_orig_codes_made_run_250602/*`
- `bad codes/*`

Use with caution:
- `extract_and_save_behav_only_by_cond_from_dataTabFULL_251201.m`
  - Header comment says it may not produce meaningful results.
  - Currently points OUTDIR to `condition_packs_behavior_only_from_dataTabFULL_v2` (not your main active folder).

## 5) Path Hygiene Checklist (Before Running Anything)

Most scripts hardcode Windows absolute paths. Before each run, verify:

1. `OUTDIR` points to the intended pack folder.
2. Input files exist:
- `C:\Users\plattlab\Tim\Stock_market_experiment\MATLAB_Copy\dataTabFULL.xlsx`
- Recovery spreadsheet paths if running extraction scripts.
3. `EXCLUDE_SESSIONS` matches the exact analysis cohort you want.

## 6) Naming Convention Note

For script names like `*_260212.m`, suffix is typically `yymmdd` (here: 2026-02-12).  
In general, prefer the latest date-stamped script within a topic unless a newer file says otherwise in its header.

## 7) Safe Working Practice Going Forward

1. Keep one canonical pack folder per analysis branch.
2. At the top of each script, keep:
- purpose
- required inputs
- output files
- last validated pack directory
3. Save new outputs into dated subfolders to avoid overwrite confusion.

