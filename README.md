# Monkey Stock Market (MSM) Analysis

Last updated: 2026-02-13

## Project Goal
This project analyzes behavior and neural activity in a monkey stock-market task to test:
- how social context changes choice behavior (`Buy` / `Hold` / `Sell`)
- whether oxytocin (OT) modulates social effects
- whether recorded neurons (mSTS-focused in prior claims) encode social context, choice, and outcome-related signals

## Current Status
This repository contains multiple generations of scripts. The active branch of work is behavior-first analysis using condition packs from:
- `condition_packs_behavior_only_from_dataTabFULL`

Recent analysis emphasis (Feb 2026):
- strategy benchmarks (dynamic programming vs fixed-fundamental)
- training-related shifts in action rates and inventory drift
- theoretical reward decomposition (instead of unreliable spreadsheet reward magnitudes)
- opponent/repetition and social-context follow-ups

## Core Data and Artifacts
Main sources:
- trial-level spreadsheet (market, choice, portfolio, condition/drug, timing fields)
- spike-sorted `.nex` files for neural reconstruction
- `MSM_neuronSortMASTER_recoveryFile.xlsx` (use as clue trail only; not source of truth)
- prior draft/poster artifacts (hypothesis references, not authoritative analysis code)

## Key Field Semantics (important)
- `option`: monkey action code (`1=Buy`, `2=Hold`, `3=Sell`)
- `SessionOfDay`: day-level grouping
- `session`: market instance within day
- recommended market ID:
  - `(monkey, year, month, day, SessionOfDay, session)`

Do not assume poster naming matches spreadsheet semantics.

## Known Confounds / Blockers
- Spreadsheet reward magnitudes (`fb1`/`fb2`) appear unreliable for value magnitude analyses.
- AI condition includes mixed opponent policies (not always iid random).
- Condition rollout is time-confounded with learning/training.
- Behavioral differences across monkeys are non-trivial.
- Historical condition workbook likely has labeling inconsistencies.
- Neural sampling is imbalanced across monkeys/roles in prior selections.

## What Has Been Implemented (to date)
- recovered session/market semantics from spreadsheet
- confirmed primary action coding (`option`)
- identified non-random AI behavior requiring subtype analysis
- shifted reward analyses to theoretical reward/value from task variables
- flagged time/experience confounding in condition comparisons
- implemented normative benchmarks:
  - state-dependent dynamic programming (DP) policy
  - fixed-fundamental heuristic policy
- found monkey behavior does not match either benchmark in initial tests
- observed training trend: `Buy` up, `Hold` down, slight `Sell` down, inventory drift up
- theoretical reward decomposition: dividend-related reward increases with training while trade-related reward decreases
- price sensitivity checks: strong training effect on `P(Buy)`, weak/inconsistent `priceBuy` effect, no robust `priceBuy x training` interaction
- initial inventory-setpoint checks did not support simple buy-until-K rule
- Live-condition role split (subject-first vs opponent-second) was analyzed in recent scripts (`260209+`): average `P(Buy)` is driven mostly by the first mover, while the second mover tends to follow.
- Live-condition asymmetry: overall `P(Buy)` is higher when M1 goes first and lower when M2 goes first.

## Active Scripts (Current Analysis Suite)
Primary recent scripts include:
- `find_best_strategy_per_market_260209.m`
- `optimal_choice_analysis_260212.m`
- `pBuy_vs_priceBuy_x_training_live_split_roles_260212.m`
- `pBuy_vs_portfolio_size_260212.m`
- `pBuy_control_trialNum_marketOrig_260212.m`
- `theoretical_tot_rew1_vs_rew2_vs_training_260212.m`
- `more_buy_or_less_sell_260213.m`
- `repeat_opp_analysis_251214.m`

Reference workflow map:
- `WORKFLOW_MAP.md`

## Recommended Run Order (Behavior-First)
If condition packs already exist:
1. `plot_condition_by_date_251204.m` (cohort/time sanity check)
2. `plot_choice_distributions_251205.m` and `GLM_behav_only_251210.m` (baseline behavior)
3. `RL_basline_and_saline_separated_251217.m` (RL comparisons)
4. run the active strategy suite listed above

## Path and Reproducibility Checklist
Before running scripts, verify hardcoded paths and cohort settings:
- `OUTDIR` points to intended condition-pack folder
- required input files exist (main spreadsheet and any recovery files)
- exclusion lists match intended cohort
- outputs are written to dated folders to avoid overwrite ambiguity

## Immediate Next Steps
1. Build/validate clean market index and per-market consistency checks.
2. Reconstruct opponent behavior features and split AI by policy subtype.
3. Continue value/reward analyses using theoretical definitions (timing markers only from `fb1`/`fb2` fields).
4. Fit behavior models with condition + OT + social + AI-policy + market-template + experience controls and monkey factor.
5. Test inventory dependence controlling for `trialNum` and market template.
6. Reconstruct/rebalance neural dataset from `.nex` and test social/context/choice/OT effects with covariate controls.

## Notes on Legacy Folders
Treat as historical unless needed for provenance:
- `AM_orig_codes_made_run_250602/`
- `bad codes/`

Use caution with older extraction scripts that target alternate output folders or contain outdated assumptions.
