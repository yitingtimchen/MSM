# Pipeline Runbook

## Working Rule
Use only the canonical behavior extractor and keep neural extraction as a separate branch.
All active scripts are in `scripts/active/`.

## A. Build/Refresh Packs (when needed)
1. `scripts/active/extract_and_save_behav_only_by_cond_and_split_251205.m`
2. `scripts/active/extract_and_save_spike_and_event_by_cond_251118.m` (only when neural/event packs are required)

## B. Baseline Behavior and RL
1. `scripts/active/plot_condition_by_date_251204.m`
2. `scripts/active/plot_choice_distributions_251205.m`
3. `scripts/active/GLM_behav_only_251210.m`
4. `scripts/active/RL_subject_251126.m`
5. `scripts/active/RL_basline_and_saline_separated_251217.m`

## C. Current Strategy/Training Suite
1. `scripts/active/find_best_strategy_per_market_260209.m`
2. `scripts/active/optimal_choice_analysis_260212.m`
3. `scripts/active/pBuy_vs_priceBuy_x_training_live_split_roles_260212.m`
4. `scripts/active/pBuy_vs_portfolio_size_260212.m`
5. `scripts/active/pBuy_control_trialNum_marketOrig_260212.m`
6. `scripts/active/theoretical_tot_rew1_vs_rew2_vs_training_260212.m`
7. `scripts/active/more_buy_or_less_sell_260213.m`
8. `scripts/active/repeat_opp_analysis_251214.m`

## D. Neural/Peri-Event Branch
- `scripts/active/GLM_with_neur_251123.m`
- `scripts/active/Plotting_PSTH_251124.m`
- `scripts/active/group_psth_normrange_by_option_CI*.m`
- `scripts/active/plot_m1_reward_onset_firing_rate_baseline_spread_across_cond_250904.m`

## E. Pre-Run Checklist Per Script
1. Confirm local input paths
2. Confirm output directory
3. Confirm cohort exclusions
4. Confirm condition pack source and date
