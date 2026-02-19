# Pipeline Runbook

## Working Rule
Use only the canonical behavior extractor and keep neural extraction as a separate branch.
All active scripts are in `scripts/active/`.

## A. Build/Refresh Packs (when needed)
1. `scripts/active/d20251205_extract_behavior_condition_packs_with_split.m`
2. `scripts/active/d20251118_extract_spike_event_condition_packs.m` (only when neural/event packs are required)

## B. Baseline Behavior and RL
1. `scripts/active/d20251204_plot_condition_by_date.m`
2. `scripts/active/d20251205_plot_choice_distributions.m`
3. `scripts/active/d20251210_fit_glm_behavior_only.m`
4. `scripts/active/d20251126_fit_rl_by_subject_and_condition.m`
5. `scripts/active/d20251217_fit_rl_baseline_vs_saline_separate.m`

## C. Current Strategy/Training Suite
1. `scripts/active/d20260209_find_best_strategy_per_market.m`
2. `scripts/active/d20260212_analyze_optimal_choice.m`
3. `scripts/active/d20260212_analyze_pbuy_vs_price_over_training_by_live_role.m`
4. `scripts/active/d20260212_analyze_pbuy_vs_portfolio_size.m`
5. `scripts/active/d20260212_analyze_pbuy_with_trialnum_market_control.m`
6. `scripts/active/d20260212_analyze_theoretical_total_reward_over_training.m`
7. `scripts/active/d20260213_analyze_buy_vs_sell_action_trends.m`
8. `scripts/active/d20251214_analyze_repeat_and_opponent_effects.m`

## D. Neural/Peri-Event Branch
- `scripts/active/d20251123_fit_glme_with_neural_features.m`
- `scripts/active/d20251124_plot_psth_multi_event.m`
- `scripts/active/group_psth_normrange_by_option_CI*.m`
- `scripts/active/d20250904_plot_m1_reward_onset_firing_rate_baseline_spread_by_condition.m`

## E. Pre-Run Checklist Per Script
1. Confirm local input paths
2. Confirm output directory
3. Confirm cohort exclusions
4. Confirm condition pack source and date
