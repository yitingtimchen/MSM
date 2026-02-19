% File: d20260213_plot_psth_choice_reward_sustained_subset.m



cd('C:\Users\plattlab\MSM\outputs_local\condition_packs_250922')

load("Live_condition_pack.mat")
%%%%%% for Choice–Reward Sustained (76 neurons)
sessions1 = [3 3 3 3  3  5 5 5 5  6 6 6 6 6  6  7 7 7 7 7  7  8 8 8 8 8 8  8  8  8  8  8  8  8  9 9 9 9 9 9 9  9  9  9  9  9  9  9  9  9  9  9  9  9  9  9  9  9  10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];  
neurons1  = [2 3 8 12 16 7 8 9 11 4 5 8 9 12 14 2 4 5 6 10 17 3 4 6 8 9 13 15 18 20 21 22 23 29 1 2 3 4 5 7 11 15 16 17 19 22 23 24 27 30 31 32 34 38 39 40 46 47 1  2  3  5  6  10 13 14 15 17 18 21 24 26 29 31 33 36];

sessions2 = [1 1 1 3 3 3 3  6  6  6  6  7  7  7  7  8 8  9 9 9  9  9  9  9  9  10 10 10 10 10 10 10 10 10]; 
neurons2 =  [1 2 3 1 4 5 24 10 16 17 18 12 13 14 15 7 19 8 9 18 20 21 25 28 43 8  11 20 22 23 25 30 38 39];

sessions = [sessions1 sessions2];
neurons = [neurons1 neurons2];

% if you want to align it to other events but keep the same baseline time window measurement
d20260213_group_psth_normalized_by_option_ci_alignable(C, sessions, neurons, 'm1rew1on',  20, 100, [-1500 1500], [-1000 -700], 'ci', 0.90, 'm1rew1on');
d20260213_group_psth_normalized_by_option_ci(C, sessions, neurons, 20, 100, [-1500 1500], [-1000 -700], 'ci', 0.90);




