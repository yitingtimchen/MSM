% File: d20250929_plot_psth_perievent_with_anova.m
% oink

close all; clc; clear;

PLOT = 0;
DO_ANOVA = 1;

BIN_SIZE = 20;     % ms
SMOOTH_WIN = 100;  % ms
BASELINE_WIN = [-1000, -700];
ERR_TYPE = 'ci';
CI_ALPHIA = 0.9;

% Define multiple TARGET events and a corresponding WIN per TARGET (rows match targets)
TARGETS = ["m1choice","m1rew1on","m1rew1off","m1rew2on","m1rew2off", ...
    "m2choice","m2rew1on","m2rew2on","m2rew2off"];
WINS    = repmat([-250 250], [length(TARGETS), 1]);

BASELINE_ANCHOR = "m1rew1on";
conds = {'AI', 'Replay', 'Decoy', 'Live'};

root_dir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_250922';

nT = numel(TARGETS);
nC = numel(conds);

% ===== One combined figure: 4 (conditions) x 5 (targets) =====
if PLOT
    f = figure('Color','w','Name','Mean normalized PSTH by Choice (4x5)');
    tl = tiledlayout(f, nC, nT, 'Padding','compact','TileSpacing','compact');
    sgtitle(tl, 'Mean normalized PSTH by Choice — rows=conditions, cols=events');
    
    for c = 1:nC
        S = load([root_dir '\' conds{c} '_condition_pack.mat']);
        C = S.C;
    
        fprintf("\nPlotting %s\n", conds{c});
    
        for t = 1:nT
            ax = nexttile(tl, (c-1)*nT + t);
            group_psth_normrange_by_option_CI_250925(ax, conds{c}, TARGETS(t), BASELINE_ANCHOR, C, ...
                BIN_SIZE, SMOOTH_WIN, WINS(t,:), BASELINE_WIN, ERR_TYPE, CI_ALPHIA);
    
            % Keep legend only on the first tile to reduce clutter
            if ~(c==1 && t==1)
                legend(ax, 'off');
            end
        end
    end

end

%% === One-way repeated-measures ANOVA per tile (BUY, HOLD, SELL) + posthoc ===

if DO_ANOVA
    
    ALPHA = 0.05; % per-comparison alpha (two-sided in ANOVA; pairwise are Tukey/Bonferroni two-sided)
    do_fdr = true; % apply Benjamini–Hochberg across all pairwise tests
    
    Results = struct([]);
    all_pair_p = []; all_pair_keys = strings(0,1);
    
    for c = 1:nC
    S = load([root_dir '\' conds{c} '_condition_pack.mat']); C = S.C;
    for t = 1:nT
    [bin_centers, N_buy_all, N_hold_all, N_sell_all] = collect_norm_timecourses_by_option_250929( ...
    C, TARGETS(t), BASELINE_ANCHOR, BIN_SIZE, WINS(t,:), BASELINE_WIN);
    
        yBuy  = mean(N_buy_all , 2, 'omitnan');
        yHold = mean(N_hold_all, 2, 'omitnan');
        ySell = mean(N_sell_all, 2, 'omitnan');
    
        ok = isfinite(yBuy) & isfinite(yHold) & isfinite(ySell);
        Y = [yBuy(ok), yHold(ok), ySell(ok)];
        nNeurons = size(Y,1);
    
        p_rm = NaN; F_rm = NaN; df1 = NaN; df2 = NaN; epsGG = NaN;
        p_BH = NaN; d_BH = NaN; p_BS = NaN; d_BS = NaN;
        
        if nNeurons >= 2
        T = array2table(Y, 'VariableNames', {'Buy','Hold','Sell'});
        within = table(categorical({'Buy';'Hold';'Sell'}), 'VariableNames', {'Condition'});
        rm = fitrm(T,'Buy-Sell ~ 1','WithinDesign',within);
        rtab = ranova(rm,'WithinModel','Condition');
        
        if ismember('pValueGG', rtab.Properties.VariableNames)
            p_rm = rtab.pValueGG(1);
        else
            p_rm = rtab.pValue(1);
        end
        F_rm = rtab.F(1);
        k = 3; df1 = k-1; df2 = (k-1)*(nNeurons-1);  % uncorrected dfs (robust & explicit)
        epsGG = NaN;  % optional: [epsGG,~,~] = epsilon(rm); epsGG = epsGG(1);
        
        ph = multcompare(rm,'Condition','ComparisonType','bonferroni');
        
        % ensure a generic 'Difference' column exists
        if ~ismember('Difference', ph.Properties.VariableNames)
            if ismember('MeanDifference', ph.Properties.VariableNames)
                ph.Difference = ph.MeanDifference;
            elseif ismember('EstDiff', ph.Properties.VariableNames)
                ph.Difference = ph.EstDiff;
            elseif ismember('Estimate', ph.Properties.VariableNames)
                ph.Difference = ph.Estimate;
            else
                error('No difference column in multcompare output.');
            end
        end
        
        % Buy vs Hold (sign so positive => Buy > Hold)
        iBH = (ph.Condition_1=="Buy" & ph.Condition_2=="Hold");
        if any(iBH)
            row = ph(iBH,:);
            d_BH = (row.Condition_1=="Buy").*row.Difference + (row.Condition_1=="Hold").*(-row.Difference);
            p_BH = row.pValue;
        end
        
        % Buy vs Sell (sign so positive => Buy > Sell)
        iBS = (ph.Condition_1=="Buy" & ph.Condition_2=="Sell"); 
        if any(iBS)
            row = ph(iBS,:);
            d_BS = (row.Condition_1=="Buy").*row.Difference + (row.Condition_1=="Sell").*(-row.Difference);
            p_BS = row.pValue;
        end
        
        
        end
    
        Results(c,t).cond   = string(conds{c});
        Results(c,t).target = string(TARGETS(t));
        Results(c,t).nNeurons = nNeurons;
        Results(c,t).p_rm   = p_rm;
        Results(c,t).F_rm   = F_rm;
        Results(c,t).df1    = df1;
        Results(c,t).df2    = df2;
        Results(c,t).epsGG  = epsGG;
        Results(c,t).p_Buy_vs_Hold = p_BH;
        Results(c,t).d_Buy_minus_Hold = d_BH;
        Results(c,t).p_Buy_vs_Sell = p_BS;
        Results(c,t).d_Buy_minus_Sell = d_BS;
    
        if isfinite(p_BH), all_pair_p(end+1,1) = p_BH; all_pair_keys(end+1,1) = sprintf('%s|%s|Buy>Hold',conds{c},TARGETS(t)); end
        if isfinite(p_BS), all_pair_p(end+1,1) = p_BS; all_pair_keys(end+1,1) = sprintf('%s|%s|Buy>Sell',conds{c},TARGETS(t)); end
    end
    
    
    end
    
    % FDR across ALL pairwise tests (optional)
    qmap = containers.Map;
    if do_fdr && ~isempty(all_pair_p)
    qvals = fdr_bh_local_250929(all_pair_p, 0.05);
    for i=1:numel(all_pair_p), qmap(all_pair_keys(i)) = qvals(i); end
    end
    
    % Summarize to table
    rows = [];
    for c = 1:nC
    for t = 1:nT
    R = Results(c,t);
    keyBH = sprintf('%s|%s|Buy>Hold',R.cond,R.target);
    keyBS = sprintf('%s|%s|Buy>Sell',R.cond,R.target);
    qBH = NaN; qBS = NaN;
    if do_fdr && isKey(qmap,keyBH), qBH = qmap(keyBH); end
    if do_fdr && isKey(qmap,keyBS), qBS = qmap(keyBS); end
    rows = [rows; {R.cond, R.target, R.nNeurons, R.p_rm, R.F_rm, R.df1, R.df2, ...
    R.p_Buy_vs_Hold, qBH, R.d_Buy_minus_Hold, ...
    R.p_Buy_vs_Sell, qBS, R.d_Buy_minus_Sell}]; %#ok<AGROW>
    end
    end
    SummaryTbl = cell2table(rows, 'VariableNames', ...
    {'Cond','Target','nNeurons','p_rm','F_rm','df1','df2', ...
    'p_Buy_vs_Hold','q_Buy_vs_Hold','effect_BuyMinusHold', ...
    'p_Buy_vs_Sell','q_Buy_vs_Sell','effect_BuyMinusSell'});
    
    filteredTbl = SummaryTbl(SummaryTbl.q_Buy_vs_Hold < 0.05 & SummaryTbl.q_Buy_vs_Sell < 0.05 & SummaryTbl.effect_BuyMinusHold > 0.01 & SummaryTbl.effect_BuyMinusSell > 0.01, :);
    filteredTbl = sortrows(filteredTbl, ["effect_BuyMinusHold", "effect_BuyMinusSell"], 'descend');
    disp(filteredTbl);
    % You can flag significance, e.g. q < 0.05 & effect > 0 meaning BUY > other.
end

% ======================== FUNCTIONS ========================

function [bin_centers, N_buy_all, N_hold_all, N_sell_all] = collect_norm_timecourses_by_option_250929( ...
C, target_evt, baseline_evt, bin_ms, win, baseline_win)

edges = win(1):bin_ms:win(2);
bin_centers = edges(1:end-1) + bin_ms/2;
nBins = numel(edges) - 1;

baseline_edges = baseline_win(1):bin_ms:baseline_win(2);
baseline_bin_centers = baseline_edges(1:end-1) + bin_ms/2;
nBins_baseline = numel(baseline_edges) - 1;

% Locate event table
if isfield(C, 'eventTables')
E = C.eventTables;
elseif isfield(C, 'evenTables')
E = C.evenTables;
else
error('Neither C.eventTables nor C.evenTables found.');
end

N_buy_all = [];
N_hold_all = [];
N_sell_all = [];

base_mask = (baseline_bin_centers >= baseline_win(1)) & (baseline_bin_centers < baseline_win(2));
if ~any(base_mask), error('Baseline mask empty. Check baseline_win vs bin_ms.'); end

bin_sec = bin_ms / 1000;

for s = 1:numel(C.files)
if ~iscell(E) || size(E,1) < s || ~istable(E{s,1}), continue; end
varnames = string(E{s,1}.Properties.VariableNames);
if sum(ismember(["option", string(target_evt), string(baseline_evt)], varnames)) ~= 3, continue; end

% Target event times
evtCell = E{s,1}.(target_evt){1,1};
if istable(evtCell), evt_ms = table2array(evtCell); else, evt_ms = evtCell; end
evt_ms_vec = double(evt_ms(:));

% Baseline anchor event times
refCell = E{s,1}.(baseline_evt){1,1};
if istable(refCell), baseline_ms = table2array(refCell); else, baseline_ms = refCell; end
baseline_ms_vec = double(baseline_ms(:));

% Choice codes
optCell = E{s,1}.option{1,1};
if istable(optCell), opt = table2array(optCell); else, opt = optCell; end
opt_vec = double(opt(:));

% Align & clean
n = min([numel(evt_ms_vec), numel(opt_vec), numel(baseline_ms_vec)]);
evt_ms_vec      = evt_ms_vec(1:n);
opt_vec         = opt_vec(1:n);
baseline_ms_vec = baseline_ms_vec(1:n);
valid = ~isnan(evt_ms_vec) & ~isnan(opt_vec) & ~isnan(baseline_ms_vec);
evt_ms_vec      = evt_ms_vec(valid);
opt_vec         = opt_vec(valid);
baseline_ms_vec = baseline_ms_vec(valid);
nTrials = numel(evt_ms_vec);
if nTrials == 0, continue; end

neuType = C.neuronType{s,1};
sustained_neu_idx = (neuType == 1) | (neuType == 2);
rawSess_sustained = C.rawSpikes{s,1}(sustained_neu_idx);
if isempty(rawSess_sustained), continue; end

for neu = 1:numel(rawSess_sustained)
    spikes_sec = rawSess_sustained{neu};
    spikes_ms = double(spikes_sec) * 1000;

    % Target event bins
    counts = zeros(nTrials, nBins);
    for tt = 1:nTrials
        rel_ms = spikes_ms - evt_ms_vec(tt);
        rel_win = rel_ms(rel_ms >= win(1) & rel_ms < win(2));
        counts(tt,:) = histcounts(rel_win, edges);
    end
    fr_trials = counts / bin_sec;

    % Baseline bins
    baseline_counts = zeros(nTrials, nBins_baseline);
    for tt = 1:nTrials
        rel_ms_b = spikes_ms - baseline_ms_vec(tt);
        rel_win_b = rel_ms_b(rel_ms_b >= baseline_win(1) & rel_ms_b < baseline_win(2));
        baseline_counts(tt,:) = histcounts(rel_win_b, baseline_edges);
    end
    baseline_fr_trials = baseline_counts / bin_sec;

    % Baseline subtraction (scalar per neuron)
    mu_base_neuron = mean(baseline_fr_trials(:, base_mask), 'all', 'omitnan');
    fr_bs = fr_trials - mu_base_neuron;

    % Global min–max per neuron (across all trials x bins)
    vec = fr_bs(:); mn = min(vec); mx = max(vec); denom = mx - mn; if ~isfinite(denom) || denom <= 0, denom = eps; end
    norm_trials = reshape((vec - mn) / denom, size(fr_bs));

    idx_buy  = (opt_vec == 1);
    idx_hold = (opt_vec == 2);
    idx_sell = (opt_vec == 3);

    if any(idx_buy),  N_buy_all  = [N_buy_all ; mean(norm_trials(idx_buy ,:), 1, 'omitnan')]; end
    if any(idx_hold), N_hold_all = [N_hold_all; mean(norm_trials(idx_hold,:), 1, 'omitnan')]; end
    if any(idx_sell), N_sell_all = [N_sell_all; mean(norm_trials(idx_sell,:), 1, 'omitnan')]; end
end


end
end

function q = fdr_bh_local_250929(p, qlevel)
    % Benjamini–Hochberg (independent/positive dependency)
    p = p(:); n = numel(p);
    [ps, idx] = sort(p);
    th = (1:n)' * (qlevel / n);
    rej = ps <= th;
    k = find(rej, 1, 'last');
    q = NaN(size(p));
    if ~isempty(k)
    cutoff = ps(k);
    q(idx) = min(1, ps * n ./ (1:n)');
    % Monotone nonincreasing adjustment
    q(idx) = flipud(cummin(flipud(q(idx))));
    else
    q(:) = 1;
    end
end

function group_psth_normrange_by_option_CI_250925(ax, cond, target_evt, baseline_evt, C, bin_ms, smooth_ms, win, baseline_win, error_type, ci_alpha)
% Baseline-subtract per neuron (pooled across trials), then global min–max
% normalize via range across all trials x time per neuron.
% Average per condition and across neurons, and plot Buy/Hold/Sell with SEM or CI bands into provided axes.

% ---- Binning setup ----
bin_sec = bin_ms / 1000;

edges = win(1):bin_ms:win(2);
bin_centers = edges(1:end-1) + bin_ms/2;
nBins = numel(edges) - 1;

baseline_edges = baseline_win(1):bin_ms:baseline_win(2);
baseline_bin_centers = baseline_edges(1:end-1) + bin_ms/2;
nBins_baseline = numel(baseline_edges) - 1;

% ---- Locate event table ----
if isfield(C, 'eventTables')
    E = C.eventTables;
elseif isfield(C, 'evenTables')
    E = C.evenTables;
else
    error('Neither C.eventTables nor C.evenTables found.');
end

% ---- Storage for per-neuron condition means (normalized) ----
N_buy_all  = [];
N_hold_all = [];
N_sell_all = [];

% track total trials per condition across neurons (for legend)
tot_buy = 0; tot_hold = 0; tot_sell = 0;

% ---- Baseline mask ----
base_mask = (baseline_bin_centers >= baseline_win(1)) & (baseline_bin_centers < baseline_win(2));
if ~any(base_mask)
    error('Baseline mask empty. Check baseline_win vs bin_ms.');
end

for s = 1:numel(C.files)
    % check required columns exist
    if ~iscell(E) || size(E,1) < s || ~istable(E{s,1})
        continue
    end
    varnames = string(E{s,1}.Properties.VariableNames);
    if sum(ismember(["option", string(target_evt), string(baseline_evt)], varnames)) ~= 3
        continue
    end

    % ----- Target event times (ms) -----
    evtCell = E{s,1}.(target_evt){1,1};
    if istable(evtCell) || isa(evtCell,'table'), evt_ms = table2array(evtCell);
    else,                                        evt_ms = evtCell;
    end
    evt_ms_vec = double(evt_ms(:));

    % ----- BASELINE_ANCHOR event times (ms) -----
    refCell = E{s,1}.(baseline_evt){1,1};
    if istable(refCell) || isa(refCell,'table'), baseline_ms = table2array(refCell);
    else,                                        baseline_ms = refCell;
    end
    baseline_ms_vec = double(baseline_ms(:));

    % ----- Choice codes (1=Buy,2=Hold,3=Sell) -----
    optCell = E{s,1}.option{1,1};
    if istable(optCell) || isa(optCell,'table'), opt = table2array(optCell);
    else,                                        opt = optCell;
    end
    opt_vec = double(opt(:));

    % ----- Align & clean -----
    n = min([numel(evt_ms_vec), numel(opt_vec), numel(baseline_ms_vec)]);
    evt_ms_vec      = evt_ms_vec(1:n);
    opt_vec         = opt_vec(1:n);
    baseline_ms_vec = baseline_ms_vec(1:n);
    valid = ~isnan(evt_ms_vec) & ~isnan(opt_vec) & ~isnan(baseline_ms_vec);
    evt_ms_vec      = evt_ms_vec(valid);
    opt_vec         = opt_vec(valid);
    baseline_ms_vec = baseline_ms_vec(valid);
    nTrials = numel(evt_ms_vec);
    if nTrials == 0
        continue
    end

    % ----- Spikes (sec -> ms) -----
    neuType = C.neuronType{s,1};
    sustained_neu_idx  = (neuType == 1) | (neuType == 2);
    % suppressed_neu_idx = (neuType == 3); % reserved if needed

    rawSess_sustained = C.rawSpikes{s,1}(sustained_neu_idx);
    if isempty(rawSess_sustained)
        continue
    end

    for neu = 1:numel(rawSess_sustained)
        spikes_sec = rawSess_sustained{neu};
        spikes_ms = double(spikes_sec) * 1000;

        % ----- Target event: Bin spikes per trial -> FR (Hz) -----
        counts = zeros(nTrials, nBins);
        for t = 1:nTrials
            rel_ms  = spikes_ms - evt_ms_vec(t);
            rel_win = rel_ms(rel_ms >= win(1) & rel_ms < win(2));
            counts(t,:) = histcounts(rel_win, edges);
        end
        fr_trials = counts / bin_sec;  % [trials x bins], Hz

        % ----- BASELINE_ANCHOR event: Bin spikes per trial -> FR (Hz) -----
        baseline_counts = zeros(nTrials, nBins_baseline);
        for t = 1:nTrials
            baseline_rel_ms  = spikes_ms - baseline_ms_vec(t);
            baseline_rel_win = baseline_rel_ms(baseline_rel_ms >= baseline_win(1) & baseline_rel_ms < baseline_win(2));
            baseline_counts(t,:) = histcounts(baseline_rel_win, baseline_edges);
        end
        baseline_fr_trials = baseline_counts / bin_sec;  % [trials x bins], Hz

        % ----- Pooled baseline subtraction (per neuron) -----
        mu_base_neuron = mean(baseline_fr_trials(:, base_mask), 'all', 'omitnan');  % scalar Hz
        fr_bs = fr_trials - mu_base_neuron;  % baseline-subtracted FR

        % ----- Global min–max normalization across ALL trials x bins (per neuron) -----
        vec = fr_bs(:);
        mn = min(vec); mx = max(vec); denom = mx - mn;
        if ~isfinite(denom) || denom <= 0, denom = eps; end
        vec = (vec - mn) / denom;
        norm_trials = reshape(vec, size(fr_bs));  % [trials x bins]

        % ----- Per-condition mean (per neuron) -----
        idx_buy  = (opt_vec == 1);
        idx_hold = (opt_vec == 2);
        idx_sell = (opt_vec == 3);

        if any(idx_buy),  N_buy_all  = [N_buy_all;  mean(norm_trials(idx_buy ,:), 1, 'omitnan')];  tot_buy  = tot_buy  + sum(idx_buy);  end
        if any(idx_hold), N_hold_all = [N_hold_all; mean(norm_trials(idx_hold,:), 1, 'omitnan')];  tot_hold = tot_hold + sum(idx_hold); end
        if any(idx_sell), N_sell_all = [N_sell_all; mean(norm_trials(idx_sell,:), 1, 'omitnan')];  tot_sell = tot_sell + sum(idx_sell); end
    end
end

% ----- Across-neuron means -----
N_buy_mean  = mean(N_buy_all , 1, 'omitnan');
N_hold_mean = mean(N_hold_all, 1, 'omitnan');
N_sell_mean = mean(N_sell_all, 1, 'omitnan');

% ----- Error bands across neurons (SEM or CI) -----
[buy_lo,  buy_hi ] = compute_bands(N_buy_all , error_type, ci_alpha);
[hold_lo, hold_hi] = compute_bands(N_hold_all, error_type, ci_alpha);
[sell_lo, sell_hi] = compute_bands(N_sell_all, error_type, ci_alpha);

% ----- Smoothing -----
smooth_bins = max(3, round(smooth_ms / bin_ms));
if exist('smoothdata','file') == 2
    smoothf = @(x) smoothdata(x,'gaussian',smooth_bins);
else
    w = ones(1, smooth_bins) / max(1, smooth_bins);
    smoothf = @(x) conv(x, w, 'same');
end
N_buy_s   = smoothf(N_buy_mean );
N_hold_s  = smoothf(N_hold_mean);
N_sell_s  = smoothf(N_sell_mean);
buy_lo_s  = smoothf(buy_lo );  buy_hi_s  = smoothf(buy_hi );
hold_lo_s = smoothf(hold_lo);  hold_hi_s = smoothf(hold_hi);
sell_lo_s = smoothf(sell_lo);  sell_hi_s = smoothf(sell_hi);

% ----- Plot into provided axes -----
axes(ax); cla(ax); hold(ax,'on');

% shaded error (if requested)
if ~strcmpi(error_type,'none')
    shaded(ax, bin_centers, buy_lo_s , buy_hi_s , [0 0.6 0]);   % green
    shaded(ax, bin_centers, hold_lo_s, hold_hi_s, [0 0.4 1]);   % blue
    shaded(ax, bin_centers, sell_lo_s, sell_hi_s, [0.85 0 0]);  % red
end

% mean lines
h = gobjects(0);
if any(isfinite(N_buy_s)),  h(end+1) = plot(ax, bin_centers, N_buy_s , 'g', 'LineWidth', 2.0); end
if any(isfinite(N_hold_s)), h(end+1) = plot(ax, bin_centers, N_hold_s, 'b', 'LineWidth', 2.0); end
if any(isfinite(N_sell_s)), h(end+1) = plot(ax, bin_centers, N_sell_s, 'r', 'LineWidth', 2.0); end

% Y limits from all lines and bands
y_all = [N_buy_s, N_hold_s, N_sell_s, buy_lo_s, buy_hi_s, hold_lo_s, hold_hi_s, sell_lo_s, sell_hi_s];
y_all = y_all(isfinite(y_all));
if isempty(y_all)
    ymin = 0; ymax = 1;
else
    ymin = min(y_all); ymax = max(y_all);
end
pad  = 0.05 * max(1e-6, ymax - ymin);
% ylim(ax, [ymin - pad, ymax + pad]);
ylim(ax, [0.05, 0.2]);

% X limits, zero line, labels
xlim(ax, win);
if exist('xline','file') == 2
    xline(ax, 0, '--', 'LineWidth', 1);
else
    yl = ylim(ax); plot(ax, [0 0], yl, '--', 'LineWidth', 1);
end
grid(ax, 'on');
% xlabel(ax, sprintf('Time from %s (ms) - smoothing approx %d ms', target_evt, smooth_ms));
% ylabel(ax, 'Mean normalized activity (0-1), baseline-subtracted');
ylabel(ax, cond, FontSize=10, FontWeight="bold")
if string(cond) == "AI"
    title(ax, sprintf('Align: %s', target_evt), FontSize=12, FontWeight="bold");
end
if ~isempty(h)
    legend(ax, h, {sprintf('Buy (trials=%d)',tot_buy), ...
                   sprintf('Hold (trials=%d)',tot_hold), ...
                   sprintf('Sell (trials=%d)',tot_sell)}, 'Location','best');
end
end

function [lo, hi] = compute_bands(M, error_type, ci_alpha)
% M is nNeurons x nBins of per-neuron means (in [0,1])
if nargin < 2 || isempty(error_type), error_type = 'ci'; end
if nargin < 3 || isempty(ci_alpha),   ci_alpha = 0.95; end
mu  = mean(M, 1, 'omitnan');
sd  = std(M, 0, 1, 'omitnan');
n   = sum(isfinite(M), 1);
sem = sd ./ max(1, sqrt(n));
switch lower(error_type)
    case 'sem'
        lo = mu - sem;
        hi = mu + sem;
    case 'ci'
        alpha2 = (1 - ci_alpha)/2;
        lo = nan(size(mu)); hi = nan(size(mu));
        for j = 1:numel(mu)
            if n(j) >= 2 && isfinite(sem(j))
                tcrit = tinv(1 - alpha2, n(j) - 1);
                lo(j) = mu(j) - tcrit * sem(j);
                hi(j) = mu(j) + tcrit * sem(j);
            end
        end
    otherwise
        lo = nan(size(mu)); hi = nan(size(mu));
end
end

function shaded(ax, t, lo, hi, col)
mask = isfinite(t) & isfinite(lo) & isfinite(hi);
if ~any(mask), return; end
idx = find(mask);
starts = [idx(1), idx(find(diff(idx) > 1) + 1)];
ends   = [idx(find(diff(idx) > 1)), idx(end)];
for s = 1:numel(starts)
    ii = starts(s):ends(s);
    xx = t(ii); yy1 = lo(ii); yy2 = hi(ii);
    fill(ax, [xx, fliplr(xx)], [yy1, fliplr(yy2)], col, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility','off');
end
end
