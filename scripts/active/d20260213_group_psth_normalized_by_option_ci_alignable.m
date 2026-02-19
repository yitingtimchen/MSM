% File: d20260213_group_psth_normalized_by_option_ci_alignable.m
function d20260213_group_psth_normalized_by_option_ci_alignable(C, sessionIdxs, neuronIdxs, align_event_name, bin_ms, smooth_ms, win, baseline_win, error_type, ci_alpha, baseline_event_name)
% Align PSTHs to an arbitrary event (align_event_name), while computing the
% baseline from baseline_win relative to baseline_event_name (default: m1rew1on).
% Then min–max normalize per neuron, average by option (Buy/Hold/Sell),
% smooth, and plot with SEM/CI shading.
%
% Example:
%   sessions = [1 1 5]; neurons = [2 4 5];
%   d20260213_group_psth_normalized_by_option_ci_alignable(C, sessions, neurons, 'm1rew2on', 20, 100);
%   d20260213_group_psth_normalized_by_option_ci_alignable(C, sessions, neurons, 'offerOn', 20, 100, [-800 1200]);
%   % If your baseline event has a different name:
%   d20260213_group_psth_normalized_by_option_ci_alignable(C, sessions, neurons, 'm1rew2on', 20, 100, [], [], 'ci', 0.95, 'm1rew1on');

% ---- Defaults ----
if nargin < 4 || isempty(align_event_name), align_event_name = 'm1rew2on'; end
if nargin < 5 || isempty(bin_ms),           bin_ms = 5;                  end
if nargin < 6 || isempty(smooth_ms),        smooth_ms = 200;             end
if nargin < 7 || isempty(win),              win = [-1000, 1500];         end % ms, relative to align_event
if nargin < 8 || isempty(baseline_win),     baseline_win = [-1000, -700]; end % ms, relative to baseline_event_name
if nargin < 9 || isempty(error_type),       error_type = 'ci';           end % 'none'|'sem'|'ci'
if nargin < 10 || isempty(ci_alpha),        ci_alpha = 0.95;             end
if nargin < 11 || isempty(baseline_event_name), baseline_event_name = 'm1rew1on'; end

% ---- Binning setup (main, aligned to chosen event) ----
edges       = win(1):bin_ms:win(2);
bin_centers = edges(1:end-1) + bin_ms/2;
bin_sec     = bin_ms / 1000;
nBins       = numel(edges) - 1;

% ---- Baseline binning (relative to baseline_event_name) ----
base_edges = baseline_win(1):bin_ms:baseline_win(2);
if numel(base_edges) < 2
    error('Baseline window too narrow for given bin size.');
end

% ---- Locate event table ----
if isfield(C, 'eventTables')
    E = C.eventTables;
elseif isfield(C, 'evenTables')
    E = C.evenTables;
else
    error('Neither C.eventTables nor C.evenTables found.');
end

% ---- Storage for per-neuron condition means (normalized) ----
N_buy_all  = nan(numel(neuronIdxs), nBins);
N_hold_all = nan(numel(neuronIdxs), nBins);
N_sell_all = nan(numel(neuronIdxs), nBins);
tot_buy = 0; tot_hold = 0; tot_sell = 0;

for k = 1:numel(neuronIdxs)
    sess = sessionIdxs(k);
    neu  = neuronIdxs(k);

    % ----- Spikes (sec -> ms) -----
    rawSess = C.rawSpikes{sess,1};
    if iscell(rawSess), spikes_sec = rawSess{neu};
    else,               spikes_sec = rawSess(:,neu);
    end
    spikes_ms = double(spikes_sec) * 1000;

    % ----- Pull events robustly from struct OR table -----
    E_sess = E{sess,1};

    [baseline_ms, ok_base] = get_event_vec(E_sess, baseline_event_name);
    if ~ok_base
        warning('Session %d missing baseline event "%s"; skipping neuron %d.', sess, baseline_event_name, neu);
        continue
    end

    [align_ms, ok_align] = get_event_vec(E_sess, align_event_name);
    if ~ok_align
        warning('Session %d missing alignment event "%s"; skipping neuron %d.', sess, align_event_name, neu);
        continue
    end

    [opt_vec, ok_opt] = get_event_vec(E_sess, 'option');
    if ~ok_opt
        warning('Session %d missing option codes; skipping neuron %d.', sess, neu);
        continue
    end

    % ----- Align & clean (require baseline_event, align_event, and option) -----
    n = min([numel(baseline_ms), numel(align_ms), numel(opt_vec)]);
    baseline_ms = baseline_ms(1:n);
    align_ms    = align_ms(1:n);
    opt_vec     = opt_vec(1:n);

    valid = isfinite(baseline_ms) & isfinite(align_ms) & isfinite(opt_vec);
    baseline_ms = baseline_ms(valid);
    align_ms    = align_ms(valid);
    opt_vec     = opt_vec(valid);

    nTrials = numel(baseline_ms);
    if nTrials == 0
        warning('No valid trials for session %d neuron %d; skipping.', sess, neu);
        continue
    end

    % ----- Main PSTH: align to chosen event -----
    counts_main = zeros(nTrials, nBins);
    for t = 1:nTrials
        rel_ms  = spikes_ms - align_ms(t);
        rel_win = rel_ms(rel_ms >= win(1) & rel_ms < win(2));
        counts_main(t,:) = histcounts(rel_win, edges);
    end
    fr_main = counts_main / bin_sec;  % Hz

    % ----- Baseline: pooled mean over baseline_win relative to baseline_event_name -----
    base_counts = zeros(nTrials, numel(base_edges)-1);
    for t = 1:nTrials
        rel_ms_b  = spikes_ms - baseline_ms(t);
        rel_win_b = rel_ms_b(rel_ms_b >= baseline_win(1) & rel_ms_b < baseline_win(2));
        base_counts(t,:) = histcounts(rel_win_b, base_edges);
    end
    mu_base_neuron = mean(base_counts, 'all') / bin_sec;   % scalar Hz

    % ----- Baseline subtraction on aligned FR -----
    fr_bs = fr_main - mu_base_neuron;

    % ----- Global min–max normalization across ALL trials x bins (per neuron) -----
    vec = fr_bs(:);
    if exist('normalize','file') == 2
        vec = normalize(vec, 'scale');  % 0..1
    else
        mn = min(vec); mx = max(vec); denom = mx - mn;
        if ~isfinite(denom) || denom <= 0, denom = eps; end
        vec = (vec - mn) / denom;
    end
    norm_trials = reshape(vec, size(fr_bs));

    % ----- Per-condition mean (per neuron) -----
    idx_buy  = (opt_vec == 1);
    idx_hold = (opt_vec == 2);
    idx_sell = (opt_vec == 3);

    if any(idx_buy),  N_buy_all(k,:)  = mean(norm_trials(idx_buy ,:), 1, 'omitnan');  tot_buy  = tot_buy  + sum(idx_buy);  end
    if any(idx_hold), N_hold_all(k,:) = mean(norm_trials(idx_hold,:), 1, 'omitnan');  tot_hold = tot_hold + sum(idx_hold); end
    if any(idx_sell), N_sell_all(k,:) = mean(norm_trials(idx_sell,:), 1, 'omitnan');  tot_sell = tot_sell + sum(idx_sell); end
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

% ----- Plot -----
figure('Color','w'); hold on
h = [];

% shaded error (if requested)
if ~strcmpi(error_type,'none')
    shaded(bin_centers, buy_lo_s , buy_hi_s , [0 0.6 0]);   % green
    shaded(bin_centers, hold_lo_s, hold_hi_s, [0 0.4 1]);   % blue
    shaded(bin_centers, sell_lo_s, sell_hi_s, [0.85 0 0]);  % red
end

% mean lines
if any(isfinite(N_buy_s)),  h(end+1) = plot(bin_centers, N_buy_s , 'Color', [0 0.6 0], 'LineWidth', 3.5); end
if any(isfinite(N_hold_s)), h(end+1) = plot(bin_centers, N_hold_s, 'Color', [0 0 0.8], 'LineWidth', 3.5); end
if any(isfinite(N_sell_s)), h(end+1) = plot(bin_centers, N_sell_s, 'r', 'LineWidth', 3.5); end

% Y limits from all lines and bands
y_all = [N_buy_s, N_hold_s, N_sell_s, buy_lo_s, buy_hi_s, hold_lo_s, hold_hi_s, sell_lo_s, sell_hi_s];
y_all = y_all(isfinite(y_all));
if isempty(y_all), ymin = 0; ymax = 1;
else,              ymin = min(y_all); ymax = max(y_all);
end
pad  = 0.05 * max(1e-6, ymax - ymin);
ylim([ymin - pad, ymax + pad]);

% X limits, zero line, labels
xlim(win);
if exist('xline','file') == 2
    xline(0,'k--','LineWidth',1);
else
    yl = ylim; plot([0 0], yl, 'k--', 'LineWidth', 1);
end
grid on
xlabel(sprintf('Time from %s (ms) - smoothing approx %d ms', align_event_name, smooth_ms));
ylabel('Mean normalized activity (0-1), baseline-subtracted');
title(sprintf('Mean normalized PSTH by Choice (aligned to %s; baseline from %s)', align_event_name, baseline_event_name));
if ~isempty(h)
    legend(h, {sprintf('Buy (trials=%d)',tot_buy), ...
               sprintf('Hold (trials=%d)',tot_hold), ...
               sprintf('Sell (trials=%d)',tot_sell)}, 'Location','best');
end
end

% ======== Helpers ========
function [vec, ok] = get_event_vec(ES, name)
% Extract event vector (ms) from ES (table or struct) robustly.
ok  = false; 
vec = [];

try
    % 1) Pick the container and pull the column/field
    if istable(ES)
        if any(strcmp(name, ES.Properties.VariableNames))
            x = ES.(name);
        else
            return
        end
    elseif isstruct(ES)
        if isfield(ES, name)
            x = ES.(name);
        else
            return
        end
    else
        return
    end

    % 2) Repeatedly unwrap cells
    while iscell(x)
        if isempty(x)
            return
        end
        x = x{1};
    end

    % 3) If it's a table, extract a usable column
    if istable(x)
        w = size(x, 2);
        got = false;
        for i = 1:w
            col = x{:, i};  % extract raw column
            % unwrap cells again if needed
            while iscell(col) && ~isempty(col)
                col = col{1};
            end
            if isnumeric(col) || islogical(col) || iscategorical(col) || isstring(col) || ischar(col)
                x = col;
                got = true;
                break
            end
        end
        if ~got
            return
        end
    end

    % 4) Final unwrapping in case we got a cell from the table column
    while iscell(x)
        if isempty(x)
            return
        end
        x = x{1};
    end

    % 5) Type coercion to numeric vector
    if iscategorical(x)
        x = double(x);
    elseif isstring(x) || ischar(x)
        x = str2double(string(x));
    elseif islogical(x)
        x = double(x);
    elseif ~isnumeric(x)
        return
    end

    vec = double(x(:));
    ok  = ~isempty(vec) && any(isfinite(vec));
catch
    ok = false; 
    vec = [];
end
end


function [lo, hi] = compute_bands(M, error_type, ci_alpha)
% M is nNeurons x nBins of per-neuron means (already in [0,1])
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
        lo = mu; hi = mu;
        for j = 1:numel(mu)
            if n(j) >= 2 && isfinite(sem(j))
                tcrit = tinv(1 - alpha2, n(j) - 1);
                lo(j) = mu(j) - tcrit * sem(j);
                hi(j) = mu(j) + tcrit * sem(j);
            else
                lo(j) = NaN; hi(j) = NaN;
            end
        end
    otherwise
        lo = nan(size(mu)); hi = nan(size(mu));
end
end

function shaded(t, lo, hi, col)
mask = isfinite(t) & isfinite(lo) & isfinite(hi);
if ~any(mask), return; end
idx = find(mask);
starts = [idx(1), idx(find(diff(idx) > 1) + 1)];
ends   = [idx(find(diff(idx) > 1)), idx(end)];
for s = 1:numel(starts)
    ii = starts(s):ends(s);
    xx = t(ii); yy1 = lo(ii); yy2 = hi(ii);
    fill([xx, fliplr(xx)], [yy1, fliplr(yy2)], col, ...
        'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility','off');
end
end
