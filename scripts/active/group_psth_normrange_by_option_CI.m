function group_psth_normrange_by_option_CI(C, sessionIdxs, neuronIdxs, bin_ms, smooth_ms, win, baseline_win, error_type, ci_alpha)
% Baseline-subtract per neuron (pooled across trials), then global min–max
% normalize via normalize(vec,'range') across all trials x time per neuron.
% Finally, average per condition and across neurons, and plot the three
% conditions overlaid (Buy/Hold/Sell), with optional SEM or CI shading.
%
% Example:
%   sessions = [1 1 5]; neurons = [2 4 5];
%   group_psth_normrange_by_option_CI(C, sessions, neurons, 20, 100);                 % 95% CI (default)
%   group_psth_normrange_by_option_CI(C, sessions, neurons, 20, 100, [], [], 'sem');  % SEM shading
%   group_psth_normrange_by_option_CI(C, sessions, neurons, 20, 100, [], [], 'ci', 0.90); % 90% CI


% ---- Defaults ----
if nargin < 4 || isempty(bin_ms),      bin_ms = 5;                  end
if nargin < 5 || isempty(smooth_ms),   smooth_ms = 200;             end
if nargin < 6 || isempty(win),         win = [-1000, 1500];         end % ms
if nargin < 7 || isempty(baseline_win), baseline_win = [-1000, -700]; end
if nargin < 8 || isempty(error_type),  error_type = 'ci';           end % 'none'|'sem'|'ci'
if nargin < 9 || isempty(ci_alpha),    ci_alpha = 0.95;             end

% ---- Binning setup ----
edges = win(1):bin_ms:win(2);
bin_centers = edges(1:end-1) + bin_ms/2;
bin_sec = bin_ms / 1000;
nBins = numel(edges) - 1;

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

% (optional) track total trials per condition across neurons (for legend)
tot_buy = 0; tot_hold = 0; tot_sell = 0;

% ---- Baseline mask ----
base_mask = (bin_centers >= baseline_win(1)) & (bin_centers < baseline_win(2));
if ~any(base_mask)
    error('Baseline mask empty. Check baseline_win vs win and bin_ms.');
end

for k = 1:numel(neuronIdxs)
    sess = sessionIdxs(k);
    neu  = neuronIdxs(k);

    % ----- Spikes (sec -> ms) -----
    rawSess = C.rawSpikes{sess,1};
    if iscell(rawSess), spikes_sec = rawSess{neu};
    else,               spikes_sec = rawSess(:,neu);
    end
    spikes_ms = double(spikes_sec) * 1000;

    % ----- Reward-on times (ms) -----
    evtCell = E{sess,1}.m1rew1on{1,1};
    if istable(evtCell) || isa(evtCell,'table'), evt_ms = table2array(evtCell);
    else,                                        evt_ms = evtCell;
    end
    evt_ms_vec = double(evt_ms(:));

    % ----- Choice codes (1=Buy,2=Hold,3=Sell) -----
    optCell = E{sess,1}.option{1,1};
    if istable(optCell) || isa(optCell,'table'), opt = table2array(optCell);
    else,                                        opt = optCell;
    end
    opt_vec = double(opt(:));

    % ----- Align & clean -----
    n = min(numel(evt_ms_vec), numel(opt_vec));
    evt_ms_vec = evt_ms_vec(1:n);
    opt_vec    = opt_vec(1:n);
    valid = ~isnan(evt_ms_vec) & ~isnan(opt_vec);
    evt_ms_vec = evt_ms_vec(valid);
    opt_vec    = opt_vec(valid);
    nTrials = numel(evt_ms_vec);
    if nTrials == 0
        warning('No valid trials for session %d neuron %d; skipping.', sess, neu);
        continue
    end

    % ----- Bin spikes per trial -> FR (Hz) -----
    counts = zeros(nTrials, nBins);
    for t = 1:nTrials
        rel_ms  = spikes_ms - evt_ms_vec(t);
        rel_win = rel_ms(rel_ms >= win(1) & rel_ms < win(2));
        counts(t,:) = histcounts(rel_win, edges);
    end
    fr_trials = counts / bin_sec;           % [trials x bins], Hz

    % ----- Pooled baseline subtraction (per neuron) -----
    mu_base_neuron = mean(fr_trials(:, base_mask), 'all', 'omitnan');  % scalar Hz
    fr_bs = fr_trials - mu_base_neuron;                                 % baseline-subtracted FR

    % ----- Global min–max normalization across ALL trials x bins (per neuron) -----
    vec = fr_bs(:);
    if exist('normalize','file') == 2
        vec = normalize(vec, 'scale');      % min->0, max->1 over ALL elements
    else
        mn = min(vec); mx = max(vec); denom = mx - mn;
        if ~isfinite(denom) || denom <= 0, denom = eps; end
        vec = (vec - mn) / denom;
    end
    norm_trials = reshape(vec, size(fr_bs));  % back to [trials x bins]

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

% Y limits from all lines and bands (allow <0 and >1)
y_all = [N_buy_s, N_hold_s, N_sell_s, ...
         buy_lo_s, buy_hi_s, hold_lo_s, hold_hi_s, sell_lo_s, sell_hi_s];
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
xlabel(sprintf('Time from m1rew1on (ms) - smoothing approx %d ms', smooth_ms));
ylabel('Mean normalized activity (0-1), baseline-subtracted');
title('Mean normalized PSTH by Choice (overlay across selected neurons)');
if ~isempty(h)
    legend(h, {sprintf('Buy (trials=%d)',tot_buy), ...
               sprintf('Hold (trials=%d)',tot_hold), ...
               sprintf('Sell (trials=%d)',tot_sell)}, 'Location','best');
end
end

% ===== Helpers =====
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
