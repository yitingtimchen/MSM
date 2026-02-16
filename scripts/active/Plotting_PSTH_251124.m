
%% Batch runner: variable-interval trimmed PSTHs (arrange figures 2x4: scPrefixes x condNames)
clc; close all; clear;
warning('off');

cd('C:\Users\plattlab\MSM\outputs_local\condition_packs_251118')

% ---------- Configuration ----------
align_sets = {
    ["m1trialStart","m1rew1on","m1rew1off","m1rew2on","m1rew2off",...
     "m2trialStart","m2rew1on","m2rew1off","m2rew2on","m2rew2off"]
    };

win_sets = {
    [ 
      -1000 0;     % m1trialStart
      -500 500;    % m1rew1on
      -500 500;    % m1rew1off
      -500 1000;   % m1rew2on
      -700 1000;   % m1rew2off/m1trialStop
      -700 100;    % m2trialStart
      -500 500;    % m2rew1on
      -500 500;    % m2rew1off
      -500 1000;   % m2rew2on
      -500 1500]   % m2rew2off/m2trialStop
    };

types_to_use = [1 2];
scPrefixes   = {'OT ', 'Saline '};          % rows (2)
condNames    = {'AI','Decoy','Replay','Live'};  % cols (4)
bin_ms       = 10;
smooth_ms    = 50;
baseline_win = [-1500 -1000];
error_type   = 'ci';
ci_alpha     = 0.90;
baseline_event_name = "m1rew1on";
combine_not_buy     = 0;

% ---------- Figure tiling config ----------
nRows = numel(scPrefixes);
nCols = numel(condNames);
scr   = get(0,'ScreenSize');                % [left bottom width height]
figW  = scr(3)/nCols;
figH  = scr(4)/nRows;

% ---------- Batch loop ----------
for p = 1:numel(scPrefixes)
    prefix = scPrefixes{p};
    for c = 1:numel(condNames)
        cond = condNames{c};

        fprintf("\nProcessing " + string(cond) + " condition...\n")
        fname = sprintf('%s%s_condition_pack.mat', prefix, cond);
        if ~isfile(fname), warning('Missing file: %s', fname); continue; end

        S = load(fname);
        if ~isfield(S,'C'), warning('No struct C in %s', fname); continue; end
        C = S.C;
        C.cond = upper(cond);

        C = filterNeuronAndAddType(C, prefix);
        [sessions, neurons] = collectSessionsNeurons(C.neuronType, types_to_use);
        if isempty(sessions)
            warning('No neurons of requested types for %s%s.', prefix, cond);
            continue
        end

        for a = 1:numel(align_sets)
            align_events = align_sets{a};
            plot_wins    = win_sets{a};

            group_psth_trimmed_overlay_and_mean_perEventWin( ...
                C, sessions, neurons, ...
                align_events, plot_wins, ...
                bin_ms, smooth_ms, ...
                baseline_win, error_type, ci_alpha, ...
                baseline_event_name, ...
                sprintf('%s %s', strtrim(prefix), cond), ...
                combine_not_buy ...
                );

            % ---- Tile this figure into 2 x 4 grid (scPrefixes x condNames) ----
            hFig = gcf;
            set(hFig,'Units','pixels');
            left   = scr(1) + (c-1)*figW;
            bottom = scr(2) + (nRows-p)*figH;   % p=1 on top row
            set(hFig,'OuterPosition',[left bottom figW figH]);
        end
    end
end

% ---------- Optional utilities after plots ----------
% % Equalize Y-limits across all open figures (all axes except colorbars)
% % ax = findall(groot,'Type','axes','-not','Tag','Colorbar');
% % if ~isempty(ax)
% %     ylc = get(ax,'YLim');
% %     if iscell(ylc), ymat = vertcat(ylc{:}); else, ymat = ylc; end
% %     yl = [min(ymat(:,1)) max(ymat(:,2))];
% %     set(ax, 'YLim', yl);
% % end
% % Link Y for synchronized zoom/pan
% % if ~isempty(ax), linkaxes(ax,'y'); end
% % Stack figures vertically on primary monitor (tweak margins/gap as desired)
% % stackFiguresVertical(1, [300 100 100 100], 20);

%% ============================= Helpers & Core =============================

function [sessions, neurons] = collectSessionsNeurons(neuronType, types_to_use)
% Convert neuronType cell array into parallel (sessions, neurons) vectors for selected types.
sessions = []; neurons = [];
for s = 1:numel(neuronType)
    lbl = neuronType{s};
    if isempty(lbl), continue; end
    keep = find(ismember(lbl, uint8(types_to_use)));
    if ~isempty(keep)
        sessions = [sessions, repmat(s, 1, numel(keep))];
        neurons  = [neurons, reshape(keep,1,[])];
    end
end
end

function group_psth_trimmed_overlay_and_mean_perEventWin(C, sessionIdxs, neuronIdxs, align_event_names, plot_wins, bin_ms, smooth_ms, baseline_win, error_type, ci_alpha, baseline_event_name, figTitle, combine_not_buy)
% Version without trial overlays: only grouped PSTHs (one row of subplots).
%
% Tim Chen / Platt Lab — per-event windows + proportional subplot widths + negative-duration exclusion.

% ---- Normalize inputs ----
if ischar(align_event_names),   align_event_names = string(align_event_names); end
if iscellstr(align_event_names),align_event_names = string(align_event_names); end
if isstring(align_event_names) && isscalar(align_event_names), align_event_names = align_event_names(:)'; end
nAlign = numel(align_event_names);
if nargin < 13 || isempty(combine_not_buy), combine_not_buy = false; end

% ---- Per-event plotting windows -> nAlign x 2 matrix ----
wins_mat = parse_plot_wins(plot_wins, nAlign);   % [tmin tmax] per event

% ---- Per-event bin edges and centers ----
edges_list       = cell(1,nAlign);
bin_centers_list = cell(1,nAlign);
nBins_list       = zeros(1,nAlign);
for a = 1:nAlign
    tmin = wins_mat(a,1); tmax = wins_mat(a,2);
    edges_list{a}       = tmin:bin_ms:tmax;
    bin_centers_list{a} = edges_list{a}(1:end-1) + bin_ms/2;
    nBins_list(a)       = numel(edges_list{a}) - 1;
end
bin_sec     = bin_ms / 1000;
smooth_bins = max(3, round(smooth_ms / bin_ms));

% ---- Baseline binning (relative to baseline_event_name) ----
base_edges  = baseline_win(1):bin_ms:baseline_win(2);
if numel(base_edges) < 2, error('Baseline window too narrow for given bin size.'); end

% ---- Locate event table ----
if isfield(C, 'eventTables'), E = C.eventTables;
elseif isfield(C, 'evenTables'), E = C.evenTables;
else, error('Neither C.eventTables nor C.evenTables found.');
end

% =================== Interval summary with invalid flagging ===================
dt_pool = cell(1, max(nAlign-1,0));
dt_neg_counts = zeros(1, max(nAlign-1,0));
dt_tot_counts = zeros(1, max(nAlign-1,0));
sess_unique = unique(sessionIdxs(:))';
for ss = sess_unique
    E_sess = E{ss,1};
    ev = cell(1,nAlign); okv = false(1,nAlign);
    for a = 1:nAlign
        [ev{a}, okv(a)] = get_event_vec(E_sess, align_event_names(a));
    end
    if ~all(okv), continue; end
    m = min(cellfun(@numel, ev));
    for a = 1:nAlign-1
        d = double(ev{a+1}(1:m) - ev{a}(1:m));   % ms
        good_mask = isfinite(d) & (d > 0);
        bad_mask  = isfinite(d) & ~good_mask;    % negative or zero
        dt_tot_counts(a) = dt_tot_counts(a) + sum(isfinite(d));
        dt_neg_counts(a) = dt_neg_counts(a) + sum(bad_mask);
        if any(good_mask)
            dt_pool{a} = [dt_pool{a}; d(good_mask)]; %#ok<AGROW>
        end
    end
end
dt_mean_ms = cellfun(@(x) mean(x,'omitnan'), dt_pool);
dt_std_ms  = cellfun(@(x) std(x,0,'omitnan'), dt_pool);
dt_n       = cellfun(@(x) numel(x), dt_pool);
pair_labels = strings(1, max(nAlign-1,0));
for a = 1:numel(pair_labels)
    pair_labels(a) = sprintf('%s→%s', align_event_names(a), align_event_names(a+1));
end
if ~isempty(pair_labels)
    fprintf('Consecutive-event intervals (negative/zero durations excluded):\n');
    for a = 1:numel(pair_labels)
        fprintf('  %s: mean=%.0f ms, std=%.0f, n=%d | excluded=%d of %d\n', ...
            pair_labels(a), dt_mean_ms(a), dt_std_ms(a), dt_n(a), dt_neg_counts(a), dt_tot_counts(a));
    end
end

% ---- Group containers ----
buy_all      = cell(1, nAlign);
hold_all     = cell(1, nAlign);    % only when ~combine_not_buy
sell_all     = cell(1, nAlign);    % only when ~combine_not_buy
notbuy_all   = cell(1, nAlign);    % only when combine_not_buy

tot_buy     = zeros(1, nAlign);
tot_hold    = zeros(1, nAlign);
tot_sell    = zeros(1, nAlign);
tot_notbuy  = zeros(1, nAlign);

% Keep counters for how many trials get dropped per panel due to invalid intervals
drop_counts = zeros(1, nAlign);

% ========================= Main accumulation loop =========================
for k = 1:numel(neuronIdxs)
    sess = sessionIdxs(k);
    neu  = neuronIdxs(k);

    % -- Spikes (sec -> ms) --
    rawSess = C.rawSpikes{sess,1};
    if iscell(rawSess), spikes_sec = rawSess{neu};
    else,               spikes_sec = rawSess(:,neu);
    end
    spikes_ms = double(spikes_sec) * 1000;

    % -- Events for this session --
    E_sess = E{sess,1};
    [baseline_ms, ok_base] = get_event_vec(E_sess, baseline_event_name);
    if ~ok_base, warning('Session %d missing baseline event "%s"; skipping neuron %d.', sess, baseline_event_name, neu); continue; end

    % Collect all requested event time vectors once (for prev/next trimming)
    align_times = cell(1, nAlign); ok_vec = false(1, nAlign);
    for a = 1:nAlign
        [v, ok] = get_event_vec(E_sess, align_event_names(a));
        align_times{a} = v; ok_vec(a) = ok;
    end
    if ~all(ok_vec)
        warning('Session %d missing some alignment events; skipping neuron %d.', sess, neu);
        continue
    end

    % Option codes for this session
    [opt_vec_full, ok_opt_all] = get_event_vec(E_sess, 'option');
    if ~ok_opt_all, warning('Session %d missing option codes; skipping neuron %d.', sess, neu); continue; end

    % -- For each event: build trimmed trial matrix and aggregate --
    for a = 1:nAlign
        align_ms    = align_times{a};
        edges       = edges_list{a};
        bin_centers = bin_centers_list{a};
        nBins       = nBins_list(a);
        tmin        = wins_mat(a,1);
        tmax        = wins_mat(a,2);

        % Neighbor events (if exist)
        if a > 1,      prev_ms = align_times{a-1}; else, prev_ms = []; end
        if a < nAlign, next_ms = align_times{a+1}; else, next_ms = []; end

        % Truncate to common trial count
        n = min([numel(baseline_ms), numel(align_ms), numel(opt_vec_full)]);
        if a > 1,      n = min(n, numel(prev_ms)); end
        if a < nAlign, n = min(n, numel(next_ms)); end
        if n == 0, continue; end

        baseline_ms_n = baseline_ms(1:n);
        align_ms_n    = align_ms(1:n);
        opt_vec_n     = opt_vec_full(1:n);
        if ~isempty(prev_ms), prev_ms_n = prev_ms(1:n); else, prev_ms_n = []; end
        if ~isempty(next_ms), next_ms_n = next_ms(1:n); else, next_ms_n = []; end

        % Keep only valid trials (finite timestamps)
        valid = isfinite(baseline_ms_n) & isfinite(align_ms_n) & isfinite(opt_vec_n);
        if ~isempty(prev_ms_n), valid = valid & isfinite(prev_ms_n); end
        if ~isempty(next_ms_n), valid = valid & isfinite(next_ms_n); end

        % ---- EXCLUDE negative/zero consecutive durations for this panel ----
        if ~isempty(prev_ms_n)
            valid = valid & (prev_ms_n < align_ms_n);
        end
        if ~isempty(next_ms_n)
            valid = valid & (align_ms_n < next_ms_n);
        end

        dropped_here = sum(~valid);
        drop_counts(a) = drop_counts(a) + dropped_here;

        baseline_ms_n = baseline_ms_n(valid);
        align_ms_n    = align_ms_n(valid);
        opt_vec_n     = opt_vec_n(valid);
        if ~isempty(prev_ms_n), prev_ms_n = prev_ms_n(valid); end
        if ~isempty(next_ms_n), next_ms_n = next_ms_n(valid); end
        nTrials = numel(align_ms_n);
        if nTrials == 0, continue; end

        % -- Histogram counts aligned to chosen event, then trim by prev/next interval --
        counts_main = nan(nTrials, nBins);
        for t = 1:nTrials
            rel_ms  = spikes_ms - align_ms_n(t);
            in_win  = rel_ms(rel_ms >= tmin & rel_ms < tmax);
            c = histcounts(in_win, edges);

            LB = tmin;  UB = tmax;
            if ~isempty(prev_ms_n)
                gap_prev = align_ms_n(t) - prev_ms_n(t);
                if isfinite(gap_prev), LB = max(LB, -gap_prev); end
            end
            if ~isempty(next_ms_n)
                gap_next = next_ms_n(t) - align_ms_n(t);
                if isfinite(gap_next), UB = min(UB,  gap_next); end
            end

            valid_bins = (bin_centers > LB) & (bin_centers < UB);

            row = nan(1, nBins);
            row(valid_bins) = c(valid_bins);
            counts_main(t,:) = row;
        end
        fr_main = counts_main / bin_sec;  % Hz; NaNs preserved

        % -- Baseline (relative to baseline_event_name): mean + std from BASELINE ONLY --
        base_counts = zeros(nTrials, numel(base_edges)-1);
        for t = 1:nTrials
            rel_ms_b  = spikes_ms - baseline_ms_n(t);
            rel_win_b = rel_ms_b(rel_ms_b >= baseline_win(1) & rel_ms_b < baseline_win(2));
            base_counts(t,:) = histcounts(rel_win_b, base_edges);
        end
        base_fr = base_counts / bin_sec;
        mu_base_neuron  = mean(base_fr(:), 'omitnan');
        sig_base_neuron = std(base_fr(:), 0, 'omitnan');
        if ~isfinite(sig_base_neuron) || sig_base_neuron <= 0
            sig_base_neuron = eps;
        end

        fr_z = (fr_main - mu_base_neuron) / sig_base_neuron;
        if ~any(isfinite(fr_z(:)))
            continue;
        end
        norm_trials = fr_z;
        norm_trials_s = nan_smooth_rows(norm_trials, smooth_bins);

        % -- Per-neuron means by group --
        if combine_not_buy
            idx_buy    = (opt_vec_n == 1);
            idx_notbuy = (opt_vec_n == 2) | (opt_vec_n == 3);
            if any(idx_buy)
                b = mean(norm_trials_s(idx_buy ,:), 1, 'omitnan');
                buy_all{a} = [buy_all{a}; b]; %#ok<AGROW>
                tot_buy(a) = tot_buy(a) + sum(idx_buy);
            end
            if any(idx_notbuy)
                n = mean(norm_trials_s(idx_notbuy,:), 1, 'omitnan');
                notbuy_all{a} = [notbuy_all{a}; n]; %#ok<AGROW>
                tot_notbuy(a) = tot_notbuy(a) + sum(idx_notbuy);
            end
        else
            idx_buy  = (opt_vec_n == 1);
            idx_hold = (opt_vec_n == 2);
            idx_sell = (opt_vec_n == 3);
            if any(idx_buy)
                b = mean(norm_trials_s(idx_buy ,:), 1, 'omitnan');
                buy_all{a} = [buy_all{a}; b]; %#ok<AGROW>
                tot_buy(a) = tot_buy(a) + sum(idx_buy);
            end
            if any(idx_hold)
                h = mean(norm_trials_s(idx_hold,:), 1, 'omitnan');
                hold_all{a} = [hold_all{a}; h]; %#ok<AGROW>
                tot_hold(a) = tot_hold(a) + sum(idx_hold);
            end
            if any(idx_sell)
                s = mean(norm_trials_s(idx_sell,:), 1, 'omitnan');
                sell_all{a} = [sell_all{a}; s]; %#ok<AGROW>
                tot_sell(a) = tot_sell(a) + sum(idx_sell);
            end
        end
    end
end

% Report how many trials were dropped per panel due to invalid intervals
if any(drop_counts>0)
    fprintf('Dropped invalid trials per panel (non-positive prev/next interval):\n');
    for a = 1:nAlign
        fprintf('  %s: %d dropped\n', string(align_event_names(a)), drop_counts(a));
    end
end

% ========================= Across-neuron summaries ========================
if combine_not_buy
    N_buy_mean     = cell(1,nAlign);
    N_notbuy_mean  = cell(1,nAlign);
    buy_lo = cell(1,nAlign); buy_hi = cell(1,nAlign);
    nbuy_lo = cell(1,nAlign); nbuy_hi = cell(1,nAlign);
else
    N_buy_mean  = cell(1,nAlign); N_hold_mean = cell(1,nAlign); N_sell_mean = cell(1,nAlign);
    buy_lo = cell(1,nAlign); buy_hi = cell(1,nAlign);
    hold_lo = cell(1,nAlign); hold_hi = cell(1,nAlign);
    sell_lo = cell(1,nAlign); sell_hi = cell(1,nAlign);
end

for a = 1:nAlign
    nBins = nBins_list(a);

    if combine_not_buy
        if isempty(buy_all{a}),    buy_all{a}    = nan(0, nBins); end
        if isempty(notbuy_all{a}), notbuy_all{a} = nan(0, nBins); end

        N_buy_mean{a}    = mean(buy_all{a}   , 1, 'omitnan');
        N_notbuy_mean{a} = mean(notbuy_all{a}, 1, 'omitnan');

        [buy_lo{a},  buy_hi{a} ]  = compute_bands(buy_all{a}   , error_type, ci_alpha);
        [nbuy_lo{a}, nbuy_hi{a}]  = compute_bands(notbuy_all{a}, error_type, ci_alpha);

        N_buy_mean{a}    = nan_smooth_row(N_buy_mean{a},     smooth_bins);
        N_notbuy_mean{a} = nan_smooth_row(N_notbuy_mean{a},  smooth_bins);
        buy_lo{a}  = nan_smooth_row(buy_lo{a},   smooth_bins);  buy_hi{a}  = nan_smooth_row(buy_hi{a},   smooth_bins);
        nbuy_lo{a} = nan_smooth_row(nbuy_lo{a},  smooth_bins);  nbuy_hi{a} = nan_smooth_row(nbuy_hi{a},  smooth_bins);
    else
        if isempty(buy_all{a}),  buy_all{a}  = nan(0, nBins); end
        if isempty(hold_all{a}), hold_all{a} = nan(0, nBins); end
        if isempty(sell_all{a}), sell_all{a} = nan(0, nBins); end

        N_buy_mean{a}  = mean(buy_all{a} , 1, 'omitnan');
        N_hold_mean{a} = mean(hold_all{a}, 1, 'omitnan');
        N_sell_mean{a} = mean(sell_all{a}, 1, 'omitnan');

        [buy_lo{a},  buy_hi{a} ] = compute_bands(buy_all{a} , error_type, ci_alpha);
        [hold_lo{a}, hold_hi{a}] = compute_bands(hold_all{a}, error_type, ci_alpha);
        [sell_lo{a}, sell_hi{a}] = compute_bands(sell_all{a}, error_type, ci_alpha);

        N_buy_mean{a}  = nan_smooth_row(N_buy_mean{a},  smooth_bins);
        N_hold_mean{a} = nan_smooth_row(N_hold_mean{a}, smooth_bins);
        N_sell_mean{a} = nan_smooth_row(N_sell_mean{a}, smooth_bins);
        buy_lo{a}  = nan_smooth_row(buy_lo{a},  smooth_bins);  buy_hi{a}  = nan_smooth_row(buy_hi{a},  smooth_bins);
        hold_lo{a} = nan_smooth_row(hold_lo{a}, smooth_bins);  hold_hi{a} = nan_smooth_row(hold_hi{a}, smooth_bins);
        sell_lo{a} = nan_smooth_row(sell_lo{a}, smooth_bins);  sell_hi{a} = nan_smooth_row(sell_hi{a}, smooth_bins);
    end
end

% ====================== Proportional subplot layout =======================
W    = wins_mat(:,2) - wins_mat(:,1);
frac = (W / sum(W)).';   % 1 x nAlign

figure('Color','w','Name',figTitle);

L=0.08; R=0.02; B=0.12; T=0.2; colGap=0.02;
totW = 1 - L - R - colGap*(nAlign-1);
totH = 1 - T - B;

x = L;  ax = gobjects(1,nAlign);
for a = 1:nAlign
    w = frac(a) * totW;
    ax(a) = axes('Position',[x, B, w, totH]);  hold(ax(a),'on');
    x = x + w + colGap;
end

% ============================== Plotting =================================
for a = 1:nAlign
    tmin        = wins_mat(a,1);
    tmax        = wins_mat(a,2);
    bin_centers = bin_centers_list{a};

    hleg = [];
    if combine_not_buy
        if exist('nbuy_lo','var') && ~strcmpi(error_type,'none')
            shaded(ax(a), bin_centers, buy_lo{a},  buy_hi{a},  [0 0.6 0]);   % Buy band
            shaded(ax(a), bin_centers, nbuy_lo{a}, nbuy_hi{a}, [0 0.4 1]);   % Not-Buy band
        end
        if any(isfinite(N_buy_mean{a})),    hleg(end+1) = plot(ax(a), bin_centers, N_buy_mean{a},    'Color',[0 0.6 0], 'LineWidth', 2.5); end %#ok<AGROW>
        if any(isfinite(N_notbuy_mean{a})), hleg(end+1) = plot(ax(a), bin_centers, N_notbuy_mean{a}, 'Color',[0 0 0],   'LineWidth', 2.5); end %#ok<AGROW>
        xline(ax(a), 0, 'k--', 'LineWidth', 1);
        grid(ax(a), 'on'); ylim(ax(a), [-0.1 1]); xlim(ax(a), [tmin tmax]);
        if a == 1, ylabel(ax(a), 'Mean normalized activity'); end
        xlabel(ax(a), {'Time','(ms)'});

        name  = char(align_event_names(a));  % e.g. 'm1trialStart'
        % Split at:
        %  - letter/digit runs starting lowercase (m1, trial, on, off, ...)
        %  - camel-case chunks starting with uppercase (Start, Stop, ...)
        parts = regexp(name, '(?:[a-z]+[0-9]*|[A-Z][a-z0-9]*)', 'match');        
        multi = strjoin(parts, sprintf('\n'));  % "m1\ntrial\nStart"
        title(ax(a), multi, 'Interpreter','none');

        if ~isempty(hleg) && a == nAlign
            legend(ax(a), hleg, {sprintf('Buy (trials=%d)', tot_buy(a)), sprintf('Not-Buy (trials=%d)', tot_notbuy(a))}, 'Location','northeast');
        end
    else
        if ~strcmpi(error_type,'none')
            shaded(ax(a), bin_centers, hold_lo{a}, hold_hi{a}, [0 0.4 1]);    % Hold band
            shaded(ax(a), bin_centers, buy_lo{a} , buy_hi{a} , [0 0.6 0]);    % Buy band
            shaded(ax(a), bin_centers, sell_lo{a}, sell_hi{a}, [0.85 0 0]);   % Sell band
        end
        if any(isfinite(N_buy_mean{a})),  hleg(end+1) = plot(ax(a), bin_centers, N_buy_mean{a},  'Color',[0 0.6 0], 'LineWidth', 2.5); end %#ok<AGROW>
        if any(isfinite(N_hold_mean{a})), hleg(end+1) = plot(ax(a), bin_centers, N_hold_mean{a}, 'Color',[0 0.4 1], 'LineWidth', 2.5); end %#ok<AGROW>
        if any(isfinite(N_sell_mean{a})), hleg(end+1) = plot(ax(a), bin_centers, N_sell_mean{a}, 'Color',[0.85 0 0], 'LineWidth', 2.5); end %#ok<AGROW>
        xline(ax(a), 0, 'k--', 'LineWidth', 1);
        grid(ax(a), 'on'); ylim(ax(a), [-0.1 1]); xlim(ax(a), [tmin tmax]);
        if a == 1, ylabel(ax(a), 'Mean normalized activity'); end
        xlabel(ax(a), {'Time','(ms)'});

        name  = char(align_event_names(a));  % e.g. 'm1trialStart'
        % Split at:
        %  - letter/digit runs starting lowercase (m1, trial, on, off, ...)
        %  - camel-case chunks starting with uppercase (Start, Stop, ...)
        parts = regexp(name, '(?:[a-z]+[0-9]*|[A-Z][a-z0-9]*)', 'match');        
        multi = strjoin(parts, sprintf('\n'));  % "m1\ntrial\nStart"
        title(ax(a), multi, 'Interpreter','none');

        if ~isempty(hleg) && a == nAlign
            legend(ax(a), hleg, {sprintf('Buy (trials=%d)', tot_buy(a)), sprintf('Hold (trials=%d)', tot_hold(a)), sprintf('Sell (trials=%d)', tot_sell(a))}, 'Location','northeast');
        end
    end
end

% Link y-limits across all panels; hide all but leftmost y-axis
linkaxes(ax, 'y');
for i = 2:numel(ax)
    ax(i).YAxis.Visible = 'off';
end
ax(1).YAxis.Visible = 'on';
set(ax, 'XGrid','off','YGrid','on');

sgtitle(sprintf('%s | Binning = %d ms | Smoothing ~%d ms', figTitle, bin_ms, smooth_ms));

FONT_AX   = 10;
FONT_LAB  = 10;
FONT_TIT  = 10;
FONT_SGT  = 12;
FONT_LEG  = 10;

all_ax = findall(gcf,'Type','axes');
for a = all_ax.'
    a.FontSize = FONT_AX;
    a.LabelFontSizeMultiplier = 1;
    a.TitleFontSizeMultiplier = 1;
    if ~isempty(a.XLabel), a.XLabel.FontSize = FONT_LAB; end
    if ~isempty(a.YLabel), a.YLabel.FontSize = FONT_LAB; end
    if ~isempty(a.Title),  a.Title.FontSize  = FONT_TIT; end
end

all_leg = findall(gcf,'Type','legend');
set(all_leg,'FontSize',FONT_LEG);

sg = findall(gcf,'Type','axes','Tag','suptitle');
if ~isempty(sg), sg.FontSize = FONT_SGT; end

end

% ------------------------------ Local helpers ------------------------------
function wins_mat = parse_plot_wins(plot_wins, nAlign)
% Normalize plot_wins into an nAlign x 2 numeric matrix of [tmin tmax] per event.
if isnumeric(plot_wins) && isequal(size(plot_wins), [1 2])
    wins_mat = repmat(plot_wins, nAlign, 1);
elseif isnumeric(plot_wins) && isequal(size(plot_wins), [nAlign 2])
    wins_mat = plot_wins;
elseif iscell(plot_wins) && numel(plot_wins) == nAlign
    wins_mat = nan(nAlign,2);
    for a = 1:nAlign
        w = plot_wins{a};
        if isnumeric(w) && numel(w)==2
            wins_mat(a,:) = w(:)';
        else
            error('plot_wins{%d} must be a numeric [tmin tmax].', a);
        end
    end
else
    error('plot_wins must be [1x2], [nAlignx2], or 1xN cell of [1x2].');
end
end

function out = ternary(cond, a, b)
% Simple ternary helper for strings
if cond, out = a; else, out = b; end
end

function drawAlphaLine(ax,x,y,color,alpha,lw)
% NaN-aware polyline with per-edge alpha using a patch trick.
% Splits at NaNs and draws contiguous segments with specified EdgeAlpha.
x = x(:); y = y(:);
m = isfinite(x) & isfinite(y);
if ~any(m), return; end
idx = find(m);
% Break into contiguous runs
breaks = [1; find(diff(idx) > 1)+1; numel(idx)+1];
for b = 1:numel(breaks)-1
    ii = idx(breaks(b):breaks(b+1)-1);
    if numel(ii) < 2, continue; end
    V = [x(ii) y(ii)];
    F = [(1:numel(ii)-1)' (2:numel(ii))'];  % open polyline (no closure)
    patch('Faces',F,'Vertices',V, ...
        'FaceColor','none','EdgeColor',color, ...
        'EdgeAlpha',alpha,'LineWidth',lw, ...
        'Parent',ax,'HandleVisibility','off');
end
end

function Y = nan_smooth_rows(X, k)
% NaN-aware moving average smoothing along columns for a 2D matrix.
% Each row is convolved with a length-k boxcar, ignoring NaNs via mask renormalization.
if k <= 1 || isempty(X), Y = X; return; end
w = ones(1, k) / k;                     % Boxcar kernel
Y = nan(size(X));
for i = 1:size(X,1)
    xi = X(i,:);
    m  = isfinite(xi);
    xi(~m) = 0;
    num = conv(xi, w, 'same');
    den = conv(double(m), w, 'same');
    yi = num ./ max(den, eps);
    yi(den==0) = NaN;
    Y(i,:) = yi;
end
end

function y = nan_smooth_row(x, k)
% NaN-aware moving average for a single row vector.
if k <= 1 || isempty(x), y = x; return; end
w = ones(1, k) / k;
m = isfinite(x);
x0 = x; x0(~m) = 0;
num = conv(x0, w, 'same');
den = conv(double(m), w, 'same');
y = num ./ max(den, eps);
y(den==0) = NaN;
end

function [vec, ok] = get_event_vec(ES, name)
% Extract event vector (ms or labels) from table/struct; ok=true on success.
ok  = false;
vec = [];
try
    if istable(ES)
        if any(strcmp(name, ES.Properties.VariableNames)), x = ES.(name); else, return, end
    elseif isstruct(ES)
        if isfield(ES, name), x = ES.(name); else, return, end
    else
        return
    end
    while iscell(x)
        if isempty(x), return, end
        x = x{1};
    end
    if istable(x)
        w = size(x,2); got = false;
        for i = 1:w
            col = x{:,i};
            while iscell(col) && ~isempty(col), col = col{1}; end
            if isnumeric(col) || islogical(col) || iscategorical(col) || isstring(col) || ischar(col)
                x = col; got = true; break
            end
        end
        if ~got, return, end
    end
    while iscell(x)
        if isempty(x), return, end
        x = x{1};
    end
    if iscategorical(x), x = double(x);
    elseif isstring(x) || ischar(x), x = str2double(string(x));
    elseif islogical(x), x = double(x);
    elseif ~isnumeric(x), return
    end
    vec = double(x(:));
    ok  = ~isempty(vec) && any(isfinite(vec));
catch
    ok = false; vec = [];
end
end

function [lo, hi] = compute_bands(M, error_type, ci_alpha)
% Compute SEM/CI bands across rows of M (observations x time); NaN-robust.
if nargin < 2 || isempty(error_type), error_type = 'ci'; end
if nargin < 3 || isempty(ci_alpha),   ci_alpha   = 0.95; end
mu  = mean(M, 1, 'omitnan');
sd  = std(M, 0, 1, 'omitnan');
n   = sum(isfinite(M), 1);
sem = sd ./ max(1, sqrt(n));
switch lower(error_type)
    case 'sem'
        lo = mu - sem; hi = mu + sem;
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

function shaded(ax, t, lo, hi, col)
% Draw a shaded band between lo and hi at times t on axes ax.
mask = isfinite(t) & isfinite(lo) & isfinite(hi);
if ~any(mask), return; end %#ok<*NOPTS>
idx = find(mask);
starts = [idx(1), idx(find(diff(idx) > 1) + 1)];
ends   = [idx(find(diff(idx) > 1)), idx(end)];
for s = 1:numel(starts)
    ii = starts(s):ends(s);
    xx = t(ii); yy1 = lo(ii); yy2 = hi(ii);
    fill(ax, [xx, fliplr(xx)], [yy1, fliplr(yy2)], col, ...
        'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility','off');
end
end

function C = filterNeuronAndAddType(C, prefix)
% Build C.neuronType from built-in neuron lists for C.cond.
% Types: 1=Choice–Reward Sustained, 2=Choice-Related Transient, 3=Choice-Suppressed

nSess = numel(C.rawSpikes);
C.neuronType = cell(nSess,1);
for s = 1:nSess
    C.neuronType{s} = zeros(numel(C.rawSpikes{s}),1,'uint8');  % default 0 = unlabeled
end

% ---------- DEFINE LABEL LISTS ----------
if strcmp(prefix, '')
    switch upper(string(C.cond))
        case "LIVE"
            sessions1 = [3 3 3 3 3 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
            neurons1  = [2 3 8 12 16 7 8 9 11 4 5 8 9 12 14 2 4 5 6 10 17 3 4 6 8 9 13 15 18 20 21 22 23 29 1 2 3 4 5 7 11 15 16 17 19 22 23 24 27 30 31 32 34 38 39 40 46 47 1 2 3 5 6 10 13 14 15 17 18 21 24 26 29 31 33 36];
            sessions2 = [1 1 1 3 3 3 3 6 6 6 6 7 7 7 7 8 8 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10];
            neurons2  = [1 2 3 1 4 5 24 10 16 17 18 12 13 14 15 7 19 8 9 18 20 21 25 28 43 8 11 20 22 23 25 30 38 39];
            sessions3 = [2 2 6 7 7 7 8 9 9];
            neurons3  = [3 4 3 3 8 9 28 36 45];
        case "REPLAY"
            sessions1 = [3 4 4 4 4 4 4 4 4 4 4 4 4 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
            neurons1  = [2 1 2 3 4 5 8 9 10 11 13 15 16 2 3 5 7 15 16 3 4 5 7 9 10 12 15 16 17 18 19 20 21 22 2 3 4 5 6 8 10 11 13 14 15 19 20 21 22 24 25 28 29 30 31 32 33 41 42];
            sessions2 = [5 5 5 5 8 8 8 10 10 10 10 10 10 10];
            neurons2  = [4 8 9 8 1 11 12 7 9 16 27 35 36 38];
            sessions3 = [1 1 5 5 5 5 5 5 5 8 8 10];
            neurons3  = [1 2 1 3 5 6 7 10 11 4 10 40];
        case "AI"
            sessions1 = [2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 5 5 5 5 5 5 7 7 7 7 7 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
            neurons1  = [1 3 4 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 5 6 9 12 16 17 3 10 15 17 18 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 25 3 4 5 8 9 10 11 12 13 14 16 17 18 19 20 23 24 27 28 33 34];
            sessions2 = [2 2 2 5 5 5 5 7 7 9 10 10 10 10];
            neurons2  = [6 7 8 4 8 11 13 12 14 27 6 7 22 31];
            sessions3 = [2 2 7 7 7 7 9 10];
            neurons3  = [5 9 5 6 8 13 24 32];
        case "DECOY"
            sessions1 = [1 1 1 1 2 2 2 2 5 5 5 6 6 6 6 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 11 11 11 11 11 11 12 12 12 12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14];
            neurons1  = [1 2 4 6 1 2 3 5 2 3 5 1 2 3 4 2 5 1 2 3 4 5 7 8 10 5 6 14 15 17 3 5 14 16 17 18 1 2 3 5 6 7 8 9 10 11 12 13 14 15 17 19 2 3 4 5 6 13 15 16 17 19 20 22 24 25 30 31 32 34 35 36 37 39 43 44 3 6 7 10 11 12 13 14 15 19 20 21 22 25 26 27 28 29 31];
            sessions2 = [2 4 8 9 9 9 9 11 11 11 11 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14];
            neurons2  = [9 1 6 13 16 19 20 6 8 12 13 8 9 10 12 18 21 26 27 40 8 9 16 17 18 23 30];
            sessions3 = [1 1 1 2 2 2 4 9 11 13 13];
            neurons3  = [7 9 11 6 8 10 2 11 10 33 42];
        otherwise
            sessions1 = []; neurons1 = [];
            sessions2 = []; neurons2 = [];
            sessions3 = []; neurons3 = [];
    end

elseif strcmp(prefix, "OT ")
    % Default: no groups defined yet for OT -> leave lists empty
    sessions1 = []; neurons1 = [];
    sessions2 = []; neurons2 = [];
    sessions3 = []; neurons3 = [];

    switch upper(string(C.cond))
        case "LIVE"
        case "REPLAY"
        case "AI"
        case "DECOY"
    end

elseif strcmp(prefix, "Saline ")
    % Default: no groups defined yet for Saline -> leave lists empty
    sessions1 = []; neurons1 = [];
    sessions2 = []; neurons2 = [];
    sessions3 = []; neurons3 = [];

    switch upper(string(C.cond))
        case "LIVE"
        case "REPLAY"
        case "AI"
        case "DECOY"
    end

else
    sessions1 = []; neurons1 = [];
    sessions2 = []; neurons2 = [];
    sessions3 = []; neurons3 = [];
end

% ---------- APPLY LABELS WHERE DEFINED ----------
sessLists = {sessions1, sessions2, sessions3};
neurLists = {neurons1,  neurons2,  neurons3};
for t = 1:3
    ss = sessLists{t}; nn = neurLists{t};
    for i = 1:numel(ss)
        s = ss(i); n = nn(i);
        if s>=1 && s<=nSess && n>=1 && n<=numel(C.rawSpikes{s})
            C.neuronType{s}(n) = uint8(t);
        end
    end
end

% ---------- DEFAULT: IF NO GROUPS, USE ALL NEURONS AS TYPE 1 ----------
anyLabeled = false;
for s = 1:nSess
    if any(C.neuronType{s} ~= 0)
        anyLabeled = true;
        break;
    end
end

if ~anyLabeled
    for s = 1:nSess
        C.neuronType{s}(:) = uint8(1);   % everyone is group 1
    end
end

end

function stackFiguresVertical(monitorIdx, margin, gap)
% Stack all open figure windows vertically on a selected monitor.
if nargin < 1 || isempty(monitorIdx), monitorIdx = 1; end
if nargin < 2 || isempty(margin),     margin     = [300 100 100 100]; end % [left right bottom top] px
if nargin < 3 || isempty(gap),        gap        = 20; end
figs = findall(groot,'Type','figure','Visible','on');
if isempty(figs), return; end
[~, ix] = sort([figs.Number]); figs = figs(ix);
set(figs, 'WindowState','normal', 'Units','pixels');
mp = get(0,'MonitorPositions');                      % [x y w h] per monitor
monitorIdx = max(1, min(monitorIdx, size(mp,1)));
M = mp(monitorIdx,:);
left   = M(1) + margin(1);
right  = M(1) + M(3) - margin(2);
bottom = M(2) + margin(3);
top    = M(2) + M(4) - margin(4);
usableW = max(1, right - left);
usableH = max(1, top - bottom);
n = numel(figs);
h = (usableH - gap*(n-1)) / n;  % height per figure
w = usableW;
for i = 1:n
    y = bottom + (n - i) * (h + gap);               % stack from top to bottom
    set(figs(i), 'Position', [left, y, w, h]);
end
end
