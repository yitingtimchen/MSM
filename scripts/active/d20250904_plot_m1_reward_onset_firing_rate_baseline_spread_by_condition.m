% File: d20250904_plot_m1_reward_onset_firing_rate_baseline_spread_by_condition.m
% =========================================================================
% Overlay & compare PSTHs across multiple conditions
% AND quantify baseline "spread" (% of post-event firing) at block & trial levels.
%
% NEW (this version):
% • PSTH figure shows ONLY the condition-level re-normalized curves ("C"):
%   - For each normalization row and condition:
%       1) Re-baseline the condition's mean PSTH using BASELINE_WIN.
%       2) Z-score using μ/σ computed over the ENTIRE window (not just baseline).
%          Option: σ pooled across conditions (default) or per condition.
%   - SEM is scaled by the same σ used for that condition (or pooled σ).
%
% • Multiple event variant sets supported; one PSTH figure per set.
% • Baseline spread analysis (with gating) unchanged and still available.
%
% IMPORTANT:
% - All trial-level smoothing occurs BEFORE any stats/normalization.
% - Normalization B inside compute_psth_for_pack uses full-window μ/σ on
%   the baseline-subtracted signal with scope-matched σ (trial/block/neuron)
%   and optional variance shrinkage (kept for completeness).
%
% UPDATE (multi-condition sets):
% - Define multiple condition sets (cond_sets) to compare.
% - For each event set × condition set, produce separate figures/outputs.
% - Output directory is namespaced by event tag and condition-set tag.
% =========================================================================

close all;
clc;
clear;

%% ----------------------------- TOGGLES -----------------------------------
PLOT_BASELINE_SPREADS = 1;  % combined figure across conditions (per event set & cond-set)
SORT_BY_TRIAL_SPREAD  = 1;  % sort neurons (within each condition segment) by trial-level % range (desc)

DRAW_CHOICE_BOX   = 1;                  % highlight choice window on PSTH
CHOICE_BOX_WIN    = [-0.25, -0.05];     % seconds (start, end)
CHOICE_BOX_COLOR  = [1.00, 0.92, 0.50]; % light yellow
CHOICE_BOX_ALPHA  = 0.18;               % transparency

% Condition-level re-normalization (column C)
DO_CONDITION_LEVEL_RENORM = 1;
% 'pooled' = single μ/σ across all conditions (best cross-condition comparability)
% 'per-condition' = each condition gets its own μ/σ
COND_RENORM_MODE = 'pooled';

% Gating thresholds (used only for baseline-spread analysis)
RATE_FLOOR_HZ        = 1.0; % minimum acceptable POST mean FR to consider trial/block/session
DELTA_FLOOR_HZ       = 5.0; % absolute minimum post−baseline increase (Hz) to count as "evoked"
K_SIGMA              = 1.0; % require Δ_t ≥ K_SIGMA * σ_BASE_t
MIN_TRIALS_PER_BLOCK = 8;   % min kept trials in a block to compute trial-level %
MIN_BLOCKS_WITH_DATA = 3;   % min blocks with kept trials to compute block-level %

is_mac = 0;

% Figure sizing / typography
FIG_SIZE_SPREADS = [80, 80, 1500, 400]; % [x y width height] for spread plots (1×2 grid)
FZ_AX    = 12; % axis tick/label size
FZ_TITLE = 12; % panel title size
FZ_SGT   = 14; % super-title size
FZ_LEGEND= 10;

%% ----------------------------- Inputs ------------------------------------
if is_mac
    packDir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs';
else
    packDir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs';
end

% === Multiple CONDITION SETS to compare (EDIT HERE) =======================
% Each cell is a set; inside each set, list the condition names to include.
% Example: three sets — all, task-only, control-only
cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    % {'OT AI', 'OT Replay', 'OT Decoy', 'OT Live'}; ...
    % {'Saline AI', 'Saline Replay', 'Saline Decoy', 'Saline Live'} ...
};

% >>> Multiple event variant sets (edit here) <<<
event_variants = { ...
    "m1rew1on", ...
    % "m2rew1on"
};

% Analysis windows (seconds, relative to event)
win          = [-0.5, 1.5];   % master window for extraction (must cover both sub-windows)
BASELINE_WIN = [-0.5, -0.4];  % baseline window for re-baseline step (and spread)
POST_WIN     = [ 0.0,  1.5];  % post-event window (used in spread analysis only)

% PSTH normalization strategies (rows)
% norm_methods = {'none','baseline_sub_trial','baseline_sub_neuron','baseline_sub_chunk','zscore_neuron'};
norm_methods = {'baseline_sub_trial'};
scale_to_hz        = true;   % convert 1 ms binary → Hz (divide by dt)
smooth_sigma_sec   = 0.02;   % Gaussian sigma (sec); 0 disables smoothing (PRE-analysis everywhere)
show_sem_shading   = true;   % shaded SEM per condition in PSTH overlays

% Base output directory (per event set × condition-set will get subfolders)
BASE_OUT_DIR = fullfile(packDir, 'overlay_outputs_multiEvent');
if ~exist(BASE_OUT_DIR, 'dir'), mkdir(BASE_OUT_DIR); end

% Task structure
TR = 15; % trials per block
CH = 6;  % number of blocks (chunks)

% --------- Z-score stability controls (used inside PSTH function for B) ---
use_robust_sigma        = false; % if true, use MAD*1.4826; else std
sigma_shrink_lambda     = 0.0;   % keep 0 to preserve scope differences cleanly
sigma_ridge             = 1e-6;  % small ridge added to variance (Hz^2)
trial_sigma_min_samples = 2;     % fallback if too few baseline bins (baseline stats)
sigma_floor             = 0.5;   % σ floor (Hz or arb. units after subtraction)

% =========================================================================
% Storage for aggregated spike matrices
% Organized as AllSpikes.(condset).(event).(condition).{raw, base, z}
% =========================================================================
AllSpikes = struct();

%% --------------------------- Prep time axis -------------------------------
% Pull dt from the first available pack across ANY condition in ANY set.
dt = [];
all_conds_flat = unique([cond_sets{:}]); % flatten cell-of-cell to 1×N cellstr
for ci = 1:numel(all_conds_flat)
    fp = fullfile(packDir, sprintf('%s_condition_pack.mat', all_conds_flat{ci}));
    if exist(fp, 'file')
        S = load(fp, 'C');
        dt = S.C.dt;
        break;
    end
end
if isempty(dt), error('No condition packs found in %s (checked across all cond_sets).', packDir); end

% Build time axis + masks
t        = (win(1):dt:win(2)); t = t(:)'; % enforce row
baseMask = (t >= BASELINE_WIN(1)) & (t <= BASELINE_WIN(2));
postMask = (t >= POST_WIN(1))     & (t <= POST_WIN(2));

%% ===================== LOOP OVER EVENT & CONDITION SETS ===================
for es = 1:numel(event_variants)
    event_variant = event_variants{es};
    ev_tag = event_variant; % short label for titles/files

    for cs = 1:numel(cond_sets)
        conds = cond_sets{cs};
        cond_tag = cond_set_label(conds); % label like "AI-Replay-Decoy-Live"
        spike_mat = cell(1, length(conds)); % initialize storage for spike matrices

        % Namespace outputs by event tag and condition-set tag
        outDir = fullfile(BASE_OUT_DIR, sprintf('%s__%s', ev_tag, cond_tag));
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        %% ------------ Combined Baseline Spread Figure (all conditions) -------
        if PLOT_BASELINE_SPREADS
            % Prepare colors and accumulators
            cond_cols = lines(numel(conds));
            X_block_all  = []; Y_block_all  = []; C_block_all  = []; % scatter
            X_trial_all  = []; Y_trial_all  = []; C_trial_all  = []; % box data
            seg_edges    = []; seg_labels   = {};                   % condition segment edges and labels

            x_offset = 0; % running x offset (concatenate condition segments)
            for ci = 1:numel(conds)
                cond = conds{ci};
                packPath = fullfile(packDir, sprintf('%s_condition_pack.mat', cond));
                if ~exist(packPath, 'file')
                    warning('Pack not found for spread analysis: %s (skipping)', packPath);
                    continue;
                end

                [block_pct_by_neuron, trial_pct_by_neuron] = compute_baseline_spreads_for_pack( ...
                    packPath, t, baseMask, postMask, event_variant, TR, CH, scale_to_hz, smooth_sigma_sec, ...
                    use_robust_sigma, trial_sigma_min_samples, RATE_FLOOR_HZ, DELTA_FLOOR_HZ, K_SIGMA, ...
                    MIN_TRIALS_PER_BLOCK, MIN_BLOCKS_WITH_DATA);

                % Sorting WITHIN CONDITION by trial-level % range
                Nneur = numel(trial_pct_by_neuron);
                trial_spread_range = nan(Nneur,1);
                for n = 1:Nneur
                    v = trial_pct_by_neuron{n};
                    v = v(isfinite(v));
                    if ~isempty(v)
                        trial_spread_range(n) = max(v) - min(v);
                    end
                end
                if SORT_BY_TRIAL_SPREAD
                    tmp = trial_spread_range; tmp(~isfinite(tmp)) = -Inf;
                    [~, sort_idx_desc] = sort(tmp, 'descend');
                    posMap = zeros(Nneur,1); posMap(sort_idx_desc) = 1:numel(sort_idx_desc);
                else
                    posMap = (1:Nneur).';
                end

                % Block-level (ONE per neuron)
                for n = 1:Nneur
                    vals = block_pct_by_neuron{n};
                    vals = vals(isfinite(vals));
                    if isempty(vals), continue; end
                    X_block_all  = [X_block_all;  x_offset + posMap(n)]; %#ok<AGROW>
                    Y_block_all  = [Y_block_all;  vals(1)];              %#ok<AGROW>
                    C_block_all  = [C_block_all;  ci];                   %#ok<AGROW>
                end

                % Trial-level (≤6 per neuron; one per block)
                for n = 1:Nneur
                    vals = trial_pct_by_neuron{n};
                    vals = vals(isfinite(vals));
                    if isempty(vals), continue; end
                    X_trial_all  = [X_trial_all;  (x_offset + posMap(n)) .* ones(numel(vals),1)]; %#ok<AGROW>
                    Y_trial_all  = [Y_trial_all;  vals(:)];                                      %#ok<AGROW>
                    C_trial_all  = [C_trial_all;  ci .* ones(numel(vals),1)];                    %#ok<AGROW>
                end

                % Segment bookkeeping
                if Nneur > 0
                    seg_edges = [seg_edges; x_offset + 0.5, x_offset + Nneur + 0.5]; %#ok<AGROW>
                    seg_labels{end+1} = cond; %#ok<AGROW>
                    x_offset = x_offset + Nneur + 2; % add small gap between condition segments
                end
            end

            % ---- Plot the combined figure (two panels) ----
            figSP = figure('Name', sprintf('Baseline Spread %% — %s — %s', ev_tag, cond_tag), ...
                           'Color','w', 'Units','pixels', 'Position', FIG_SIZE_SPREADS);
            tl = tiledlayout(figSP, 1, 2, 'TileSpacing','compact', 'Padding','compact');

            infoStr = sprintf('BL %d..%d ms\nPOST %d..%d ms\n\\sigma = %d ms', ...
                round(1000*BASELINE_WIN(1)), round(1000*BASELINE_WIN(2)), ...
                round(1000*POST_WIN(1)),     round(1000*POST_WIN(2)), ...
                round(1000*smooth_sigma_sec));

            % (Left) Block-level scatter (one value/neuron)
            axL = nexttile(tl, 1);  hold(axL, 'on');
            if ~isempty(X_block_all)
                for ci = 1:numel(conds)
                    pick = (C_block_all == ci);
                    if any(pick)
                        jitter = 0.05 .* randn(sum(pick),1);
                        scatter(axL, X_block_all(pick) + jitter, Y_block_all(pick), 24, ...
                                'filled', 'MarkerFaceAlpha',0.8, 'MarkerEdgeAlpha',0.6, ...
                                'MarkerFaceColor', cond_cols(ci,:), 'DisplayName', conds{ci});
                    end
                end
            end
            xlabel(axL, 'Neuron index (concatenated by condition)');
            ylabel(axL, 'Spread (% of post-event mean)');
            title(axL, sprintf('[Left] Block-level spread%% per neuron\n%s', infoStr), 'Interpreter','tex');
            if ~isempty(seg_edges), draw_cond_segments(axL, seg_edges, seg_labels); end
            grid(axL, 'on'); set(axL,'FontSize',FZ_AX); axL.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
            legend(axL, 'Location','bestoutside');

            % (Right) Trial-level box plots (≤6 values/neuron; one per block), colored by condition
            axR = nexttile(tl, 2);  hold(axR, 'on');
            usedBoxchart = false;
            try
                for ci = 1:numel(conds)
                    pick = (C_trial_all == ci);
                    if any(pick)
                        bc = boxchart(axR, X_trial_all(pick), Y_trial_all(pick), 'BoxWidth', 0.6);
                        bc.BoxFaceColor = cond_cols(ci,:);
                        bc.BoxFaceAlpha = 0.35;
                        bc.MarkerColor  = cond_cols(ci,:);
                        bc.WhiskerLineColor = cond_cols(ci,:);
                    end
                end
                usedBoxchart = true;
            catch
                % Fallback (no color grouping support here)
                plot_box_groups(axR, X_trial_all, Y_trial_all, 0.38);
            end
            xlabel(axR, 'Neuron index (concatenated by condition)');
            ylabel(axR, 'Spread (% of post-event mean)');
            title(axR, sprintf('[Right] Trial-level spread%% per neuron\n%s', infoStr), 'Interpreter','tex');
            if ~isempty(seg_edges), draw_cond_segments(axR, seg_edges, seg_labels); end
            grid(axR, 'on'); set(axR,'FontSize',FZ_AX); axR.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
            if usedBoxchart, legend(axR, conds, 'Location','bestoutside'); end

            sg = sgtitle(tl, sprintf('Baseline Spread as %% of Post-Event FR — %s — %s', ev_tag, cond_tag), 'Interpreter','none');
            set(sg,'FontSize',FZ_SGT);

            % Save
            savefig_smart(figSP, fullfile(outDir, sprintf('baseline_spread_percent_%s_CONDSET_%s_gated.png', ev_tag, cond_tag)));
        end

    end
end

%% ============================== FUNCTIONS ================================
function label = cond_set_label(conds)
% Build a label from a list of condition names (cellstr), joined with '-'
% and sanitized to be filesystem-friendly.
if isempty(conds)
    label = 'conds';
    return;
end
if isstring(conds), conds = cellstr(conds); end
label = strjoin(conds, '-');
label = regexprep(label, '[^\w-]', '');
end

function draw_cond_segments(ax, seg_edges, seg_labels)
% (kept for spread plots; not used in PSTH C-only figure)
yl = ylim(ax);
for k = 1:size(seg_edges,1)
    x0 = seg_edges(k,1); x1 = seg_edges(k,2);
    plot(ax, [x0, x0], yl, '-', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.5, 'HandleVisibility','off');
    plot(ax, [x1, x1], yl, '-', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.5, 'HandleVisibility','off');
    xc = (x0 + x1) ./ 2.0;
    text(ax, xc, yl(2), [' ' seg_labels{k} ' '], 'HorizontalAlignment','center', ...
        'VerticalAlignment','top', 'Color',[0.2,0.2,0.2], 'FontWeight','bold', 'Rotation',0, ...
        'BackgroundColor',[1,1,1,0.6], 'Margin',2, 'Clipping','on');
end
uistack(findobj(ax,'Type','text'), 'top');
end

function [block_pct_by_neuron, trial_pct_by_neuron] = compute_baseline_spreads_for_pack( ...
    packPath, t, baseMask, postMask, event_variant, TR, CH, ...
    scale_to_hz, smooth_sigma_sec, ...
    use_robust_sigma, trial_sigma_min_samples, RATE_FLOOR_HZ, DELTA_FLOOR_HZ, K_SIGMA, ...
    MIN_TRIALS_PER_BLOCK, MIN_BLOCKS_WITH_DATA)
% (unchanged from your gated version)
S = load(packPath,'C');
C = S.C;
dt = C.dt;
wRelBins = round(t ./ dt);
block_pct_by_neuron = {};
trial_pct_by_neuron = {};
for f = 1:numel(C.files)
    if isempty(C.eventTables{f}), continue; end
    props = string(C.eventTables{f}.Properties.VariableNames);
    pick  = find(ismember(event_variant, props), 1, 'first');
    if isempty(pick), continue; end
    varName = char(event_variant(pick));
    idxMat = C.eventTables{f}.(varName){1};
    if ~ismatrix(idxMat), continue; end
    Sbin     = C.S{f};
    [Nf, Tf] = size(Sbin);

    row_centers  = [];
    row_trialIdx = [];
    row_blockIdx = [];
    for m = 1:CH
        for tr = 1:TR
            center = idxMat(tr, m);
            if ~isnan(center) && center >= 1 && center <= Tf
                row_centers(end+1,1)  = center; %#ok<AGROW>
                row_trialIdx(end+1,1) = tr;     %#ok<AGROW>
                row_blockIdx(end+1,1) = m;      %#ok<AGROW>
            end
        end
    end
    if isempty(row_centers), continue; end

    for n = 1:Nf
        R    = numel(row_centers);
        Tmat = nan(R, numel(t));
        for r = 1:R
            segIdx = row_centers(r) + wRelBins;
            valid  = (segIdx >= 1) & (segIdx <= Tf);
            if any(valid)
                Tmat(r, valid) = double(Sbin(n, segIdx(valid)));
            end
        end
        if scale_to_hz, Tmat = Tmat ./ dt; end
        if smooth_sigma_sec > 0
            sigma_bins = smooth_sigma_sec ./ dt;
            for rr = 1:R
                Tmat(rr,:) = conv_gauss_reflect(Tmat(rr,:), sigma_bins);
            end
        end
        b_rows_mu    = mean(Tmat(:, baseMask), 2, 'omitnan');
        post_rows_mu = mean(Tmat(:, postMask), 2, 'omitnan');

        b_rows_sd    = nan(R,1);
        n_base_bins  = sum(isfinite(Tmat(:, baseMask)), 2);
        for r = 1:R
            if n_base_bins(r) >= trial_sigma_min_samples
                b_rows_sd(r) = robust_or_std(Tmat(r, baseMask), use_robust_sigma);
            else
                b_rows_sd(r) = NaN;
            end
        end
        xb_all  = Tmat(:, baseMask); xb_all = xb_all(isfinite(xb_all));
        sig_neu = iff(isempty(xb_all), NaN, robust_or_std(xb_all, use_robust_sigma));

        trial_pct_vals        = nan(0,1);
        block_baseline_means  = nan(1, CH);
        kept_post_all_session = [];

        for m = 1:CH
            rows_m = find(row_blockIdx == m);
            if isempty(rows_m), continue; end
            xb_block = Tmat(rows_m, baseMask); xb_block = xb_block(isfinite(xb_block));
            sig_blk  = iff(isempty(xb_block), NaN, robust_or_std(xb_block, use_robust_sigma));

            keep_mask = false(numel(rows_m),1);
            for k = 1:numel(rows_m)
                rIdx  = rows_m(k);
                sig_r = b_rows_sd(rIdx);
                if ~isfinite(sig_r), sig_r = sig_blk; end
                if ~isfinite(sig_r), sig_r = sig_neu; end
                if ~isfinite(sig_r), sig_r = NaN; end

                delta  = post_rows_mu(rIdx) - b_rows_mu(rIdx);
                thresh = max(DELTA_FLOOR_HZ, iff(isfinite(sig_r), K_SIGMA .* sig_r, 0));
                keep_mask(k) = (post_rows_mu(rIdx) >= RATE_FLOOR_HZ) && (delta >= thresh);
            end

            kept_idx_block = rows_m(keep_mask);
            if numel(kept_idx_block) >= MIN_TRIALS_PER_BLOCK
                bvals_kept   = b_rows_mu(kept_idx_block);
                post_kept    = post_rows_mu(kept_idx_block);
                spread_block = max(bvals_kept) - min(bvals_kept);

                denom_block_post = mean(post_kept, 'omitnan');
                if isfinite(denom_block_post) && (denom_block_post >= RATE_FLOOR_HZ)
                    pct_trial_level_block = 100 .* spread_block ./ denom_block_post;
                else
                    pct_trial_level_block = NaN;
                end
                trial_pct_vals(end+1,1) = pct_trial_level_block; %#ok<AGROW>
                block_baseline_means(m) = mean(bvals_kept, 'omitnan');
                kept_post_all_session   = [kept_post_all_session; post_kept(:)]; %#ok<AGROW>
            else
                trial_pct_vals(end+1,1) = NaN; %#ok<AGROW>
                block_baseline_means(m) = NaN;
            end
        end

        bmeans_finite = block_baseline_means(isfinite(block_baseline_means));
        if numel(bmeans_finite) >= MIN_BLOCKS_WITH_DATA
            denom_session_post = mean(kept_post_all_session, 'omitnan');
            if isfinite(denom_session_post) && (denom_session_post >= RATE_FLOOR_HZ)
                spread_blocks   = max(bmeans_finite) - min(bmeans_finite);
                pct_block_level = 100 .* spread_blocks ./ denom_session_post;
            else
                pct_block_level = NaN;
            end
        else
            pct_block_level = NaN;
        end

        block_pct_by_neuron{end+1,1} = pct_block_level; %#ok<AGROW>
        trial_pct_by_neuron{end+1,1} = trial_pct_vals;  %#ok<AGROW>
    end
end
end

function [muA, semA, muB, semB, Nused, nfiles, NpbB, sdB, AllSpikes] = compute_psth_for_pack( ...
    packPath, t, baseMask, event_variant, TR, CH, norm_method, ...
    scale_to_hz, smooth_sigma_sec, ...
    use_robust_sigma, sigma_shrink_lambda, sigma_ridge, ...
    trial_sigma_min_samples, sigma_floor, AllSpikes, conds, ci)
    % COMPUTE_PSTH_FOR_PACK
    % - Smooth each trial BEFORE stats.
    % - Build Tmat, then compute:
    %   A: baseline-sub (scope determined by norm_method)
    %   B: A + scope-matched z-score using μ/σ over the ENTIRE window of A
    %      (trial/block/neuron scope). Optional variance shrinkage knob.
    
    S = load(packPath,'C');
    C = S.C;
    dt = C.dt;
    
    wRelBins = round(t ./ dt);
    all_raw = [];
    all_A = [];
    all_B = [];
    nfiles = 0;
    
    for f = 1:numel(C.files)
        if isempty(C.eventTables{f}), continue; end
        props = string(C.eventTables{f}.Properties.VariableNames);
        pick  = find(ismember(event_variant, props), 1, 'first');
        if isempty(pick), continue; end
        varName = char(event_variant(pick));
    
        idxMat = C.eventTables{f}.(varName){1};
        if ~ismatrix(idxMat), continue; end
    
        Sbin = C.S{f};
        [Nf, Tf] = size(Sbin);
    
        row_centers  = [];
        row_blockIdx = [];
        for m = 1:CH
            for tr = 1:TR
                center = idxMat(tr, m);
                if ~isnan(center) && center >= 1 && center <= Tf
                    row_centers(end+1,1)  = center; %#ok<AGROW>
                    row_blockIdx(end+1,1) = m;      %#ok<AGROW>
                end
            end
        end
        if isempty(row_centers), continue; end
        nfiles = nfiles + 1;
    
        for n = 1:Nf
            R    = numel(row_centers);
            Tmat = nan(R, numel(t));
            for r = 1:R
                segIdx = row_centers(r) + wRelBins;
                valid  = (segIdx >= 1) & (segIdx <= Tf);
                if any(valid)
                    Tmat(r, valid) = double(Sbin(n, segIdx(valid)));
                end
            end
    
            if scale_to_hz, Tmat = Tmat ./ dt; end
            if smooth_sigma_sec > 0
                sigma_bins = smooth_sigma_sec ./ dt;
                for rr = 1:R
                    Tmat(rr,:) = conv_gauss_reflect(Tmat(rr,:), sigma_bins);
                end
            end
    
            % ---- Baseline stats for subtraction scopes (use BASELINE_WIN) ----
            b_rows_mu = mean(Tmat(:, baseMask), 2, 'omitnan');  % trial-level μ
            xb        = Tmat(:, baseMask); xb = xb(isfinite(xb));
            mu_neu    = iff(isempty(xb), 0, mean(xb));
            mu_block  = nan(1,CH);
            for m = 1:CH
                rows_m    = find(row_blockIdx == m);
                xm        = Tmat(rows_m, baseMask);
                xm        = xm(isfinite(xm));
                mu_block(m) = iff(isempty(xm), NaN, mean(xm));
            end
    
            % ------------------- A: baseline-sub -------------------------------
            Tmat_A = Tmat;
            switch char(norm_tag(norm_method))
                case 'none'
                    % no-op
                case 'baseline_sub_trial'
                    Tmat_A = Tmat_A - b_rows_mu;
                case 'baseline_sub_neuron'
                    Tmat_A = Tmat_A - mu_neu;
                case 'baseline_sub_chunk'
                    for m = 1:CH
                        rows_m = find(row_blockIdx == m);
                        bm     = mu_block(m);
                        if ~isfinite(bm), bm = 0; end
                        Tmat_A(rows_m,:) = Tmat_A(rows_m,:) - bm;
                    end
                case 'zscore_neuron'
                    % do a classic neuron z just for A in this mode
                    s_neu = robust_or_std(xb, use_robust_sigma);
                    if ~isfinite(s_neu) || s_neu == 0, s_neu = 1; end
                    Tmat_A = (Tmat - mu_neu) ./ max(s_neu, 1e-9);
                otherwise
                    error('Unknown norm_method: %s', norm_method);
            end

            % ------------------- B: A + scope-matched z (full window) ----------
            % Compute μ/σ over the ENTIRE window of Tmat_A at the chosen scope.
            Tmat_B = Tmat_A;
    
            switch char(norm_tag(norm_method))
                case {'baseline_sub_trial','baseline_sub_neuron','baseline_sub_chunk','none','zscore_neuron'}
                    xaA    = Tmat_A(:); xaA = xaA(isfinite(xaA));
                    muAall = mean(xaA, 'all', 'omitnan');
                    sdAall = std (xaA, [], 'all', 'omitnan');
                    if ~isfinite(sdAall) || sdAall <= 0, sdAall = sigma_floor; end
                    % shrink-to-self == no-op except ridge
                    if sigma_shrink_lambda > 0
                        var_shrunk = sdAall.^2 + sigma_ridge;
                        sdAall = max(sqrt(var_shrunk), sigma_floor);
                    end
                    Tmat_B = (Tmat_A - muAall) ./ sdAall;
            end
    
            % ---------------------- Per-neuron PSTHs ---------------------------
            psth_raw = mean(Tmat, 1, 'omitnan');
            psth_A = mean(Tmat_A, 1, 'omitnan');
            psth_B = mean(Tmat_B, 1, 'omitnan');
            all_raw(end+1,:) = psth_raw; %#ok<AGROW>
            all_A(end+1,:) = psth_A; %#ok<AGROW>
            all_B(end+1,:) = psth_B; %#ok<AGROW>
        end
    end
    
    spike_mat = all_raw;
    spike_mat_bs = all_A;
    spike_mat_z = all_B;

    % --------------------------------------------------------------
    % Store in organized container
    % --------------------------------------------------------------
    condset_label = replace(cell2mat(conds), ' ', '_');   % e.g. 'AI_Replay_Decoy_Live'
    event_label   = event_variant;          % e.g. 'stimOn'
    cond_label    = replace(conds{ci}, ' ', '_');      % e.g. 'AI'
    
    % Ensure nested fields exist
    if ~isfield(AllSpikes, condset_label)
        AllSpikes.(condset_label) = struct();
    end
    if ~isfield(AllSpikes.(condset_label), event_label)
        AllSpikes.(condset_label).(event_label) = struct();
    end
    if ~isfield(AllSpikes.(condset_label).(event_label), cond_label)
        AllSpikes.(condset_label).(event_label).(cond_label) = struct('raw',[],'base',[],'z',[]);
    end
    
    % Append along 3rd dimension = neurons
    AllSpikes.(condset_label).(event_label).(cond_label).raw  = cat(3, AllSpikes.(condset_label).(event_label).(cond_label).raw,  spike_mat);
    AllSpikes.(condset_label).(event_label).(cond_label).base = cat(3, AllSpikes.(condset_label).(event_label).(cond_label).base, spike_mat_bs);
    AllSpikes.(condset_label).(event_label).(cond_label).z    = cat(3, AllSpikes.(condset_label).(event_label).(cond_label).z,   spike_mat_z);

    % --- OPTIONAL NEURON FILTER ---
    % Loads once per session; does nothing if file/keys are missing.
    persistent KEEP_STATE;
    neuron_validation_mat = 'C:\Users\plattlab\MSM\outputs_local\neuron_validation\keep_indices.mat';
    if isempty(KEEP_STATE)
        if exist(neuron_validation_mat,'file')
            Skeep = load(neuron_validation_mat','keep_indices','variant');
            KEEP_STATE = Skeep;          % fields: keep_indices.(condset).(group).(variant) = [idx...]
        else
            KEEP_STATE = struct();       % no-op
        end
    end
    
    fprintf("\n" + string(cond_label) + "\n")
    fprintf("\tBefore: N = " + size(all_A, 1) + "\n");

    % Labels here already match tmp.m keys:
    %   condset_label = replace(cell2mat(conds), ' ', '_');   % e.g. 'AIReplayDecoyLive' or 'OT_AIOT_ReplayOT_DecoyOT_Live'
    %   cond_label    = replace(conds{ci}, ' ', '_');         % e.g. 'AI', 'OT_AI', 'Saline_Live'
    if isfield(KEEP_STATE,'keep_indices') && isfield(KEEP_STATE,'variant')
        setname = condset_label;
        grpname = cond_label;
        varname = KEEP_STATE.variant;    % 'raw'\n'base'\n'z' (must match which matrix you use below)
        if isfield(KEEP_STATE.keep_indices, setname) && ...
           isfield(KEEP_STATE.keep_indices.(setname), grpname) && ...
           isfield(KEEP_STATE.keep_indices.(setname).(grpname), varname)
            keep_j = KEEP_STATE.keep_indices.(setname).(grpname).(varname);
            % Apply to rows (neurons) of the corresponding matrices
            if ~isempty(keep_j)
                all_raw = all_raw(keep_j, :);
                all_A   = all_A(keep_j, :);
                all_B   = all_B(keep_j, :);
            end
        end
    end

    fprintf("\tAfter: N = " + size(all_A, 1) + "\n");

    [muA, semA, NpbA, sdA] = pop_mean_sem(all_A);
    [muB, semB, NpbB, sdB] = pop_mean_sem(all_B);
    Nused = size(all_A,1);
end

%% ------------------------------ Utils ------------------------------------
function out = norm_tag(x)
x = string(x);
if     strcmpi(x,'none'),               out = "none";
elseif strcmpi(x,'baseline_sub_trial'), out = "baseline_sub_trial";
elseif strcmpi(x,'baseline_sub_neuron'),out = "baseline_sub_neuron";
elseif strcmpi(x,'baseline_sub_chunk'), out = "baseline_sub_chunk";
elseif strcmpi(x,'zscore_neuron'),      out = "zscore_neuron";
else, error('Unknown norm_method: %s', x);
end
end

function s = robust_or_std(x, use_robust)
x = x(:); x = x(isfinite(x));
if isempty(x), s = NaN; return; end
if use_robust, s = 1.4826 .* mad(x,1);
else,          s = std(x,0);
end
end

function y = iff(cond, a, b)
if cond, y = a; else, y = b; end
end

function [mu, sem, Npb, sd] = pop_mean_sem(M)
% M: (#neurons x #time)
mu  = mean(M, 1, 'omitnan');
sd  = std (M, 0, 1, 'omitnan');
Npb = sum(~isnan(M), 1);
sem = sd ./ max(1, sqrt(Npb));
mu  = mu(:)'; sem = sem(:)'; Npb = Npb(:)'; sd = sd(:)';
end

function k = gauss1d(sigma_bins)
if sigma_bins <= 0, k = 1; return; end
hw = ceil(3.0 .* sigma_bins);
x  = -hw:hw;
k  = exp(-0.5 .* (x ./ sigma_bins).^2);
k  = k ./ sum(k);
end

function y = conv_gauss_reflect(x, sigma_bins)
if sigma_bins <= 0, y = x; return; end
x = x(:)'; k = gauss1d(sigma_bins);
hw = (numel(k) - 1) ./ 2;

m  = double(isfinite(x));
x0 = x; x0(~isfinite(x0)) = 0;

if hw < 1
    num = conv(x0, k, 'same');
    den = conv(m,  k, 'same');
    y   = num ./ den; y(den < eps) = NaN; return;
end

leftPadX = fliplr(x0(1:hw));
rightPadX= fliplr(x0(end-hw+1:end));
leftPadM = fliplr(m(1:hw));
rightPadM= fliplr(m(end-hw+1:end));

xp = [leftPadX, x0, rightPadX];
mp = [leftPadM, m,  rightPadM];

yp_num = conv(xp, k, 'same');
yp_den = conv(mp, k, 'same');

core_num = yp_num(hw+1 : hw+numel(x));
core_den = yp_den(hw+1 : hw+numel(x));

y = core_num ./ core_den;
y(core_den < eps) = NaN;
end

function savefig_smart(figH, fname)
try
    exportgraphics(figH, fname, 'Resolution', 300);
catch
    [p, n, ~] = fileparts(fname);
    if ~exist(p,'dir'), mkdir(p); end
    saveas(figH, fullfile(p, [n, '.png']));
end
end

function plot_box_groups(ax, X, Y, halfWidth)
% (kept for spread fallback; not used in PSTH C-only figure)
hold(ax, 'on');
groups = unique(X(:)');
for g = groups
    yi = Y(X==g); yi = yi(isfinite(yi));
    if numel(yi) < 2
        scatter(ax, g, yi, 12, 'filled', 'MarkerFaceAlpha',0.8);
        continue;
    end
    q  = quantile(yi, [0.25, 0.5, 0.75]);
    iqr= q(3) - q(1);
    wL = q(1) - 1.5 .* iqr;
    wU = q(3) + 1.5 .* iqr;
    yi_sorted  = sort(yi);
    whisk_low  = yi_sorted(find(yi_sorted >= wL, 1, 'first'));  if isempty(whisk_low),  whisk_low  = min(yi_sorted); end
    whisk_high = yi_sorted(find(yi_sorted <= wU, 1, 'last'));   if isempty(whisk_high), whisk_high = max(yi_sorted); end
    patch(ax, [g-halfWidth, g+halfWidth, g+halfWidth, g-halfWidth], [q(1), q(1), q(3), q(3)], ...
          [0.3, 0.5, 0.85], 'FaceAlpha',0.35, 'EdgeColor','k');
    plot(ax, [g-halfWidth, g+halfWidth], [q(2), q(2)], 'k-', 'LineWidth',1.4);
    plot(ax, [g, g], [q(3), whisk_high], 'k-');
    plot(ax, [g, g], [q(1), whisk_low ], 'k-');
    plot(ax, [g-0.15, g+0.15], [whisk_high, whisk_high], 'k-');
    plot(ax, [g-0.15, g+0.15], [whisk_low,  whisk_low ], 'k-');
end
hold(ax, 'off');
end

