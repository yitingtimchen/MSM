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
% =========================================================================

close all;
clc;
clear;

%% ----------------------------- TOGGLES -----------------------------------
PLOT_BASELINE_SPREADS = 0;  % combined figure across conditions (per event set)
PLOT_PSTH_OVERLAYS    = 1;  % PSTH grid (per event set), shows ONLY column C
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
FIG_SIZE_SPREADS = [80, 80, 1850, 850]; % [x y width height] for spread plots (1×2 grid)
FZ_AX    = 14; % axis tick/label size
FZ_TITLE = 14; % panel title size
FZ_SGT   = 16; % super-title size
FZ_LEGEND= 12;

%% ----------------------------- Inputs ------------------------------------
if is_mac
    packDir = '/Users/yitingchen/Desktop/Lab/Stock_market/Tim/condition_packs';
else
    packDir = 'C:\Users\plattlab\Tim\Stock_market_experiment\Tim\condition_packs';
end

conds = {'AI','Replay','Decoy','Live'}; % conditions to compare (order used for segments)

% >>> Multiple event variant sets (edit here) <<<
event_variant_sets = { ...
    ["m1rew1on"] ...
};

% Analysis windows (seconds, relative to event)
win          = [-0.5, 1];   % master window for extraction (must cover both sub-windows)
BASELINE_WIN = [-0.5, -0.4];  % baseline window for re-baseline step (and spread)
POST_WIN     = [ 0.0,  1.5];  % post-event window (used in spread analysis only)

% PSTH normalization strategies (rows)
% norm_methods = {'none','baseline_sub_trial','baseline_sub_neuron','baseline_sub_chunk','zscore_neuron'};
norm_methods = {'baseline_sub_trial'};
scale_to_hz        = true;   % convert 1 ms binary → Hz (divide by dt)
smooth_sigma_sec   = 0.02;   % Gaussian sigma (sec); 0 disables smoothing (PRE-analysis everywhere)
show_sem_shading   = true;   % shaded SEM per condition in PSTH overlays

outDir = fullfile(packDir, 'overlay_outputs_multiEvent');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Task structure
TR = 15; % trials per block
CH = 6;  % number of blocks (chunks)

% --------- Z-score stability controls (used inside PSTH function for B) ---
use_robust_sigma        = false; % if true, use MAD*1.4826; else std
sigma_shrink_lambda     = 0.0;   % keep 0 to preserve scope differences cleanly
sigma_ridge             = 1e-6;  % small ridge added to variance (Hz^2)
trial_sigma_min_samples = 2;     % fallback if too few baseline bins (baseline stats)
sigma_floor             = 0.5;   % σ floor (Hz or arb. units after subtraction)

%% --------------------------- Prep time axis -------------------------------
% Pull dt from the first available pack on disk.
dt = [];
for ci = 1:numel(conds)
    fp = fullfile(packDir, sprintf('%s_condition_pack.mat', conds{ci}));
    if exist(fp, 'file')
        S = load(fp, 'C');
        dt = S.C.dt;
        break;
    end
end
if isempty(dt), error('No condition packs found in %s.', packDir); end

% Build time axis + masks
t        = (win(1):dt:win(2)); t = t(:)'; % enforce row
baseMask = (t >= BASELINE_WIN(1)) & (t <= BASELINE_WIN(2));
postMask = (t >= POST_WIN(1))     & (t <= POST_WIN(2));

%% ===================== LOOP OVER EVENT VARIANT SETS =======================
for es = 1:numel(event_variant_sets)
    event_variants = event_variant_sets{es};
    ev_tag = event_set_label(event_variants); % short label for titles/files

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
                packPath, t, baseMask, postMask, event_variants, TR, CH, scale_to_hz, smooth_sigma_sec, ...
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
        figSP = figure('Name', sprintf('Baseline Spread %% — %s — All Conditions', ev_tag), ...
                       'Color','w', 'Units','pixels', 'Position', FIG_SIZE_SPREADS);
        tl = tiledlayout(figSP, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    
        infoStr = sprintf('BL %d..%d ms | POST %d..%d ms | \\sigma=%d ms', ...
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
        title(axL, sprintf('[Left] Block-level spread%% per neuron — %s', infoStr), 'Interpreter','tex');
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
        title(axR, sprintf('[Right] Trial-level spread%% per neuron — %s', infoStr), 'Interpreter','tex');
        if ~isempty(seg_edges), draw_cond_segments(axR, seg_edges, seg_labels); end
        grid(axR, 'on'); set(axR,'FontSize',FZ_AX); axR.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
        if usedBoxchart, legend(axR, conds, 'Location','bestoutside'); end
    
        sg = sgtitle(tl, sprintf('Baseline Spread as %% of Post-Event FR — %s — All Conditions', ev_tag), 'Interpreter','none');
        set(sg,'FontSize',FZ_SGT);
    
        % Save
        savefig_smart(figSP, fullfile(outDir, sprintf('baseline_spread_percent_ALLCONDS_%s_gated.png', ev_tag)));
    end

    %% -------- PSTH overlays: show [A] [B] and [C] side-by-side (per event set) ---
    if PLOT_PSTH_OVERLAYS
        rows = numel(norm_methods);
        W_in = 25;  H_in = 4.2 .* rows + 2;  dpi = 300;

        figP = figure('Name', sprintf('PSTH Overlays — %s — A & C only', ev_tag), ...
            'Color','w', 'Units','inches', 'Position',[0.5, 0.5, W_in, H_in], ...
            'PaperUnits','inches', 'PaperSize',[W_in, H_in], ...
            'PaperPosition',[0, 0, W_in, H_in], 'PaperPositionMode','manual', ...
            'Resize','off', 'WindowState','normal', 'WindowStyle','normal', ...
            'Renderer','painters', 'GraphicsSmoothing','on');
        tlP = tiledlayout(figP, rows, 3, 'TileSpacing','loose', 'Padding','loose');

        for nm = 1:rows
            norm_method = norm_methods{nm};

            % Collect per-condition PSTH stats under this normalization
            stats = struct('cond',[],'muA',[],'semA',[],'muB',[],'semB',[],'N',[],'NpbB',[],'sdB',[],'files',[]);
            for ci = 1:numel(conds)
                cond     = conds{ci};
                packPath = fullfile(packDir, sprintf('%s_condition_pack.mat', cond));
                if ~exist(packPath,'file')
                    warning('Pack not found: %s (skipping)', packPath);
                    continue;
                end

                [muA, semA, muB, semB, Nused, nfiles, NpbB, sdB] = compute_psth_for_pack( ...
                    packPath, t, baseMask, event_variants, TR, CH, norm_method, ...
                    scale_to_hz, smooth_sigma_sec, ...
                    use_robust_sigma, sigma_shrink_lambda, sigma_ridge, ...
                    trial_sigma_min_samples, sigma_floor);

                stats(end+1) = struct('cond',cond, 'muA',muA, 'semA',semA, ...
                                      'muB',muB, 'semB',semB, 'N',Nused, 'NpbB', NpbB,'sdB', sdB, 'files',nfiles); %#ok<SAGROW>
            end
            stats = stats(2:end);
            if isempty(stats)
                warning('No usable conditions for norm=%s in %s', norm_method, ev_tag);
                continue;
            end

            cols = lines(numel(stats));

            % ---------- Column 1: [A] Baseline-sub (per-condition) ----------
            axA = nexttile(tlP, (nm-1).*2 + 1); hold(axA, 'on');
            for i = 1:numel(stats)
                if show_sem_shading
                    shade_mean_sem(t, stats(i).muA, stats(i).semA, cols(i,:).*0.6 + 0.4);
                end
                plot(axA, t, stats(i).muA, 'Color', cols(i,:), 'LineWidth', 1.8, ...
                     'DisplayName', sprintf('%s (N=%d, files=%d)', stats(i).cond, stats(i).N, stats(i).files));
            end
            if DRAW_CHOICE_BOX
                add_xbox(axA, CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA);
            end
            xline(axA, 0, '--', 'Color', [0.2, 0.2, 0.2], 'HandleVisibility','off');
            xlabel(axA, 'Time from event (s)');
            ylabel(axA, ylabel_for(norm_method, scale_to_hz));
            title (axA, sprintf('[%dA] Baseline-sub | norm=%s | %s', nm, norm_method, ev_tag), 'Interpreter','none');
            grid(axA, 'on'); xlim(axA, [t(1), t(end)]);
            set(axA,'FontSize',FZ_AX); axA.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
            legend(axA, 'Location','best');

            % ---------- Column 2: [B] Baseline-sub (z-scored) (per-condition) ----------
            axB = nexttile(tlP, (nm-1).*2 + 2); hold(axB, 'on');
            for i = 1:numel(stats)
                if show_sem_shading
                    shade_mean_sem(t, stats(i).muB, stats(i).semB, cols(i,:).*0.6 + 0.4);
                end
                plot(axB, t, stats(i).muB, 'Color', cols(i,:), 'LineWidth', 1.8, ...
                     'DisplayName', sprintf('%s (N=%d, files=%d)', stats(i).cond, stats(i).N, stats(i).files));
            end
            if DRAW_CHOICE_BOX
                add_xbox(axB, CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA);
            end
            xline(axB, 0, '--', 'Color', [0.2, 0.2, 0.2], 'HandleVisibility','off');
            xlabel(axB, 'Time from event (s)');
            ylabel(axB, ylabel_for('zscore_neuron', scale_to_hz));
            title (axB, sprintf('[%dB] Baseline-sub z-scored per condition | norm=%s | %s', nm, norm_method, ev_tag), 'Interpreter','none');
            grid(axB, 'on'); xlim(axB, [t(1), t(end)]);
            set(axB,'FontSize',FZ_AX); axB.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
            legend(axB, 'Location','best');

            % ---------- Column 3: [C] pooled re-z across all conditions ----------
            % 1) Re-baseline each condition's z-scored B by BASELINE_WIN
            muC0 = cell(1,numel(stats)); % baseline-sub'd B (per condition)
            seC0 = cell(1,numel(stats)); % SEM carried along (will be scaled by pooled σ)
            wC = cell(1,numel(stats)); % per-bin weights (e.g., N per bin)
            for i = 1:numel(stats)
                muB_i = stats(i).muB;
                seB_i = stats(i).semB;
                muC0{i} = muB_i - mean(muB_i(baseMask), 'omitnan'); % baseline subtraction on B
                seC0{i} = seB_i; % will rescale after pooled σ
                wC{i} = stats(i).NpbB; % use per-bin counts as weights
            end
            
            % 2) Pool ALL conditions and ALL time bins to get global μ/σ (weights = per-bin counts)
            all_vals = []; all_w = [];
            for i = 1:numel(stats)
                msk = isfinite(muC0{i}) & isfinite(wC{i}) & (wC{i} > 0);
                if any(msk)
                    all_vals = [all_vals, muC0{i}(msk)]; %#ok<AGROW>
                    all_w = [all_w, wC{i}(msk) ]; %#ok<AGROW>
                end
            end
            if isempty(all_vals)
                warning('No finite values for pooled z in Column C. Skipping.');
                continue;
            end
            all_w = all_w / max(eps, sum(all_w)); % normalize weights
            global_mu = sum(all_w .* all_vals);
            global_sd = sqrt( sum(all_w .* (all_vals - global_mu).^2) );
            if ~isfinite(global_sd) || global_sd < sigma_floor, global_sd = sigma_floor; end
            
            % 3) Re-z each condition using pooled (global) μ/σ
            axC = nexttile(tlP, (nm-1).*3 + 3); hold(axC,'on');
            cols = lines(numel(stats)); % reuse condition colors
            for i = 1:numel(stats)
                muC = (muC0{i} - global_mu) ./ global_sd;
                seC = seC0{i} ./ global_sd;
                if show_sem_shading
                    shade_mean_sem(t, muC, seC, cols(i,:).*0.6 + 0.4);
                end
                plot(axC, t, muC, 'Color', cols(i,:), 'LineWidth', 1.8, ...
                'DisplayName', sprintf('%s (N=%d, files=%d)', stats(i).cond, stats(i).N, stats(i).files));
            end
            if DRAW_CHOICE_BOX, add_xbox(axC, CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA); end
            xline(axC, 0, '--', 'Color', [0.2, 0.2, 0.2], 'HandleVisibility','off');
            xlabel(axC, 'Time from event (s)');
            ylabel(axC, 'Pooled z (B baseline-sub, all conds)');
            title (axC, sprintf('[%dC] Pooled re-z across conditions | norm=%s | %s', nm, norm_method, ev_tag), 'Interpreter','none');
            grid(axC, 'on'); xlim(axC, [t(1), t(end)]); set(axC,'FontSize',FZ_AX); axC.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
            legend(axC, 'Location','best');
        end

        sg = sgtitle(tlP, sprintf('Population PSTHs — %s — Columns: A (baseline-sub) & C (cross-cond z)', ev_tag), 'Interpreter','none');
        set(sg,'FontSize',FZ_SGT);
        if ~exist(outDir,'dir'); mkdir(outDir); end
        outPNG = fullfile(outDir, sprintf('overlay_PSTH_A_and_C_%s.png', ev_tag));
        print(figP, outPNG, '-dpng', sprintf('-r%d', dpi));
    end

end

%% ============================== FUNCTIONS ================================
function label = event_set_label(evset)
% Build a compact label from an event variant set (string array).
% Example: ["m1rew1ON","m1rew1On","m1rew1on"] -> "m1rew1"
evset = string(evset);
if isempty(evset)
    label = 'events';
    return;
end
s = char(evset(1));
digits = regexp(s, '^[^A-Za-z]*([A-Za-z0-9]+)', 'tokens', 'once');
if ~isempty(digits)
    label = digits{1};
else
    label = s;
end
label = regexprep(label, '(ON|On|on)$', '');
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
    packPath, t, baseMask, postMask, event_variants, TR, CH, ...
    scale_to_hz, smooth_sigma_sec, ...
    use_robust_sigma, trial_sigma_min_samples, RATE_FLOOR_HZ, DELTA_FLOOR_HZ, K_SIGMA, ...
    MIN_TRIALS_PER_BLOCK, MIN_BLOCKS_WITH_DATA)
% (unchanged from your gated version; omitted here for brevity)
% Return values are used only if PLOT_BASELINE_SPREADS==1.
S = load(packPath,'C');
C = S.C;
dt = C.dt;
wRelBins = round(t ./ dt);
block_pct_by_neuron = {};
trial_pct_by_neuron = {};
for f = 1:numel(C.files)
    if isempty(C.eventTables{f}), continue; end
    props = string(C.eventTables{f}.Properties.VariableNames);
    pick  = find(ismember(event_variants, props), 1, 'first');
    if isempty(pick), continue; end
    varName = char(event_variants(pick));
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

function [muA, semA, muB, semB, Nused, nfiles, NpbB, sdB] = compute_psth_for_pack( ...
    packPath, t, baseMask, event_variants, TR, CH, norm_method, ...
    scale_to_hz, smooth_sigma_sec, ...
    use_robust_sigma, sigma_shrink_lambda, sigma_ridge, ...
    trial_sigma_min_samples, sigma_floor)
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
all_A = [];
all_B = [];
nfiles = 0;

for f = 1:numel(C.files)
% for f = 1 % test
    if isempty(C.eventTables{f}), continue; end
    props = string(C.eventTables{f}.Properties.VariableNames);
    pick  = find(ismember(event_variants, props), 1, 'first');
    if isempty(pick), continue; end
    varName = char(event_variants(pick));

    idxMat = C.eventTables{f}.(varName){1};
    if ~ismatrix(idxMat), continue; end

    Sbin = C.S{f};
    % Sbin = Sbin(1, :); % test
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
        % option: take neuron mean first
        psth_A = mean(Tmat_A, 1, 'omitnan');
        psth_B = mean(Tmat_B, 1, 'omitnan');
        all_A(end+1,:) = psth_A; %#ok<AGROW>
        all_B(end+1,:) = psth_B; %#ok<AGROW>

        % % b option: concatenate each trial
        % all_A = [all_A; Tmat_A]; %#ok<AGROW>
        % all_B = [all_B; Tmat_B]; %#ok<AGROW>
    end
end

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

function y = ylabel_for(norm_method, scale_to_hz)
nm = char(norm_tag(norm_method));
switch nm
    case 'zscore_neuron', y = 'Z-scored (baseline)';
    case {'baseline_sub_trial','baseline_sub_neuron','baseline_sub_chunk'}
        if scale_to_hz, y = '\Delta rate from baseline (Hz)'; else, y = '\Delta (prob.)'; end
    otherwise
        if scale_to_hz, y = 'Firing rate (Hz)'; else, y = '1 ms spike prob.'; end
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

function shade_mean_sem(x, mu, sem, rgb)
x = x(:)'; mu = mu(:)'; sem = sem(:)';
upper = mu + sem; lower = mu - sem;
if numel(x) == numel(mu) && all(isfinite(upper)) && all(isfinite(lower))
    Xpoly = [x, fliplr(x)];
    Ypoly = [upper, fliplr(lower)];
    hFill = fill(Xpoly, Ypoly, rgb, 'EdgeColor','none', 'FaceAlpha',0.18);
    set(hFill, 'HandleVisibility','off');
    if isprop(hFill,'Annotation')
        hFill.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
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

function add_xbox(ax, xwin, color, alpha)
yl = ylim(ax);
h = patch(ax, [xwin(1), xwin(2), xwin(2), xwin(1)], [yl(1), yl(1), yl(2), yl(2)], ...
          color, 'EdgeColor','none', 'FaceAlpha', alpha, 'HandleVisibility','off');
uistack(h, 'bottom');
end
