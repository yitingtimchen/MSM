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
PLOT_PSTH_OVERLAYS    = 1;  % PSTH grid (per event set & cond-set), shows columns A, B, C

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
    {'OT AI', 'OT Replay', 'OT Decoy', 'OT Live'}; ...
    {'Saline AI', 'Saline Replay', 'Saline Decoy', 'Saline Live'} ...
};

% >>> Multiple event variant sets (edit here) <<<
event_variants = { ...
    "m1rew1on", ...
    "m2rew1on"
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

        %% -------- PSTH overlays: show [A] [B] and [C] side-by-side (per event & cond set) ---
        if PLOT_PSTH_OVERLAYS
            rows = numel(norm_methods);
            W_in = 15;  H_in = 2.5 .* rows + 2;  dpi = 300;

            figP = figure('Name', sprintf('PSTH Overlays — %s — %s — A & C only', ev_tag, cond_tag), ...
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

                    [muA, semA, muB, semB, Nused, nfiles, NpbB, sdB, AllSpikes] = compute_psth_for_pack( ...
                        packPath, t, baseMask, event_variant, TR, CH, norm_method, ...
                        scale_to_hz, smooth_sigma_sec, ...
                        use_robust_sigma, sigma_shrink_lambda, sigma_ridge, ...
                        trial_sigma_min_samples, sigma_floor, AllSpikes, conds, ci);

                    stats(end+1) = struct('cond',cond, 'muA',muA, 'semA',semA, ...
                                          'muB',muB, 'semB',semB, 'N',Nused, 'NpbB', NpbB,'sdB', sdB, 'files',nfiles); %#ok<SAGROW>
                end
                stats = stats(2:end);
                if isempty(stats)
                    warning('No usable conditions for norm=%s in %s / %s', norm_method, ev_tag, cond_tag);
                    continue;
                end

                cols = lines(numel(stats));

                % ---------- Column 1: [A] Baseline-sub (per-condition) ----------
                axA = nexttile(tlP, (nm-1).*3 + 1); hold(axA, 'on');
                for i = 1:numel(stats)
                    if show_sem_shading
                        shade_mean_sem(t, stats(i).muA, stats(i).semA, cols(i,:).*0.6 + 0.4);
                    end
                    plot(axA, t, stats(i).muA, 'Color', cols(i,:), 'LineWidth', 1.8, ...
                         'DisplayName', sprintf('%s (N=%d, files=%d)', stats(i).cond, stats(i).N, stats(i).files));
                end
                if DRAW_CHOICE_BOX && string(event_variant) == "m1rew1on"
                    add_xbox(axA, CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA);
                end
                xline(axA, 0, '--', 'Color', [0.2, 0.2, 0.2], 'HandleVisibility','off');
                xlabel(axA, 'Time from event (s)');
                ylabel(axA, ylabel_for(norm_method, scale_to_hz));
                title (axA, sprintf('[%dA] Baseline-sub\nnorm = %s', nm, norm_method), 'Interpreter','none');
                grid(axA, 'on'); xlim(axA, [t(1), t(end)]);
                set(axA,'FontSize',FZ_AX); axA.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
                legend(axA, 'Location','best');

                % ---------- Column 2: [B] Baseline-sub (z-scored) (per-condition) ----------
                axB = nexttile(tlP, (nm-1).*3 + 2); hold(axB, 'on');
                for i = 1:numel(stats)
                    if show_sem_shading
                        shade_mean_sem(t, stats(i).muB, stats(i).semB, cols(i,:).*0.6 + 0.4);
                    end
                    plot(axB, t, stats(i).muB, 'Color', cols(i,:), 'LineWidth', 1.8, ...
                         'DisplayName', sprintf('%s (N=%d, files=%d)', stats(i).cond, stats(i).N, stats(i).files));
                end
                if DRAW_CHOICE_BOX && string(event_variant) == "m1rew1on"
                    add_xbox(axB, CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA);
                end
                xline(axB, 0, '--', 'Color', [0.2, 0.2, 0.2], 'HandleVisibility','off');
                xlabel(axB, 'Time from event (s)');
                ylabel(axB, ylabel_for('zscore_neuron', scale_to_hz));
                title (axB, sprintf('[%dB] Baseline-sub z-scored per condition\nnorm = %s', nm, norm_method), 'Interpreter','none');
                grid(axB, 'on'); xlim(axB, [t(1), t(end)]);
                set(axB,'FontSize',FZ_AX); axB.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
                legend(axB, 'Location','best');

                % ---------- Column 3: [C] pooled re-z across conditions ----------
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
                if DRAW_CHOICE_BOX && string(event_variant) == "m1rew1on"
                    add_xbox(axC, CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA); 
                end
                xline(axC, 0, '--', 'Color', [0.2, 0.2, 0.2], 'HandleVisibility','off');
                xlabel(axC, 'Time from event (s)');
                ylabel(axC, 'Pooled z (B baseline-sub, all conds)');
                title (axC, sprintf('[%dC] Pooled re-z across conditions\nnorm = %s', nm, norm_method), 'Interpreter','none');
                grid(axC, 'on'); xlim(axC, [t(1), t(end)]); set(axC,'FontSize',FZ_AX); axC.TitleFontSizeMultiplier = FZ_TITLE/FZ_AX;
                legend(axC, 'Location','best');
            end

            sg = sgtitle(tlP, sprintf('Population PSTHs — %s — %s', ev_tag, cond_tag), 'Interpreter','none');
            set(sg,'FontSize',FZ_SGT);
            if ~exist(outDir,'dir'); mkdir(outDir); end
            outPNG = fullfile(outDir, sprintf('overlay_PSTH_A_and_C_%s__CONDSET_%s.png', ev_tag, cond_tag));
            print(figP, outPNG, '-dpng', sprintf('-r%d', dpi));
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

function add_xbox(ax, xwin, color, alpha)
yl = ylim(ax);
h = patch(ax, [xwin(1), xwin(2), xwin(2), xwin(1)], [yl(1), yl(1), yl(2), yl(2)], ...
          color, 'EdgeColor','none', 'FaceAlpha', alpha, 'HandleVisibility','off');
uistack(h, 'bottom');
end
