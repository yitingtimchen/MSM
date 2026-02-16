% =========================================================================
% Unified PSTH Overlays — Event × Condition-Set × Condition, over Choices
% Columns: [A] baseline-sub, [B] z-scored (full-window), [C] pooled re-z across choices
%
% Based on:
% - plot_m1_reward_onset_firing_rate_across_cond_250903.m

% - plot_m1_reward_onset_firing_rate_across_choice_250903.m

% =========================================================================
close all; clc; clear;

%% ----------------------------- TOGGLES -----------------------------------
PLOT_PSTH_OVERLAYS = 0; % generate A/B/C overlays per (event × cond-set × condition)
DO_VENN_HIST = 0; % toggle for venn+hist figure per condition set
PLOT_PSTH_OVERLAP = 1; % set 1 to make the second PSTH version for the Venn overlap

VENN_DIFF_THRESH = 10; % threshold for ">" (e.g., 0 for strictly greater)
HIST_NUM_BINS = 100; % number of bins (symmetric around 0)
SAVE_VENN_HIST = 1; % save figures

SC_WIN = [0, 1.5]; % seconds post-reward (500–1500 ms)
DRAW_CHOICE_BOX = 1; % highlight c hoice window on PSTH
CHOICE_BOX_WIN = [-0.25, -0.05]; % seconds (start, end)
CHOICE_BOX_COLOR = [1.00, 0.92, 0.50]; % light yellow
CHOICE_BOX_ALPHA = 0.18; % transparency

% Figure sizing / typography
FZ_AX = 8; % axis tick/label size
FZ_TITLE = 8; % panel title size
FZ_SGT = 10; % super-title size

% ----------------------------- Inputs -------------------------------------
is_mac = 0;
if is_mac
packDir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs';
else
packDir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs';
end

% === CONDITION SETS =======================================================
% cond_sets = { ...
% {'AI','Replay','Decoy','Live'}; ...
% {'OT AI', 'OT Replay', 'OT Decoy', 'OT Live'}; ...
% {'Saline AI', 'Saline Replay', 'Saline Decoy', 'Saline Live'} ...
% };
cond_sets = { {'AI','Replay','Decoy','Live'} }; % test

% === Multiple EVENT VARIANTS ==============================================
% event_variants = { "m1rew1on", "m2rew1on" };
event_variants = { "m1rew1on" }; % test

win          = [-0.5, 1.5];   % master window for extraction (must cover both sub-windows)
BASELINE_WIN = [-0.5, -0.4]; % baseline window for re-baseline step
% POST_WIN = [ 0.0, 1.5]; % reserved if you add spread analysis

% PSTH normalization strategies (rows)
norm_methods = {'baseline_sub_trial'};
scale_to_hz = true; % convert 1 ms binary → Hz (divide by dt)
smooth_sigma_sec = 0.05; % Gaussian sigma (sec); 0 disables smoothing (PRE-analysis everywhere)
show_sem_shading = true; % shaded SEM per choice in PSTH overlays

% Base output directory (per event × condition-set × condition)
BASE_OUT_DIR = fullfile(packDir, 'overlay_outputs_multiEvent_choiceUnified');
if ~exist(BASE_OUT_DIR, 'dir'), mkdir(BASE_OUT_DIR); end

% Task structure
TR = 15; % trials per block
CH = 6; % number of blocks (chunks)

% --------- Z-score stability controls (used inside PSTH function for B) ---
use_robust_sigma = false; % if true, use MAD*1.4826; else std
sigma_shrink_lambda = 0.0; % keep 0 to preserve scope differences cleanly
sigma_ridge = 1e-6; % small ridge added to variance (Hz^2)
trial_sigma_min_samples = 2; % fallback if too few baseline bins (baseline stats)
sigma_floor = 0.5; % σ floor (Hz or arb. units after subtraction)

% Choice labels/colors (Buy/Hold/Sell)
CHOICE_LABELS = {'Buy (1)','Hold (2)','Sell (3)'};
CHOICE_VALS = [1, 2, 3];
CHOICE_COLS = [0, 0.6, 0; % Buy - green
0, 0.4, 1; % Hold - blue
0.9, 0.2, 0.2]; % Sell - red

%% --------------------------- Prep time axis -------------------------------
% Pull dt from the first available pack across ANY condition in ANY set.
dt = [];
all_conds_flat = unique([cond_sets{:}]); % flatten cell-of-cell to 1×N cellstr
for ci = 1:numel(all_conds_flat)
fp = fullfile(packDir, sprintf('%s_condition_pack.mat', all_conds_flat{ci}));
if exist(fp, 'file')
S = load(fp, 'C'); dt = S.C.dt; break;
end
end
if isempty(dt), error('No condition packs found in %s (checked across all cond_sets).', packDir); end

% Build time axis + masks
t = (win(1):dt:win(2)); t = t(:)'; % enforce row
baseMask = (t >= BASELINE_WIN(1)) & (t <= BASELINE_WIN(2));

%% ===================== LOOP OVER EVENT & CONDITION SETS ===================
% ================= FIGURE: 4 conditions as columns × 7 methods (A–G) as rows =================
% Assumes all variables and helper functions from your script are in scope.
if PLOT_PSTH_OVERLAYS
    for es = 1:numel(event_variants)
    event_variant = event_variants{es};
    ev_tag = event_variant;
    
    for cs = 1:numel(cond_sets)
        conds = cond_sets{cs};
        condset_tag = cond_set_label(conds); % e.g., "AI-Replay-Decoy-Live"
    
        % ---------- Compute stats for each condition in this set ----------
        condStats = struct('cond',[],'stats',[]);
        for ci = 1:numel(conds)
            cond_name = conds{ci};
            packPath = fullfile(packDir, sprintf('%s_condition_pack.mat', cond_name));
            if ~exist(packPath, 'file')
                warning('Pack not found: %s (skipping)', packPath);
                continue;
            end
    
            [muRaw_all, semRaw_all, muA_all, semA_all, muB_all, semB_all, ...
             muArowz_all, semArowz_all, muBrowz_all, semBrowz_all, ...
             muAmatz_all, semAmatz_all, muBmatz_all, semBmatz_all, ...
             N_all, NpbB_all, sdB_all] = ...
                compute_psth_choice_for_pack( ...
                    packPath, t, baseMask, event_variant, TR, CH, norm_methods{1}, ...
                    scale_to_hz, smooth_sigma_sec, ...
                    use_robust_sigma, sigma_shrink_lambda, sigma_ridge, ...
                    trial_sigma_min_samples, sigma_floor, CHOICE_VALS, conds, ci, ev_tag);
    
            % Pack per-choice stats
            stats = struct('choice',[],'muRaw',[],'semRaw',[], ...
                           'muA',[],'semA',[],'muB',[],'semB',[], ...
                           'muArowz',[],'semArowz',[],'muBrowz',[],'semBrowz',[], ...
                           'muAmatz',[],'semAmatz',[],'muBmatz',[],'semBmatz',[], ...
                           'N',[],'NpbB',[],'sdB',[],'files',[]);
            for k = 1:numel(CHOICE_VALS)
                stats(k).choice   = CHOICE_VALS(k);
                stats(k).muRaw    = muRaw_all{k};   stats(k).semRaw    = semRaw_all{k};
                stats(k).muA      = muA_all{k};     stats(k).semA      = semA_all{k};
                stats(k).muB      = muB_all{k};     stats(k).semB      = semB_all{k};
                stats(k).muArowz  = muArowz_all{k}; stats(k).semArowz  = semArowz_all{k};
                stats(k).muBrowz  = muBrowz_all{k}; stats(k).semBrowz  = semBrowz_all{k};
                stats(k).muAmatz  = muAmatz_all{k}; stats(k).semAmatz  = semAmatz_all{k};
                stats(k).muBmatz  = muBmatz_all{k}; stats(k).semBmatz  = semBmatz_all{k};
                stats(k).N        = N_all(k);
                stats(k).NpbB     = NpbB_all{k};
                stats(k).sdB      = sdB_all{k};
                stats(k).files    = 1; % not meaningful across combined input; adjust if you track per-cond files
            end
            condStats(ci).cond  = cond_name;
            condStats(ci).stats = stats;
        end
    
        % Remove empties (in case of missing packs)
        condStats = condStats(~arrayfun(@isempty, {condStats.cond}));
        if isempty(condStats), continue; end
    
        % ---------- Build figure: rows = 7 methods (A–G), cols = conditions ----------
        rows = 7; cols = numel(condStats);
        W_in = 1.5 * cols;  H_in = 1.3 * rows + 1.2;  dpi = 450;
    
        figG = figure('Name', sprintf('PSTH Grid — %s — %s — A–G × Conditions', char(ev_tag), condset_tag), ...
            'Color','w','Units','inches','Position',[0.5,0.5,W_in,H_in], ...
            'PaperUnits','inches','PaperSize',[W_in,H_in], ...
            'PaperPosition',[0,0,W_in,H_in],'PaperPositionMode','manual', ...
            'Resize','off','Renderer','painters','GraphicsSmoothing','on');
    
        tlG = tiledlayout(figG, rows, cols, 'TileSpacing','loose', 'Padding','compact');
    
        % Row labels and y-axis labels
        rowTitles = { ...
            "Raw Firing Rate", ...
            "Base-sub trial", ...
            "Base-sub neuron", ...
            ["Base-sub trial"; "z trial"], ...
            ["Base-sub neuron"; "z trial"], ...
            ["Base-sub trial"; "z neuron"], ...
            ["Base-sub neuron"; "z neuron"]};
        rowYlabs = { ...
            "Firing rate (Hz)", ...
            "Δ from baseline (Hz)", ...
            "Δ from baseline (Hz)", ...
            "z-score", ...
            "z-score", ...
            "z-score", ...
            "z-score"};
    
        AX = gobjects(rows, cols); % place before the r/c loops
        for r = 1:rows
            for c = 1:cols
                ax = nexttile(tlG, (r-1)*cols + c); hold(ax, 'on');
                AX(r,c) = ax;
                st = condStats(c).stats;
    
                % plot choices
                for i = 1:numel(st)
                    switch r
                        case 1, mu = st(i).muRaw;    se = st(i).semRaw;
                        case 2, mu = st(i).muA;      se = st(i).semA;
                        case 3, mu = st(i).muB;      se = st(i).semB;
                        case 4, mu = st(i).muArowz;  se = st(i).semArowz;
                        case 5, mu = st(i).muBrowz;  se = st(i).semBrowz;
                        case 6, mu = st(i).muAmatz;  se = st(i).semAmatz;
                        case 7, mu = st(i).muBmatz;  se = st(i).semBmatz;
                    end
                    if isempty(mu), continue; end
                    if show_sem_shading
                        shade_mean_sem(t, mu, se, CHOICE_COLS(i,:).*0.6 + 0.4);
                    end
                    plot(ax, t, mu, 'Color', CHOICE_COLS(i,:), 'LineWidth', 1, ...
                         'DisplayName', sprintf('%s (N=%d)', CHOICE_LABELS{i}, st(i).N));
                end
                xline(ax, 0, '--', 'Color', [0.2,0.2,0.2], 'HandleVisibility','off');
                grid(ax, 'on'); xlim(ax, [t(1), t(end)]);
                set(ax,'FontSize',FZ_AX);
    
                % Axes labels/titles
                if c == 1
                    ylabel(ax, rowYlabs{r}, 'Interpreter','none');
                end
                if r == 1
                    title(ax, sprintf('%s\n\n%s', condStats(c).cond, rowTitles{r}), 'Interpreter','none');
                else
                    title(ax, rowTitles{r}, 'Interpreter','none');
                end
                if r == rows
                    xlabel(ax, 'Time from event (s)');
                end
    
                % Only one legend (top-left tile)
                if r == 1
                    legend(ax, 'Location','southeast', 'Box','off','FontSize',6,'IconColumnWidth',4);
                else
                    legend(ax,'off');
                end
            end
        end
    
        for r = 1:rows
            ra = AX(r, isgraphics(AX(r,:))); % valid axes in this row
            if isempty(ra), continue; end
            yl = get(ra, 'YLim'); % cell array when multiple handles
            if iscell(yl), yl = cell2mat(yl); end % N×2
            lo = min(yl(:,1)); hi = max(yl(:,2)); % union of extents
            set(ra, 'YLim', [lo hi]); % apply to all axes in row
            linkaxes(ra, 'y'); % keep synced on zoom/pan
    
            for c = 1:cols
                if DRAW_CHOICE_BOX && string(event_variant) == "m1rew1on"
                    add_xbox(AX(r, c), CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA);
                end
            end
        end
    
        sg = sgtitle(tlG, sprintf('Population PSTHs — %s — %s (rows A–G; cols = conditions)', ...
                     char(ev_tag), condset_tag), 'Interpreter','none');
        set(sg,'FontSize',FZ_SGT);
    
        % Save
        outDir = fullfile(BASE_OUT_DIR, sprintf('%s_%s_grid_rowsAtoG_colsConds', char(ev_tag), condset_tag));
        if ~exist(outDir,'dir'); mkdir(outDir); end
        outPNG = fullfile(outDir, sprintf('PSTH_grid_rowsAtoG_colsConds__%s__%s.png', char(ev_tag), condset_tag));
        print(figG, outPNG, '-dpng', sprintf('-r%d', dpi));
    end
    
    
    end
end

%%
if DO_VENN_HIST
    for es = 1:numel(event_variants)
    event_variant = event_variants{es}; ev_tag = event_variant;
        for cs = 1:numel(cond_sets)
        conds = cond_sets{cs};
        condset_tag = cond_set_label(conds);
        outDir = fullfile(BASE_OUT_DIR, sprintf('%s_%s_venn_hist', char(ev_tag), condset_tag));
        if ~exist(outDir,'dir'), mkdir(outDir); end
        
            % Precompute diffs per condition, gather ranges for common histogram edges
            Cinfo = struct('name',[],'diffs',[],'N',[]);
            all_vals = [];
            for ci = 1:numel(conds)
                cond_name = conds{ci};
                packPath = fullfile(packDir, sprintf('%s_condition_pack.mat', cond_name));
                if ~exist(packPath,'file')
                    warning('Pack not found: %s (skipping)', packPath); continue;
                end
                [~, diffs, ~, info] = compute_spikecount_diffs_for_pack( ...
                    packPath, event_variant, TR, CH, [0.5,1.5], CHOICE_VALS, conds, ci, ev_tag);
                Cinfo(ci).name = cond_name;
                Cinfo(ci).diffs = diffs;           % columns: [Buy−Sell, Buy−Hold]
                Cinfo(ci).N = size(diffs,1);
                if ~isempty(diffs), all_vals = [all_vals; diffs(:)]; end %#ok<AGROW>
            end
    
            % Common histogram edges across conditions
            if isempty(all_vals) || ~any(isfinite(all_vals))
                warning('No spike-count diffs available for set %s.', condset_tag);
                continue;
            end
            maxabsv = max(abs(all_vals(isfinite(all_vals))));
            if maxabsv == 0, maxabsv = 1; end
            edges = linspace(-maxabsv, maxabsv, HIST_NUM_BINS+1);
    
            % Colors for differences (blend of choice colors)
            col_BH = mean([CHOICE_COLS(1,:); CHOICE_COLS(2,:)],1); % Buy vs Hold
            col_BS = mean([CHOICE_COLS(1,:); CHOICE_COLS(3,:)],1); % Buy vs Sell
    
            % Figure layout: rows=2 (top: Venn, bottom: Hist), cols=#conditions
            rows = 2; cols = numel(conds);
            W_in = 4 * cols;  H_in = 6;  dpi = 450;
            figVH = figure('Name', sprintf('Venn+Hist — %s — %s', char(ev_tag), condset_tag), ...
                'Color','w','Units','inches','Position',[0.5,0.5,W_in,H_in], ...
                'PaperUnits','inches','PaperSize',[W_in,H_in], 'Resize','off', ...
                'Renderer','painters','GraphicsSmoothing','on');
            tl = tiledlayout(figVH, rows, cols, 'TileSpacing','tight','Padding','tight');
    
            for ci = 1:numel(conds)
                if isempty(Cinfo(ci).diffs)
                    % Create empty placeholders to keep grid shape
                    nexttile(tl, ci); axis off; title(sprintf('%s\n(no data)', Cinfo(ci).name),'Interpreter','none');
                    nexttile(tl, cols+ci); axis off;
                    continue;
                end
                D = Cinfo(ci).diffs; % [Buy−Sell, Buy−Hold]
                BS = D(:,1); BH = D(:,2);
    
                % --- Top row: Venn (A: BUY>HOLD, B: BUY>SELL) ---
                A = BH > VENN_DIFF_THRESH;     % BUY>HOLD
                B = BS > VENN_DIFF_THRESH;     % BUY>SELL
                nA = sum(A); nB = sum(B); nAB = sum(A & B); N = numel(A);
    
                axV = nexttile(tl, ci); cla(axV); axis(axV,'off'); hold(axV,'on');
                draw_venn_two(axV, nA, nB, nAB, N, {'BUY>HOLD','BUY>SELL'});
                title(axV, sprintf('%s — Venn', Cinfo(ci).name), 'Interpreter','none');
    
                % --- Bottom row: Histogram overlay of diffs ---
                axH = nexttile(tl, cols+ci); hold(axH,'on'); box(axH,'on'); 
                histogram(axH, BH, 'BinEdges', edges, 'DisplayStyle','bar', 'FaceAlpha',0.45, 'EdgeColor','none', 'FaceColor', col_BH);
                histogram(axH, BS, 'BinEdges', edges, 'DisplayStyle','bar', 'FaceAlpha',0.45, 'EdgeColor','none', 'FaceColor', col_BS);
                ylims = ylim(axH); plot(axH, [0 0], ylims, '--', 'Color', [0.2 0.2 0.2]); ylim(axH, ylims);
                xlabel(axH, 'Spike-count difference'); ylabel(axH, 'Neurons');
                legend(axH, {'Buy−Hold','Buy−Sell'}, 'Location','northeast', 'Box','off');
                title(axH, sprintf('%s — Diff Distributions', Cinfo(ci).name), 'Interpreter','none');
            end
    
            sg = sgtitle(tl, sprintf('Condition Set: %s | Event: %s\n', condset_tag, char(ev_tag)), 'Interpreter','none', 'Fontsize', 14);

            if SAVE_VENN_HIST
                outPNG = fullfile(outDir, sprintf('VennHist__%s__%s.png', char(ev_tag), condset_tag));
                print(figVH, outPNG, '-dpng', sprintf('-r%d', dpi));
            end
        end
    end

end

%%
if PLOT_PSTH_OVERLAP
for es = 1:numel(event_variants)
event_variant = event_variants{es}; ev_tag = event_variant;
for cs = 1:numel(cond_sets)
conds = cond_sets{cs};
condset_tag = cond_set_label(conds);

        % ---------- Compute stats for each condition in this set (OVERLAP ONLY) ----------
        condStats = struct('cond',[],'stats',[]);
        for ci = 1:numel(conds)
            cond_name = conds{ci};
            packPath = fullfile(packDir, sprintf('%s_condition_pack.mat', cond_name));
            if ~exist(packPath, 'file'), warning('Pack not found: %s (skipping)', packPath); continue; end

            % overlap mask (logical per-neuron, AFTER standard keep/has_any filtering)
            overlap_mask = get_overlap_mask_for_pack( ...
                packPath, event_variant, TR, CH, SC_WIN, CHOICE_VALS, conds, ci, ev_tag, VENN_DIFF_THRESH);

            [muRaw_all, semRaw_all, muA_all, semA_all, muB_all, semB_all, ...
             muArowz_all, semArowz_all, muBrowz_all, semBrowz_all, ...
             muAmatz_all, semAmatz_all, muBmatz_all, semBmatz_all, ...
             N_all, NpbB_all, sdB_all] = ...
                compute_psth_choice_for_pack( ...
                    packPath, t, baseMask, event_variant, TR, CH, norm_methods{1}, ...
                    scale_to_hz, smooth_sigma_sec, ...
                    use_robust_sigma, sigma_shrink_lambda, sigma_ridge, ...
                    trial_sigma_min_samples, sigma_floor, CHOICE_VALS, conds, ci, ev_tag, overlap_mask);

            stats = struct('choice',[],'muRaw',[],'semRaw',[], ...
                           'muA',[],'semA',[],'muB',[],'semB',[], ...
                           'muArowz',[],'semArowz',[],'muBrowz',[],'semBrowz',[], ...
                           'muAmatz',[],'semAmatz',[],'muBmatz',[],'semBmatz',[], ...
                           'N',[],'NpbB',[],'sdB',[],'files',[]);
            for k = 1:numel(CHOICE_VALS)
                stats(k).choice   = CHOICE_VALS(k);
                stats(k).muRaw    = muRaw_all{k};   stats(k).semRaw    = semRaw_all{k};
                stats(k).muA      = muA_all{k};     stats(k).semA      = semA_all{k};
                stats(k).muB      = muB_all{k};     stats(k).semB      = semB_all{k};
                stats(k).muArowz  = muArowz_all{k}; stats(k).semArowz  = semArowz_all{k};
                stats(k).muBrowz  = muBrowz_all{k}; stats(k).semBrowz  = semBrowz_all{k};
                stats(k).muAmatz  = muAmatz_all{k}; stats(k).semAmatz  = semAmatz_all{k};
                stats(k).muBmatz  = muBmatz_all{k}; stats(k).semBmatz  = semBmatz_all{k};
                stats(k).N        = N_all(k);
                stats(k).NpbB     = NpbB_all{k};
                stats(k).sdB      = sdB_all{k};
                stats(k).files    = 1;
            end
            condStats(ci).cond  = sprintf('%s (overlap)', cond_name);
            condStats(ci).stats = stats;
        end

        % Remove empties
        condStats = condStats(~arrayfun(@isempty, {condStats.cond}));
        if isempty(condStats), continue; end

        % ---------- Plot grid (same style as your main PSTH figure) ----------
        rows = 7; cols = numel(condStats);
        W_in = 1.5 * cols;  H_in = 1.3 * rows + 1.2;  dpi = 450;
        figG = figure('Name', sprintf('PSTH Grid (OVERLAP) — %s — %s', char(ev_tag), condset_tag), ...
            'Color','w','Units','inches','Position',[0.5,0.5,W_in,H_in], ...
            'PaperUnits','inches','PaperSize',[W_in,H_in], 'Resize','off', ...
            'Renderer','painters','GraphicsSmoothing','on');
        tlG = tiledlayout(figG, rows, cols, 'TileSpacing','loose', 'Padding','compact');

        rowTitles = { ...
            "Raw Firing Rate", "Base-sub trial", "Base-sub neuron", ...
            ["Base-sub trial"; "z trial"], ["Base-sub neuron"; "z trial"], ...
            ["Base-sub trial"; "z neuron"], ["Base-sub neuron"; "z neuron"]};
        rowYlabs = {"Firing rate (Hz)","Δ from baseline (Hz)","Δ from baseline (Hz)", ...
                    "z-score","z-score","z-score","z-score"};

        AX = gobjects(rows, cols);
        for r = 1:rows
            for c = 1:cols
                ax = nexttile(tlG, (r-1)*cols + c); hold(ax,'on'); AX(r,c)=ax;
                st = condStats(c).stats;
                for i = 1:numel(st)
                    switch r
                        case 1, mu = st(i).muRaw;    se = st(i).semRaw;
                        case 2, mu = st(i).muA;      se = st(i).semA;
                        case 3, mu = st(i).muB;      se = st(i).semB;
                        case 4, mu = st(i).muArowz;  se = st(i).semArowz;
                        case 5, mu = st(i).muBrowz;  se = st(i).semBrowz;
                        case 6, mu = st(i).muAmatz;  se = st(i).semAmatz;
                        case 7, mu = st(i).muBmatz;  se = st(i).semBmatz;
                    end
                    if isempty(mu), continue; end
                    if show_sem_shading, shade_mean_sem(t, mu, se, CHOICE_COLS(i,:).*0.6 + 0.4); end
                    plot(ax, t, mu, 'Color', CHOICE_COLS(i,:), 'LineWidth', 1, ...
                         'DisplayName', sprintf('%s (N=%d)', CHOICE_LABELS{i}, st(i).N));
                end
                xline(ax, 0, '--', 'Color',[0.2,0.2,0.2], 'HandleVisibility','off'); grid(ax,'on');
                xlim(ax, [t(1), t(end)]); set(ax,'FontSize',FZ_AX);
                if c==1, ylabel(ax, rowYlabs{r}, 'Interpreter','none'); end
                if r==1
                    title(ax, sprintf('%s\n\n%s', condStats(c).cond, rowTitles{r}), 'Interpreter','none');
                else
                    title(ax, rowTitles{r}, 'Interpreter','none');
                end
                if r==rows, xlabel(ax, 'Time from event (s)'); end
                if r==1, legend(ax,'Location','southeast','Box','off','FontSize',6,'IconColumnWidth',4); else, legend(ax,'off'); end
            end
        end
        for r=1:rows
            ra = AX(r, isgraphics(AX(r,:))); if isempty(ra), continue; end
            yl = get(ra,'YLim'); if iscell(yl), yl = cell2mat(yl); end
            set(ra,'YLim',[min(yl(:,1)), max(yl(:,2))]); linkaxes(ra,'y');
            for c=1:cols
                if DRAW_CHOICE_BOX && string(event_variant)=="m1rew1on"
                    add_xbox(AX(r,c), CHOICE_BOX_WIN, CHOICE_BOX_COLOR, CHOICE_BOX_ALPHA);
                end
            end
        end
        sg = sgtitle(tlG, sprintf('Population PSTHs — OVERLAP (BUY>HOLD ∩ BUY>SELL) — %s — %s', ...
                    char(ev_tag), condset_tag), 'Interpreter','none');
        set(sg,'FontSize',FZ_SGT);

        outDir = fullfile(BASE_OUT_DIR, sprintf('%s_%s_grid_rowsAtoG_colsConds_OVERLAP', char(ev_tag), condset_tag));
        if ~exist(outDir,'dir'); mkdir(outDir); end
        outPNG = fullfile(outDir, sprintf('PSTH_grid_OVERLAP__%s__%s.png', char(ev_tag), condset_tag));
        print(figG, outPNG, '-dpng', sprintf('-r%d', dpi));
    end
end


end

%% ============================== FUNCTIONS ================================
function mask = get_overlap_mask_for_pack(packPath, event_variant, TR, CH, SC_WIN, CHOICE_VALS, conds, ci, ev_tag, thr)
[~, diffs, ~, ~] = compute_spikecount_diffs_for_pack( ...
packPath, event_variant, TR, CH, SC_WIN, CHOICE_VALS, conds, ci, ev_tag);
if isempty(diffs)
mask = false(0,1);
else
BS = diffs(:,1); % Buy−Sell
BH = diffs(:,2); % Buy−Hold
mask = (BH > thr) & (BS > thr);
end
end

function draw_venn_two(ax, nA, nB, nAB, N, labels)
% Quick two-set Venn (schematic, not area-accurate). Displays counts.
% labels: {labelA, labelB}
if nargin < 6 || isempty(labels), labels = {'A','B'}; end
cla(ax); axis(ax,'equal'); axis(ax,[0 10 0 6]); axis(ax,'off'); hold(ax,'on');

% Layout params
r = 2.2; cA = [3.5, 3.0]; cB = [6.0, 3.0]; % centers
tcol = [0.1 0.1 0.1];
faceA = [0.6 0.8 1.0]; faceB = [1.0 0.8 0.6];

th = linspace(0, 2*pi, 400);
plot(ax, cA(1)+r*cos(th), cA(2)+r*sin(th), 'Color',[0.3 0.5 0.9], 'LineWidth',1.2);
plot(ax, cB(1)+r*cos(th), cB(2)+r*sin(th), 'Color',[0.9 0.6 0.3], 'LineWidth',1.2);
patch(ax, cA(1)+r*cos(th), cA(2)+r*sin(th), faceA, 'FaceAlpha',0.20, 'EdgeColor','none');
patch(ax, cB(1)+r*cos(th), cB(2)+r*sin(th), faceB, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------------- Counts & positions ---------------------------
onlyA = max(nA - nAB, 0);
onlyB = max(nB - nAB, 0);
neither = max(N - (onlyA + onlyB + nAB), 0);

% Labels
text(ax, cA(1)-r, cA(2)+r+0.6, sprintf('%s', labels{1}), 'FontWeight','bold','Color',[0.2 0.35 0.8]);
text(ax, cB(1)+r-0.2, cB(2)+r+0.6, sprintf('%s', labels{2}), 'FontWeight','bold','Color',[0.8 0.45 0.2], 'HorizontalAlignment','right');

% Counts in regions
text(ax, cA(1)-0.9, cA(2), sprintf('%d', onlyA), 'Color', tcol, 'FontSize', 10, 'HorizontalAlignment','center');
text(ax, (cA(1)+cB(1))/2, cA(2), sprintf('%d', nAB), 'Color', tcol, 'FontSize', 10, 'HorizontalAlignment','center', 'FontWeight','bold');
text(ax, cB(1)+0.9, cB(2), sprintf('%d', onlyB), 'Color', tcol, 'FontSize', 10, 'HorizontalAlignment','center');

% Neither (bottom center)
text(ax, 5, 0.45, sprintf('Neither: %d | Total: %d', neither, N), ...
'Color', [0.25 0.25 0.25], 'HorizontalAlignment','center');
end

function [meanCounts, diffs, nTrials, info] = compute_spikecount_diffs_for_pack( ...
packPath, event_variant, TR, CH, SC_WIN, CHOICE_VALS, conds, ci, ev_tag)

S = load(packPath,'C'); C = S.C; dt = C.dt;

K = numel(CHOICE_VALS);
meanCounts = []; % (#neurons x K)
nTrials = []; % (#neurons x K)
rows_added = 0;

sc_bins = round(SC_WIN./dt); % [start end] in bins, event-aligned
ibeg = sc_bins(1); iend = sc_bins(2)-1; % [start, end) style

for f = 1:numel(C.files)
T = C.eventTables{f};
if isempty(T) || height(T)~=1, continue; end

props = string(T.Properties.VariableNames);
evIdx = find(ismember(event_variant, props), 1, 'first');
if isempty(evIdx), continue; end
evName = char(event_variant(evIdx));
if ~iscell(T.(evName)) || isempty(T.(evName){1}), continue; end
idxMat = T.(evName){1};                % 15x6 centers
Sbin   = C.S{f}; [Nf,Tf] = size(Sbin);

hasOption = ismember('option', T.Properties.VariableNames) && iscell(T.option) && ~isempty(T.option{1});
if hasOption
    optMat = T.option{1};
    if ~ismatrix(optMat), hasOption = false; end
end

centers = nan(TR*CH,1); optLabels = nan(TR*CH,1);
k = 0;
for m = 1:CH
    for tr = 1:TR
        k = k + 1;
        centers(k) = idxMat(tr,m);
        if hasOption, optLabels(k) = optMat(tr,m); end
    end
end
valid = isfinite(centers) & centers>=1 & centers<=Tf;
centers = centers(valid);
if hasOption, optLabels = optLabels(valid); end
if isempty(centers), continue; end

% per-neuron summaries in this file
mCounts_f = nan(Nf, K);
nTrials_f = zeros(Nf, K);

for n = 1:Nf
    sc_counts = nan(numel(centers), 1);
    for r = 1:numel(centers)
        rng = (centers(r)+ibeg) : (centers(r)+iend-1);
        v = (rng>=1) & (rng<=Tf);
        if any(v)
            sc_counts(r) = sum(Sbin(n, rng(v)));
        end
    end

    for ii = 1:K
        if hasOption
            pick = ismember(optLabels, CHOICE_VALS(ii));
            if any(pick)
                mCounts_f(n, ii) = mean(sc_counts(pick), 'omitnan');
                nTrials_f(n, ii) = sum(pick);
            end
        end
    end
end

meanCounts = [meanCounts; mCounts_f]; %#ok<AGROW>
nTrials    = [nTrials;    nTrials_f]; %#ok<AGROW>
rows_added = rows_added + Nf;


end

% === neuron filter (match keep_indices logic used elsewhere) ===
if ~isempty(meanCounts)
neurons_before = size(meanCounts,1);
else
neurons_before = 0;
end

persistent KEEP_STATE;
neuron_validation_mat = 'C:\Users\plattlab\MSM\outputs_local\neuron_validation\keep_indices.mat';
if isempty(KEEP_STATE)
if exist(neuron_validation_mat,'file')
Skeep = load(neuron_validation_mat,'keep_indices','variant');
KEEP_STATE = Skeep;
else
KEEP_STATE = struct();
end
end

condset_label = replace(cell2mat(conds), ' ', '');
cond_label = replace(conds{ci}, ' ', '');
varname = 'z'; if isfield(KEEP_STATE,'variant'); varname = KEEP_STATE.variant; end

keep_j = 1:neurons_before;
if isfield(KEEP_STATE,'keep_indices') && ...
isfield(KEEP_STATE.keep_indices, condset_label) && ...
isfield(KEEP_STATE.keep_indices.(condset_label), cond_label) && ...
isfield(KEEP_STATE.keep_indices.(condset_label).(cond_label), varname)
% keep_j = KEEP_STATE.keep_indices.(condset_label).(cond_label).(varname);
if ~isempty(meanCounts), meanCounts = meanCounts(keep_j, :); end
if ~isempty(nTrials), nTrials = nTrials(keep_j, :); end
end

if ~isempty(meanCounts)
neurons_after_keep = size(meanCounts,1);
else
neurons_after_keep = 0;
end

% drop rows with no choice data
if neurons_after_keep > 0
has_any = any(isfinite(meanCounts), 2);
else
has_any = [];
end
additional_removed_no_option = sum(~has_any);
meanCounts = meanCounts(has_any, :);
nTrials = nTrials(has_any, :);
neurons_after_final = size(meanCounts,1);

% diffs: [Buy−Sell, Buy−Hold]
bIdx = find(CHOICE_VALS==1,1); hIdx = find(CHOICE_VALS==2,1); sIdx = find(CHOICE_VALS==3,1);
diffs = [meanCounts(:, bIdx) - meanCounts(:, sIdx), ...
meanCounts(:, bIdx) - meanCounts(:, hIdx)];

info = struct('neurons_before', neurons_before, ...
'neurons_after_keep', neurons_after_keep, ...
'additional_removed_no_option', additional_removed_no_option, ...
'neurons_after_final', neurons_after_final, ...
'keep_j', keep_j);
end

function label = cond_set_label(conds)
if isempty(conds), label = 'conds'; return; end
if isstring(conds), conds = cellstr(conds); end
label = strjoin(conds, '-');
label = regexprep(label, '[^\w-]', '');
end

function s = sanitize(x)
s = regexprep(x, '[^\w-]', '');
end

function [muRaw_all, semRaw_all, muA_all, semA_all, muB_all, semB_all, ...
muArowz_all, semArowz_all, muBrowz_all, semBrowz_all, ...
muAmatz_all, semAmatz_all, muBmatz_all, semBmatz_all, ...
N_all, NpbB_all, sdB_all, ...
neurons_before, neurons_after_keep, additional_removed_no_option, neurons_after_final, keep_j, ...
nfiles] = ...
compute_psth_choice_for_pack( ...
packPath, t, baseMask, event_variant, TR, CH, norm_method, ...
scale_to_hz, smooth_sigma_sec, ...
use_robust_sigma, sigma_shrink_lambda, sigma_ridge, ...
trial_sigma_min_samples, sigma_floor, CHOICE_VALS, conds, ci, ev_tag, restrict_rows)
% Mirrors compute_psth_for_pack, but splits by choice option,
% and applies the SAME neuron filter logic as the cond version.
% IMPORTANT: Keeps per-neuron row alignment across choices even if 'option'
% is missing or a neuron has no trials of a given choice (fills NaNs).

if nargin < 19, restrict_rows = []; end

S = load(packPath,'C'); C = S.C; dt = C.dt;
relBins = round(t ./ dt);

% Per-choice accumulators (rows = neurons, aligned across choices)
K = numel(CHOICE_VALS);
M_raw = cell(1,K);
M_A = cell(1,K);
M_B = cell(1,K);
M_Arowz = cell(1,K);
M_Browz = cell(1,K);
M_Amatz = cell(1,K);
M_Bmatz = cell(1,K);

for k = 1:K
M_raw{k} = [];
M_A{k} = [];
M_B{k} = [];
M_Arowz{k} = [];
M_Browz{k} = [];
M_Amatz{k} = [];
M_Bmatz{k} = [];
end

nfiles = 0;

for f = 1:numel(C.files)
T = C.eventTables{f};
if isempty(T) || height(T)~=1, continue; end

% Event presence check (same behavior as cond script)
props = string(T.Properties.VariableNames);
evIdx = find(ismember(event_variant, props), 1, 'first');
if isempty(evIdx), continue; end
evName = char(event_variant(evIdx));
if ~iscell(T.(evName)) || isempty(T.(evName){1}), continue; end
idxMat = T.(evName){1};               % 15×6 indices
if ~ismatrix(idxMat), continue; end

Sbin = C.S{f}; [Nf, Tf] = size(Sbin);

% If 'option' missing/empty, we will fill NaN rows per neuron for EACH choice.
hasOption = ismember('option', T.Properties.VariableNames) && iscell(T.option) && ~isempty(T.option{1});
if hasOption
    optMat = T.option{1};             % 15×6 (1..3)
    if ~ismatrix(optMat), hasOption = false; end
end

% Flatten valid trial centers (+ options if available)
centers = nan(TR*CH,1); optLabels = nan(TR*CH,1);
k = 0;
for m = 1:CH
    for tr = 1:TR
        k = k + 1;
        centers(k) = idxMat(tr,m);
        if hasOption, optLabels(k) = optMat(tr,m); end
    end
end
valid = isfinite(centers) & centers>=1 & centers<=Tf;
centers = centers(valid);
if hasOption, optLabels = optLabels(valid); end
if isempty(centers), continue; end

segIdxAll = centers + relBins;                 % (#rows × #time)
validSeg  = segIdxAll>=1 & segIdxAll<=Tf;

% ---------- Per neuron ----------
for n = 1:size(Sbin,1)
    % Build trials×time matrix (NaN outside valid span)
    X = nan(numel(centers), numel(t));
    for r = 1:numel(centers)
        vmask = validSeg(r,:);
        if any(vmask)
            X(r, vmask) = double(Sbin(n, segIdxAll(r, vmask)));
        end
    end

    % Binary→Hz and smoothing
    if scale_to_hz, X = X ./ dt; end
    if smooth_sigma_sec > 0
        sigma_bins = smooth_sigma_sec ./ dt;
        for rr = 1:size(X,1)
            X(rr,:) = conv_gauss_reflect(X(rr,:), sigma_bins);
        end
    end

    % ---- A: baseline-sub per trial ----
    b_trial_mu = mean(X(:, baseMask), 2, 'omitnan');
    XA = X - b_trial_mu;

    % ---- B: baseline-sub per neuron ----
    XB = X - mean(b_trial_mu, 'omitnan');

    % ---- Arowz: baseline-sub per trial, then z-score per trial ----
    muArow = mean(XA(:, baseMask), 2, 'omitnan');
    sdArow = std (XA(:, baseMask), [], 2, 'omitnan');
    XArowz = (XA - muArow) ./ sdArow;

    % ---- Browz: baseline-sub per neuron, then z-score per trial ----
    muBrow = mean(XB(:, baseMask), 2, 'omitnan');
    sdBrow = std (XB(:, baseMask), [], 2, 'omitnan');
    XBrowz = (XB - muBrow) ./ sdBrow;

    % ---- Amatz: baseline-sub per trial, then z-score per neuron ----
    muAmat = mean(XA(:, baseMask), 'all', 'omitnan');
    sdAmat = std (XA(:, baseMask), [], 'all', 'omitnan');
    XAmatz = (XA - muAmat) ./ sdAmat;

    % ---- Bmatz: baseline-sub per neuron, then z-score per neuron ----
    muBmat = mean(XB(:, baseMask), 'all', 'omitnan');
    sdBmat = std (XB(:, baseMask), [], 'all', 'omitnan');
    XBmatz = (XB - muBmat) ./ sdBmat;

    % ---- Append one row per CHOICE for THIS neuron (keep row alignment) ----
    for ii = 1:K
        if hasOption
            pick = ismember(optLabels, CHOICE_VALS(ii));
            if any(pick)
                % xr = mean(X (pick,:), 1, 'omitnan');
                % xa = mean(XA(pick,:), 1, 'omitnan');
                % xb = mean(XB(pick,:), 1, 'omitnan');
                % xarowz = mean(XArowz(pick,:), 1, 'omitnan');
                % xbrowz = mean(XBrowz(pick,:), 1, 'omitnan');
                % xamatz = mean(XAmatz(pick,:), 1, 'omitnan');
                % xbmatz = mean(XBmatz(pick,:), 1, 'omitnan');
                xr = X (pick,:);
                xa = XA(pick,:);
                xb = XB(pick,:);
                xarowz = (XA(pick,:) - mean(XA(pick,:), 2)) ./ std(XA(pick,:), [], 2);
                xbrowz = (XB(pick,:) - mean(XB(pick,:), 2)) ./ std(XB(pick,:), [], 2);
                xamatz = (XA(pick,:) - mean(XA, 'all')) ./ std(XA, [], 'all');
                xbmatz = (XB(pick,:) - mean(XB, 'all')) ./ std(XB, [], 'all');
            else
                xr = nan(1, numel(t));
                xa = xr;
                xb = xr;
                xarowz = xr;
                xbrowz = xr;
                xamatz = xr;
                xbmatz = xr;
            end
        else
            % 'option' missing → preserve neuron row with NaNs (no contribution to means)
            xr = nan(1, numel(t));
            xa = xr;
            xb = xr;
            xarowz = xr;
            xbrowz = xr;
            xamatz = xr;
            xbmatz = xr;
        end
        % M_raw{ii}   = [M_raw{ii};   xr]; %#ok<AGROW>
        % M_A{ii}     = [M_A{ii};     xa]; %#ok<AGROW>
        % M_B{ii}     = [M_B{ii};     xb]; %#ok<AGROW>
        % M_Arowz{ii} = [M_Arowz{ii}; xarowz]; %#ok<AGROW>
        % M_Browz{ii} = [M_Browz{ii}; xbrowz]; %#ok<AGROW>
        % M_Amatz{ii} = [M_Amatz{ii}; xamatz]; %#ok<AGROW>
        % M_Bmatz{ii} = [M_Bmatz{ii}; xbmatz]; %#ok<AGROW>
        M_raw{ii}   = [M_raw{ii};   mean(xr, 1)]; %#ok<AGROW>
        M_A{ii}     = [M_A{ii};     mean(xa, 1)]; %#ok<AGROW>
        M_B{ii}     = [M_B{ii};     mean(xb, 1)]; %#ok<AGROW>
        M_Arowz{ii} = [M_Arowz{ii}; mean(xarowz, 1)]; %#ok<AGROW>
        M_Browz{ii} = [M_Browz{ii}; mean(xbrowz, 1)]; %#ok<AGROW>
        M_Amatz{ii} = [M_Amatz{ii}; mean(xamatz, 1)]; %#ok<AGROW>
        M_Bmatz{ii} = [M_Bmatz{ii}; mean(xbmatz, 1)]; %#ok<AGROW>
    end
end

nfiles = nfiles + 1;


end

% --- A) capture BEFORE FILTERING (place right BEFORE the OPTIONAL NEURON FILTER block) ---
if ~isempty(M_A) && ~isempty(M_A{1})
neurons_before = size(M_A{1}, 1);
else
neurons_before = 0;
end

% ====================== OPTIONAL NEURON FILTER (MATCH COND VERSION) ======================
persistent KEEP_STATE;
neuron_validation_mat = 'C:\Users\plattlab\MSM\outputs_local\neuron_validation\keep_indices.mat';
if isempty(KEEP_STATE)
if exist(neuron_validation_mat,'file')
Skeep = load(neuron_validation_mat,'keep_indices','variant'); % fix stray quote
KEEP_STATE = Skeep; % keep_indices.(setname).(grpname).(variant) = [idx...]
else
KEEP_STATE = struct();
end
end

condset_label = replace(cell2mat(conds), ' ', ''); % e.g., 'AI_Replay_Decoy_Live'
cond_label = replace(conds{ci}, ' ', ''); % e.g., 'AI'
varname = 'z';
if isfield(KEEP_STATE,'variant'); varname = KEEP_STATE.variant; end

keep_j = 1:neurons_before; % default (no filtering if matrix missing)
if isfield(KEEP_STATE,'keep_indices') && ...
isfield(KEEP_STATE.keep_indices, condset_label) && ...
isfield(KEEP_STATE.keep_indices.(condset_label), cond_label) && ...
isfield(KEEP_STATE.keep_indices.(condset_label).(cond_label), varname)
% keep_j = KEEP_STATE.keep_indices.(condset_label).(cond_label).(varname);
for ii = 1:K
if ~isempty(M_raw{ii}), M_raw{ii} = M_raw{ii}(keep_j, :); end
if ~isempty(M_A{ii}), M_A{ii} = M_A{ii}(keep_j, :); end
if ~isempty(M_B{ii}), M_B{ii} = M_B{ii}(keep_j, :); end
if ~isempty(M_Arowz{ii}), M_Arowz{ii} = M_Arowz{ii}(keep_j, :); end
if ~isempty(M_Browz{ii}), M_Browz{ii} = M_Browz{ii}(keep_j, :); end
if ~isempty(M_Amatz{ii}), M_Amatz{ii} = M_Amatz{ii}(keep_j, :); end
if ~isempty(M_Bmatz{ii}), M_Bmatz{ii} = M_Bmatz{ii}(keep_j, :); end
end
end

% plot_FR_dist_by_choice(M_raw, cond_label);
% plot_FR_dist_by_choice(M_A, cond_label);
% plot_FR_dist_by_choice(M_B, cond_label);

% --- B) counts/removals (place RIGHT AFTER applying keep_j, and BEFORE "%% Aggregate per-choice") ---
if ~isempty(M_A) && ~isempty(M_A{1})
neurons_after_keep = size(M_A{1}, 1);
else
neurons_after_keep = 0;
end

% rows that have ANY finite data in ANY choice (NaN-only rows came from missing/no 'option')
if neurons_after_keep > 0
has_any = false(neurons_after_keep,1);
for ii = 1:numel(M_A)
if ~isempty(M_A{ii})
has_any = has_any | any(isfinite(M_A{ii}), 2);
end
end
else
has_any = [];
end

additional_removed_no_option = sum(~has_any);
neurons_after_final = neurons_after_keep - additional_removed_no_option;

% drop NaN-only rows so aggregation reflects the final counts
if neurons_after_keep > 0 && any(~has_any)
for ii = 1:numel(M_A)
if ~isempty(M_A{ii}), M_A{ii} = M_A{ii}(has_any, :); end
if ~isempty(M_B{ii}), M_B{ii} = M_B{ii}(has_any, :); end
if ~isempty(M_raw{ii}), M_raw{ii} = M_raw{ii}(has_any, :); end
if ~isempty(M_Arowz{ii}), M_Arowz{ii} = M_Arowz{ii}(has_any,:); end
if ~isempty(M_Browz{ii}), M_Browz{ii} = M_Browz{ii}(has_any,:); end
if ~isempty(M_Amatz{ii}), M_Amatz{ii} = M_Amatz{ii}(has_any,:); end
if ~isempty(M_Bmatz{ii}), M_Bmatz{ii} = M_Bmatz{ii}(has_any,:); end
end
end

fprintf('\n[%s | %s]\n', char(ev_tag), cond_label);
fprintf('\t1) Neurons before filtering: %d\n', neurons_before);
fprintf('\t2) Neurons after filtering: %d\n', neurons_after_final);
fprintf('\t\t Removed due to no match in Annamarie''s matrix: %d\n', max(neurons_before-length(keep_j),0));
fprintf('\t\t Removed due to no option: %d\n', additional_removed_no_option);

if ~isempty(restrict_rows)
if islogical(restrict_rows), idx = find(restrict_rows); else, idx = restrict_rows(:); end
for ii = 1:numel(M_A)
if ~isempty(M_A{ii}), M_A{ii} = M_A{ii}(idx, :); end
if ~isempty(M_B{ii}), M_B{ii} = M_B{ii}(idx, :); end
if ~isempty(M_raw{ii}), M_raw{ii} = M_raw{ii}(idx, :); end
if ~isempty(M_Arowz{ii}), M_Arowz{ii} = M_Arowz{ii}(idx, :); end
if ~isempty(M_Browz{ii}), M_Browz{ii} = M_Browz{ii}(idx, :); end
if ~isempty(M_Amatz{ii}), M_Amatz{ii} = M_Amatz{ii}(idx, :); end
if ~isempty(M_Bmatz{ii}), M_Bmatz{ii} = M_Bmatz{ii}(idx, :); end
end
neurons_after_final = numel(idx); % update count
end

% =========================================================================================

% ---- Aggregate per-choice ----
muRaw_all = cell(1,K); semRaw_all = cell(1,K);
muA_all = cell(1,K); semA_all = cell(1,K);
muB_all = cell(1,K); semB_all = cell(1,K);

% Completed/returned variants:
muArowz_all = cell(1,K); semArowz_all = cell(1,K);
muBrowz_all = cell(1,K); semBrowz_all = cell(1,K);
muAmatz_all = cell(1,K); semAmatz_all = cell(1,K);
muBmatz_all = cell(1,K); semBmatz_all = cell(1,K);

NpbB_all = cell(1,K); sdB_all = cell(1,K);
N_all = zeros(1,K);

for ii = 1:K
if isempty(M_raw{ii})
muRaw_all{ii} = nan(1,numel(t)); semRaw_all{ii} = nan(1,numel(t));
muA_all{ii} = nan(1,numel(t)); semA_all{ii} = nan(1,numel(t));
muB_all{ii} = nan(1,numel(t)); semB_all{ii} = nan(1,numel(t));

    muArowz_all{ii} = nan(1,numel(t)); semArowz_all{ii} = nan(1,numel(t));
    muBrowz_all{ii} = nan(1,numel(t)); semBrowz_all{ii} = nan(1,numel(t));
    muAmatz_all{ii} = nan(1,numel(t)); semAmatz_all{ii} = nan(1,numel(t));
    muBmatz_all{ii} = nan(1,numel(t)); semBmatz_all{ii} = nan(1,numel(t));

    NpbB_all{ii}= zeros(1,numel(t)); sdB_all{ii} = nan(1,numel(t));
    N_all(ii) = 0;
else
    % N counts neurons with ANY finite value after filtering (ignores NaN-only placeholders)
    row_has_data = any(isfinite(M_raw{ii}), 2);
    N_all(ii) = sum(row_has_data);

    [muRaw, semRaw]                   = pop_mean_sem(M_raw{ii});
    [muA, semA]                       = pop_mean_sem(M_A{ii});
    [muB, semB, NpbB, sdB]            = pop_mean_sem(M_B{ii});
    [muArowz, semArowz]               = pop_mean_sem(M_Arowz{ii});
    [muBrowz, semBrowz]               = pop_mean_sem(M_Browz{ii});
    [muAmatz, semAmatz]               = pop_mean_sem(M_Amatz{ii});
    [muBmatz, semBmatz]               = pop_mean_sem(M_Bmatz{ii});

    muRaw_all{ii} = muRaw; semRaw_all{ii} = semRaw;
    muA_all{ii}   = muA;   semA_all{ii}   = semA;
    muB_all{ii}   = muB;   semB_all{ii}   = semB;

    muArowz_all{ii} = muArowz; semArowz_all{ii} = semArowz;
    muBrowz_all{ii} = muBrowz; semBrowz_all{ii} = semBrowz;
    muAmatz_all{ii} = muAmatz; semAmatz_all{ii} = semAmatz;
    muBmatz_all{ii} = muBmatz; semBmatz_all{ii} = semBmatz;

    NpbB_all{ii} = NpbB; sdB_all{ii} = sdB;
end


end
end

function [mu, sem, Npb, sd] = pop_mean_sem(M)
% M: (#neurons x #time)
mu = mean(M, 1, 'omitnan');
sd = std (M, 0, 1, 'omitnan');
Npb = sum(~isnan(M), 1);
sem = sd ./ max(1, sqrt(Npb));
mu = mu(:)'; sem = sem(:)'; Npb = Npb(:)'; sd = sd(:)';
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
x = -hw:hw;
k = exp(-0.5 .* (x ./ sigma_bins).^2);
k = k ./ sum(k);
end

function y = conv_gauss_reflect(x, sigma_bins)
if sigma_bins <= 0, y = x; return; end
x = x(:)'; k = gauss1d(sigma_bins);
hw = (numel(k) - 1) ./ 2;

m = double(isfinite(x));
x0 = x; x0(~isfinite(x0)) = 0;

if hw < 1
num = conv(x0, k, 'same');
den = conv(m, k, 'same');
y = num ./ den; y(den < eps) = NaN; return;
end

leftPadX = fliplr(x0(1:hw));
rightPadX= fliplr(x0(end-hw+1:end));
leftPadM = fliplr(m(1:hw));
rightPadM= fliplr(m(end-hw+1:end));

xp = [leftPadX, x0, rightPadX];
mp = [leftPadM, m, rightPadM];

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

function plot_FR_dist_by_choice(M, cond_label)
testBuy = mean(M{1}(:, 500:1500), 2);
testHold = mean(M{2}(:, 500:1500), 2);
testSell = mean(M{3}(:, 500:1500), 2);

figure

subplot(3, 1, 1)
histogram(testBuy, 0:10:200)
ylim([0 50])
ylabel("Neurons")
xlabel("Firing rate after reward onset (Hz)")
title("Buy")
subplot(3, 1, 2)
histogram(testHold, 0:10:200)
ylim([0 50])
ylabel("Neurons")
xlabel("Firing rate after reward onset (Hz)")
title("Hold")
subplot(3, 1, 3)
histogram(testSell, 0:10:200)
ylim([0 50])
ylabel("Neurons")
xlabel("Firing rate after reward onset (Hz)")
title("Sell")

sgtitle(string(cond_label) + " (N = " + size(testBuy, 1) + ")")
end