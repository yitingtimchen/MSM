% File: d20251214_analyze_timestamped_total_reward_vs_training.m
%% Total Reward over Training Sessions per Monkey
% This script computes, for each session:
%   - sumFb1      = sum over trials of fb1
%   - sumFb2      = sum over trials of fb2
%   - totalReward = sum over trials of (fb1 + fb2)
%   - nTrials     = number of valid trials with non-NaN (fb1+fb2)
% It then:
%   1) Builds a per-session summary table across all treatments.
%   2) Assigns a training index (sessionIdx) within each monkey
%      based on chronological session order.
%   3) Plots a chosen reward metric vs sessionIdx for each monkey with a
%      single linear trendline (pooling all treatments),
%      shades the 95% CI around the fit,
%      and flags outlier sessions (|standardized residual| > 3).
%
% REQUIREMENTS:
%   - Condition pack .mat files built from dataTabFULL, located in OUTDIR.
%   - Each pack contains a struct C with field C.eventTables, where each
%     entry is a table "etab" having variables:
%       fb1, fb2, year, month, day, SessionOfDay, monkey
%
% NOTE:
%   - "Treatment" is defined by which group of condition names the session
%     came from: baseline, OT, or Saline.
%   - This script only uses reward, date/session info, and monkey ID;
%     it ignores choices and other behavioral covariates.

clear; clc; 
% close all;

%% ---------- CONFIG ----------
% Path to condition packs (behavior_only_from_dataTabFULL)
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';

% Condition sets grouped by treatment
% (These must match the naming used when building the packs.)
cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames = {'baseline','OT','Saline'};

% TOGGLE: which metric to plot?  'fb1', 'fb2', or 'total'
plotMetric = 'fb1';  % <<< CHANGE THIS TO 'fb1', 'fb2', OR 'total' AS NEEDED

%% ---------- MASTER SESSION SUMMARY ACROSS ALL TREATMENTS ----------
% One row per session/file:
%   Treatment, ConditionRaw, Monkey, Session (yyyy-mm-dd-ss),
%   SessionDate (datetime), sumFb1, sumFb2, totalReward, meanReward, nTrials
sessionSummaryAll = table();

%% ================== MAIN LOOP OVER TREATMENT SETS ==================
for s = 1:numel(cond_sets)
    conds     = cond_sets{s};
    treatName = setNames{s};

    fprintf('\n=== Treatment: %s ===\n', treatName);

    % Per-treatment session summary table
    sessLocal = table();

    % ---------- LOOP OVER CONDITIONS IN THIS TREATMENT ----------
    for c = 1:numel(conds)
        cname   = conds{c};  % e.g. 'AI', 'OT Live', 'Saline Decoy'
        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));

        if ~exist(matFile,'file')
            warning('Missing condition pack: %s', matFile);
            continue;
        end

        tmp = load(matFile, 'C');
        C   = tmp.C;

        % ---------- LOOP OVER FILES (SESSIONS) IN THIS CONDITION ----------
        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            % We need reward columns fb1 and fb2
            if ~ismember('fb1', etab.Properties.VariableNames) || ...
               ~ismember('fb2', etab.Properties.VariableNames)
                warning('Missing fb1/fb2 in condition %s, file %d; skipping.', cname, f);
                continue;
            end

            if isempty(etab.fb1{1}) || isempty(etab.fb2{1})
                continue;
            end

            % ---- Per-trial reward arrays ----
            fb1 = etab.fb1{1};
            fb2 = etab.fb2{1};

            % Flatten to vectors (TR x MKT -> TR*MKT)
            fb1_vec = fb1(:);
            fb2_vec = fb2(:);

            % Use fb1+fb2 to define valid trials (both present)
            rewardVec = fb1_vec + fb2_vec;
            valid     = ~isnan(rewardVec);
            if ~any(valid)
                continue;
            end

            % Session-level sums
            sumFb1      = sum(fb1_vec(valid));
            sumFb2      = sum(fb2_vec(valid));
            totalReward = sumFb1 + sumFb2;
            nTrials     = sum(valid);
            meanReward  = totalReward / nTrials;

            % ---- True session ID (yyyy-mm-dd-ss) ----
            neededDateVars = {'year','month','day','SessionOfDay'};
            if ~all(ismember(neededDateVars, etab.Properties.VariableNames))
                warning('Missing date/session fields; skipping a file in %s.', cname);
                continue;
            end

            y   = etab.year{1}(:);
            m   = etab.month{1}(:);
            d   = etab.day{1}(:);
            sod = etab.SessionOfDay{1}(:);

            if isempty(y) || isempty(m) || isempty(d) || isempty(sod)
                continue;
            end

            y0   = y(1);
            m0   = m(1);
            d0   = d(1);
            sod0 = sod(1);

            sessionStr  = string(sprintf('%04d-%02d-%02d-%02d', y0, m0, d0, sod0));
            sessionDate = datetime(y0, m0, d0);

            % ---- Monkey ID (1 or 2) ----
            if ~ismember('monkey', etab.Properties.VariableNames) || isempty(etab.monkey{1})
                warning('Missing monkey field; skipping a file in %s.', cname);
                continue;
            end
            monkeyVal = etab.monkey{1}(1);
            monkeyStr = string(monkeyVal);

            % ---- Append one row per session ----
            newRow = table( ...
                string(treatName), ...   % Treatment
                string(cname),   ...     % ConditionRaw
                monkeyStr,       ...     % Monkey
                sessionStr,      ...     % Session (yyyy-mm-dd-ss)
                sessionDate,     ...     % SessionDate (datetime, date only)
                sumFb1,          ...     % sumFb1 for this session
                sumFb2,          ...     % sumFb2 for this session
                totalReward,     ...     % totalReward for this session
                meanReward,      ...     % meanReward per valid trial
                nTrials,         ...     % number of valid trials
                'VariableNames', { ...
                    'Treatment','ConditionRaw','Monkey','Session','SessionDate', ...
                    'sumFb1','sumFb2','totalReward','meanReward','nTrials' ...
                } ...
            );

            sessLocal = [sessLocal; newRow]; %#ok<AGROW>
        end
    end

    if isempty(sessLocal)
        warning('No sessions for treatment %s; skipping.', treatName);
        continue;
    end

    % Append to global summary
    sessionSummaryAll = [sessionSummaryAll; sessLocal]; %#ok<AGROW>
end

%% ---------- CHECK & BUILD TRAINING SESSION INDEX ----------
if isempty(sessionSummaryAll)
    warning('No sessions found in any treatment. Nothing to analyze.');
    return;
end

% Convert to categoricals for indexing
sessionSummaryAll.Monkey  = categorical(sessionSummaryAll.Monkey);
% sessionSummaryAll.Session = categorical(sessionSummaryAll.Session);

% For each monkey, sort sessions chronologically by Session string and
% assign an integer sessionIdx = 1,2,3,... in that order.
sessTab = sessionSummaryAll(:, {'Monkey','Session'});
[sessTabSorted, sortIdx] = sortrows(sessTab, {'Monkey','Session'});

sessionIdxSorted = zeros(height(sessTabSorted),1);
monkeys = categories(sessTabSorted.Monkey);

for im = 1:numel(monkeys)
    mID  = monkeys{im};
    mask = sessTabSorted.Monkey == mID;
    nS   = nnz(mask);
    sessionIdxSorted(mask) = (1:nS).';
end

sessionIdx = zeros(height(sessionSummaryAll),1);
sessionIdx(sortIdx) = sessionIdxSorted;
sessionSummaryAll.sessionIdx = sessionIdx;

%% ---------- DERIVE CONDITION BASE (for color coding) ----------
% ConditionBase is the "base" condition (AI, Replay, Decoy, Live),
% stripping OT/Saline prefixes, etc.
condBase = strings(height(sessionSummaryAll),1);
for i = 1:height(sessionSummaryAll)
    parts = split(string(sessionSummaryAll.ConditionRaw(i)));
    condBase(i) = parts(end);
end
sessionSummaryAll.ConditionBase = categorical(condBase);

%% ---------- PLOT SETTINGS (markers for treatment, colors for condition base) ----------
% Treatments: keep original setNames order in legends/loops
treatAllCats = categories(categorical(sessionSummaryAll.Treatment));
treatCats    = setNames(ismember(setNames, treatAllCats));  % preserves baseline, OT, Saline order

markerList   = {'o','s','^','d','v','>'};  % different shape per treatment (use first N)

% Condition bases: keep base order as in cond_sets: AI, Replay, Decoy, Live
baseOrder    = {'AI','Replay','Decoy','Live'};
allBaseCats  = categories(sessionSummaryAll.ConditionBase);
baseCats     = baseOrder(ismember(baseOrder, allBaseCats));  % preserves AI, Replay, Decoy, Live order
baseColors   = lines(numel(baseCats));                       % different color per base

% ---------- SELECT FIELD / LABELS BASED ON plotMetric ----------
switch lower(plotMetric)
    case 'fb1'
        yField   = 'sumFb1';
        yLabel   = 'Total fb1 (seconds) per session';
        fig1Name = 'Total choice reward vs training (by session order)';
        fig2Name = 'Total choice reward vs training (by date)';
        sg1Title = 'Total choice reward over training sessions (sessionIdx, all treatments)';
        sg2Title = 'Total choice reward over training sessions (session date, all treatments)';
    case 'fb2'
        yField   = 'sumFb2';
        yLabel   = 'Total fb2 (seconds) per session';
        fig1Name = 'Total portfolio reward vs training (by session order)';
        fig2Name = 'Total portfolio reward vs training (by date)';
        sg1Title = 'Total portfolio reward over training sessions (sessionIdx, all treatments)';
        sg2Title = 'Total portfolio reward over training sessions (session date, all treatments)';
    otherwise
        yField   = 'totalReward';
        yLabel   = 'Total reward (seconds) per session';
        fig1Name = 'Total (choice + portfolio) reward vs training (by session order)';
        fig2Name = 'Total (choice + portfolio) reward vs training (by date)';
        sg1Title = 'Total (choice + portfolio) reward over training sessions (sessionIdx, all treatments)';
        sg2Title = 'Total (choice + portfolio) reward over training sessions (session date, all treatments)';
end

%% ---------- OPTIONAL: SAVE SUMMARY TABLE ----------
% Uncomment to save the per-session summary to disk for downstream analysis.
% save('sessionSummary_totalReward.mat', 'sessionSummaryAll');

%% ---------- PLOTTING (refactored) ----------
fprintf('\nPlotting %s vs training session and date per monkey...\n', yField);

% Plot vs sessionIdx
figure('Name', fig1Name, 'Color','w');
plot_metric_vs_x(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, ...
    markerList, yField, yLabel, 'idx', sg1Title);

% Plot vs session date
figure('Name', fig2Name, 'Color','w');
plot_metric_vs_x(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, ...
    markerList, yField, yLabel, 'date', sg2Title);


%% ---------- LOCAL FUNCTION (put this at end of file) ----------
function plot_metric_vs_x(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, ...
                          markerList, yField, yLabel, xMode, sgTitle)

nMonk = numel(monkeys);

for im = 1:nMonk
    mID   = monkeys{im};
    rowsM = sessionSummaryAll.Monkey == mID;

    subplot(1, nMonk, im); hold on;

    x_all = [];
    y_all = [];
    idx_all = [];

    % For legend: one handle per condition base (color), in baseCats order
    baseLegendHandles = gobjects(numel(baseCats),1);

    % Scatter points by treatment (marker shape) and condition base (color)
    for it = 1:numel(treatCats)
        tName  = treatCats{it};
        rowsT  = categorical(sessionSummaryAll.Treatment) == tName;
        marker = markerList{min(it, numel(markerList))};

        for ib = 1:numel(baseCats)
            bName = baseCats{ib};
            rowsB = sessionSummaryAll.ConditionBase == bName;

            rows = rowsM & rowsT & rowsB;
            if ~any(rows)
                continue;
            end

            switch lower(xMode)
                case 'idx'
                    x = double(sessionSummaryAll.sessionIdx(rows));
                case 'date'
                    x = datenum(sessionSummaryAll.SessionDate(rows));
                otherwise
                    error('Unknown xMode: %s', xMode);
            end
            y = sessionSummaryAll.(yField)(rows);

            h = scatter(x, y, 50, ...
                        'Marker', marker, ...
                        'MarkerFaceColor', baseColors(ib,:), ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerFaceAlpha', 0.7);

            x_all   = [x_all; x(:)];
            y_all   = [y_all; y(:)];
            idx_all = [idx_all; find(rows)];   % map back to table rows

            % Store one handle per base for legend (in baseCats order)
            if ~isgraphics(baseLegendHandles(ib))
                baseLegendHandles(ib) = h;
            end
        end
    end

    % Fit, shade 95% CI, plot linear trendline, and flag + label outliers
    if numel(x_all) > 1
        xvec = x_all(:);
        yvec = y_all(:);
        mdl  = fitlm(xvec, yvec);

        slope = mdl.Coefficients.Estimate(2);
        R2    = mdl.Rsquared.Ordinary;
        pval  = mdl.Coefficients.pValue(2);

        xfit = linspace(min(xvec), max(xvec), 100)';

        % 95% CI around regression line
        [yfit, yCI] = predict(mdl, xfit, 'Alpha', 0.05);
        fill([xfit; flipud(xfit)], [yCI(:,1); flipud(yCI(:,2))], ...
             [0.7 0.7 0.7], 'EdgeColor','none', 'FaceAlpha', 0.3);

        % Regression line
        plot(xfit, yfit, 'k-', 'LineWidth', 1.5);

        % Outlier detection via standardized residuals
        residStd    = mdl.Residuals.Standardized;
        outlierMask = abs(residStd) > 2.5;  % |z|>2.5 ~ outlier

        if any(outlierMask)
            x_out   = xvec(outlierMask);
            y_out   = yvec(outlierMask);
            idx_out = idx_all(outlierMask);

            plot(x_out, y_out, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, ...
                 'MarkerFaceColor', 'none');  % highlight outliers

            % Label outliers with sessionIdx or date
            for k = 1:numel(x_out)
                switch lower(xMode)
                    case 'idx'
                        lbl = sprintf('%d', sessionSummaryAll.sessionIdx(idx_out(k)));
                    case 'date'
                        thisDate = sessionSummaryAll.SessionDate(idx_out(k));
                        % lbl = datestr(thisDate, 'yyyy-mm-dd');
                        lbl = char(sessionSummaryAll.Session(idx_out(k)));
                end
                text(x_out(k), y_out(k), ['  ' lbl], ...
                     'HorizontalAlignment','left', ...
                     'VerticalAlignment','bottom', ...
                     'FontSize', 8, ...
                     'Color', 'r');
            end
        end

        % Annotate slope, R^2, and p-value on the plot
        txt = sprintf('slope = %.2f\nR^2 = %.3f\np = %.3g', slope, R2, pval);
        text(0.05, 0.7, txt, ...
             'Units','normalized', ...
             'HorizontalAlignment','left', ...
             'VerticalAlignment','top', ...
             'FontSize', 10);
    end

    % Axis labels / formatting
    switch lower(xMode)
        case 'idx'
            xlabel('sessionIdx (within monkey, all treatments)');
        case 'date'
            xlabel('Session date');

            if ~isempty(x_all)
                ax = gca;
                ax.XLim = [min(x_all) max(x_all)];

                uniqX   = unique(x_all);
                nTicks  = min(8, numel(uniqX));
                if nTicks > 1
                    ax.XTick = linspace(ax.XLim(1), ax.XLim(2), nTicks);
                else
                    ax.XTick = ax.XLim(1);
                end

                xt    = ax.XTick;
                xt_dt = datetime(xt, 'ConvertFrom', 'datenum');
                ax.XTickLabel = cellstr(datestr(xt_dt, 'yyyy-mm-dd'));
                ax.XTickLabelRotation = 45;
            end
    end

    ylabel(yLabel);
    title(sprintf('Monkey %s', char(mID)), 'Interpreter','none');

    % Legend: colors (condition base) + marker shapes (treatment)
    validMask = isgraphics(baseLegendHandles);

    % Dummy handles for marker shapes (treatments)
    treatLegendHandles = gobjects(numel(treatCats),1);
    for it = 1:numel(treatCats)
        mk = markerList{min(it, numel(markerList))};
        treatLegendHandles(it) = scatter(nan, nan, 50, ...
            'Marker', mk, ...
            'MarkerFaceColor', [0.85 0.85 0.85], ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 1.0);
    end

    % Build labels as string -> then to cellstr to keep dimensions consistent
    baseLabels  = "Condition: "  + string(baseCats(validMask));
    treatLabels = "Treatment: " + string(treatCats(:));

    legHandles  = [baseLegendHandles(validMask); treatLegendHandles(:)];
    legLabels   = cellstr([baseLabels(:); treatLabels(:)]);

    legend(legHandles, legLabels, 'Location','northwest', 'Box','off','NumColumns',2);
end

if exist('sgtitle','file')
    sgtitle(sgTitle, 'Interpreter','none');
end

end