% File: d20251214_analyze_repeat_and_opponent_effects.m
%% Repeat Opponent Choice Analysis (Decontaminated + Live-Only Opponent + Comparisons)
% -------------------------------------------------------------------------
% Main questions:
%   1) Live vs non-live effects within each Treatment
%   2) Difference-in-differences across Treatment (LiveFlag*Treatment)
%   3) Do the two monkeys show the same patterns / trends?
%
% Plot encodings:
%   - Marker SHAPE  = Treatment group (baseline / OT / Saline)
%   - Marker FILL   = Base condition (AI / Replay / Decoy / Live)
%   - Columns       = Monkey 1 / Monkey 2 / Both (pooled)
%
% Notes:
%   - "LiveFlag" means base condition == Live (real opponent).
%   - Some opponent metrics are only meaningful in Live sessions.
%   - Session-level models treat each session as binomial counts.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ---------- CONFIG ----------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';

% cond_sets   { ...
%     {'AI','Replay','Decoy','Live'}; ...
%     {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
%     {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
% };
% setNames     = {'baseline','OT','Saline'};
% baseCondList = {'AI','Replay','Decoy','Live'};
cond_sets = { ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames     = {'OT','Saline'};
baseCondList = {'AI','Replay','Decoy','Live'};

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
EXCLUDE_SESSIONS = { ...
    % '2018-07-07-01', ... % M1, baseline AI
    % '2018-07-20-02', ... % M1, baseline Replay
    % '2018-07-21-01', ... % M1, baseline Decoy
    % '2018-07-06-02', ... % M2, baseline AI
    % '2018-07-06-01', ... % M2, baseline AI
    % '2018-08-01-03', ... % M2, baseline Decoy
    % '2018-08-22-02', ... % M2, OT Decoy
    % '2018-06-07-01', ... % M1, baseline Live
    % '2018-06-07-02', ... % M1, baseline Decoy
    % '2018-06-08-01', ... % M2, baseline Replay
    % '2018-06-08-02', ... % M2, baseline Live
    % '2018-06-08-03', ... % M2, baseline Decoy
};
EXCLUDE_SESSIONS = string(EXCLUDE_SESSIONS(:));
EXCLUDE_SESSIONS = EXCLUDE_SESSIONS(strlength(EXCLUDE_SESSIONS) > 0);
excludedSessSeen = strings(0,1);

% ---------- TOGGLES ----------
PLOT_BAR_COND_MODE   = "Collapsed";  % "Base4" or "Collapsed"  (Collapsed = NonLive vs Live)
PLOT_2to5_TREAT_MODE = "Separate";   % "Combined" or "Separate"
MONKEY_EFFECT_MODE   = "Fixed";      % "Fixed" or "Random"  (Random uses (1|Monkey) when >1 monkey)
PLOT_SHOW_SIG        = true;         % add significance brackets + stars on bar plots

% Treatment -> marker (shape coding)
treatOrder  = {'baseline','OT','Saline'};   % intended reference order
markerList  = {'o','s','^','d','v','>'};
treatMarker = containers.Map(treatOrder, markerList(1:3));

% Base condition -> face color (fill coding)
condFaceColors = containers.Map();
condFaceColors('AI')     = [0.0000 0.4470 0.7410];  % default color 1
condFaceColors('Replay') = [0.8500 0.3250 0.0980];  % default color 2
condFaceColors('Decoy')  = [0.9290 0.6940 0.1250];  % default color 3
condFaceColors('Live')   = [0.4940 0.1840 0.5560];  % default color 4

%% ---------- BUILD MASTER SESSION SUMMARY ----------
sessionSummaryAll = table();

fprintf('\n=== BUILDING sessionSummaryAll ===\n');

for s = 1:numel(cond_sets)
    conds     = cond_sets{s};
    treatName = setNames{s};
    fprintf('Loading treatment: %s\n', treatName);

    for c = 1:numel(conds)
        cname   = conds{c};
        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));
        if ~exist(matFile,'file')
            warning('Missing condition pack: %s', matFile);
            continue;
        end

        tmp = load(matFile, 'C');
        C   = tmp.C;

        tok = regexp(cname, '(AI|Replay|Decoy|Live)$', 'tokens', 'once');
        if isempty(tok)
            warning('Could not parse base condition from name: %s', cname);
            baseCond = string(cname);
        else
            baseCond = string(tok{1});
        end

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            if ~ismember('option', etab.Properties.VariableNames) || ...
               ~ismember('optionOpp', etab.Properties.VariableNames)
                continue;
            end
            if isempty(etab.option{1}) || isempty(etab.optionOpp{1})
                continue;
            end

            needed = {'year','month','day','SessionOfDay'};
            if ~all(ismember(needed, etab.Properties.VariableNames))
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

            sessionStr = string(sprintf('%04d-%02d-%02d-%02d', y(1), m(1), d(1), sod(1)));

            if ~isempty(EXCLUDE_SESSIONS) && any(sessionStr == EXCLUDE_SESSIONS)
                excludedSessSeen(end+1,1) = sessionStr; %#ok<AGROW>
                continue;
            end

            if ~ismember('monkey', etab.Properties.VariableNames) || isempty(etab.monkey{1})
                warning('Missing monkey field; skipping a file in %s.', cname);
                continue;
            end
            monkeyStr = string(etab.monkey{1}(1));

            subjChoices = mapChoicesToLabels(etab.option{1});
            oppRaw      = mapChoicesToLabels(etab.optionOpp{1});

            subjChoices = subjChoices(:);
            oppRaw      = oppRaw(:);

            nT = min(numel(subjChoices), numel(oppRaw));
            subjChoices = subjChoices(1:nT);
            oppRaw      = oppRaw(1:nT);

            oppPrev  = [missing; oppRaw(1:end-1)];
            subjPrev = [missing; subjChoices(1:end-1)];

            % ---------- METRICS ----------
            % RAW subject match to opponent previous (diagnostic; can be contaminated)
            validRaw = ~ismissing(subjChoices) & ~ismissing(oppPrev);
            if ~any(validRaw), continue; end
            isRaw   = (subjChoices == oppPrev);
            pRepeat = mean(isRaw(validRaw));
            nRaw    = sum(validRaw);

            % Subject self-repeat (stickiness)
            validSelf = ~ismissing(subjChoices) & ~ismissing(subjPrev);
            isSelf = (subjChoices == subjPrev);
            if any(validSelf)
                pSelfRepeat = mean(isSelf(validSelf));
                nSelf       = sum(validSelf);
            else
                pSelfRepeat = NaN; nSelf = 0;
            end

            % Opponent followed subject on previous trial (diagnostic; mainly meaningful Live)
            validPrev = ~ismissing(oppPrev) & ~ismissing(subjPrev);
            isOppPrevMatch = (oppPrev == subjPrev);
            if any(validPrev)
                pOppPrevMatchSubjPrev = mean(isOppPrevMatch(validPrev));
                nPrevMatch            = sum(validPrev);
            else
                pOppPrevMatchSubjPrev = NaN; nPrevMatch = 0;
            end

            % Primary: subject repeats opponent previous choice, excluding "copy" trials
            validCond = validRaw & ~ismissing(subjPrev) & (oppPrev ~= subjPrev);
            if any(validCond)
                pRepeatCond = mean(isRaw(validCond));
                nCond       = sum(validCond);
            else
                pRepeatCond = NaN; nCond = 0;
            end

            % Live-only: opponent matches subject on same trial
            pRepeatOpp = NaN;
            nOpp       = 0;
            if baseCond == "Live"
                validOpp = ~ismissing(subjChoices) & ~ismissing(oppRaw);
                if any(validOpp)
                    isOpp = (oppRaw == subjChoices);
                    pRepeatOpp = mean(isOpp(validOpp));
                    nOpp       = sum(validOpp);
                end
            end

            newRow = table( ...
                string(treatName), string(cname), baseCond, ...
                monkeyStr, sessionStr, ...
                pRepeat, nRaw, ...
                pRepeatCond, nCond, ...
                pSelfRepeat, nSelf, ...
                pOppPrevMatchSubjPrev, nPrevMatch, ...
                pRepeatOpp, nOpp, ...
                'VariableNames', { ...
                    'Treatment','ConditionRaw','ConditionBase', ...
                    'Monkey','Session', ...
                    'pRepeat','nTrialsRaw', ...
                    'pRepeatCond','nTrialsCond', ...
                    'pSelfRepeat','nTrialsSelf', ...
                    'pOppPrevMatchSubjPrev','nTrialsPrevMatch', ...
                    'pRepeatOpp','nTrialsOpp' ...
                } ...
            );

            sessionSummaryAll = [sessionSummaryAll; newRow]; %#ok<AGROW>
        end
    end
end

%% ---------- EXCLUSION REPORT ----------
if ~isempty(EXCLUDE_SESSIONS)
    excludedSessSeenU = unique(excludedSessSeen);
    notFound = setdiff(EXCLUDE_SESSIONS, excludedSessSeenU);

    fprintf('\n=== Exclusion report ===\n');
    fprintf('Requested exclusions: %d\n', numel(EXCLUDE_SESSIONS));
    fprintf('Excluded (found+skipped): %d\n', numel(excludedSessSeenU));
    if ~isempty(excludedSessSeenU)
        disp('Excluded sessions actually seen:');
        disp(excludedSessSeenU);
    end
    if ~isempty(notFound)
        disp('Exclusion sessions NOT found in loaded data (check strings):');
        disp(notFound);
    end
end

%% ---------- GLOBAL SESSION INDEX + FACTORS ----------
if isempty(sessionSummaryAll)
    warning('No sessions found. Nothing to plot/model.');
    return;
end

sessionSummaryAll.Monkey = categorical(string(sessionSummaryAll.Monkey));

% sessionIdx: chronological order within each monkey
sessTab = sessionSummaryAll(:, {'Monkey','Session'});
[sessTabSorted, sortIdx] = sortrows(sessTab, {'Monkey','Session'});

sessionIdxSorted = zeros(height(sessTabSorted),1);
monkeys = categories(sessTabSorted.Monkey);
for im = 1:numel(monkeys)
    mask = sessTabSorted.Monkey == monkeys{im};
    sessionIdxSorted(mask) = (1:nnz(mask)).';
end

sessionIdx = zeros(height(sessionSummaryAll),1);
sessionIdx(sortIdx) = sessionIdxSorted;
sessionSummaryAll.sessionIdx = sessionIdx;

% factors
sessionSummaryAll.Treatment = categorical(string(sessionSummaryAll.Treatment));
sessionSummaryAll.Treatment = removecats(sessionSummaryAll.Treatment);
sessionSummaryAll.Treatment = reordercats(sessionSummaryAll.Treatment, ...
    treatOrder(ismember(treatOrder, categories(sessionSummaryAll.Treatment))));

sessionSummaryAll.ConditionBase = categorical(string(sessionSummaryAll.ConditionBase), baseCondList);
sessionSummaryAll.LiveFlag = categorical( ...
    string(sessionSummaryAll.ConditionBase) == "Live", ...
    [false true], {'NonLive','Live'} ...
);

%% ========================================================================
%% ================================ PLOTS =================================
%% ========================================================================

%% (1) PER-TREATMENT BAR PLOTS (PRIMARY metric = pRepeatCond)
fprintf('\nPlotting: Subject repeats opponent previous choice (filtered), by condition group...\n');
plotPrimaryBarsByTreatment_3col(sessionSummaryAll, treatOrder, baseCondList, ...
    PLOT_BAR_COND_MODE, PLOT_SHOW_SIG);

%% (2)-(5) Learning-style scatter plots
if PLOT_2to5_TREAT_MODE == "Combined"
    fprintf('\nPlotting learning-style scatter plots (all treatments together)...\n');
    plotNonLiveVsLiveScatter_3col(sessionSummaryAll, 'pRepeatCond', treatOrder, treatMarker, condFaceColors);
    plotNonLiveVsLiveScatter_3col(sessionSummaryAll, 'pSelfRepeat', treatOrder, treatMarker, condFaceColors);
    plotNonLiveVsLiveScatter_3col(sessionSummaryAll, 'pOppPrevMatchSubjPrev', treatOrder, treatMarker, condFaceColors);
    plotLiveOnlyOpponentMetrics_3col(sessionSummaryAll, treatOrder, treatMarker, condFaceColors);
elseif PLOT_2to5_TREAT_MODE == "Separate"
    for it = 1:numel(treatOrder)
        tName = treatOrder{it};
        ST = sessionSummaryAll(string(sessionSummaryAll.Treatment) == string(tName), :);
        if isempty(ST), continue; end
        fprintf('\n[Scatter] Treatment=%s\n', tName);

        plotNonLiveVsLiveScatter_3col(ST, 'pRepeatCond', treatOrder, treatMarker, condFaceColors, "Treatment: " + string(tName));
        plotNonLiveVsLiveScatter_3col(ST, 'pSelfRepeat', treatOrder, treatMarker, condFaceColors, "Treatment: " + string(tName));
        plotNonLiveVsLiveScatter_3col(ST, 'pOppPrevMatchSubjPrev', treatOrder, treatMarker, condFaceColors, "Treatment: " + string(tName));
        plotLiveOnlyOpponentMetrics_3col(ST, treatOrder, treatMarker, condFaceColors, "Treatment: " + string(tName));
    end
else
    error('Unknown PLOT_2to5_TREAT_MODE: %s', PLOT_2to5_TREAT_MODE);
end

%% (6) DIRECT MEAN COMPARISONS (bars)
fprintf('\nPlotting mean comparisons (bars): Live vs non-live, and within Live by treatment...\n');

% pooled across all treatments
plotMeanBars_LiveVsNonLive_3col(sessionSummaryAll, PLOT_SHOW_SIG, "All treatments pooled");

% per treatment (one figure per treatment)
treatLevels = treatOrder(ismember(treatOrder, categories(sessionSummaryAll.Treatment)));
for i = 1:numel(treatLevels)
    Ti = sessionSummaryAll(sessionSummaryAll.Treatment == treatLevels{i}, :);
    if isempty(Ti), continue; end
    plotMeanBars_LiveVsNonLive_3col(Ti, PLOT_SHOW_SIG, "Treatment: " + string(treatLevels{i}));
end

% within-Live: compare treatments (already per-treatment on x-axis)
plotMeanBars_WithinLiveByTreatment_3col(sessionSummaryAll, treatOrder, PLOT_SHOW_SIG);

%% ========================================================================
%% ============================== MODELING ================================
%% ========================================================================

fprintf('\n=== FITTING GLMMs (session-level binomial) ===\n');

ResultsKey = table();

metricFields = { ...
    'pRepeatCond','nTrialsCond', true; ...
    'pSelfRepeat','nTrialsSelf', true; ...
    'pOppPrevMatchSubjPrev','nTrialsPrevMatch', true; ...
    'pRepeatOpp','nTrialsOpp', false; ...
    'pRepeat','nTrialsRaw', true ...
};

for i = 1:size(metricFields,1)
    pField = metricFields{i,1};
    nField = metricFields{i,2};
    doLVN  = metricFields{i,3};
    info   = getMetricInfo(pField);

    ResultsKey = [ResultsKey; runMetricGLMMs(sessionSummaryAll, ...
        pField, nField, info.modelLabel, doLVN, MONKEY_EFFECT_MODE)]; %#ok<AGROW>
end

disp('=== Key hypothesis tests (p-values + signed effect sizes) ===');
disp(ResultsKey);

%% ---------- Per-monkey fits (safe; no forced levels) ----------
disp('=== Key hypothesis tests by monkey (separate fits) ===');
ResultsKey_ByMonkey = table();
monkCats = categories(sessionSummaryAll.Monkey);

for im = 1:numel(monkCats)
    Sm = sessionSummaryAll(sessionSummaryAll.Monkey == monkCats{im}, :);
    if isempty(Sm), continue; end

    for i = 1:size(metricFields,1)
        pField = metricFields{i,1};
        nField = metricFields{i,2};
        doLVN  = metricFields{i,3};
        info   = getMetricInfo(pField);

        ResultsKey_ByMonkey = [ResultsKey_ByMonkey; ...
            runMetricGLMMs(Sm, pField, nField, info.modelLabel + " | Monkey " + string(monkCats{im}), doLVN, "None")]; %#ok<AGROW>
    end
end
disp(ResultsKey_ByMonkey);

%% ---------- Monkey similarity tests (do they share patterns / trends?) ----------
fprintf('\n=== Monkey similarity tests (fixed-effect LRTs) ===\n');
for i = 1:size(metricFields,1)
    pField = metricFields{i,1};
    nField = metricFields{i,2};
    doLVN  = metricFields{i,3};
    info   = getMetricInfo(pField);
    testMonkeySimilarity(sessionSummaryAll, pField, nField, info.modelLabel, doLVN);
end

%% ========================================================================
%% ============================== HELPERS ================================
%% ========================================================================

function info = getMetricInfo(metricField)
% Descriptive names for plots and modeling outputs.
    info = struct();

    switch string(metricField)
        case "pRepeatCond"
            info.short = "Subject repeats opponent (conditional)";
            info.long  = "Probability of the subject matches the opponent's previous choice";
            info.yLab  = "P(subj(t) == opp(t-1) | opp(t-1) != subj(t-1))";
            info.modelLabel = "Subject repeats opponent prev choice (conditional)";
        case "pSelfRepeat"
            info.short = "Subject repeats self (stickiness)";
            info.long  = "Probability the subject repeats its own previous choice";
            info.yLab  = "P(subj(t) == subj(t-1)";
            info.modelLabel = "Subject choice stickiness";
        case "pOppPrevMatchSubjPrev"
            info.short = "Opponent repeats subject (prev trial)";
            info.long  = "P(opp(t-1) == subj(t-1)";
            info.yLab  = "Probability opponent repeats subject (previous trial)";
            info.modelLabel = "Opponent followed on previous trial";
        case "pRepeatOpp"
            info.short = "Opponent repeats subject (Live only)";
            info.long  = "P(opp(t) == subj(t) (only meaningful in Live sessions)";
            info.yLab  = "Probability opponent repeats subject (same trial)";
            info.modelLabel = "Opponent same-trial matching (Live only)";
        case "pRepeat"
            info.short = "Subject matches opponent (no prior, diagnostic)";
            info.long  = "Raw probability the subject matches the opponent's previous choice; can be inflated if the opponent copies the subject";
            info.yLab  = "P(subj(t) == opp(t-1))";
            info.modelLabel = "Subject matches opponent prev choice (no prior)";
        otherwise
            info.short = char(metricField);
            info.long  = char(metricField);
            info.yLab  = "Probability";
            info.modelLabel = string(metricField);
    end
end

function [masks, titles] = monkeyColumnMasks(monkeyCat)
% Always return 3 columns: Monkey 1, Monkey 2, Both (pooled).
% If only 1 monkey is present, column 2 is empty.
    monkeys = categories(monkeyCat);
    titles = cell(1,3);
    masks  = cell(1,3);

    if numel(monkeys) >= 1
        titles{1} = sprintf('Monkey %s', char(monkeys{1}));
        masks{1}  = (monkeyCat == monkeys{1});
    else
        titles{1} = 'Monkey 1';
        masks{1}  = false(size(monkeyCat));
    end

    if numel(monkeys) >= 2
        titles{2} = sprintf('Monkey %s', char(monkeys{2}));
        masks{2}  = (monkeyCat == monkeys{2});
    else
        titles{2} = 'Monkey 2 (no data)';
        masks{2}  = false(size(monkeyCat));
    end

    titles{3} = 'Both monkeys (pooled)';
    masks{3}  = true(size(monkeyCat));
end

function plotPrimaryBarsByTreatment_3col(S, treatOrder, baseCondList, plotMode, doSig)
    info = getMetricInfo("pRepeatCond");
    [condGroups, condNames] = getCondGroups(plotMode);

    for itrt = 1:numel(treatOrder)
        tName = treatOrder{itrt};
        T = S(string(S.Treatment) == string(tName), :);
        if isempty(T), continue; end

        T.Monkey = removecats(categorical(string(T.Monkey)));
        T.ConditionBase = categorical(string(T.ConditionBase), baseCondList);

        [masks, titles] = monkeyColumnMasks(T.Monkey);

        figure('Name', sprintf('Bars | %s | %s', info.short, tName), 'Color','w');

        for col = 1:3
            Tc = T(masks{col}, :);

            subplot(1,3,col); hold on;

            if isempty(Tc)
                title(sprintf('%s\n(no data)', titles{col}), 'Interpreter','none');
                axis off;
                continue;
            end

            nCond = numel(condGroups);
            means = NaN(nCond,1);
            sems  = NaN(nCond,1);
            nSess = zeros(nCond,1);

            for ic = 1:nCond
                rowsC = ismember(string(Tc.ConditionBase), string(condGroups{ic}));
                nSess(ic) = nnz(rowsC);
                vals = Tc.pRepeatCond(rowsC);
                [means(ic), sems(ic)] = meanSEM(vals);
            end

            xtickLbl = cell(nCond,1);
            for ic = 1:nCond
                xtickLbl{ic} = sprintf('%s (%d)', condNames{ic}, nSess(ic));
            end

            bar(1:nCond, means, 0.7);
            errorbar(1:nCond, means, sems, 'k', 'linestyle','none', 'LineWidth',1);

            set(gca,'XTick',1:nCond,'XTickLabel',xtickLbl);
            xtickangle(0);
            ylim([0.25 0.75]);

            xlabel('Session type');
            ylabel(info.yLab);

            title(sprintf('%s\n%s', titles{col}, info.short), 'Interpreter','none');

            % Significance
            if doSig
                ax = gca;
                if plotMode == "Collapsed" && nCond == 2
                    p = pBinomTwoGroup(Tc, 'pRepeatCond','nTrialsCond', 'LiveFlag', 'NonLive', 'Live');
                    if ~isnan(p)
                        y = max(means + sems) + 0.06*diff(ylim(ax));
                        addSigBracket(ax, 1, 2, y, p);
                    end
                else
                    pOmni = pBinomOmnibus(Tc, 'pRepeatCond','nTrialsCond', 'ConditionBase');
                    if ~isnan(pOmni)
                        yl = ylim(ax);
                        text(mean(1:nCond), yl(2)-0.06*diff(yl), sprintf('omnibus p=%.3g', pOmni), ...
                            'HorizontalAlignment','center','FontSize',9);
                    end
                end
            end
        end

        if exist('sgtitle','file')
            if plotMode == "Collapsed"
                sgtitle(sprintf('%s\nTreatment: %s | Comparing non-live vs live sessions', info.long, tName), 'Interpreter','none');
            else
                sgtitle(sprintf('%s\nTreatment: %s | Comparing AI / Replay / Decoy / Live', info.long, tName), 'Interpreter','none');
            end
        end
    end
end

function [condGroups, condNames] = getCondGroups(plotMode)
    if plotMode == "Base4"
        condGroups = { {'AI'}, {'Replay'}, {'Decoy'}, {'Live'} };
        condNames  = { 'AI', 'Replay', 'Decoy', 'Live' };
    elseif plotMode == "Collapsed"
        condGroups = { {'AI','Replay','Decoy'}, {'Live'} };
        condNames  = { 'Non-live', 'Live' };
    else
        error('Unknown PLOT_BAR_COND_MODE: %s', plotMode);
    end
end

function plotNonLiveVsLiveScatter_3col(S, metricField, treatOrder, treatMarker, condFaceColors, extraTitle)
% 2 rows (NonLive, Live) x 3 columns (Monkey1, Monkey2, Both)
    if nargin < 6, extraTitle = ""; end

    info = getMetricInfo(metricField);

    S.Monkey = removecats(categorical(string(S.Monkey)));
    [masks, titles] = monkeyColumnMasks(S.Monkey);

    treatPresent = categories(removecats(categorical(string(S.Treatment))));
    treatLegend  = treatOrder(ismember(treatOrder, treatPresent));

    condPresent = categories(removecats(categorical(string(S.ConditionBase))));
    condOrder = {'AI','Replay','Decoy','Live'};
    condLegend = condOrder(ismember(condOrder, condPresent));

    figure('Name', sprintf('Scatter | %s', info.short), 'Color','w');

    rowDefs = { ...
        'NonLive', "Non-live sessions (AI / Replay / Decoy)"; ...
        'Live',    "Live sessions (real opponent)" ...
    };

    for r = 1:2
        rowFlag = rowDefs{r,1};
        rowTitle = rowDefs{r,2};

        for col = 1:3
            subplot(2,3,(r-1)*3+col); hold on;

            rows = masks{col} & (S.LiveFlag == rowFlag) & ~isnan(S.(metricField)) & ~isnan(S.sessionIdx);
            x = double(S.sessionIdx(rows));
            y = S.(metricField)(rows);
            trt = S.Treatment(rows);
            cnd = S.ConditionBase(rows);

            if isempty(x)
                title(sprintf('%s\n%s\n(no data)', titles{col}, rowTitle), 'Interpreter','none');
                ylim([0 1]);
                continue;
            end

            scatterByTreatmentAndCondition(x, y, trt, cnd, treatLegend, treatMarker, condFaceColors, 40);

            addTrendlineWithStats(x, y);

            ylim([0 1]);
            ylabel(info.yLab);

            if r == 2
                xlabel('Session index (chronological within monkey)');
            end

            title(sprintf('%s\n%s', titles{col}, rowTitle), 'Interpreter','none');

            % Single legend only (top-left subplot)
            if r == 1 && col == 1
                addScatterLegend(treatLegend, treatMarker, condLegend, condFaceColors);
            end
        end
    end

    if exist('sgtitle','file')
        if strlength(extraTitle) > 0
            sgtitle(sprintf('%s\n%s', info.long, extraTitle), 'Interpreter','none');
        else
            sgtitle(info.long, 'Interpreter','none');
        end
    end
end

function plotLiveOnlyOpponentMetrics_3col(S, treatOrder, treatMarker, condFaceColors, extraTitle)
% 2 rows (metrics) x 3 columns (Monkey1, Monkey2, Both)
    if nargin < 5, extraTitle = ""; end

    SL = S(S.LiveFlag == 'Live', :);
    if isempty(SL)
        warning('No Live sessions found; skipping live-only opponent scatter plots.');
        return;
    end

    SL.Monkey = removecats(categorical(string(SL.Monkey)));
    [masks, titles] = monkeyColumnMasks(SL.Monkey);

    treatPresent = categories(removecats(categorical(string(SL.Treatment))));
    treatLegend  = treatOrder(ismember(treatOrder, treatPresent));

    condLegend = {'Live'}; % only Live is present here

    metricList = { 'pRepeatOpp'; 'pOppPrevMatchSubjPrev' };
    figure('Name','Scatter | Opponent coupling (Live only)', 'Color','w');

    for r = 1:2
        mf = metricList{r};
        info = getMetricInfo(mf);

        for col = 1:3
            subplot(2,3,(r-1)*3+col); hold on;

            rows = masks{col} & ~isnan(SL.(mf)) & ~isnan(SL.sessionIdx);
            x = double(SL.sessionIdx(rows));
            y = SL.(mf)(rows);
            trt = SL.Treatment(rows);
            cnd = SL.ConditionBase(rows); % Live only

            if isempty(x)
                title(sprintf('%s\n%s\n(no data)', titles{col}, info.short), 'Interpreter','none');
                ylim([0 1]);
                continue;
            end

            scatterByTreatmentAndCondition(x, y, trt, cnd, treatLegend, treatMarker, condFaceColors, 40);
            addTrendlineWithStats(x, y);

            ylim([0 1]);
            ylabel(info.yLab);
            xlabel('Session index (chronological within monkey)');

            title(sprintf('%s\n%s', titles{col}, info.short), 'Interpreter','none');

            if r == 1 && col == 1
                addScatterLegend(treatLegend, treatMarker, condLegend, condFaceColors);
            end
        end
    end

    if exist('sgtitle','file')
        baseTitle = "Opponent coupling in Live sessions (only Live base condition is present)";
        if strlength(extraTitle) > 0
            sgtitle(sprintf('%s\n%s', baseTitle, extraTitle), 'Interpreter','none');
        else
            sgtitle(baseTitle, 'Interpreter','none');
        end
    end
end

function scatterByTreatmentAndCondition(x, y, trt, cnd, treatLegend, treatMarker, condFaceColors, sz)
% Marker shape = treatment, marker face color = base condition
    x = x(:); y = y(:);
    trtS = string(trt(:));
    cndS = string(cnd(:));

    for it = 1:numel(treatLegend)
        tName = string(treatLegend{it});
        mk = treatMarker(char(tName));

        rowsT = (trtS == tName);
        if ~any(rowsT), continue; end

        condsHere = unique(cndS(rowsT));
        for ic = 1:numel(condsHere)
            cName = string(condsHere(ic));
            rows = rowsT & (cndS == cName);

            xi = x(rows);
            yi = y(rows);

            good = ~isnan(xi) & ~isnan(yi);
            xi = xi(good); yi = yi(good);
            if isempty(xi), continue; end

            fc = [0.7 0.7 0.7];
            if isKey(condFaceColors, char(cName))
                fc = condFaceColors(char(cName));
            end

            scatter(xi, yi, sz, ...
                'Marker', mk, ...
                'MarkerEdgeColor','k', ...
                'MarkerFaceColor', fc, ...
                'LineWidth', 1.0);
        end
    end
end

function addScatterLegend(treatLegend, treatMarker, condLegend, condFaceColors)
% Combined legend: (1) treatment marker shapes and (2) condition face colors.
    hold on;

    h = gobjects(0,1);
    labels = strings(0,1);

    % Treatment (shape)
    for it = 1:numel(treatLegend)
        tName = string(treatLegend{it});
        mk = treatMarker(char(tName));

        hh = scatter(nan, nan, 50, ...
            'Marker', mk, ...
            'MarkerEdgeColor','k', ...
            'MarkerFaceColor','w', ...
            'LineWidth', 1.0);
        h(end+1,1) = hh; %#ok<AGROW>
        labels(end+1,1) = "Treatment: " + tName; %#ok<AGROW>
    end

    % Condition (face color)
    for ic = 1:numel(condLegend)
        cName = string(condLegend{ic});
        fc = [0.7 0.7 0.7];
        if isKey(condFaceColors, char(cName))
            fc = condFaceColors(char(cName));
        end

        hh = scatter(nan, nan, 50, ...
            'Marker', 'o', ...
            'MarkerEdgeColor','k', ...
            'MarkerFaceColor', fc, ...
            'LineWidth', 1.0);
        h(end+1,1) = hh; %#ok<AGROW>
        labels(end+1,1) = "" + cName; %#ok<AGROW>
    end

    legend(h, cellstr(labels), 'Location','northeast', 'Box','off');
end

function plotMeanBars_LiveVsNonLive_3col(S, doSig, treatLabel)
% rows = metrics, cols = Monkey1/Monkey2/Both
% treatLabel: string for figure title (e.g., "baseline", "OT", "Saline", "pooled")
    if nargin < 3 || strlength(string(treatLabel)) == 0
        treatLabel = "pooled over treatments";
    end

    metricDefs = { ...
        'pRepeatCond','nTrialsCond'; ...
        'pSelfRepeat','nTrialsSelf'; ...
        'pOppPrevMatchSubjPrev','nTrialsPrevMatch'; ...
        'pRepeat','nTrialsRaw' ...
    };

    [masks, titles] = monkeyColumnMasks(removecats(categorical(string(S.Monkey))));

    figure('Name', sprintf('Bars | Non-live vs live | %s', treatLabel), 'Color','w');

    for r = 1:size(metricDefs,1)
        pField = metricDefs{r,1};
        nField = metricDefs{r,2};
        info   = getMetricInfo(pField);

        for col = 1:3
            subplot(size(metricDefs,1), 3, (r-1)*3 + col); hold on;

            rowsM = masks{col};

            rowsNL = rowsM & (S.LiveFlag=='NonLive') & ~isnan(S.(pField));
            rowsL  = rowsM & (S.LiveFlag=='Live')    & ~isnan(S.(pField));

            vNL = S.(pField)(rowsNL);
            vL  = S.(pField)(rowsL);

            [mNL, seNL] = meanSEM(vNL);
            [mL,  seL]  = meanSEM(vL);

            bar(1, mNL, 0.7);
            bar(2, mL,  0.7);

            errorbar(1, mNL, seNL, 'k', 'linestyle','none', 'LineWidth',1);
            errorbar(2, mL,  seL,  'k', 'linestyle','none', 'LineWidth',1);

            set(gca,'XTick',[1 2],'XTickLabel',{'Non-live','Live'});
            xtickangle(0);
            ylim([0.25 0.75]);

            if r == size(metricDefs,1)
                xlabel('Session type');
            end
            ylabel('Mean probability');

            title(sprintf('%s\n%s', titles{col}, info.short), 'Interpreter','none');

            if doSig
                p = pBinomTwoGroup(S(rowsM,:), pField, nField, 'LiveFlag','NonLive','Live');
                if ~isnan(p)
                    ax = gca;
                    y  = max([mNL+seNL, mL+seL]) + 0.06*diff(ylim(ax));
                    addSigBracket(ax, 1, 2, y, p);
                end
            end
        end
    end

    if exist('sgtitle','file')
        sgtitle(sprintf('Mean probabilities: non-live vs live (%s)', treatLabel), 'Interpreter','none');
    end
end

function plotMeanBars_WithinLiveByTreatment_3col(S, treatOrder, doSig)
% rows = metrics, cols = Monkey1/Monkey2/Both; Live only
    SL = S(S.LiveFlag=='Live', :);
    if isempty(SL)
        warning('No Live sessions found; skipping within-Live mean-by-treatment bars.');
        return;
    end

    metricDefs = { ...
        'pRepeatCond','nTrialsCond'; ...
        'pRepeatOpp','nTrialsOpp'; ...
        'pSelfRepeat','nTrialsSelf'; ...
        % 'pOppPrevMatchSubjPrev','nTrialsPrevMatch' ...
    };

    [masks, titles] = monkeyColumnMasks(removecats(categorical(string(SL.Monkey))));

    figure('Name','Bars | Within Live: mean probability by treatment', 'Color','w');

    for r = 1:size(metricDefs,1)
        pField = metricDefs{r,1};
        nField = metricDefs{r,2};
        info = getMetricInfo(pField);

        for col = 1:3
            subplot(size(metricDefs,1), 3, (r-1)*3 + col); hold on;

            rowsM = masks{col};
            TL = SL(rowsM, :);

            means = NaN(numel(treatOrder),1);
            sems  = NaN(numel(treatOrder),1);
            for it = 1:numel(treatOrder)
                tname = treatOrder{it};
                rowsT = (string(TL.Treatment) == string(tname)) & ~isnan(TL.(pField));
                v = TL.(pField)(rowsT);
                [means(it), sems(it)] = meanSEM(v);
            end

            bar(1:numel(treatOrder), means, 0.7);
            errorbar(1:numel(treatOrder), means, sems, 'k', 'linestyle','none', 'LineWidth',1);

            set(gca,'XTick',1:numel(treatOrder),'XTickLabel',treatOrder);
            xtickangle(0);
            ylim([0.25 0.75]);

            if r == size(metricDefs,1)
                xlabel('Treatment group');
            end
            ylabel('Mean probability');

            title(sprintf('%s\n%s', titles{col}, info.short), 'Interpreter','none');

            if doSig
                ax = gca;
                p12 = pBinomTwoGroup(TL, pField, nField, 'Treatment', treatOrder{1}, treatOrder{2});
                p13 = pBinomTwoGroup(TL, pField, nField, 'Treatment', treatOrder{1}, treatOrder{3});
                p23 = pBinomTwoGroup(TL, pField, nField, 'Treatment', treatOrder{2}, treatOrder{3});

                yBase = max(means + sems) + 0.06*diff(ylim(ax));
                if ~isnan(p12), addSigBracket(ax, 1, 2, yBase + 0.00*diff(ylim(ax)), p12); end
                if ~isnan(p13), addSigBracket(ax, 1, 3, yBase + 0.06*diff(ylim(ax)), p13); end
                if ~isnan(p23), addSigBracket(ax, 2, 3, yBase + 0.12*diff(ylim(ax)), p23); end
            end
        end
    end

    if exist('sgtitle','file')
        sgtitle('Mean probabilities in Live sessions: comparing treatment groups', 'Interpreter','none');
    end
end

function row = makeResultRow(metricName, testName, pval, estLogOdds, oddsRatio)
    if nargin < 4 || isempty(estLogOdds), estLogOdds = NaN; end
    if nargin < 5 || isempty(oddsRatio),  oddsRatio  = NaN; end
    row = table(string(metricName), string(testName), pval, estLogOdds, oddsRatio, ...
        'VariableNames', {'Metric','Test','pValue','EstLogOdds','OddsRatio'});
end

function [est,se,or] = contrastEffect(mdl, L)
    beta = fixedEffects(mdl);

    CovB = [];
    try
        if isprop(mdl,'CoefficientCovariance') && ~isempty(mdl.CoefficientCovariance)
            CovB = mdl.CoefficientCovariance;
        end
    catch
    end
    if isempty(CovB)
        try
            CovB = mdl.CovarianceMatrix;
        catch
        end
    end
    if isempty(CovB)
        try
            [~, CovB] = fixedEffects(mdl);
        catch
            CovB = NaN(numel(beta));
        end
    end

    est = L * beta;
    se  = sqrt(max(0, L * CovB * L.'));
    or  = exp(est);
end

function monkTerm = makeMonkeyTerm(T, monkeyMode)
    monkTerm = "";
    if monkeyMode == "None"
        return;
    end
    if ~ismember('Monkey', T.Properties.VariableNames)
        return;
    end
    if numel(categories(T.Monkey)) < 2
        return;
    end
    if monkeyMode == "Fixed"
        monkTerm = " + Monkey";
    elseif monkeyMode == "Random"
        monkTerm = " + (1|Monkey)";
    else
        error('Unknown MONKEY_EFFECT_MODE: %s', monkeyMode);
    end
end

function outTbl = runMetricGLMMs(S, pField, nField, metricName, doLiveVsNonLive, monkeyMode)
    outTbl = table();

    T = S(~isnan(S.(pField)) & S.(nField) > 0 & ~isnan(S.sessionIdx), :);
    if isempty(T)
        outTbl = [outTbl; makeResultRow(metricName, "NO DATA", NaN)];
        return;
    end

    T.nSuccess = round(T.(pField) .* T.(nField));
    T.nSuccess = max(0, min(T.(nField), T.nSuccess));

    % Categorical (do not force absent levels)
    T.Monkey   = removecats(categorical(string(T.Monkey)));
    T.Treatment= removecats(categorical(string(T.Treatment)));
    T.LiveFlag = removecats(categorical(string(T.LiveFlag)));

    ordT = {'baseline','OT','Saline'};
    T.Treatment = reordercats(T.Treatment, ordT(ismember(ordT, categories(T.Treatment))));
    ordL = {'NonLive','Live'};
    T.LiveFlag  = reordercats(T.LiveFlag,  ordL(ismember(ordL, categories(T.LiveFlag))));

    % Guards
    if numel(categories(T.Treatment)) < 2
        outTbl = [outTbl; makeResultRow(metricName, "Insufficient Treatment levels", NaN)];
        return;
    end
    if doLiveVsNonLive && numel(categories(T.LiveFlag)) < 2
        outTbl = [outTbl; makeResultRow(metricName, "Insufficient LiveFlag levels", NaN)];
        doLiveVsNonLive = false;
    end

    monkTerm = makeMonkeyTerm(T, monkeyMode);

    % --- Live vs NonLive + diff-in-diff (LiveFlag*Treatment) ---
    if doLiveVsNonLive
        mdl0 = fitglme(T, ...
            sprintf('nSuccess ~ Treatment + sessionIdx%s', monkTerm), ...
            'Distribution','Binomial','Link','logit', ...
            'BinomialSize', T.(nField));

        mdl1 = fitglme(T, ...
            sprintf('nSuccess ~ LiveFlag*Treatment + sessionIdx%s', monkTerm), ...
            'Distribution','Binomial','Link','logit', ...
            'BinomialSize', T.(nField));

        cmp = compare(mdl0, mdl1);
        outTbl = [outTbl; makeResultRow(metricName, "Does Live-vs-non-live differ by treatment? (LiveFlag*Treatment omnibus)", cmp.pValue(2))];

        cn = mdl1.CoefficientNames;

        % Live vs NonLive within baseline (reference)
        if any(strcmp(cn,'LiveFlag_Live'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'LiveFlag_Live')) = 1;
            p = coefTest(mdl1, L, 0);
            [est,~,or] = contrastEffect(mdl1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Live vs non-live (baseline)", p, est, or)];
        end

        % Live vs NonLive within OT
        if any(strcmp(cn,'LiveFlag_Live')) && any(strcmp(cn,'LiveFlag_Live:Treatment_OT'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'LiveFlag_Live')) = 1;
            L(strcmp(cn,'LiveFlag_Live:Treatment_OT')) = 1;
            p = coefTest(mdl1, L, 0);
            [est,~,or] = contrastEffect(mdl1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Live vs non-live (OT)", p, est, or)];
        end

        % Live vs NonLive within Saline
        if any(strcmp(cn,'LiveFlag_Live')) && any(strcmp(cn,'LiveFlag_Live:Treatment_Saline'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'LiveFlag_Live')) = 1;
            L(strcmp(cn,'LiveFlag_Live:Treatment_Saline')) = 1;
            p = coefTest(mdl1, L, 0);
            [est,~,or] = contrastEffect(mdl1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Live vs non-live (Saline)", p, est, or)];
        end

        % Diff-in-diff: (Live-NonLive) OT vs baseline
        if any(strcmp(cn,'LiveFlag_Live:Treatment_OT'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'LiveFlag_Live:Treatment_OT')) = 1;
            p = coefTest(mdl1, L, 0);
            [est,~,or] = contrastEffect(mdl1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Difference-in-differences: OT vs baseline", p, est, or)];
        end

        % Diff-in-diff: (Live-NonLive) Saline vs baseline
        if any(strcmp(cn,'LiveFlag_Live:Treatment_Saline'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'LiveFlag_Live:Treatment_Saline')) = 1;
            p = coefTest(mdl1, L, 0);
            [est,~,or] = contrastEffect(mdl1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Difference-in-differences: Saline vs baseline", p, est, or)];
        end

        % Diff-in-diff: OT vs Saline
        if any(strcmp(cn,'LiveFlag_Live:Treatment_OT')) && any(strcmp(cn,'LiveFlag_Live:Treatment_Saline'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'LiveFlag_Live:Treatment_OT'))     =  1;
            L(strcmp(cn,'LiveFlag_Live:Treatment_Saline')) = -1;
            p = coefTest(mdl1, L, 0);
            [est,~,or] = contrastEffect(mdl1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Difference-in-differences: OT vs Saline", p, est, or)];
        end
    end

    % --- Within Live: Treatment contrasts (optional but useful) ---
    TL = T(string(T.LiveFlag) == "Live", :);
    if isempty(TL)
        outTbl = [outTbl; makeResultRow(metricName, "Within Live: NO DATA", NaN)];
        return;
    end

    monkTermL = makeMonkeyTerm(TL, monkeyMode);

    if numel(categories(TL.Treatment)) >= 2
        mdlL0 = fitglme(TL, ...
            sprintf('nSuccess ~ sessionIdx%s', monkTermL), ...
            'Distribution','Binomial','Link','logit', ...
            'BinomialSize', TL.(nField));

        mdlL1 = fitglme(TL, ...
            sprintf('nSuccess ~ Treatment + sessionIdx%s', monkTermL), ...
            'Distribution','Binomial','Link','logit', ...
            'BinomialSize', TL.(nField));

        cmpL = compare(mdlL0, mdlL1);
        outTbl = [outTbl; makeResultRow(metricName, "Within Live: any treatment effect? (Treatment omnibus)", cmpL.pValue(2))];

        cn = mdlL1.CoefficientNames;

        if any(strcmp(cn,'Treatment_OT'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'Treatment_OT')) = 1;
            p = coefTest(mdlL1, L, 0);
            [est,~,or] = contrastEffect(mdlL1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Within Live: OT vs baseline", p, est, or)];
        end

        if any(strcmp(cn,'Treatment_Saline'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'Treatment_Saline')) = 1;
            p = coefTest(mdlL1, L, 0);
            [est,~,or] = contrastEffect(mdlL1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Within Live: Saline vs baseline", p, est, or)];
        end

        if any(strcmp(cn,'Treatment_OT')) && any(strcmp(cn,'Treatment_Saline'))
            L = zeros(1,numel(cn));
            L(strcmp(cn,'Treatment_OT'))     =  1;
            L(strcmp(cn,'Treatment_Saline')) = -1;
            p = coefTest(mdlL1, L, 0);
            [est,~,or] = contrastEffect(mdlL1, L);
            outTbl = [outTbl; makeResultRow(metricName, "Within Live: OT vs Saline", p, est, or)];
        end
    end
end

function testMonkeySimilarity(S, pField, nField, label, doLiveVsNonLive)
    T = S(~isnan(S.(pField)) & S.(nField)>0 & ~isnan(S.sessionIdx), :);
    if isempty(T)
        fprintf('%s: NO DATA\n', label); return;
    end

    T.nSuccess = round(T.(pField).*T.(nField));
    T.nSuccess = max(0, min(T.(nField), T.nSuccess));

    T.Monkey   = removecats(categorical(string(T.Monkey)));
    T.Treatment= removecats(categorical(string(T.Treatment)));
    T.LiveFlag = removecats(categorical(string(T.LiveFlag)));

    ordT = {'baseline','OT','Saline'};
    T.Treatment = reordercats(T.Treatment, ordT(ismember(ordT, categories(T.Treatment))));
    ordL = {'NonLive','Live'};
    T.LiveFlag  = reordercats(T.LiveFlag,  ordL(ismember(ordL, categories(T.LiveFlag))));

    if numel(categories(T.Monkey)) < 2
        fprintf('%s: only one monkey present\n', label); return;
    end

    if doLiveVsNonLive
        if numel(categories(T.LiveFlag)) < 2 || numel(categories(T.Treatment)) < 2
            fprintf('%s: insufficient LiveFlag/Treatment levels for similarity tests\n', label);
            return;
        end

        % A) Do monkeys share the same session trend?
        mdl_sharedTrend = fitglme(T, ...
            'nSuccess ~ Monkey + LiveFlag*Treatment + sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', T.(nField));

        mdl_diffTrend = fitglme(T, ...
            'nSuccess ~ Monkey + LiveFlag*Treatment + Monkey*sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', T.(nField));

        cmpA = compare(mdl_sharedTrend, mdl_diffTrend);
        fprintf('%s | Do monkeys differ in session trend? (Monkey*sessionIdx) p = %.3g\n', label, cmpA.pValue(2));

        % B) Do monkeys share the same Live-vs-non-live pattern / diff-in-diff?
        mdl_sharedDD = fitglme(T, ...
            'nSuccess ~ Monkey + LiveFlag*Treatment + sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', T.(nField));

        mdl_diffDD = fitglme(T, ...
            'nSuccess ~ Monkey*LiveFlag*Treatment + sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', T.(nField));

        cmpB = compare(mdl_sharedDD, mdl_diffDD);
        fprintf('%s | Do monkeys differ in Live-vs-non-live by treatment? (Monkey*LiveFlag*Treatment) p = %.3g\n', label, cmpB.pValue(2));

    else
        % Live-only metric: test within Live
        TL = T(string(T.LiveFlag)=="Live", :);
        if isempty(TL) || numel(categories(TL.Treatment)) < 2
            fprintf('%s: insufficient Live/Treatment data for similarity tests\n', label);
            return;
        end

        mdl_sharedTrend = fitglme(TL, ...
            'nSuccess ~ Monkey + Treatment + sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', TL.(nField));

        mdl_diffTrend = fitglme(TL, ...
            'nSuccess ~ Monkey + Treatment + Monkey*sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', TL.(nField));

        cmpA = compare(mdl_sharedTrend, mdl_diffTrend);
        fprintf('%s | (Live-only) Do monkeys differ in session trend? (Monkey*sessionIdx) p = %.3g\n', label, cmpA.pValue(2));

        mdl_sharedTrt = fitglme(TL, ...
            'nSuccess ~ Monkey + Treatment + sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', TL.(nField));

        mdl_diffTrt = fitglme(TL, ...
            'nSuccess ~ Monkey*Treatment + sessionIdx', ...
            'Distribution','Binomial','Link','logit', 'BinomialSize', TL.(nField));

        cmpB = compare(mdl_sharedTrt, mdl_diffTrt);
        fprintf('%s | (Live-only) Do monkeys differ in treatment effects? (Monkey*Treatment) p = %.3g\n', label, cmpB.pValue(2));
    end
end

function [m, sem] = meanSEM(vals)
    vals = vals(~isnan(vals));
    n = numel(vals);
    if n == 0
        m = NaN; sem = NaN;
    else
        m = mean(vals);
        if n > 1
            sem = std(vals) / sqrt(n);
        else
            sem = NaN;
        end
    end
end

function addTrendlineWithStats(x_all, y_all)
% Adds pooled linear trendline + annotation: n, slope, p-value, R^2 (fitlm).
    x_all = x_all(:); y_all = y_all(:);
    good = ~isnan(x_all) & ~isnan(y_all);
    x_all = x_all(good); y_all = y_all(good);

    n = numel(x_all);
    if n < 3
        text(0.05, 0.95, sprintf('n=%d', n), 'Units','normalized', ...
            'HorizontalAlignment','left','VerticalAlignment','top', 'FontSize',9);
        return;
    end

    if std(x_all) < eps || std(y_all) < eps
        slope = 0; pval = 1; R2 = 0;
        xfit = linspace(min(x_all), max(x_all)+1, 100);
        yfit = mean(y_all) + 0*xfit;
        plot(xfit, yfit, 'k-', 'LineWidth', 1.5);
    else
        mdl = fitlm(x_all, y_all);
        slope = mdl.Coefficients.Estimate(2);
        pval  = mdl.Coefficients.pValue(2);
        R2    = mdl.Rsquared.Ordinary;

        xfit = linspace(min(x_all), max(x_all), 100);
        yfit = mdl.Coefficients.Estimate(1) + slope*xfit;
        plot(xfit, yfit, 'k-', 'LineWidth', 1.5);
    end

    text(0.05, 0.95, sprintf('n=%d\nslope=%.3f\np=%.3g\nR^2=%.2f', n, slope, pval, R2), ...
        'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 9);
end

function out = mapChoicesToLabels(x)
    if iscell(x); x = x(:); end
    if isnumeric(x)
        out = strings(numel(x),1);
        out(:) = missing;
        out(x==1) = "BUY";
        out(x==2) = "HOLD";
        out(x==3) = "SELL";
    elseif isstring(x)
        out = upper(x(:));
    elseif iscellstr(x)
        out = upper(string(x(:)));
    else
        out = strings(numel(x),1);
        out(:) = missing;
    end
end

%% ============================== SIGNIFICANCE HELPERS ==============================

function p = pBinomTwoGroup(T, pField, nField, groupVar, levelA, levelB)
% Binomial GLM on session-level counts: nSuccess ~ group (2-level)
    p = NaN;

    if isempty(T) || ~ismember(pField, T.Properties.VariableNames) || ~ismember(nField, T.Properties.VariableNames)
        return;
    end
    if ~ismember(groupVar, T.Properties.VariableNames)
        return;
    end

    rows = ~isnan(T.(pField)) & (T.(nField) > 0) & ~ismissing(string(T.(groupVar)));
    T = T(rows,:);
    if isempty(T), return; end

    g = string(T.(groupVar));
    keep = (g==string(levelA)) | (g==string(levelB));
    T = T(keep,:);
    if height(T) < 2, return; end

    g = categorical(string(T.(groupVar)));
    g = removecats(g);
    if numel(categories(g)) ~= 2
        return;
    end
    T.(groupVar) = g;

    T.nSuccess = round(T.(pField).*T.(nField));
    T.nSuccess = max(0, min(T.(nField), T.nSuccess));

    try
        mdl = fitglm(T, sprintf('nSuccess ~ %s', groupVar), ...
            'Distribution','binomial', 'BinomialSize', T.(nField));
        cn = mdl.CoefficientNames;
        idx = find(startsWith(cn, groupVar + "_"), 1);
        if isempty(idx), return; end
        L = zeros(1,numel(cn)); L(idx) = 1;
        p = coefTest(mdl, L, 0);
    catch
        p = NaN;
    end
end

function p = pBinomOmnibus(T, pField, nField, factorVar)
% Deviance test: factor vs intercept-only
    p = NaN;

    if isempty(T) || ~ismember(pField, T.Properties.VariableNames) || ~ismember(nField, T.Properties.VariableNames)
        return;
    end
    if ~ismember(factorVar, T.Properties.VariableNames)
        return;
    end

    rows = ~isnan(T.(pField)) & (T.(nField) > 0) & ~ismissing(string(T.(factorVar)));
    T = T(rows,:);
    if isempty(T), return; end

    f = categorical(string(T.(factorVar)));
    f = removecats(f);
    if numel(categories(f)) < 2
        return;
    end
    T.(factorVar) = f;

    T.nSuccess = round(T.(pField).*T.(nField));
    T.nSuccess = max(0, min(T.(nField), T.nSuccess));

    try
        mdl0 = fitglm(T, 'nSuccess ~ 1', ...
            'Distribution','binomial', 'BinomialSize', T.(nField));
        mdl1 = fitglm(T, sprintf('nSuccess ~ %s', factorVar), ...
            'Distribution','binomial', 'BinomialSize', T.(nField));

        d0 = mdl0.Deviance;
        d1 = mdl1.Deviance;
        df = mdl1.NumCoefficients - mdl0.NumCoefficients;
        if df <= 0, return; end
        p = 1 - chi2cdf(d0 - d1, df);
    catch
        p = NaN;
    end
end

function addSigBracket(ax, x1, x2, y, p)
    if isnan(p) || ~isfinite(p), return; end
    if isempty(ax) || ~isgraphics(ax), ax = gca; end

    stars = pToStars(p);
    if ~contains(stars, '*'), return; end

    yl = ylim(ax);
    h = 0.02*diff(yl);

    yNeed = y + 2*h;
    if yNeed > yl(2)
        ylim(ax, [yl(1) min(1, yNeed + 0.04*diff(yl))]);
        yl = ylim(ax);
        h = 0.02*diff(yl);
    end

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth',1);
    text(ax, mean([x1 x2]), y+h, stars, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',10, 'FontWeight','bold', ...
        'Color','r');
end

function s = pToStars(p)
    if p < 1e-3
        s = '***';
    elseif p < 1e-2
        s = '**';
    elseif p < 5e-2
        s = '*';
    else
        s = 'n.s.';
    end
end
