%% Buy Bias vs Reward Learning Analysis: Phase 0 + Phase 1 + Phase 1b (pBuy–Reward–Training)
%
% PHASE 0:
%   Build a session-level summary table across all treatments, with:
%       - pBuy       : fraction of BUY choices per session
%       - meanFB1    : mean choice reward (fb1) per session
%       - meanFB2    : mean portfolio reward (fb2) per session
%       - meanReward : mean total reward (fb1+fb2) per session
%       - meanRewFB1_* : per-choice mean fb1 (BUY/HOLD/SELL)
%       - meanRewFB2_* : per-choice mean fb2
%       - meanRewTot_* : per-choice mean (fb1+fb2)
%
% PHASE 1b (direct relationship analysis):
%   For each Treatment, make an N×4 grid (N = number of ConditionGroups present):
%       row: ConditionGroup (not_Live and/or Live)
%       col1: pBuy vs sessionIdx (+fit)  [naive training effect]
%       col2: meanReward vs sessionIdx (+fit)  [reward learning with training]
%       col3: pBuy vs meanReward (+fit)  [reward drives buy bias]
%       col4: partial training effect: residual(pBuy|meanReward) vs residual(sessionIdx|meanReward)
%
%   Additionally, for each Treatment create a separate coefficient summary figure with 4 subplots:
%       (A) pBuy ~ sessionIdx
%       (B) meanReward ~ sessionIdx
%       (C) pBuy ~ meanReward
%       (D) pBuy ~ meanReward + sessionIdx
%   where each subplot shows coefficient(s) with 95% CI for each ConditionGroup present.
%
%   Also run session-level LMs (pooled + per monkey):
%       (A) pBuy ~ sessionIdx
%       (B) meanReward ~ sessionIdx
%       (C) pBuy ~ meanReward
%       (D) pBuy ~ meanReward + sessionIdx
%   Diagnostic: if meanReward is significant while sessionIdx is not (Model D), the apparent
%   pBuy increase with training is explained by reward rather than training per se.
%
% PHASE 1 (original; reward-by-choice + fb1/fb2 relationships):
%   For each Treatment, make an N×6 grid (N = number of ConditionGroups present):
%       row: ConditionGroup (not_Live and/or Live)
%       col1: reward by choice  (fb1 & fb2)
%       col2: pBuy vs sessionIdx (+fit)
%       col3: meanFB1 vs sessionIdx (+fit)
%       col4: meanFB2 vs sessionIdx (+fit)
%       col5: meanFB1 vs pBuy (+fit)
%       col6: meanFB2 vs pBuy (+fit)
%
% This script is stand-alone: set OUTDIR correctly.

clear; clc; close all;

%% ---------- CONFIG ----------
% Path to condition packs built from dataTabFULL
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
EXCLUDE_SESSIONS = { ...
    '2018-07-07-01', ... % M1, baseline AI
    '2018-07-20-02', ... % M1, baseline Replay
    '2018-07-21-01', ... % M1, baseline Decoy
    '2018-07-06-02', ... % M2, baseline AI
    '2018-07-06-01', ... % M2, baseline AI
    '2018-08-01-03', ... % M2, baseline Decoy
    '2018-08-22-02', ... % M2, OT Decoy
    '2018-06-07-01', ... % M1, baseline Live
    '2018-06-07-02', ... % M1, baseline Decoy
    '2018-06-08-01', ... % M2, baseline Replay
    '2018-06-08-02', ... % M2, baseline Live
    '2018-06-08-03', ... % M2, baseline Decoy
};

% Condition sets grouped by pharmacological treatment
cond_sets = { ...
    % {'AI','Replay','Decoy','Live'}; ...
    % {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames = {
    % 'baseline', ...
    % 'OT', ...
    'Saline' ...
    };   % treatment labels (same order as cond_sets)

% Canonical base conditions (used for categorical ordering)
baseCondList = {'AI','Replay','Decoy','Live'};

% Choice labels used throughout
choiceLabels = {'BUY','HOLD','SELL'};

% Row order for plots (stacked in one figure)
groupOrder = {'not_Live','Live'};  % row1 = not_Live, row2 = Live

%% ---------- MASTER SESSION SUMMARY ACROSS ALL TREATMENTS ----------
% One row per session / monkey / condition pack:
%   Treatment, ConditionRaw, ConditionBase, ConditionGroup,
%   Monkey, Session, SessionDate,
%   pBuy,
%   meanFB1, meanFB2, meanReward,
%   meanRewFB1_* (BUY/HOLD/SELL),
%   meanRewFB2_* (BUY/HOLD/SELL),
%   meanRewTot_* (BUY/HOLD/SELL),
%   n_* (BUY/HOLD/SELL), nTrials

sessionSummaryAll = table();

%% ================== MAIN LOOP OVER TREATMENT SETS (PHASE 0) ==================
for s = 1:numel(cond_sets)
    conds     = cond_sets{s};      % condition names for this treatment
    treatName = setNames{s};       % 'baseline', 'OT', or 'Saline'

    fprintf('\n=== Treatment: %s ===\n', treatName);

    for c = 1:numel(conds)
        cname   = conds{c};        % e.g. 'AI', 'OT Live', 'Saline Decoy'
        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));
        if ~exist(matFile,'file')
            warning('Missing condition pack: %s', matFile);
            continue;
        end

        tmp = load(matFile, 'C');
        C   = tmp.C;

        % Extract base condition label: AI / Replay / Decoy / Live
        tok = regexp(cname, '(AI|Replay|Decoy|Live)$', 'tokens', 'once');
        if isempty(tok)
            warning('Could not parse base condition from name: %s', cname);
            baseCond = string(cname);
        else
            baseCond = string(tok{1});
        end

        % Define social vs non-social grouping:
        %   Live -> "Live", everything else -> "not_Live"
        if baseCond == "Live"
            condGroup = "Live";
        else
            condGroup = "not_Live";
        end

        % ---------- LOOP OVER FILES (SESSIONS) IN THIS CONDITION ----------
        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            % Must have subject choices and reward components
            neededVars = {'option','fb1','fb2','year','month','day','SessionOfDay','monkey'};
            if any(~ismember(neededVars, etab.Properties.VariableNames))
                warning('Skipping file in %s: missing one of %s', ...
                        cname, strjoin(neededVars, ', '));
                continue;
            end
            if isempty(etab.option{1}) || isempty(etab.fb1{1}) || isempty(etab.fb2{1})
                continue;
            end

            % ---- True session ID (yyyy-mm-dd-ss, from first trial) ----
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

            % ---- Exclusion filter ----
            if ~isempty(EXCLUDE_SESSIONS)
                ex = string(EXCLUDE_SESSIONS(:));
                if any(sessionStr == ex)
                    fprintf('[Exclude] Skipping session %s (%s)\n', sessionStr, cname);
                    continue;
                end
            end

            % ---- Subject choices as string labels (BUY/HOLD/SELL) ----
            subjChoices = mapChoicesToLabels(etab.option{1});  % TR x MKT -> string array
            subjChoices = subjChoices(:);                      % flatten across markets

            % ---- Reward per trial: fb1 (choice), fb2 (portfolio), total ----
            fb1 = etab.fb1{1};                   % numeric TR x MKT
            fb2 = etab.fb2{1};                   % numeric TR x MKT
            R1  = fb1(:);                        % choice reward vector
            R2  = fb2(:);                        % portfolio reward vector
            RT  = R1 + R2;                       % total reward vector

            % Safety: align lengths (in case of minor inconsistencies)
            nT = min([numel(subjChoices), numel(R1), numel(R2), numel(RT)]);
            subjChoices = subjChoices(1:nT);
            R1 = R1(1:nT);
            R2 = R2(1:nT);
            RT = RT(1:nT);

            % Valid trials: well-defined choice and non-NaN total reward
            valid = ~ismissing(subjChoices) & ~isnan(RT);
            if ~any(valid)
                continue;
            end

            % ---- Session-level overall BUY fraction and mean rewards ----
            pBuy       = mean(subjChoices(valid) == "BUY");
            meanFB1    = mean(R1(valid));
            meanFB2    = mean(R2(valid));
            meanReward = mean(RT(valid));   % fb1 + fb2
            nValid     = sum(valid);

            % ---- Reward per choice (BUY / HOLD / SELL) within session ----
            meanFB1_perChoice  = nan(1, 3);
            meanFB2_perChoice  = nan(1, 3);
            meanTot_perChoice  = nan(1, 3);
            n_perChoice        = zeros(1, 3);

            for ic = 1:3
                choiceName = choiceLabels{ic};
                mask = valid & (subjChoices == choiceName);
                n_perChoice(ic) = nnz(mask);
                if n_perChoice(ic) > 0
                    meanFB1_perChoice(ic) = mean(R1(mask));
                    meanFB2_perChoice(ic) = mean(R2(mask));
                    meanTot_perChoice(ic) = mean(RT(mask));
                end
            end

            % ---- Monkey ID (1 or 2) ----
            if isempty(etab.monkey{1})
                warning('Missing monkey field; skipping a file in %s.', cname);
                continue;
            end
            monkeyVal = etab.monkey{1}(1);
            monkeyStr = string(monkeyVal);

            % ---- Append one row to master session summary ----
            newRow = table( ...
                string(treatName), ...        % Treatment
                string(cname),   ...          % ConditionRaw
                baseCond,        ...          % ConditionBase
                condGroup,       ...          % ConditionGroup (Live vs not_Live)
                monkeyStr,       ...          % Monkey
                sessionStr,      ...          % Session
                sessionDate,     ...          % SessionDate
                pBuy,            ...          % pBuy
                meanFB1,         ...          % meanFB1 (choice reward)
                meanFB2,         ...          % meanFB2 (portfolio reward)
                meanReward,      ...          % meanReward (fb1+fb2)
                meanFB1_perChoice(1), ...     % meanRewFB1_BUY
                meanFB1_perChoice(2), ...     % meanRewFB1_HOLD
                meanFB1_perChoice(3), ...     % meanRewFB1_SELL
                meanFB2_perChoice(1), ...     % meanRewFB2_BUY
                meanFB2_perChoice(2), ...     % meanRewFB2_HOLD
                meanFB2_perChoice(3), ...     % meanRewFB2_SELL
                meanTot_perChoice(1), ...     % meanRewTot_BUY
                meanTot_perChoice(2), ...     % meanRewTot_HOLD
                meanTot_perChoice(3), ...     % meanRewTot_SELL
                n_perChoice(1),  ...          % n_BUY
                n_perChoice(2),  ...          % n_HOLD
                n_perChoice(3),  ...          % n_SELL
                nValid,          ...          % nTrials
                'VariableNames', { ...
                    'Treatment','ConditionRaw','ConditionBase','ConditionGroup', ...
                    'Monkey','Session','SessionDate','pBuy', ...
                    'meanFB1','meanFB2','meanReward', ...
                    'meanRewFB1_BUY','meanRewFB1_HOLD','meanRewFB1_SELL', ...
                    'meanRewFB2_BUY','meanRewFB2_HOLD','meanRewFB2_SELL', ...
                    'meanRewTot_BUY','meanRewTot_HOLD','meanRewTot_SELL', ...
                    'n_BUY','n_HOLD','n_SELL','nTrials' ...
                } ...
            );

            sessionSummaryAll = [sessionSummaryAll; newRow]; %#ok<AGROW>
        end
    end
end

%% ---------- CHECK IF WE HAVE ANY DATA ----------
if isempty(sessionSummaryAll)
    error('No sessions found in any treatment. Check OUTDIR and condition packs.');
end

%% ---------- CONVERT TO CATEGORICALS ----------
sessionSummaryAll.Monkey         = categorical(sessionSummaryAll.Monkey);
sessionSummaryAll.ConditionBase  = categorical(sessionSummaryAll.ConditionBase, baseCondList);
sessionSummaryAll.ConditionGroup = categorical(sessionSummaryAll.ConditionGroup, groupOrder);
sessionSummaryAll.Treatment      = categorical(sessionSummaryAll.Treatment, setNames);

%% ---------- GLOBAL SESSION INDEX PER MONKEY ----------
% sessionIdx = training order per monkey across ALL treatments (after exclusions)

sessTab = sessionSummaryAll(:, {'Monkey','Session'});
[sessTabSorted, sortIdx] = sortrows(sessTab, {'Monkey','Session'});

sessionIdxSorted = zeros(height(sessTabSorted),1);
monkeys = categories(sessTabSorted.Monkey);  % list of monkey IDs (e.g. '1','2')

for im = 1:numel(monkeys)
    mID  = monkeys{im};
    mask = sessTabSorted.Monkey == mID;
    nS   = nnz(mask);
    sessionIdxSorted(mask) = (1:nS).';
end

sessionIdx = zeros(height(sessionSummaryAll),1);
sessionIdx(sortIdx) = sessionIdxSorted;
sessionSummaryAll.sessionIdx = sessionIdx;

%% ================== PHASE 1b: DIRECT pBuy–Reward–Training ANALYSIS ==================
colorM1      = [0.0000 0.4470 0.7410]; % Monkey 1
colorM2      = [0.8500 0.3250 0.0980]; % Monkey 2

RUN_ALL_ENV = true;  % pooled across all treatments (still split by ConditionGroup rows)

fprintf('\n\n================== PHASE 1b: pBuy–Reward–Training ==================\n');

for s = 1:numel(setNames)
    treatName = setNames{s};

    % Determine which ConditionGroup rows actually exist for this treatment
    groupsToPlot = {};
    for gg = 1:numel(groupOrder)
        gName = groupOrder{gg};
        if any(sessionSummaryAll.Treatment == treatName & sessionSummaryAll.ConditionGroup == gName)
            groupsToPlot{end+1} = gName; %#ok<AGROW>
        end
    end
    nRows = numel(groupsToPlot);
    if nRows == 0
        fprintf('\n[Info] No sessions for treatment %s (any group). Skipping Phase 1b.\n', treatName);
        continue;
    end

    % Store per-row model fits for coefficient summary figure
    statsA_list = cell(nRows,1);
    statsB_list = cell(nRows,1);
    statsC_list = cell(nRows,1);
    statsD_list = cell(nRows,1);

    figName = sprintf('Phase1b: %s (not_Live vs Live)', treatName);
    figure('Name', figName, 'Color','w');

    for g = 1:nRows
        groupName = groupsToPlot{g};

        rowsEnv = sessionSummaryAll.Treatment == treatName & ...
                  sessionSummaryAll.ConditionGroup == groupName;
        envData = sessionSummaryAll(rowsEnv, :);

        envLabel = sprintf('%s-%s', treatName, groupName);

        x_idx      = double(envData.sessionIdx);
        pBuy       = envData.pBuy;
        meanReward = envData.meanReward;

        fprintf('\n----------------------------------------------\n');
        fprintf('Environment (Phase 1b): %s\n', envLabel);
        fprintf('N sessions: %d\n', height(envData));

        % ---- Model A: pBuy ~ sessionIdx ----
        stats_pBuy_vs_idx = fitLinearModel(x_idx, pBuy, ...
            sprintf('%s: (A) p(BUY) ~ sessionIdx', envLabel), {'sessionIdx'});

        % ---- Model B: meanReward ~ sessionIdx ----
        stats_reward_vs_idx = fitLinearModel(x_idx, meanReward, ...
            sprintf('%s: (B) meanReward ~ sessionIdx', envLabel), {'sessionIdx'});

        % ---- Model C: pBuy ~ meanReward ----
        stats_pBuy_vs_reward = fitLinearModel(meanReward, pBuy, ...
            sprintf('%s: (C) p(BUY) ~ meanReward', envLabel), {'meanReward'});

        % ---- Model D: pBuy ~ meanReward + sessionIdx ----
        X_full = [meanReward(:), x_idx(:)];
        stats_pBuy_vs_reward_idx = fitLinearModel(X_full, pBuy, ...
            sprintf('%s: (D) p(BUY) ~ meanReward + sessionIdx', envLabel), {'meanReward','sessionIdx'});

        % Store fits for later coefficient summary figure
        statsA_list{g} = stats_pBuy_vs_idx;
        statsB_list{g} = stats_reward_vs_idx;
        statsC_list{g} = stats_pBuy_vs_reward;
        statsD_list{g} = stats_pBuy_vs_reward_idx;

        % ---- Quick interpretation printout ----
        printPBuyRewardTrainingInterpretation(envLabel, stats_pBuy_vs_reward_idx);

        % ---- Per-monkey versions of Model D ----
        fprintf('--- Per-monkey Model D for %s ---\n', envLabel);
        for im = 1:numel(monkeys)
            mID = monkeys{im};
            rowsM = envData.Monkey == mID;
            if sum(rowsM) < 3
                fprintf('%s Monkey %s: <3 sessions, skipping Model D.\n', envLabel, char(mID));
                continue;
            end
            Xm = [meanReward(rowsM), x_idx(rowsM)];
            ym = pBuy(rowsM);
            fitLinearModel(Xm, ym, sprintf('%s Monkey %s: (D) p(BUY) ~ meanReward + sessionIdx', envLabel, char(mID)), {'meanReward','sessionIdx'});
        end

        % ---- Partial residual plot for training after controlling meanReward ----
        [res_pBuy, res_idx] = partialResiduals(pBuy, x_idx, meanReward);

        stats_partial = fitLinearModel(res_idx, res_pBuy, ...
            sprintf('%s: partial p(BUY) vs sessionIdx | meanReward', envLabel), {'sessionIdx_resid'});

        % ---------- PLOTS FOR THIS ENVIRONMENT (nRows×4 GRID) ----------
        baseIdx = (g-1)*4;

        % col1: pBuy vs sessionIdx
        subplot(nRows,4, baseIdx + 1);
        plotScatterWithFit(x_idx, pBuy, envData, monkeys, colorM1, colorM2, ...
                           stats_pBuy_vs_idx, ...
                           'sessionIdx (within monkey, global)', ...
                           'p(BUY)', ...
                           sprintf('%s: p(BUY) vs training', groupName), ...
                           true);

        % col2: meanReward vs sessionIdx
        subplot(nRows,4, baseIdx + 2);
        plotScatterWithFit(x_idx, meanReward, envData, monkeys, colorM1, colorM2, ...
                           stats_reward_vs_idx, ...
                           'sessionIdx (within monkey, global)', ...
                           'mean reward (fb1+fb2)', ...
                           sprintf('%s: meanReward vs training', groupName), ...
                           false);

        % col3: pBuy vs meanReward
        subplot(nRows,4, baseIdx + 3);
        plotScatterWithFit(meanReward, pBuy, envData, monkeys, colorM1, colorM2, ...
                           stats_pBuy_vs_reward, ...
                           'mean reward (fb1+fb2)', ...
                           'p(BUY)', ...
                           sprintf('%s: p(BUY) vs meanReward', groupName), ...
                           true);

        % col4: partial training effect (residuals)
        subplot(nRows,4, baseIdx + 4);
        plotScatterWithFit(res_idx, res_pBuy, envData, monkeys, colorM1, colorM2, ...
                           stats_partial, ...
                           'residual(sessionIdx | meanReward)', ...
                           'residual(p(BUY) | meanReward)', ...
                           sprintf('%s: training effect after controlling reward', groupName), ...
                           false);
    end

    if exist('sgtitle','file')
        sgtitle(sprintf('Phase 1b (direct): Treatment %s', treatName), 'Interpreter','none');
    end

    % ---- Coefficient summary figure for this treatment (A, B, C, D) ----
    plotPhase1bCoefficientSummary(sprintf('Phase1b Coefs: %s', treatName), ...
                                  groupsToPlot, statsA_list, statsB_list, statsC_list, statsD_list);
end

%% ---------- (NEW) Combined environment: ALL treatments (still split by ConditionGroup rows) ----------
if RUN_ALL_ENV
    % Determine which groups exist globally
    groupsToPlot = {};
    for gg = 1:numel(groupOrder)
        gName = groupOrder{gg};
        if any(sessionSummaryAll.ConditionGroup == gName)
            groupsToPlot{end+1} = gName; %#ok<AGROW>
        end
    end
    nRows = numel(groupsToPlot);

    if nRows > 0
        % Store per-row model fits for coefficient summary figure
        statsA_list = cell(nRows,1);
        statsB_list = cell(nRows,1);
        statsC_list = cell(nRows,1);
        statsD_list = cell(nRows,1);

        figure('Name', 'Phase1b: ALL treatments (not_Live vs Live)', 'Color','w');

        for g = 1:nRows
            groupName = groupsToPlot{g};
            envData = sessionSummaryAll(sessionSummaryAll.ConditionGroup == groupName, :);
            envLabel = sprintf('ALL-%s', groupName);

            fprintf('\n----------------------------------------------\n');
            fprintf('Environment (Phase 1b): %s\n', envLabel);
            fprintf('N sessions: %d\n', height(envData));

            % Predictors / responses
            x_idx      = double(envData.sessionIdx);
            pBuy       = envData.pBuy;
            meanReward = envData.meanReward;

            % ---- Models for p(BUY) / meanReward ----
            stats_pBuy_vs_idx    = fitLinearModel(x_idx, pBuy, sprintf('%s: (A) p(BUY) ~ sessionIdx', envLabel), {'sessionIdx'});
            stats_reward_vs_idx  = fitLinearModel(x_idx, meanReward, sprintf('%s: (B) meanReward ~ sessionIdx', envLabel), {'sessionIdx'});
            stats_pBuy_vs_reward = fitLinearModel(meanReward, pBuy, sprintf('%s: (C) p(BUY) ~ meanReward', envLabel), {'meanReward'});

            X_full = [meanReward(:), x_idx(:)];
            stats_pBuy_vs_reward_idx = fitLinearModel(X_full, pBuy, sprintf('%s: (D) p(BUY) ~ meanReward + sessionIdx', envLabel), {'meanReward','sessionIdx'});

            % Store fits for later coefficient summary figure
            statsA_list{g} = stats_pBuy_vs_idx;
            statsB_list{g} = stats_reward_vs_idx;
            statsC_list{g} = stats_pBuy_vs_reward;
            statsD_list{g} = stats_pBuy_vs_reward_idx;

            % ---- Quick interpretation printout ----
            printPBuyRewardTrainingInterpretation(envLabel, stats_pBuy_vs_reward_idx);

            % ---- Per-monkey Model D ----
            fprintf('--- Per-monkey Model D for %s ---\n', envLabel);
            for im = 1:numel(monkeys)
                mID = monkeys{im};
                rowsM = envData.Monkey == mID;
                if sum(rowsM) < 3
                    fprintf('%s Monkey %s: <3 sessions, skipping Model D.\n', envLabel, char(mID));
                    continue;
                end
                Xm = [meanReward(rowsM), x_idx(rowsM)];
                ym = pBuy(rowsM);
                fitLinearModel(Xm, ym, sprintf('%s Monkey %s: (D) p(BUY) ~ meanReward + sessionIdx', envLabel, char(mID)), {'meanReward','sessionIdx'});
            end

            % ---- Partial residual plot for training after controlling meanReward ----
            [res_pBuy, res_idx] = partialResiduals(pBuy, x_idx, meanReward);
            stats_partial = fitLinearModel(res_idx, res_pBuy, sprintf('%s: partial p(BUY) vs sessionIdx | meanReward', envLabel), {'sessionIdx_resid'});

            % ---------- PLOTS (nRows×4 GRID) ----------
            baseIdx = (g-1)*4;

            subplot(nRows,4, baseIdx + 1);
            plotScatterWithFit(x_idx, pBuy, envData, monkeys, colorM1, colorM2, ...
                               stats_pBuy_vs_idx, ...
                               'sessionIdx (within monkey, global)', ...
                               'p(BUY)', ...
                               sprintf('%s: p(BUY) vs training', groupName), ...
                               true);

            subplot(nRows,4, baseIdx + 2);
            plotScatterWithFit(x_idx, meanReward, envData, monkeys, colorM1, colorM2, ...
                               stats_reward_vs_idx, ...
                               'sessionIdx (within monkey, global)', ...
                               'mean reward (fb1+fb2)', ...
                               sprintf('%s: meanReward vs training', groupName), ...
                               false);

            subplot(nRows,4, baseIdx + 3);
            plotScatterWithFit(meanReward, pBuy, envData, monkeys, colorM1, colorM2, ...
                               stats_pBuy_vs_reward, ...
                               'mean reward (fb1+fb2)', ...
                               'p(BUY)', ...
                               sprintf('%s: p(BUY) vs meanReward', groupName), ...
                               true);

            subplot(nRows,4, baseIdx + 4);
            plotScatterWithFit(res_idx, res_pBuy, envData, monkeys, colorM1, colorM2, ...
                               stats_partial, ...
                               'residual(sessionIdx | meanReward)', ...
                               'residual(p(BUY) | meanReward)', ...
                               sprintf('%s: training effect after controlling reward', groupName), ...
                               false);
        end

        if exist('sgtitle','file')
            sgtitle('Phase 1b (direct): ALL treatments (rows = not_Live / Live)', 'Interpreter','none');
        end

        % ---- Coefficient summary figure for ALL treatments (A, B, C, D) ----
        plotPhase1bCoefficientSummary('Phase1b Coefs: ALL treatments', ...
                                      groupsToPlot, statsA_list, statsB_list, statsC_list, statsD_list);
    end
end

%% ================== PHASE 1: ORIGINAL N×6 PLOTS PER TREATMENT ==================
for s = 1:numel(setNames)
    treatName = setNames{s};

    % Determine which ConditionGroup rows actually exist for this treatment
    groupsToPlot = {};
    for gg = 1:numel(groupOrder)
        gName = groupOrder{gg};
        if any(sessionSummaryAll.Treatment == treatName & sessionSummaryAll.ConditionGroup == gName)
            groupsToPlot{end+1} = gName; %#ok<AGROW>
        end
    end
    nRows = numel(groupsToPlot);
    if nRows == 0
        fprintf('\n[Info] No sessions for treatment %s (any group). Skipping Phase 1.\n', treatName);
        continue;
    end

    figure('Name', sprintf('Phase1: %s (not_Live vs Live)', treatName), 'Color','w');

    for g = 1:nRows
        groupName = groupsToPlot{g};

        rowsEnv = sessionSummaryAll.Treatment == treatName & ...
                  sessionSummaryAll.ConditionGroup == groupName;
        envData = sessionSummaryAll(rowsEnv, :);

        envLabel = sprintf('%s-%s', treatName, groupName);

        fprintf('\n==============================================\n');
        fprintf('Environment: Treatment = %s, ConditionGroup = %s\n', treatName, groupName);
        fprintf('Number of sessions: %d\n', height(envData));

        % Common handles for predictors and responses
        x_idx    = double(envData.sessionIdx);   % training order per monkey
        pBuy     = envData.pBuy;
        y_FB1    = envData.meanFB1;
        y_FB2    = envData.meanFB2;

        % ---- LMs (pooled) ----
        fprintf('\n--- POOLED LMs for %s ---\n', envLabel);

        stats_pBuy_vs_idx   = fitLinearModel(x_idx, pBuy, ...
            sprintf('%s: p(BUY) ~ sessionIdx (pooled)', envLabel));

        stats_FB1_vs_idx    = fitLinearModel(x_idx, y_FB1, ...
            sprintf('%s: meanFB1 ~ sessionIdx (pooled)', envLabel));

        stats_FB2_vs_idx    = fitLinearModel(x_idx, y_FB2, ...
            sprintf('%s: meanFB2 ~ sessionIdx (pooled)', envLabel));

        stats_FB1_vs_pBuy   = fitLinearModel(pBuy, y_FB1, ...
            sprintf('%s: meanFB1 ~ p(BUY) (pooled)', envLabel));

        stats_FB2_vs_pBuy   = fitLinearModel(pBuy, y_FB2, ...
            sprintf('%s: meanFB2 ~ p(BUY) (pooled)', envLabel));

        % ---- Per-monkey LMs ----
        fprintf('--- Per-monkey LMs for %s ---\n', envLabel);

        % pBuy ~ sessionIdx
        runPerMonkeyLM(envData, monkeys, x_idx, pBuy, envLabel, 'sessionIdx', 'p(BUY)');

        % meanFB1 ~ sessionIdx
        runPerMonkeyLM(envData, monkeys, x_idx, y_FB1, envLabel, 'sessionIdx', 'meanFB1');

        % meanFB2 ~ sessionIdx
        runPerMonkeyLM(envData, monkeys, x_idx, y_FB2, envLabel, 'sessionIdx', 'meanFB2');

        % meanFB1 ~ p(BUY)
        runPerMonkeyLM(envData, monkeys, pBuy, y_FB1, envLabel, 'p(BUY)', 'meanFB1');

        % meanFB2 ~ p(BUY)
        runPerMonkeyLM(envData, monkeys, pBuy, y_FB2, envLabel, 'p(BUY)', 'meanFB2');

        % ---- Aggregated reward by choice (for fb1 and fb2) ----
        [aggFB1, aggSemFB1] = aggregateRewardByChoice(envData, 'meanRewFB1');
        [aggFB2, aggSemFB2] = aggregateRewardByChoice(envData, 'meanRewFB2');

        % ---------- PLOTS FOR THIS ENVIRONMENT (nRows×6 GRID) ----------
        baseIdx = (g-1)*6;

        % col1: reward by choice (fb1 & fb2)
        subplot(nRows,6, baseIdx + 1); hold on;
        barData = [aggFB1(:), aggFB2(:)];   % 3×2: [fb1, fb2] for BUY/HOLD/SELL
        b = bar(1:3, barData, 'grouped');

        % Error bars for both fb1 and fb2
        numGroups = size(barData,1);
        numBars   = size(barData,2);
        xCenters  = nan(numGroups, numBars);
        for ib = 1:numBars
            xCenters(:,ib) = b(ib).XEndPoints.';
        end
        errData = [aggSemFB1(:), aggSemFB2(:)];
        errorbar(xCenters, barData, errData, 'k.', 'LineWidth', 1, ...
                 'HandleVisibility','off');

        set(gca, 'XTick', 1:3, 'XTickLabel', choiceLabels);
        ylabel('Reward');
        title(sprintf('%s: reward by choice', groupName), 'Interpreter','none');
        legend({'fb1 (choice)', 'fb2 (portfolio)'}, 'Location','northeast', 'Box','off');

        % col2: pBuy vs sessionIdx
        subplot(nRows,6, baseIdx + 2);
        plotScatterWithFit(x_idx, pBuy, envData, monkeys, colorM1, colorM2, ...
                           stats_pBuy_vs_idx, ...
                           'sessionIdx (within monkey, global)', ...
                           'p(BUY)', ...
                           sprintf('%s: p(BUY) vs training', groupName), ...
                           true);

        % col3: meanFB1 vs sessionIdx
        subplot(nRows,6, baseIdx + 3);
        plotScatterWithFit(x_idx, y_FB1, envData, monkeys, colorM1, colorM2, ...
                           stats_FB1_vs_idx, ...
                           'sessionIdx (within monkey, global)', ...
                           'mean fb1 (sec)', ...
                           sprintf('%s: mean choice reward vs training', groupName), ...
                           false);

        % col4: meanFB2 vs sessionIdx
        subplot(nRows,6, baseIdx + 4);
        plotScatterWithFit(x_idx, y_FB2, envData, monkeys, colorM1, colorM2, ...
                           stats_FB2_vs_idx, ...
                           'sessionIdx (within monkey, global)', ...
                           'mean fb2 (sec)', ...
                           sprintf('%s: mean portfolio reward vs training', groupName), ...
                           false);

        % col5: meanFB1 vs pBuy
        subplot(nRows,6, baseIdx + 5);
        plotScatterWithFit(pBuy, y_FB1, envData, monkeys, colorM1, colorM2, ...
                           stats_FB1_vs_pBuy, ...
                           'p(BUY)', ...
                           'mean fb1 (sec)', ...
                           sprintf('%s: mean choice reward vs p(BUY)', groupName), ...
                           false);

        % col6: meanFB2 vs pBuy
        subplot(nRows,6, baseIdx + 6);
        plotScatterWithFit(pBuy, y_FB2, envData, monkeys, colorM1, colorM2, ...
                           stats_FB2_vs_pBuy, ...
                           'p(BUY)', ...
                           'mean fb2 (sec)', ...
                           sprintf('%s: mean portfolio reward vs p(BUY)', groupName), ...
                           false);
    end

    if exist('sgtitle','file')
        sgtitle(sprintf('Phase 1 (original): Treatment %s (rows = not_Live / Live)', treatName), 'Interpreter','none');
    end
end

%% Optional: save session summary to disk for later analyses
% save('sessionSummary_buyBias_rewardLearning_phase01b.mat','sessionSummaryAll');

%% ================== HELPER FUNCTIONS ==================

function out = mapChoicesToLabels(x)
% mapChoicesToLabels
%   Convert choice codes to standardized string labels:
%     numeric 1/2/3  -> "BUY"/"HOLD"/"SELL"
%     strings/cellstr -> uppercased
%   Missing / unknown -> missing string.

    if iscell(x)
        x = x(:);
    end

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

function stats = fitLinearModel(x, y, labelStr, predictorNames)
% fitLinearModel
%   Wrapper around fitlm to handle:
%     - scalar predictor (vector x)
%     - 2D predictor matrix x

    stats = [];

    if nargin < 3
        labelStr = '';
    end

    % Ensure column-oriented data and drop NaNs
    if isvector(x)
        X = x(:);
    else
        X = x;
    end
    Yorig = y(:);

    valid = all(~isnan(X),2) & ~isnan(Yorig);
    X = X(valid,:);
    Y = Yorig(valid);

    if numel(Y) < 3 || all(var(X,0,1) == 0)
        fprintf('[fitLinearModel] %s: insufficient data or no variance, skipping.\n', labelStr);
        return;
    end

    try
        if nargin >= 4 && ~isempty(predictorNames)
            predictorNames = cellstr(string(predictorNames));
            if numel(predictorNames) == size(X,2)
                varNames = [predictorNames(:)', {'Y'}];
                mdl = fitlm(X, Y, 'VarNames', varNames);
            else
                mdl = fitlm(X, Y);
            end
        else
            mdl = fitlm(X, Y);
        end
    catch ME
        fprintf('[fitLinearModel] %s: fitlm failed: %s\n', labelStr, ME.message);
        return;
    end

    coeffs    = mdl.Coefficients;
    intercept = coeffs.Estimate(1);
    R2        = mdl.Rsquared.Ordinary;

    if size(X,2) == 1
        slope  = coeffs.Estimate(2);
        pSlope = coeffs.pValue(2);
        fprintf('%s: slope = %.4f, p = %.3g, R^2 = %.3f (N = %d)\n', ...
            labelStr, slope, pSlope, R2, numel(Y));
    else
        slope  = coeffs.Estimate(2:end);
        pSlope = coeffs.pValue(2:end);
        fprintf('%s: slopes = %s, p-values = %s, R^2 = %.3f (N = %d)\n', ...
            labelStr, mat2str(slope,4), mat2str(pSlope,3), R2, numel(Y));
    end

    stats = struct();
    stats.slope     = slope;
    stats.intercept = intercept;
    stats.R2        = R2;
    stats.pSlope    = pSlope;
    stats.model     = mdl;
end

function runPerMonkeyLM(envData, monkeys, xVec, yVec, envLabel, xName, yName)
% runPerMonkeyLM
%   Fit simple 1D linear models per monkey: yVec ~ xVec

    for im = 1:numel(monkeys)
        mID = monkeys{im};
        rowsM = envData.Monkey == mID;
        if sum(rowsM) < 3
            fprintf('%s Monkey %s: <3 sessions, skipping %s ~ %s.\n', ...
                    envLabel, char(mID), yName, xName);
            continue;
        end

        xM = xVec(rowsM);
        yM = yVec(rowsM);

        labelStr = sprintf('%s Monkey %s: %s ~ %s', ...
                           envLabel, char(mID), yName, xName);

        fitLinearModel(xM, yM, labelStr);
    end
end

function [aggMean, aggSem] = aggregateRewardByChoice(envData, meanFieldPrefix)
% aggregateRewardByChoice
%   Compute weighted mean reward per choice (BUY/HOLD/SELL) across sessions.

    choiceSuffixes = {'BUY','HOLD','SELL'};
    aggMean = nan(1,3);
    aggSem  = nan(1,3);

    for ic = 1:3
        meanFieldName  = sprintf('%s_%s', meanFieldPrefix, choiceSuffixes{ic});
        countFieldName = sprintf('n_%s', choiceSuffixes{ic});

        mu = envData.(meanFieldName);
        w  = envData.(countFieldName);

        valid = w > 0 & ~isnan(mu);
        if any(valid)
            aggMean(ic) = sum(w(valid) .* mu(valid)) / sum(w(valid));
            aggSem(ic)  = std(mu(valid)) / sqrt(nnz(valid));
        end
    end
end

function plotScatterWithFit(xVec, yVec, envData, monkeys, colorM1, colorM2, ...
                            stats, xLabel, yLabel, titleStr, yIsProb)
% plotScatterWithFit
%   Generic scatter + pooled linear fit + 95% CI + text annotation.
%   NOTE: Does not open new figures (safe to call from subplots).

    cla; hold on;

    % Scatter points by monkey
    for im = 1:numel(monkeys)
        mID   = monkeys{im};
        rowsM = envData.Monkey == mID;
        if ~any(rowsM), continue; end

        x = xVec(rowsM);
        y = yVec(rowsM);

        if im == 1
            scatter(x, y, 40, 'MarkerFaceColor', colorM1, ...
                'MarkerEdgeColor','k', 'DisplayName', sprintf('M%s', char(mID)));
        else
            scatter(x, y, 40, 'MarkerFaceColor', colorM2, ...
                'MarkerEdgeColor','k', 'DisplayName', sprintf('M%s', char(mID)));
        end
    end

    % Pooled trendline
    txt = '';
    if ~isempty(stats) && isfield(stats,'model') && ~isempty(stats.model) && ~any(isnan(stats.slope(:)))
        xfit = linspace(min(xVec), max(xVec), 100)';

        try
            [yfit, yCI] = predict(stats.model, xfit);
        catch
            yfit = nan(size(xfit));
            yCI  = nan(numel(xfit),2);
        end

        % 95% CI shading (only if valid)
        if all(~isnan(yCI(:)))
            ciX = [xfit; flipud(xfit)];
            ciY = [yCI(:,1); flipud(yCI(:,2))];
            patch(ciX, ciY, [0 0 0], 'FaceAlpha', 0.15, ...
                  'EdgeColor','none', 'HandleVisibility','off');
        end

        % Fit line
        if all(~isnan(yfit))
            plot(xfit, yfit, 'k-', 'LineWidth', 1.5, 'DisplayName','pooled fit');
        end

        % Text annotation (first slope/p-value if multi)
        slopeStr = stats.slope;
        pStr     = stats.pSlope;
        if numel(slopeStr) > 1
            slopeStr = slopeStr(1);
            pStr     = pStr(1);
        end
        txt = sprintf('slope=%.3f\nR^2=%.3f\np=%.3g', ...
            slopeStr, stats.R2, pStr);
    end

    xlabel(xLabel);
    ylabel(yLabel);
    title(titleStr, 'Interpreter','none');
    legend('Location','northwest','Box','off');

    if yIsProb
        ylim([0 1]);
    end

    ax  = gca;
    xl  = xlim(ax);
    yl  = ylim(ax);
    xText = xl(1) + 0.05*(xl(2)-xl(1));
    yText = yl(1) + 0.05*(yl(2)-yl(1));
    if ~isempty(txt)
        text(xText, yText, txt, 'VerticalAlignment','bottom','FontSize',10);
    end
end

function plotPhase1bCoefficientSummary(figTitleStr, groupsToPlot, statsA_list, statsB_list, statsC_list, statsD_list)
% plotPhase1bCoefficientSummary
%   Separate coefficient summary figure with 4 subplots (A,B,C,D).
%   Each subplot shows coefficient(s) with 95% CI for each ConditionGroup present.

    if nargin < 2 || isempty(groupsToPlot)
        return;
    end

    groupLabels = cellstr(string(groupsToPlot));
    nG = numel(groupLabels);

    % Extract CIs for each model / group
    A_est = nan(nG,1); A_lo = nan(nG,1); A_hi = nan(nG,1);
    B_est = nan(nG,1); B_lo = nan(nG,1); B_hi = nan(nG,1);
    C_est = nan(nG,1); C_lo = nan(nG,1); C_hi = nan(nG,1);

    Drew_est = nan(nG,1); Drew_lo = nan(nG,1); Drew_hi = nan(nG,1);
    Didx_est = nan(nG,1); Didx_lo = nan(nG,1); Didx_hi = nan(nG,1);

    for i = 1:nG
        [A_est(i), A_lo(i), A_hi(i)] = getCoefCI(statsA_list{i}, 'sessionIdx');
        [B_est(i), B_lo(i), B_hi(i)] = getCoefCI(statsB_list{i}, 'sessionIdx');
        [C_est(i), C_lo(i), C_hi(i)] = getCoefCI(statsC_list{i}, 'meanReward');

        [Drew_est(i), Drew_lo(i), Drew_hi(i)] = getCoefCI(statsD_list{i}, 'meanReward');
        [Didx_est(i), Didx_lo(i), Didx_hi(i)] = getCoefCI(statsD_list{i}, 'sessionIdx');
    end

    figure('Name', figTitleStr, 'Color','w');

    % (A)
    subplot(2,2,1);
    plotCoefAcrossGroups(gca, A_est, A_lo, A_hi, groupLabels, ...
        '(A) p(BUY) ~ sessionIdx', 'coef(sessionIdx)');

    % (B)
    subplot(2,2,2);
    plotCoefAcrossGroups(gca, B_est, B_lo, B_hi, groupLabels, ...
        '(B) meanReward ~ sessionIdx', 'coef(sessionIdx)');

    % (C)
    subplot(2,2,3);
    plotCoefAcrossGroups(gca, C_est, C_lo, C_hi, groupLabels, ...
        '(C) p(BUY) ~ meanReward', 'coef(meanReward)');

    % (D) two coefficients per group
    subplot(2,2,4);
    ax = gca;
    cla(ax); hold(ax,'on');

    if all(isnan(Drew_est)) && all(isnan(Didx_est))
        axis(ax, 'off');
        text(0.5, 0.5, 'No fit', 'Units','normalized', 'HorizontalAlignment','center');
    else
        x = 1:nG;
        barData = [Drew_est(:), Didx_est(:)]; % nG x 2
        b = bar(x, barData, 'grouped');

        % Error bars for each bar series
        numGroups = size(barData,1);
        numBars   = size(barData,2);
        xCenters  = nan(numGroups, numBars);
        for ib = 1:numBars
            xCenters(:,ib) = b(ib).XEndPoints(:);
        end

        loMat = [Drew_lo(:), Didx_lo(:)];
        hiMat = [Drew_hi(:), Didx_hi(:)];

        for ib = 1:numBars
            estVec = barData(:,ib);
            loVec  = loMat(:,ib);
            hiVec  = hiMat(:,ib);
            mask   = ~isnan(estVec) & ~isnan(loVec) & ~isnan(hiVec);
            if any(mask)
                errorbar(xCenters(mask,ib), estVec(mask), estVec(mask)-loVec(mask), hiVec(mask)-estVec(mask), ...
                         'k.', 'LineWidth', 1, 'HandleVisibility','off');
            end
        end

        if exist('yline','file')
            yline(0, 'k-', 'HandleVisibility','off');
        else
            xl = xlim(ax);
            plot(xl, [0 0], 'k-', 'HandleVisibility','off');
        end

        set(ax, 'XTick', x, 'XTickLabel', groupLabels, 'TickLabelInterpreter','none');
        title(ax, '(D) p(BUY) ~ meanReward + sessionIdx', 'Interpreter','none');
        ylabel(ax, 'Coefficient');
        legend({'meanReward','sessionIdx'}, 'Location','northeast', 'Box','off');

        % Reasonable y-limits based on CIs
        allLo = [Drew_lo(:); Didx_lo(:)];
        allHi = [Drew_hi(:); Didx_hi(:)];
        allLo = allLo(~isnan(allLo));
        allHi = allHi(~isnan(allHi));
        if ~isempty(allLo) && ~isempty(allHi)
            yMin = min(allLo);
            yMax = max(allHi);
            pad  = 0.08 * max(1e-6, (yMax - yMin));
            ylim(ax, [yMin - pad, yMax + pad]);
        end
    end

    if exist('sgtitle','file')
        sgtitle(figTitleStr, 'Interpreter','none');
    end
end

function plotCoefAcrossGroups(ax, est, lo, hi, groupLabels, titleStr, yLabelStr)
% plotCoefAcrossGroups
%   Single-coefficient bar plot across groups with 95% CI.

    axes(ax); %#ok<LAXES>
    cla(ax); hold(ax,'on');

    x = 1:numel(groupLabels);

    if all(isnan(est))
        axis(ax, 'off');
        text(0.5, 0.5, 'No fit', 'Units','normalized', 'HorizontalAlignment','center');
        return;
    end

    bar(x, est);

    errLow  = est - lo;
    errHigh = hi - est;
    mask = ~isnan(est) & ~isnan(errLow) & ~isnan(errHigh);
    if any(mask)
        errorbar(x(mask), est(mask), errLow(mask), errHigh(mask), ...
                 'k.', 'LineWidth', 1, 'HandleVisibility','off');
    end

    if exist('yline','file')
        yline(0, 'k-', 'HandleVisibility','off');
    else
        xl = xlim(ax);
        plot(xl, [0 0], 'k-', 'HandleVisibility','off');
    end

    set(ax, 'XTick', x, 'XTickLabel', groupLabels, 'TickLabelInterpreter','none');
    grid(ax, 'on');
    title(ax, titleStr, 'Interpreter','none');
    ylabel(ax, yLabelStr);

    % Reasonable y-limits based on CIs
    allLo = lo(~isnan(lo));
    allHi = hi(~isnan(hi));
    if ~isempty(allLo) && ~isempty(allHi)
        yMin = min(allLo);
        yMax = max(allHi);
        pad  = 0.08 * max(1e-6, (yMax - yMin));
        ylim(ax, [yMin - pad, yMax + pad]);
    end
end

function [est, lo, hi] = getCoefCI(stats, coefName)
% getCoefCI
%   Return estimate + 95% CI for a named coefficient in stats.model.

    est = nan; lo = nan; hi = nan;

    if isempty(stats) || ~isfield(stats,'model') || isempty(stats.model)
        return;
    end

    mdl = stats.model;

    try
        ciAll = coefCI(mdl); % Nx2
    catch
        try
            ciAll = mdl.coefCI;
        catch
            ciAll = nan(height(mdl.Coefficients), 2);
        end
    end

    namesAll = mdl.CoefficientNames(:);
    idx = find(strcmpi(namesAll, coefName), 1, 'first');
    if isempty(idx)
        return;
    end

    est = mdl.Coefficients.Estimate(idx);
    lo  = ciAll(idx,1);
    hi  = ciAll(idx,2);
end

function [resY, resX] = partialResiduals(y, x, z)
% partialResiduals
%   Frisch–Waugh–Lovell residuals:
%       resY = residuals(y ~ z)
%       resX = residuals(x ~ z)

    y = y(:);
    x = x(:);
    z = z(:);

    mdlY = fitlm(z, y);
    mdlX = fitlm(z, x);

    resY = mdlY.Residuals.Raw;
    resX = mdlX.Residuals.Raw;
end

function printPBuyRewardTrainingInterpretation(envLabel, statsFull)
% printPBuyRewardTrainingInterpretation
%   Diagnostic based on Model D: pBuy ~ meanReward + sessionIdx
%   Assumes predictor order was [meanReward, sessionIdx].

    if isempty(statsFull) || ~isfield(statsFull,'pSlope') || numel(statsFull.pSlope) < 2
        fprintf('[Interpretation] %s: Model D unavailable.\n', envLabel);
        return;
    end

    p_rew = statsFull.pSlope(1);
    p_idx = statsFull.pSlope(2);

    msg = string(sprintf('[Interpretation] %s: Model D p(meanReward)=%.3g, p(sessionIdx)=%.3g.', envLabel, p_rew, p_idx));

    if p_rew < 0.05 && p_idx >= 0.05
        msg = msg + " Reward explains buy bias; training term is not needed.";
    elseif p_rew < 0.05 && p_idx < 0.05
        msg = msg + " Both reward and training contribute (within this environment).";
    elseif p_rew >= 0.05 && p_idx < 0.05
        msg = msg + " Training relates to buy bias even after controlling reward (check robustness).";
    else
        msg = msg + " Neither predictor is significant (check sample size / variance).";
    end

    fprintf('%s\n', msg);
end
