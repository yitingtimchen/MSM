% File: d20251210_fit_glm_behavior_only.m
%% Buy Bias Analysis with Flexible Condition Grouping
% This script fits GLME models for BUY vs non-BUY and visualizes:
%   (1) Mean p(BUY) by condition & monkey (with error bars)
%   (2) Session-level p(BUY) vs training index (sessionIdx) with a single
%       regression line per monkey (no connected data points).
%
% It also prints GLME coefficients for each treatment set.
%
% CONDITION_MODE controls how conditions are collapsed:
%   1 = original 4 base conditions: AI, Replay, Decoy, Live
%   2 = social vs non-social:
%         nonSocial : AI, Replay
%         social    : Decoy, Live
%   3 = Live vs non-Live:
%         Live      : Live
%         not_Live  : AI, Replay, Decoy
%
% NOTE:
%   - Base condition (AI/Replay/Decoy/Live) is always derived from the
%     original condition label and stored in trialData.BaseCond. So you can
%     still compute Live vs non-Live subsets even when CONDITION_MODE is 1 or 2.
%
% REQUIREMENTS:
%   - condition_packs_behavior_only_from_dataTabFULL/*.mat
%   - mapChoicesToLabels() at the end of this file.

clear; clc; close all;

%% ---------------- CONFIG ----------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';

% Condition grouping mode (see header comment)
CONDITION_MODE     = 3;      % 1=4-level, 2=social vs non-social, 3=Live vs non-Live

USE_BINARY_CHOICE  = 1;      % 1 = TARGET_CHOICE vs non-TARGET_CHOICE
TARGET_CHOICE      = "BUY";  % which label is coded as 1 in the binary response
INCLUDE_TIME       = 1;      % include sessionIdx in GLME fixed effects

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
% NOTE: These are matched against trialData.Session which is built from
%       etab.year/month/day + etab.SessionOfDay (NOT etab.Session block).
EXCLUDE_SESSIONS = { ...
    '2018-07-07-01', ... % M1, baseline AI
    '2018-07-20-02', ... % M1, baseline Replay
    '2018-07-21-01', ... % M1, baseline Decoy
    '2018-07-06-02', ... % M2, baseline AI
    '2018-07-06-01', ... % M2, baseline AI
    '2018-08-01-03', ... % M2, baseline Decoy
    '2018-08-22-02', ... % M2, OT Decoy
};
excludeSessionsSet = string(EXCLUDE_SESSIONS(:));

% Condition sets grouped by treatment
cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames = {
    'baseline',...
    'OT',...
    'Saline'...
    };

%% ================== GLME LOOP: BUY vs non-BUY ==================
allResults = struct();

for s = 1:numel(cond_sets)
    conds = cond_sets{s};
    fprintf('\n=== Analyzing Set: %s ===\n', setNames{s});

    % ---------- LOAD CONDITION PACKS ----------
    Cstruct = struct();
    for c = 1:numel(conds)
        cname   = conds{c};
        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));
        if ~exist(matFile,'file')
            warning('Missing file: %s', matFile);
            continue;
        end
        tmp = load(matFile, 'C');
        Cstruct.(matlab.lang.makeValidName(cname)) = tmp.C;
    end
    condNames = fieldnames(Cstruct);
    if isempty(condNames)
        warning('No condition packs loaded for set %s, skipping.', setNames{s});
        continue;
    end

    % ---------- BUILD TRIAL-LEVEL TABLE ----------
    trialData = table();    % master table for this condition set

    for c = 1:numel(condNames)
        cname = condNames{c};
        C     = Cstruct.(cname);

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};
            if ~ismember('option', etab.Properties.VariableNames)
                continue;
            end
            if isempty(etab.option{1})
                continue;
            end

            % ---- true session ID (NOT etab.Session which is block) ----
            y   = etab.year{1}(:);
            m   = etab.month{1}(:);
            d   = etab.day{1}(:);
            sod = etab.SessionOfDay{1}(:);  % session-of-day (1,2,...)

            % yyyy-mm-dd-ss
            date_session = string(compose('%04d-%02d-%02d-%02d', y, m, d, sod));

            % ---- EXCLUDE specified sessions early (session-level filter) ----
            if ~isempty(excludeSessionsSet)
                % Most datasets have one session ID per eventTable; use ismember for safety
                if any(ismember(date_session, excludeSessionsSet))
                    continue;
                end
            end

            % ---- subject choices & binary response ----
            subjChoices = mapChoicesToLabels(etab.option{1});  % string BUY/HOLD/SELL
            valid       = ~ismissing(subjChoices);
            if ~any(valid)
                continue;
            end

            if ~USE_BINARY_CHOICE
                error('Only binary TARGET_CHOICE vs non-TARGET_CHOICE is implemented. Set USE_BINARY_CHOICE = true.');
            end
            y_bin = subjChoices(valid) == TARGET_CHOICE;   % 1 = TARGET_CHOICE (BUY), 0 = others

            % ---- condition label for this file (string like 'AI','OT Live',...) ----
            condLabel = repmat(string(cname), sum(valid), 1);

            % ---- basic covariates ----
            monkey    = string(etab.monkey{1}(:));
            trialNum  = double(etab.trialNum{1}(:));
            market    = string(etab.marketOrig{1}(:));
            bubble    = string(etab.bubbleMarket{1}(:));
            portfolio = double(etab.prePortfolio{1}(:));

            % ---- opponent choice / portfolio (time-aligned) ----
            oppRaw            = mapChoicesToLabels(etab.optionOpp{1});
            opponentChoice    = [missing; oppRaw(1:end-1)];
            opponentPortfolio = etab.postPortfolioOpp{1}(:);
            opponentPortfolio = [missing; opponentPortfolio(1:end-1)];

            % ---- previous choice / reward ----
            prevChoice       = [missing; subjChoices(1:end-1)];
            choiceReward     = etab.postPortfolio{1}(:) .* etab.divPerShare{1}(:);
            prevChoiceReward = [missing; choiceReward(1:end-1)];

            % ---- bubbleness (per-session scalar) ----
            fundamental        = 2.35 * (15 - etab.trialNum{1} + 1);
            relative_deviation = sum(etab.priceBuy{1} - fundamental) ./ sum(fundamental);
            nTrials            = numel(subjChoices);
            bubbleness         = repmat(relative_deviation, nTrials, 1);

            % ---- assemble file-level table ----
            tmpT = table( ...
                y_bin, ...
                condLabel(valid), ...
                string(monkey(valid)), ...
                date_session(valid), ...
                trialNum(valid), ...
                string(market(valid)), ...
                string(bubble(valid)), ...
                string(subjChoices(valid)), ...
                string(opponentChoice(valid)), ...
                string(prevChoice(valid)), ...
                choiceReward(valid), ...
                prevChoiceReward(valid), ...
                portfolio(valid), ...
                opponentPortfolio(valid), ...
                bubbleness(valid), ...
                'VariableNames', { ...
                    'Buy', 'Condition', 'Monkey', 'Session', 'Trial', 'Market', 'Bubble', ...
                    'subjChoice', 'opponentChoice', 'prevChoice', ...
                    'choiceReward', 'prevChoiceReward', ...
                    'portfolio', 'opponentPortfolio', ...
                    'bubbleness' ...
                });

            trialData = [trialData; tmpT]; %#ok<AGROW>
        end
    end

    if isempty(trialData)
        warning('No usable trials for set %s, skipping.', setNames{s});
        continue;
    end

    % ---------- CATEGORICALS & TRUE SESSION FACTOR ----------
    trialData.Monkey  = categorical(trialData.Monkey);
    trialData.Session = categorical(trialData.Session);
    trialData.Market  = categorical(trialData.Market);
    trialData.Bubble  = categorical(trialData.Bubble);

    % ---------- PER-MONKEY SESSION INDEX ----------
    sessTab = table(trialData.Monkey, trialData.Session, ...
                    'VariableNames', {'Monkey','Session'});
    [uniqSess, ~, grpSess] = unique(sessTab, 'rows');

    [uniqSessSorted, sortIdx] = sortrows(uniqSess, {'Monkey','Session'});
    sessionIdxSorted = zeros(height(uniqSessSorted),1);

    monkeys = categories(uniqSessSorted.Monkey);
    for ii = 1:numel(monkeys)
        mID   = monkeys{ii};
        mask  = (uniqSessSorted.Monkey == mID);
        nSess = nnz(mask);
        sessionIdxSorted(mask) = 1:nSess;
    end

    sessionIdxUniq          = zeros(height(uniqSess),1);
    sessionIdxUniq(sortIdx) = sessionIdxSorted;
    trialData.sessionIdx    = sessionIdxUniq(grpSess);

    % ---------- BASE CONDITION (AI / Replay / Decoy / Live) ----------
    % Derive base condition from the original Condition label (e.g. 'OT Live' -> 'Live')
    condOrig = string(trialData.Condition);
    baseCond = strings(size(condOrig));
    baseCond(:) = missing;
    for ii = 1:numel(condOrig)
        tok = regexp(condOrig(ii), '(AI|Replay|Decoy|Live)$', 'tokens', 'once');
        if ~isempty(tok)
            baseCond(ii) = tok{1};
        end
    end
    trialData.BaseCond = categorical(baseCond, {'AI','Replay','Decoy','Live'});
    % Drop any rows where we failed to parse a base condition (should be none)
    trialData = trialData(~ismissing(trialData.BaseCond), :);

    % ---------- CONDITION RECODING (depends on CONDITION_MODE) ----------
    % This controls the Condition factor used in GLME and plots.
    baseStr = string(trialData.BaseCond);
    condStrGroup = strings(size(baseStr));
    condStrGroup(:) = missing;

    switch CONDITION_MODE
        case 1  % original 4 base conditions
            condStrGroup = baseStr;
            condLevels   = {'AI','Replay','Decoy','Live'};

        case 2  % social vs non-social
            isNonSocial = baseStr=="AI" | baseStr=="Replay";
            isSocial    = baseStr=="Decoy" | baseStr=="Live";
            condStrGroup(isNonSocial) = "nonSocial";
            condStrGroup(isSocial)    = "social";
            condLevels = {'nonSocial','social'};

        case 3  % Live vs non-Live
            isLive      = baseStr=="Live";
            isNotLive   = baseStr=="AI" | baseStr=="Replay" | baseStr=="Decoy";
            % isNotLive   = baseStr=="AI" | baseStr=="Replay";
            condStrGroup(isLive)    = "Live";
            condStrGroup(isNotLive) = "not_Live";
            condLevels = {'not_Live','Live'};

        otherwise
            error('Unknown CONDITION_MODE=%d. Use 1, 2, or 3.', CONDITION_MODE);
    end

    trialData.Condition = categorical(condStrGroup, condLevels);
    % Drop any rows with missing Condition (safety)
    trialData = trialData(~ismissing(trialData.Condition), :);

    % ---------- COLLAPSE PREV/OPP CHOICES TO TARGET vs non-TARGET ----------
    if USE_BINARY_CHOICE
        nonName = "non_" + TARGET_CHOICE;

        % Opponent choice
        oc = repmat(nonName, height(trialData), 1);
        isTargetOC = trialData.opponentChoice == TARGET_CHOICE;
        oc(isTargetOC) = TARGET_CHOICE;
        oc(ismissing(trialData.opponentChoice)) = missing;
        trialData.opponentChoice = categorical(oc, [nonName, TARGET_CHOICE]);

        % Previous choice
        pc = repmat(nonName, height(trialData), 1);
        isTargetPC = trialData.prevChoice == TARGET_CHOICE;
        pc(isTargetPC) = TARGET_CHOICE;
        pc(ismissing(trialData.prevChoice)) = missing;
        trialData.prevChoice = categorical(pc, [nonName, TARGET_CHOICE]);
    else
        trialData.opponentChoice = setcats(categorical(trialData.opponentChoice), {'HOLD','BUY','SELL'});
        trialData.prevChoice     = setcats(categorical(trialData.prevChoice),     {'HOLD','BUY','SELL'});
    end

    % Drop any rows with NaN in response (shouldn't happen)
    trialData = trialData(~isnan(trialData.Buy), :);

    % ---------- SPLIT BY MONKEY & BY BASE Live / not_Live ----------
    trialDataM1 = trialData(trialData.Monkey == "1", :);
    trialDataM2 = trialData(trialData.Monkey == "2", :);

    % Live / not_Live are always defined using BaseCond
    trialDataLive      = trialData(trialData.BaseCond == "Live", :);
    trialDataNotLive   = trialData(trialData.BaseCond ~= "Live", :);

    trialDataM1Live    = trialData(trialData.Monkey == "1" & trialData.BaseCond == "Live", :);
    trialDataM2Live    = trialData(trialData.Monkey == "2" & trialData.BaseCond == "Live", :);

    trialDataM1NotLive = trialData(trialData.Monkey == "1" & trialData.BaseCond ~= "Live", :);
    trialDataM2NotLive = trialData(trialData.Monkey == "2" & trialData.BaseCond ~= "Live", :);

    % ================== MODEL LADDER ==================
    % 1) BASIC BETWEEN-CONDITION MODEL:
    %    Buy ~ 1 + Condition + (1|Session) + (1|Monkey)
    %
    % 2) ADD TRAINING (sessionIdx):
    %    Buy ~ 1 + Condition + sessionIdx + (1|Session) + (1|Monkey)
    %
    % 3) FULL MODEL (training + history + state):
    %    Buy ~ 1 + Condition + sessionIdx
    %          + opponentChoice + prevChoice + prevChoiceReward
    %          + portfolio + opponentPortfolio + bubbleness
    %          + (1|Session) + (1|Monkey)
    %
    % Live-only / not_Live-only models drop Condition (constant).

    % ----- ALL MONKEYS, ALL CONDITIONS -----
    formula_basic_all = 'Buy ~ 1 + Condition + (1|Session) + (1|Monkey)';
    if INCLUDE_TIME
        formula_time_all = 'Buy ~ 1 + Condition + sessionIdx + (1|Session) + (1|Monkey)';
    else
        formula_time_all = formula_basic_all;
    end
    formula_full_all  = [ ...
        'Buy ~ 1 + Condition + sessionIdx + ' ...
        'opponentChoice + prevChoice + prevChoiceReward + ' ...
        'portfolio + opponentPortfolio + bubbleness + ' ...
        '(1|Session) + (1|Monkey)' ...
    ];

    % ----- SINGLE MONKEY, ALL CONDITIONS -----
    formula_basic_M = 'Buy ~ 1 + Condition + (1|Session)';
    if INCLUDE_TIME
        formula_time_M  = 'Buy ~ 1 + Condition + sessionIdx + (1|Session)';
    else
        formula_time_M  = formula_basic_M;
    end
    formula_full_M  = [ ...
        'Buy ~ 1 + Condition + sessionIdx + ' ...
        'opponentChoice + prevChoice + prevChoiceReward + ' ...
        'portfolio + opponentPortfolio + bubbleness + ' ...
        '(1|Session)' ...
    ];

    % ----- SINGLE-CONDITION MODELS (Live or not_Live) -----
    formula_basic_single_all = 'Buy ~ 1 + (1|Session) + (1|Monkey)';
    if INCLUDE_TIME
        formula_time_single_all  = 'Buy ~ 1 + sessionIdx + (1|Session) + (1|Monkey)';
    else
        formula_time_single_all  = formula_basic_single_all;
    end
    formula_full_single_all  = [ ...
        'Buy ~ 1 + sessionIdx + ' ...
        'opponentChoice + prevChoice + prevChoiceReward + ' ...
        'portfolio + opponentPortfolio + bubbleness + ' ...
        '(1|Session) + (1|Monkey)' ...
    ];

    formula_basic_single_M = 'Buy ~ 1 + (1|Session)';
    if INCLUDE_TIME
        formula_time_single_M  = 'Buy ~ 1 + sessionIdx + (1|Session)';
    else
        formula_time_single_M  = formula_basic_single_M;
    end
    formula_full_single_M  = [ ...
        'Buy ~ 1 + sessionIdx + ' ...
        'opponentChoice + prevChoice + prevChoiceReward + ' ...
        'portfolio + opponentPortfolio + bubbleness + ' ...
        '(1|Session)' ...
    ];

    mdl        = struct();
    mdl_live   = struct();
    mdl_nlive  = struct();

    % ---------- ALL MONKEYS: ALL CONDITIONS ----------
    fprintf('\n--- %s: ALL MONKEYS, ALL CONDITIONS ---\n', setNames{s});

    mdl.basic_all = fitglme(trialData, formula_basic_all, ...
        'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    printCoefs(sprintf('%s: basic_all (all monkeys, all conds)', setNames{s}), mdl.basic_all);

    mdl.time_all = fitglme(trialData, formula_time_all, ...
        'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    printCoefs(sprintf('%s: time_all (all monkeys, all conds)', setNames{s}), mdl.time_all);

    mdl.full_all = fitglme(trialData, formula_full_all, ...
        'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    printCoefs(sprintf('%s: full_all (all monkeys, all conds)', setNames{s}), mdl.full_all);

    % ---------- ALL MONKEYS: Live ONLY (BaseCond == Live) ----------
    if ~isempty(trialDataLive)
        fprintf('\n--- %s: ALL MONKEYS, Live ONLY ---\n', setNames{s});

        mdl_live.basic_all = fitglme(trialDataLive, formula_basic_single_all, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_all (all monkeys, Live only)', setNames{s}), mdl_live.basic_all);

        mdl_live.time_all = fitglme(trialDataLive, formula_time_single_all, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_all (all monkeys, Live only)', setNames{s}), mdl_live.time_all);

        mdl_live.full_all = fitglme(trialDataLive, formula_full_single_all, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_all (all monkeys, Live only)', setNames{s}), mdl_live.full_all);
    else
        mdl_live.basic_all = [];
        mdl_live.time_all  = [];
        mdl_live.full_all  = [];
        fprintf('No Live trials for ALL monkeys in set %s.\n', setNames{s});
    end

    % ---------- ALL MONKEYS: not_Live ONLY (BaseCond ~= Live) ----------
    if ~isempty(trialDataNotLive)
        fprintf('\n--- %s: ALL MONKEYS, not_Live ONLY ---\n', setNames{s});

        mdl_nlive.basic_all = fitglme(trialDataNotLive, formula_basic_single_all, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_all (all monkeys, not_Live only)', setNames{s}), mdl_nlive.basic_all);

        mdl_nlive.time_all = fitglme(trialDataNotLive, formula_time_single_all, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_all (all monkeys, not_Live only)', setNames{s}), mdl_nlive.time_all);

        mdl_nlive.full_all = fitglme(trialDataNotLive, formula_full_single_all, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_all (all monkeys, not_Live only)', setNames{s}), mdl_nlive.full_all);
    else
        mdl_nlive.basic_all = [];
        mdl_nlive.time_all  = [];
        mdl_nlive.full_all  = [];
        fprintf('No not_Live trials for ALL monkeys in set %s.\n', setNames{s});
    end

    % ---------- MONKEY 1: ALL CONDITIONS ----------
    if ~isempty(trialDataM1)
        fprintf('\n--- %s: MONKEY 1, ALL CONDITIONS ---\n', setNames{s});

        mdl.basic_M1 = fitglme(trialDataM1, formula_basic_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_M1 (all conds)', setNames{s}), mdl.basic_M1);

        mdl.time_M1 = fitglme(trialDataM1, formula_time_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_M1 (all conds)', setNames{s}), mdl.time_M1);

        mdl.full_M1 = fitglme(trialDataM1, formula_full_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_M1 (all conds)', setNames{s}), mdl.full_M1);
    else
        mdl.basic_M1 = [];
        mdl.time_M1  = [];
        mdl.full_M1  = [];
        fprintf('No data for Monkey 1 in set %s.\n', setNames{s});
    end

    % ---------- MONKEY 1: Live ONLY ----------
    if ~isempty(trialDataM1Live)
        fprintf('\n--- %s: MONKEY 1, Live ONLY ---\n', setNames{s});

        mdl_live.basic_M1 = fitglme(trialDataM1Live, formula_basic_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_M1 (Live only)', setNames{s}), mdl_live.basic_M1);

        mdl_live.time_M1 = fitglme(trialDataM1Live, formula_time_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_M1 (Live only)', setNames{s}), mdl_live.time_M1);

        mdl_live.full_M1 = fitglme(trialDataM1Live, formula_full_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_M1 (Live only)', setNames{s}), mdl_live.full_M1);
    else
        mdl_live.basic_M1 = [];
        mdl_live.time_M1  = [];
        mdl_live.full_M1  = [];
        fprintf('No Live trials for Monkey 1 in set %s.\n', setNames{s});
    end

    % ---------- MONKEY 1: not_Live ONLY ----------
    if ~isempty(trialDataM1NotLive)
        fprintf('\n--- %s: MONKEY 1, not_Live ONLY ---\n', setNames{s});

        mdl_nlive.basic_M1 = fitglme(trialDataM1NotLive, formula_basic_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_M1 (not_Live only)', setNames{s}), mdl_nlive.basic_M1);

        mdl_nlive.time_M1 = fitglme(trialDataM1NotLive, formula_time_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_M1 (not_Live only)', setNames{s}), mdl_nlive.time_M1);

        mdl_nlive.full_M1 = fitglme(trialDataM1NotLive, formula_full_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_M1 (not_Live only)', setNames{s}), mdl_nlive.full_M1);
    else
        mdl_nlive.basic_M1 = [];
        mdl_nlive.time_M1  = [];
        mdl_nlive.full_M1  = [];
        fprintf('No not_Live trials for Monkey 1 in set %s.\n', setNames{s});
    end

    % ---------- MONKEY 2: ALL CONDITIONS ----------
    if ~isempty(trialDataM2)
        fprintf('\n--- %s: MONKEY 2, ALL CONDITIONS ---\n', setNames{s});

        mdl.basic_M2 = fitglme(trialDataM2, formula_basic_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_M2 (all conds)', setNames{s}), mdl.basic_M2);

        mdl.time_M2 = fitglme(trialDataM2, formula_time_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_M2 (all conds)', setNames{s}), mdl.time_M2);

        mdl.full_M2 = fitglme(trialDataM2, formula_full_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_M2 (all conds)', setNames{s}), mdl.full_M2);
    else
        mdl.basic_M2 = [];
        mdl.time_M2  = [];
        mdl.full_M2  = [];
        fprintf('No data for Monkey 2 in set %s.\n', setNames{s});
    end

    % ---------- MONKEY 2: Live ONLY ----------
    if ~isempty(trialDataM2Live)
        fprintf('\n--- %s: MONKEY 2, Live ONLY ---\n', setNames{s});

        mdl_live.basic_M2 = fitglme(trialDataM2Live, formula_basic_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_M2 (Live only)', setNames{s}), mdl_live.basic_M2);

        mdl_live.time_M2 = fitglme(trialDataM2Live, formula_time_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_M2 (Live only)', setNames{s}), mdl_live.time_M2);

        mdl_live.full_M2 = fitglme(trialDataM2Live, formula_full_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_M2 (Live only)', setNames{s}), mdl_live.full_M2);
    else
        mdl_live.basic_M2 = [];
        mdl_live.time_M2  = [];
        mdl_live.full_M2  = [];
        fprintf('No Live trials for Monkey 2 in set %s.\n', setNames{s});
    end

    % ---------- MONKEY 2: not_Live ONLY ----------
    if ~isempty(trialDataM2NotLive)
        fprintf('\n--- %s: MONKEY 2, not_Live ONLY ---\n', setNames{s});

        mdl_nlive.basic_M2 = fitglme(trialDataM2NotLive, formula_basic_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: basic_M2 (not_Live only)', setNames{s}), mdl_nlive.basic_M2);

        mdl_nlive.time_M2 = fitglme(trialDataM2NotLive, formula_time_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: time_M2 (not_Live only)', setNames{s}), mdl_nlive.time_M2);

        mdl_nlive.full_M2 = fitglme(trialDataM2NotLive, formula_full_single_M, ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        printCoefs(sprintf('%s: full_M2 (not_Live only)', setNames{s}), mdl_nlive.full_M2);
    else
        mdl_nlive.basic_M2 = [];
        mdl_nlive.time_M2  = [];
        mdl_nlive.full_M2  = [];
        fprintf('No not_Live trials for Monkey 2 in set %s.\n', setNames{s});
    end

    % ---------- STORE RESULTS ----------
    allResults.(setNames{s}).trialData      = trialData;
    allResults.(setNames{s}).GLME          = mdl;
    allResults.(setNames{s}).GLME_Live     = mdl_live;
    allResults.(setNames{s}).GLME_notLive  = mdl_nlive;
end

%% ============================================================
%  VISUALIZATIONS
%  Assumes:
%    - allResults.(setName).trialData     exists
%    - allResults.(setName).GLME.basic_all, time_all
%    - Condition is whatever you defined via CONDITION_MODE
%    - Buy is binary (0/1), TARGET_CHOICE = "BUY"
% ============================================================

setNamesResults = fieldnames(allResults);

for s = 1:numel(setNamesResults)
    setName  = setNamesResults{s};
    fprintf('\n=== Visualization for set: %s ===\n', setName);

    trialData = allResults.(setName).trialData;
    mdl_basic = allResults.(setName).GLME.basic_all;
    mdl_time  = allResults.(setName).GLME.time_all;

    if isempty(trialData) || isempty(mdl_basic)
        warning('Missing data or models for set %s, skipping.', setName);
        continue;
    end

    %% ---------- 1) MEAN p(BUY) BY CONDITION & MONKEY (with error bars) ----------
    % Compute session-level p(BUY) first, then mean & SEM across sessions
    gSess = groupsummary(trialData, {'Monkey','Condition','Session'}, 'mean', 'Buy');
    gSess.pBuy = gSess.mean_Buy;

    monkeys = categories(gSess.Monkey);
    nM      = numel(monkeys);
    condCats = categories(gSess.Condition);
    nC       = numel(condCats);

    figure('Name', sprintf('%s: p(BUY) by Condition & Monkey', setName), ...
           'Color','w');

    for iM = 1:nM
        subplot(1,nM,iM); hold on;
        mID = monkeys{iM};
        sub = gSess(gSess.Monkey == mID, :);

        meanVals = nan(nC,1);
        semVals  = nan(nC,1);
        for ic = 1:nC
            cond = condCats{ic};
            rows = sub.Condition == cond;
            vals = sub.pBuy(rows);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                meanVals(ic) = mean(vals);
                if numel(vals) > 1
                    semVals(ic) = std(vals) / sqrt(numel(vals));
                else
                    semVals(ic) = NaN;
                end
            end
        end

        x = 1:nC;
        bar(x, meanVals);
        errorbar(x, meanVals, semVals, ...
                 'k', 'LineStyle','none', 'LineWidth',1, 'HandleVisibility','off');

        set(gca,'XTick',x,'XTickLabel',cellstr(condCats), ...
            'TickLabelInterpreter','none');
        ylim([0 1]);
        ylabel('p(BUY)', 'Interpreter','none');
        title(sprintf('Monkey %s', char(mID)), 'Interpreter','none');
        box off;
    end
    sgtitle(sprintf('%s: Mean p(BUY) by Condition & Monkey', setName), 'Interpreter','none');


    %% ---------- 2) SESSION-LEVEL p(BUY) vs TRAINING (sessionIdx) ----------
    % Aggregate to session level: mean BUY per (Monkey, Condition, Session, sessionIdx)
    g2 = groupsummary(trialData, {'Monkey','Condition','Session','sessionIdx'}, 'mean', 'Buy');
    g2.pBuy = g2.mean_Buy;

    figure('Name', sprintf('%s: p(BUY) vs sessionIdx', setName), ...
           'Color','w');

    for iM = 1:nM
        subplot(1,nM,iM); hold on;
        mID = monkeys{iM};
        subM = g2(g2.Monkey == mID, :);

        condCatsLocal = categories(subM.Condition);
        nC_local      = numel(condCatsLocal);

        x_all = [];  % for regression
        y_all = [];

        for iC = 1:nC_local
            cond = condCatsLocal{iC};
            subMC = subM(subM.Condition == cond, :);
            if isempty(subMC), continue; end

            x = double(subMC.sessionIdx);
            y = subMC.pBuy;

            scatter(x, y, 40, 'filled', 'DisplayName', cond);
            x_all = [x_all; x(:)];
            y_all = [y_all; y(:)];
        end

        % Regression line pooling all conditions for this monkey
        if numel(x_all) > 1
            lm = fitlm(x_all, y_all);
            xfit = linspace(min(x_all), max(x_all), 100)';
            yfit = predict(lm, xfit);
            plot(xfit, yfit, 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');

            slope = lm.Coefficients.Estimate(2);
            pval  = lm.Coefficients.pValue(2);
            R2    = lm.Rsquared.Ordinary;

            txt = sprintf('slope = %.4f\nR^2 = %.3f\np = %.3g', slope, R2, pval);
            text(0.02, 0.98, txt, ...
                'Units','normalized', ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','top', ...
                'Interpreter','none');
        end

        xlabel('sessionIdx (within monkey)', 'Interpreter','none');
        ylabel('p(BUY)', 'Interpreter','none');
        ylim([0 1]);
        title(sprintf('Monkey %s', char(mID)), 'Interpreter','none');
        legend('Location','northeast', 'Interpreter','none');
        box off;
        set(gca,'TickLabelInterpreter','none');
    end
    sgtitle(sprintf('%s: Session-level p(BUY) vs Training', setName), 'Interpreter','none');


    %% ---------- 3) GLME FIXED EFFECTS: CONDITION COEFFICIENTS ----------
    % Visualize Condition coefficients from basic_all model
    coefTbl = mdl_basic.Coefficients;

    % pick out Condition_* rows (effect of each non-baseline level)
    condMask = startsWith(coefTbl.Name, 'Condition_');
    condRows = coefTbl(condMask, :);

    if ~isempty(condRows)
        condLabels = erase(condRows.Name, 'Condition_');

        est  = condRows.Estimate;
        se   = condRows.SE;
        ciLo = est - 1.96*se;
        ciHi = est + 1.96*se;

        figure('Name', sprintf('%s: GLME Condition Effects (basic_all)', setName), ...
               'Color','w');
        hold on;
        nb = numel(est);
        bar(1:nb, est);
        errorbar(1:nb, est, est - ciLo, ciHi - est, ...
                 'k', 'LineStyle','none', 'LineWidth', 1.5);
        yline(0,'k--');
        set(gca,'XTick',1:nb,'XTickLabel',cellstr(condLabels), ...
            'TickLabelInterpreter','none');
        ylabel('Log-odds effect on p(BUY)', 'Interpreter','none');
        title(sprintf('%s: Condition fixed effects (baseline = %s)', ...
              setName, coefTbl.Name{strcmp(coefTbl.Name,'(Intercept)')}), ...
              'Interpreter','none');
        box off;
    end

    %% ---------- 4) GLME TIME EFFECT: sessionIdx COEFFICIENT ----------
    coefTbl_time = mdl_time.Coefficients;
    rowTime = coefTbl_time(strcmp(coefTbl_time.Name,'sessionIdx'), :);
    if ~isempty(rowTime)
        est_t  = rowTime.Estimate;
        se_t   = rowTime.SE;
        ciLo_t = est_t - 1.96*se_t;
        ciHi_t = est_t + 1.96*se_t;

        figure('Name', sprintf('%s: GLME Time Effect (sessionIdx)', setName), ...
               'Color','w');
        hold on;
        bar(1, est_t);
        errorbar(1, est_t, est_t - ciLo_t, ciHi_t - est_t, ...
                 'k', 'LineStyle','none', 'LineWidth',1.5);
        yline(0,'k--');
        set(gca,'XTick',1,'XTickLabel',{'sessionIdx'}, 'TickLabelInterpreter','none');
        ylabel('Log-odds per unit sessionIdx', 'Interpreter','none');
        title(sprintf('%s: Training effect (basic_all + sessionIdx)', setName), ...
              'Interpreter','none');
        box off;
    end
end

%% ---------- Helpers ----------

function printCoefs(label, mdlObj)
    if isempty(mdlObj), return; end
    fprintf('\nCoefficients: %s\n', label);
    disp(mdlObj.Coefficients);
end

function out = mapChoicesToLabels(x)
% Returns a string array of same length with values "BUY","HOLD","SELL" (or missing)
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
        out = strings(numel(x),1); out(:) = missing;
    end
end
