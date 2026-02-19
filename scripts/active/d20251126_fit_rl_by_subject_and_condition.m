% File: d20251126_fit_rl_by_subject_and_condition.m
%% Buy Bias Analysis: RL-based Behavior Modeling
% =========================================================================
% High-level pipeline:
%   PART 0  : Build trial tables and neuron structs from condition packs.
%   PART 1  : Behavior-only RL model fitting (event-agnostic).
%   PART 1b : Block-structured parametric bootstrap for RL parameters.
%   PART 1c : Recompute RL latents (event-agnostic) and attach to all events.
%   PART 2  : Per-neuron GLM fits (with RL latents and task covariates).
%   PART 2a : Optional coefficient summary plots.
%   PART 2b : Coefficient-based clustering of neurons.
%
% IMPORTANT DESIGN CHOICES:
%   - TARGET_EVTS is ONLY used where event timing matters
%     (firing-rate windows, FR extraction, GLMs, clustering).
%   - RL fitting and RL latents are now event-agnostic:
%       * A single "canonical" event (TARGET_EVTS{1}) is used to build the
%         behavior table for RL.
%       * RL parameters are fit per (Set × Monkey × Condition × Session),
%         independent of which event is being analyzed.
%       * Once RL latents are computed for this canonical behavior table,
%         they are copied into every event-specific trial table for GLM.
%   - A simpler RL model 'SelfOnly' (no BUY-bias kappa) is added
%     alongside 'SelfOnlyBuyBias', 'SelfBubble', and 'SelfOpp'.
%   - Progress messages (fprintf) are sprinkled throughout long loops
%     so you can monitor run progress in the MATLAB console.
% =========================================================================

clear; clc; 
% close all;
warning('off');  % suppress noisy warnings (fitlm rank deficiency, etc.)

%% ======================= PATHS & CONDITION SETS ==========================
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';

% Condition sets to analyze (each cell-array is a "set")
cond_sets = {
    {'AI','Replay','Decoy','Live'}, ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}, ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'}
    };

setNames = {
    'baseline', ...
    'OT', ...
    'Saline'
    };

% Target events and windows for event-based analysis (ms).
% These affect only FR extraction and event-specific GLMs / clustering.
TARGET_EVTS = {
    "m1rew1on","m1rew1off","m1rew2on","m1rew2off", ...
    "m2rew1on","m2rew1off","m2rew2on","m2rew2off"
    };

TARGET_EVT_WINS = {
      [-250 250];   % m1rew1on
      [-250 250];   % m1rew1off
      [-250 250];   % m1rew2on
      [-250 250];   % m1rew2off/m1trialStop
      [-250 250];   % m2rew1on
      [-250 250];   % m2rew1off
      [-250 250];   % m2rew2on
      [-250 250]    % m2rew2off/m2trialStop
    };

% Baseline window for z-scoring firing rate (ms)
BASELINE_EVT = 'm1rew1on';
BASELINE_WIN = [-1000, -700];

%% =========================== HIGH-LEVEL TOGGLES ==========================
% Collapse all non-Live conditions into "not_Live"
LIVE_NON_LIVE_ONLY = 0;
% NOTE: "Live" is detected by substring, so 'OT Live' / 'Saline Live' etc.
%       are treated as "Live"; everything else as "not_Live" when enabled.

% Collapse HOLD/SELL into "Non-BUY" (binary choice space)
BUY_NON_BUY_ONLY = 0;

% FIT_SCOPE: 'per_monkey' or 'pooled'
%   - 'per_monkey' : one RL parameter vector per (Monkey × Condition × Session)
%   - 'pooled'     : one RL parameter vector per (Condition × Session)
%                    (Monkey is stored as 'ALL')
FIT_SCOPE = 'per_monkey';

% select which models to run
MODEL_TO_RUN = {'SelfOnly','SelfOnlyBuyBias','OppOnly','OppOnlyBuyBias','SelfOpp','SelfBubble'};

% RL model used to generate latents for Q/DeltaQ time series and
% terminal-Q analyses (and optionally any FR-based GLMs elsewhere):
%   - Fixed model name:
%       'SelfOnlyBuyBias' | 'SelfOnly' | 'SelfBubble' | 'SelfOpp'
%   - Automatic: 'best_by_AIC' | 'best_by_BIC' (per group, per set)
RL_MODEL_FOR_Q_TIMESERIES = 'SelfOnlyBuyBias';

% Bootstrap settings
N_BOOT = 0;
modelsForBootstrap = {'SelfOnly','SelfOnlyBuyBias','OppOnly','OppOnlyBuyBias','SelfOpp','SelfBubble'};

% Optimization controls
NUM_RESTARTS = 5;  % restarts per fit group per model

% Optionally fix RNG for reproducibility
% SEED = 1;
% rng(SEED);

%% =========================================================================
% PART 0: LOAD CONDITION PACKS AND BUILD TRIAL TABLE + NEURON STRUCT
% =========================================================================
%
% For each condition set (baseline / OT / Saline) and neural event:
%   - Build trialTableRaw: trial-wise behavior + covariates.
%   - Build neurons struct array linking each neuron to its FR vector
%     (z-scored relative to a baseline window) and the trial rows it uses.
%
% Behavior is pooled across blocks within a .nex/session but not across
% sessions. RL fitting later will use ONLY a canonical event table
% (TARGET_EVTS{1}); the per-event tables are here for neural purposes.
% =========================================================================

fprintf('\n=== PART 0: Build trial tables and neuron structs ===\n');

allResults = struct();

for pair = 1:numel(TARGET_EVTS)
    TARGET_EVT     = TARGET_EVTS{pair};
    TARGET_EVT_WIN = TARGET_EVT_WINS{pair};

    % fprintf('\n[PART 0] Event %s (%d/%d)\n', string(TARGET_EVT), pair, numel(TARGET_EVTS));

    for s = 1:numel(cond_sets)
        conds   = cond_sets{s};
        setName = setNames{s};

        % fprintf('  Set: %s | Event: %s\n', setName, string(TARGET_EVT));

        % ------------------------------ 0.1) Load condition packs -------
        Cstruct = struct();
        for c = 1:numel(conds)
            conditionName = conds{c};
            matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', conditionName));
            if ~exist(matFile,'file')
                warning('  [WARN] Missing file: %s', matFile);
                continue;
            end
            tmp = load(matFile, 'C');
            Cstruct.(matlab.lang.makeValidName(conditionName)) = tmp.C;
        end
        condNames = fieldnames(Cstruct);
        if isempty(condNames)
            warning('  [WARN] No condition packs loaded for set %s', setName);
            continue;
        end

        % --------------------- 0.2) Build trialTable & neurons ----------
        trialTable = table();
        neurons    = struct('Monkey',{},'Condition',{},'Session',{}, ...
                            'FileIndex',{},'NeuronIndex',{}, ...
                            'trialRowIdx',{},'FR',{});
        neuronCounter = 0;

        for c = 1:numel(condNames)
            conditionName = condNames{c};
            C = Cstruct.(conditionName);

            % Loop over .nex files (eventTables)
            for f = 1:numel(C.eventTables)
                etab = C.eventTables{f};

                % Require subject choice column
                if ~ismember('option', etab.Properties.VariableNames)
                    continue;
                end
                if isempty(etab.option{1})
                    continue;
                end

                % ------------- Choices and optional BUY collapse ----------
                subjectChoiceStr  = local_mapChoicesToLabels(etab.option{1});
                opponentChoiceStr = local_mapChoicesToLabels(etab.optionOpp{1});

                if BUY_NON_BUY_ONLY
                    subjectChoiceStr(subjectChoiceStr ~= "BUY")   = "Non-BUY";
                    opponentChoiceStr(opponentChoiceStr ~= "BUY") = "Non-BUY";
                end

                % Valid trials (we drop trials with missing subject choice)
                isValidTrial = ~ismissing(subjectChoiceStr);

                % Condition label per trial (original pack name)
                conditionLabelPerTrial = repmat(string(conditionName), numel(subjectChoiceStr), 1);

                % -------------------- Covariates from event table --------
                monkeyId    = string(etab.monkey{1}(:));   % monkey ID per trial
                sessionId   = double(etab.session{1}(:));  % session index (.nex)
                trialIndex  = double(etab.trialNum{1}(:)); % trial number within session
                marketLabel = string(etab.marketOrig{1}(:));
                bubbleLabel = string(etab.bubbleMarket{1}(:));
                y = etab.year{1}(:); m = etab.month{1}(:); d = etab.day{1}(:); sod = etab.SessionOfDay{1}(:);
                date = string(compose('%04d-%02d-%02d-%02d', y, m, d, sod));

                % Portfolio value at analysis event (pre or post)
                if string(TARGET_EVT).contains("m1rew2") || string(TARGET_EVT).contains("m2")
                    portfolioValueAtEvent = double(etab.postPortfolio{1}(:));
                else
                    portfolioValueAtEvent = double(etab.prePortfolio{1}(:));
                end

                % Opponent portfolio (per trial)
                opponentPortfolioValue = etab.prePortfolioOpp{1}(:);

                % Reward components
                choiceOutcomeReward = etab.fb1{1}(:);
                assert(sum(choiceOutcomeReward < 0) == 0);

                % Dividend-based portfolio reward
                % portfolioRewardFromDividend = etab.postPortfolio{1}(:) .* etab.divPerShare{1}(:);
                portfolioRewardFromDividend = etab.fb2{1}(:);
                assert(sum(portfolioRewardFromDividend < 0) == 0);

                % -------------------- Live vs not_Live collapse ----------
                if LIVE_NON_LIVE_ONLY
                    isLiveCond = contains(conditionLabelPerTrial, "Live", 'IgnoreCase', true);
                    conditionLabelForModeling = repmat("not_Live", numel(conditionLabelPerTrial), 1);
                    conditionLabelForModeling(isLiveCond) = "Live";
                else
                    conditionLabelForModeling = conditionLabelPerTrial;
                end

                % -------------------- Neural FR extraction ----------------
                Z_allNeurons = local_get_Z_fr_allNeurons(C, f, TARGET_EVT, TARGET_EVT_WIN, ...
                                                         BASELINE_EVT, BASELINE_WIN);
                if isempty(Z_allNeurons)
                    Z_allNeurons = nan(0, numel(subjectChoiceStr));
                end

                % -------------------- Assemble trial rows -----------------
                nBefore = height(trialTable);
                tmpT = table( ...
                    date(isValidTrial), ...
                    conditionLabelForModeling(isValidTrial), ...
                    string(monkeyId(isValidTrial)), ...
                    sessionId(isValidTrial), ...
                    trialIndex(isValidTrial), ...
                    string(marketLabel(isValidTrial)), ...
                    string(bubbleLabel(isValidTrial)), ...
                    string(subjectChoiceStr(isValidTrial)), ...
                    string(opponentChoiceStr(isValidTrial)), ...
                    choiceOutcomeReward(isValidTrial), ...
                    portfolioRewardFromDividend(isValidTrial), ...
                    portfolioValueAtEvent(isValidTrial), ...
                    opponentPortfolioValue(isValidTrial), ...
                    'VariableNames', { ...
                        'Date', ...
                        'Condition', ...
                        'Monkey', ...
                        'Session', ...
                        'Trial', ...
                        'Market', ...
                        'Bubble', ...
                        'subjectChoice', ...
                        'opponentChoice', ...
                        'choiceOutcomeReward', ...
                        'portfolioRewardFromDividend', ...
                        'portfolioValueAtEvent', ...
                        'opponentPortfolioValue'});

                trialTable = [trialTable; tmpT]; %#ok<AGROW>
                nAfter = height(trialTable);
                trialRowsThisFile = (nBefore+1:nAfter).';

                % -------------------- Create neuron entries ---------------
                nNeurons = size(Z_allNeurons, 1);
                for n = 1:nNeurons
                    neuronCounter = neuronCounter + 1;
                    neurons(neuronCounter).Monkey      = monkeyId(1);
                    neurons(neuronCounter).Condition   = conditionLabelForModeling(find(isValidTrial,1));
                    neurons(neuronCounter).Session     = sessionId(1);
                    neurons(neuronCounter).FileIndex   = f;
                    neurons(neuronCounter).NeuronIndex = n;
                    neurons(neuronCounter).trialRowIdx = trialRowsThisFile;
                    if ~isempty(Z_allNeurons)
                        neurons(neuronCounter).FR = Z_allNeurons(n, isValidTrial).';
                    else
                        neurons(neuronCounter).FR = nan(sum(isValidTrial),1);
                    end
                end
            end
        end

        if isempty(trialTable)
            warning('  [WARN] No trialTable rows for set %s, event %s', setName, TARGET_EVT);
            continue;
        end

        % -------------------- 0.3) Categorical encodings ----------------
        trialTable.RowIndex = (1:height(trialTable)).';  % for RL alignment

        trialTable.MonkeyRaw = categorical(trialTable.Monkey);
        trialTable.Monkey    = categorical(trialTable.Monkey);
        trialTable.Market    = categorical(trialTable.Market);
        trialTable.Bubble    = categorical(trialTable.Bubble);

        % Condition categories
        if LIVE_NON_LIVE_ONLY
            trialTable.Condition = categorical(trialTable.Condition);
            trialTable.Condition = setcats(trialTable.Condition, {'not_Live','Live'});
        else
            trialTable.Condition = categorical(trialTable.Condition);
            baseCats = {'AI','Replay','Decoy','Live'};
            if all(ismember(baseCats, categories(trialTable.Condition)))
                trialTable.Condition = setcats(trialTable.Condition, baseCats);
            end
        end

        % Choice categories
        if BUY_NON_BUY_ONLY
            trialTable.subjectChoice  = categorical(trialTable.subjectChoice,  {'Non-BUY','BUY'});
            trialTable.opponentChoice = categorical(trialTable.opponentChoice, {'Non-BUY','BUY'});
        else
            trialTable.subjectChoice  = categorical(trialTable.subjectChoice,  {'HOLD','BUY','SELL'});
            trialTable.opponentChoice = categorical(trialTable.opponentChoice, {'HOLD','BUY','SELL'});
        end

        % -------------------- 0.4) Block structure ----------------------
        [blockID, ~] = findgroups(trialTable(:,{'Date','Monkey','Condition','Session','Market'}));
        trialTable.BlockID = blockID;

        trialTable.TrialInBlock = nan(height(trialTable),1);
        uBlocks = unique(blockID);
        for b = uBlocks.'
            idxBlock = find(blockID == b);
            [~, sortIdx] = sort(trialTable.Trial(idxBlock));
            assert(max(sortIdx) == 15);
            trialTable.TrialInBlock(idxBlock(sortIdx)) = (1:numel(idxBlock)).';
        end

        % Bubble flag (1 if Bubble name contains 'bubble')
        bubbleCats   = categories(trialTable.Bubble);
        isBubbleCat  = contains(bubbleCats, '1', 'IgnoreCase', true);
        trialTable.BubbleFlag = zeros(height(trialTable),1);
        for i = 1:numel(isBubbleCat)
            if isBubbleCat(i)
                trialTable.BubbleFlag(trialTable.Bubble == bubbleCats{i}) = 1;
            end
        end

        % -------------------- 0.5) Pooling (Monkey vs ALL) --------------
        switch lower(FIT_SCOPE)
            case 'pooled'
                trialTable.Monkey = categorical(repmat("ALL", height(trialTable), 1));
            case 'per_monkey'
                % leave Monkey unchanged
            otherwise
                error('Unknown FIT_SCOPE: %s', FIT_SCOPE);
        end

        % Store raw table & neurons
        allResults.(setName).(TARGET_EVT).trialTableRaw = trialTable;
        allResults.(setName).(TARGET_EVT).neurons       = neurons;

        % Only print a summary for the canonical event (TARGET_EVTS{1})
        if pair == 1
            % Session ≈ unique Date per (Monkey × Condition)
            [sessGroup, sessKey] = findgroups(trialTable(:,{'Monkey','Condition'}));
            nSessionsPerGroup = splitapply(@(d) numel(unique(d)), trialTable.Date, sessGroup);

            totalSessions = sum(nSessionsPerGroup);
            trialsPerSession_nominal = 90;  % 6 blocks × 15 trials

            fprintf('    -> %d sessions (~%d trials total), %d neurons\n', ...
                    totalSessions, totalSessions * trialsPerSession_nominal, numel(neurons));

            fprintf('       Session counts by Monkey × Condition:\n');
            for g = 1:height(sessKey)
                fprintf('         Monkey %s | %s : %d sessions\n', ...
                    string(sessKey.Monkey(g)), string(sessKey.Condition(g)), nSessionsPerGroup(g));
            end
        end
    end
end

%% =========================================================================
% PART 1: BEHAVIOR-ONLY RL MODELING (EVENT-AGNOSTIC)
% =========================================================================
%
% RL is now fit using only a single "canonical" event table:
%   canonicalEvt = TARGET_EVTS{1}
%
% Because choices, rewards, and block structure are identical across
% neural events, this avoids redundant fitting while preserving all
% behavior for RL.
% =========================================================================

fprintf('\n=== PART 1: Behavior-only RL modeling (event-agnostic) ===\n');

behaviorFitSummary = table( ...
    string.empty(0,1), ...  % SetName
    string.empty(0,1), ...  % Event (canonical tag; for bookkeeping only)
    string.empty(0,1), ...  % Monkey
    string.empty(0,1), ...  % Condition
    string.empty(0,1), ...  % DateID
    string.empty(0,1), ...  % ModelName
    double.empty(0,1), ...  % NumTrials
    double.empty(0,1), ...  % NumParams
    double.empty(0,1), ...  % NegLogLik
    double.empty(0,1), ...  % AIC
    double.empty(0,1), ...  % BIC
    cell(0,1), ...          % Theta
    'VariableNames', {'SetName','Event','Monkey','Condition','DateID','ModelName', ...
                      'NumTrials','NumParams','NegLogLik','AIC','BIC','Theta'});

try
    optimOpts = optimoptions('fmincon','Display','off','Algorithm','interior-point');
catch
    optimOpts = [];
end

canonicalEvt = TARGET_EVTS{1};

for s = 1:numel(cond_sets)
    setName = setNames{s};
    fprintf('  [PART 1] Set %s (canonical event %s)\n', setName, string(canonicalEvt));

    if ~isfield(allResults, setName) || ~isfield(allResults.(setName), canonicalEvt)
        fprintf('    -> Skipping (no canonical event table).\n');
        continue;
    end

    trialTable = allResults.(setName).(canonicalEvt).trialTableRaw;

    switch lower(FIT_SCOPE)
        case 'per_monkey'
            groupVars = {'Monkey','Condition','Date'};
        case 'pooled'
            groupVars = {'Condition','Date'};
        otherwise
            error('Unknown FIT_SCOPE: %s', FIT_SCOPE);
    end

    [mcGroup, mcKey] = findgroups(trialTable(:,groupVars));
    nGroupsMC = max(mcGroup);
    fprintf('    -> %d behavior groups for RL\n', nGroupsMC);

    groupTables = cell(nGroupsMC,1);

    parfor g = 1:nGroupsMC
        idxMC   = (mcGroup == g);
        sessData = trialTable(idxMC,:);
        sessData = sortrows(sessData, {'Date','BlockID','TrialInBlock'});

        if ismember('Monkey', groupVars)
            monkeyID = char(mcKey.Monkey(g));
        else
            monkeyID = 'ALL';
        end
        conditionID = char(mcKey.Condition(g));
        sessionID   = mcKey.Date(g);

        localTbl = behaviorFitSummary([],:);  % empty template

        if height(sessData) < 10
            groupTables{g} = localTbl;
            continue;
        end

        for m = 1:numel(MODEL_TO_RUN)
            modelName = MODEL_TO_RUN{m};
            [thetaHat, negLL, nParams] = local_fit_one_model( ...
                modelName, sessData, BUY_NON_BUY_ONLY, NUM_RESTARTS, optimOpts);

            nTrials = height(sessData);
            AIC = 2*nParams + 2*negLL;
            BIC = log(nTrials)*nParams + 2*negLL;

            newRow = table( ...
                string(setName), ...
                string(canonicalEvt), ...
                string(monkeyID), ...
                string(conditionID), ...
                string(sessionID), ...   % DateID column as string
                string(modelName), ...
                nTrials, nParams, negLL, AIC, BIC, {thetaHat}, ...
                'VariableNames', behaviorFitSummary.Properties.VariableNames);

            localTbl = [localTbl; newRow]; %#ok<AGROW>
        end

        groupTables{g} = localTbl;
    end

    nonEmpty = ~cellfun(@isempty, groupTables);
    if any(nonEmpty)
        behaviorFitSummary = [behaviorFitSummary; vertcat(groupTables{nonEmpty})];
    end
end

allResults.behaviorFits = behaviorFitSummary;
fprintf('  -> behaviorFitSummary rows: %d\n', height(behaviorFitSummary));

%% =========================================================================
% PART 1b: PARAMETRIC BLOCK BOOTSTRAP (per Set)
% =========================================================================
% comment: I'm not sure if this is helpful at all. My goal was to determine
% if those coefficients in the RL model are real, but I don't know if the
% RL model fully predicts the behavior, so the CI might be very large
% anyway. TBD.-- Tim, 251126
fprintf('\n=== PART 1b: Bootstrap parameter uncertainty ===\n');

allResults.behaviorBootstrap = struct();
canonicalEvt = TARGET_EVTS{1};

for s = 1:numel(cond_sets)
    setName = setNames{s};
    fprintf('  [PART 1b] Set %s\n', setName);

    if ~isfield(allResults, setName) || ~isfield(allResults.(setName), canonicalEvt)
        fprintf('    -> Skipping (no canonical event table).\n');
        continue;
    end

    trialTable = allResults.(setName).(canonicalEvt).trialTableRaw;

    fitMaskSE = strcmp(allResults.behaviorFits.SetName, setName);
    fitTableSE = allResults.behaviorFits(fitMaskSE,:);

    behaviorBootstrap = local_rl_parametric_block_bootstrap( ...
        trialTable, fitTableSE, BUY_NON_BUY_ONLY, modelsForBootstrap, N_BOOT, FIT_SCOPE);

    allResults.behaviorBootstrap.(setName) = behaviorBootstrap;
    fprintf('    -> Bootstrap rows: %d\n', height(behaviorBootstrap));
end

%% =========================================================================
% PART 1c-1: SUMMARIZE RL PARAMETERS INTO LONG-FORM TABLE
% =========================================================================
fprintf('\n=== PART 1c-1: Summarize RL parameters ===\n');

paramNamesPerModel = struct();
paramNamesPerModel.SelfOnly = {'alpha','beta','omega'};
paramNamesPerModel.SelfOnlyBuyBias       = {'alpha','beta','omega','kappaBuy'};
paramNamesPerModel.OppOnly          = {'alpha_opp','beta'};
paramNamesPerModel.OppOnlyBuyBias   = {'alpha_opp','beta','kappaOppBuy'};
paramNamesPerModel.SelfOpp        = {'alpha_self','alpha_opp','beta','omega_self','gamma','kappaBuy'};
paramNamesPerModel.SelfBubble     = {'alpha','beta','omega_nonBubble','omega_bubble','kappa_nonBubble','kappa_bubble'};

allResults.behaviorParamTable = local_behaviorParamTable(allResults.behaviorFits, paramNamesPerModel);
fprintf('  -> behaviorParamTable rows: %d\n', height(allResults.behaviorParamTable));

%% =========================================================================
% PART 1c-2: RECOMPUTE RL LATENTS FOR TIMECOURSES / TERMINAL ANALYSIS
% =========================================================================
fprintf('\n=== PART 1c-2: Compute RL latents and attach to all events ===\n');

canonicalEvt = TARGET_EVTS{1};

% ----------------- 1) Compute RL latents once per Set --------------------
for s = 1:numel(cond_sets)
    setName = setNames{s};
    fprintf('  [LATENTS] Set %s\n', setName);

    if ~isfield(allResults, setName) || ...
       ~isfield(allResults.(setName), canonicalEvt) || ...
       ~isfield(allResults.(setName).(canonicalEvt), 'trialTableRaw')
        fprintf('    -> Skipping (no canonical trialTableRaw).\n');
        continue;
    end

    trialTable = allResults.(setName).(canonicalEvt).trialTableRaw;

    % Preallocate RL columns in the canonical table
    trialTable.RL_predictionError_self = nan(height(trialTable),1);
    trialTable.RL_chosenValue_self     = nan(height(trialTable),1);
    trialTable.RL_Q_self_buy           = nan(height(trialTable),1);
    trialTable.RL_Q_self_hold          = nan(height(trialTable),1);
    trialTable.RL_Q_self_sell          = nan(height(trialTable),1);

    trialTable.RL_predictionError_opp  = nan(height(trialTable),1);
    trialTable.RL_chosenValue_opp      = nan(height(trialTable),1);
    trialTable.RL_Q_opp_buy            = nan(height(trialTable),1);
    trialTable.RL_Q_opp_hold           = nan(height(trialTable),1);
    trialTable.RL_Q_opp_sell           = nan(height(trialTable),1);

    switch lower(FIT_SCOPE)
        case 'per_monkey'
            groupVars = {'Monkey','Condition','Date'};
        case 'pooled'
            groupVars = {'Condition','Date'};
        otherwise
            error('Unknown FIT_SCOPE: %s', FIT_SCOPE);
    end

    [mcGroup, mcKey] = findgroups(trialTable(:,groupVars));
    nGroupsMC = max(mcGroup);
    fprintf('    -> %d groups for RL latents\n', nGroupsMC);

    behaviorFits_local                = allResults.behaviorFits;
    RL_MODEL_FOR_Q_TIMESERIES_local   = RL_MODEL_FOR_Q_TIMESERIES;

    rowIdxCell    = cell(nGroupsMC,1);
    dSelfCell     = cell(nGroupsMC,1);
    vSelfCell     = cell(nGroupsMC,1);
    qBuySelfCell  = cell(nGroupsMC,1);
    qHoldSelfCell = cell(nGroupsMC,1);
    qSellSelfCell = cell(nGroupsMC,1);

    dOppCell      = cell(nGroupsMC,1);
    vOppCell      = cell(nGroupsMC,1);
    qBuyOppCell   = cell(nGroupsMC,1);
    qHoldOppCell  = cell(nGroupsMC,1);
    qSellOppCell  = cell(nGroupsMC,1);

    parfor g = 1:nGroupsMC
        rowIdxCell{g}    = [];
        dSelfCell{g}     = [];
        vSelfCell{g}     = [];
        qBuySelfCell{g}  = [];
        qHoldSelfCell{g} = [];
        qSellSelfCell{g} = [];

        dOppCell{g}      = [];
        vOppCell{g}      = [];
        qBuyOppCell{g}   = [];
        qHoldOppCell{g}  = [];
        qSellOppCell{g}  = [];

        idxMC   = (mcGroup == g);
        sessData = trialTable(idxMC,:);
        sessData = sortrows(sessData, {'Date','BlockID','TrialInBlock'});
        if isempty(sessData)
            continue;
        end

        if ismember('Monkey', groupVars)
            monkeyID = char(mcKey.Monkey(g));
        else
            monkeyID = 'ALL';
        end
        conditionID = char(mcKey.Condition(g));
        sessionID   = mcKey.Date(g);

        modelNameForGroup = RL_MODEL_FOR_Q_TIMESERIES_local;

        baseMask = strcmp(behaviorFits_local.SetName,   setName) & ...
                   strcmp(behaviorFits_local.Monkey,    monkeyID) & ...
                   strcmp(behaviorFits_local.Condition, conditionID) & ...
                   behaviorFits_local.DateID  == sessionID;

        if strcmpi(RL_MODEL_FOR_Q_TIMESERIES_local,'best_by_AIC') || strcmpi(RL_MODEL_FOR_Q_TIMESERIES_local,'best_by_BIC')
            fitsMC = behaviorFits_local(baseMask,:);
            if isempty(fitsMC)
                continue;
            end
            if strcmpi(RL_MODEL_FOR_Q_TIMESERIES_local,'best_by_AIC')
                [~,iBest] = min(fitsMC.AIC);
            else
                [~,iBest] = min(fitsMC.BIC);
            end
            modelNameForGroup = fitsMC.ModelName{iBest};
        end

        maskFit = baseMask & strcmp(behaviorFits_local.ModelName, modelNameForGroup);
        if ~any(maskFit)
            continue;
        end

        thetaHat = behaviorFits_local.Theta{find(maskFit,1,'first')};

        latentsSelf = [];
        latentsOpp  = [];
        switch modelNameForGroup
            case 'SelfOnlyBuyBias'
                [~, latentsSelf] = local_rlNegLogLik_SelfOnlyBuyBias(thetaHat, sessData, BUY_NON_BUY_ONLY);
            case 'SelfOnly'
                [~, latentsSelf] = local_rlNegLogLik_SelfOnly(thetaHat, sessData, BUY_NON_BUY_ONLY);
            case 'SelfBubble'
                [~, latentsSelf] = local_rlNegLogLik_SelfBubble(thetaHat, sessData, BUY_NON_BUY_ONLY);
            case 'SelfOpp'
                [~, latentsSelf, latentsOpp] = local_rlNegLogLik_SelfOpp(thetaHat, sessData, BUY_NON_BUY_ONLY);
            otherwise
                continue;
        end

        rowIdxCell{g}    = sessData.RowIndex(:);
        dSelfCell{g}     = latentsSelf.delta(:);
        vSelfCell{g}     = latentsSelf.value(:);
        qBuySelfCell{g}  = latentsSelf.Q(:,1);

        if size(latentsSelf.Q,2) >= 2
            qHoldSelfCell{g} = latentsSelf.Q(:,2);
        end
        if size(latentsSelf.Q,2) >= 3
            qSellSelfCell{g} = latentsSelf.Q(:,3);
        end

        if ~isempty(latentsOpp) && strcmp(modelNameForGroup,'SelfOpp')
            dOppCell{g}    = latentsOpp.delta(:);
            vOppCell{g}    = latentsOpp.value(:);
            qBuyOppCell{g} = latentsOpp.Q(:,1);
            if size(latentsOpp.Q,2) >= 2
                qHoldOppCell{g} = latentsOpp.Q(:,2);
            end
            if size(latentsOpp.Q,2) >= 3
                qSellOppCell{g} = latentsOpp.Q(:,3);
            end
        end
    end

    % Serial write-back
    for g = 1:nGroupsMC
        idx = rowIdxCell{g};
        if isempty(idx)
            continue;
        end

        trialTable.RL_predictionError_self(idx) = dSelfCell{g};
        trialTable.RL_chosenValue_self(idx)     = vSelfCell{g};
        trialTable.RL_Q_self_buy(idx)           = qBuySelfCell{g};

        if ~isempty(qHoldSelfCell{g})
            trialTable.RL_Q_self_hold(idx)      = qHoldSelfCell{g};
        end
        if ~isempty(qSellSelfCell{g})
            trialTable.RL_Q_self_sell(idx)      = qSellSelfCell{g};
        end

        if ~isempty(dOppCell{g})
            trialTable.RL_predictionError_opp(idx) = dOppCell{g};
            trialTable.RL_chosenValue_opp(idx)     = vOppCell{g};
            trialTable.RL_Q_opp_buy(idx)           = qBuyOppCell{g};
            if ~isempty(qHoldOppCell{g})
                trialTable.RL_Q_opp_hold(idx)      = qHoldOppCell{g};
            end
            if ~isempty(qSellOppCell{g})
                trialTable.RL_Q_opp_sell(idx)      = qSellOppCell{g};
            end
        end
    end

    % Alias for GLM: use self prediction error
    trialTable.RL_predictionError = trialTable.RL_predictionError_self;

    % Lag opponent portfolio by 1 within block
    trialTable.opponentPortfolioValue_lag1 = nan(height(trialTable),1);
    uBlocks = unique(trialTable.BlockID);
    for b = uBlocks.'
        idxBlock = find(trialTable.BlockID == b);
        blockData = trialTable(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);
        oppNow  = trialTable.opponentPortfolioValue(ordIdx);
        oppLag1 = [NaN; oppNow(1:end-1)];
        trialTable.opponentPortfolioValue_lag1(ordIdx) = oppLag1;
    end

    % Choice dummies for GLM
    if BUY_NON_BUY_ONLY
        trialTable.isBuy    = double(trialTable.subjectChoice == 'BUY');
        trialTable.isNonBuy = double(trialTable.subjectChoice == 'Non-BUY');
    else
        trialTable.isBuy  = double(trialTable.subjectChoice == 'BUY');
        trialTable.isHold = double(trialTable.subjectChoice == 'HOLD');
        trialTable.isSell = double(trialTable.subjectChoice == 'SELL');
    end

    % Store canonical table with RL & behavioral covariates
    allResults.(setName).behaviorTrialTableWithRL = trialTable;
    fprintf('    -> RL latents attached on canonical table (%d trials).\n', height(trialTable));
end

% ----------------- 2) Attach RL + covariates to each event ----------------
for pair = 1:numel(TARGET_EVTS)
    TARGET_EVT = TARGET_EVTS{pair};

    for s = 1:numel(cond_sets)
        setName = setNames{s};

        if ~isfield(allResults, setName) || ...
           ~isfield(allResults.(setName), TARGET_EVT) || ...
           ~isfield(allResults.(setName), 'behaviorTrialTableWithRL')
            continue;
        end

        trialTable_evt = allResults.(setName).(TARGET_EVT).trialTableRaw;
        behTable       = allResults.(setName).behaviorTrialTableWithRL;

        if height(trialTable_evt) ~= height(behTable)
            error('Row mismatch between event %s and canonical behavior table for set %s.', ...
                  TARGET_EVT, setName);
        end

        % Copy RL and behavioral covariates
        trialTable_evt.RL_predictionError_self = behTable.RL_predictionError_self;
        trialTable_evt.RL_chosenValue_self     = behTable.RL_chosenValue_self;
        trialTable_evt.RL_Q_self_buy           = behTable.RL_Q_self_buy;
        trialTable_evt.RL_Q_self_hold          = behTable.RL_Q_self_hold;
        trialTable_evt.RL_Q_self_sell          = behTable.RL_Q_self_sell;

        trialTable_evt.RL_predictionError_opp  = behTable.RL_predictionError_opp;
        trialTable_evt.RL_chosenValue_opp      = behTable.RL_chosenValue_opp;
        trialTable_evt.RL_Q_opp_buy            = behTable.RL_Q_opp_buy;
        trialTable_evt.RL_Q_opp_hold           = behTable.RL_Q_opp_hold;
        trialTable_evt.RL_Q_opp_sell           = behTable.RL_Q_opp_sell;

        trialTable_evt.RL_predictionError          = behTable.RL_predictionError;
        trialTable_evt.opponentPortfolioValue_lag1 = behTable.opponentPortfolioValue_lag1;

        if BUY_NON_BUY_ONLY
            trialTable_evt.isBuy    = behTable.isBuy;
            trialTable_evt.isNonBuy = behTable.isNonBuy;
        else
            trialTable_evt.isBuy  = behTable.isBuy;
            trialTable_evt.isHold = behTable.isHold;
            trialTable_evt.isSell = behTable.isSell;
        end

        allResults.(setName).(TARGET_EVT).trialTableWithRL = trialTable_evt;
    end
end

% ----------------- 3) Collect terminal Q-values at TrialInBlock==15 -------
allResults.terminalQTable = local_buildTerminalQTable(allResults, setNames, canonicalEvt, BUY_NON_BUY_ONLY);

fprintf('\n=== Analysis complete. Results are stored in allResults. ===\n');

%%
fprintf('\n=== Advanced RL analysis: starting (DeltaQ, kappa, Q timecourses) ===\n');

%% ====================== PART A: ADD DeltaQ LATENT =======================
% DeltaQ is defined here as:
%   DeltaQ(t) = Q(chosen action, t) - max Q(other available actions, t)
%
% We compute this for each trial in the canonical behavior table for each
% set, then copy the same column into all event-specific trial tables.

fprintf('\n[PART A] Computing RL_DeltaQ (decision advantage) per trial...\n');

allResults = local_addDeltaQ(allResults, setNames, TARGET_EVTS, BUY_NON_BUY_ONLY);

fprintf('  -> RL_DeltaQ and related columns attached to all behavior/event tables.\n');

%% ================== PART B: KAPPA / MODEL COMPARISON ====================

fprintf('\n[PART B] Pairwise RL model comparisons (BIC / BF)...\n');
% -------------------------------------------------------------------------
% Parameter counts / hierarchy
%   OppOnly          : [alpha_opp, beta]                      -> 2 params
%   SelfOnly         : [alpha, beta, omega]                   -> 3 params
%   OppOnlyBuyBias   : [alpha_opp, beta, kappaOppBuy]         -> 3 params
%   SelfOnlyBuyBias  : [alpha, beta, omega, kappaBuy]         -> 4 params
%   SelfBubble       : [alpha, beta, omega_nonB, omega_B,
%                       kappa_nonB, kappa_B]                  -> 6 params
%   SelfOpp          : [alpha_self, alpha_opp, beta,
%                       omega_self, gamma, kappaBuy]          -> 6 params
%
% In each comparison below, ModelA is the simpler reference model,
% ModelB is the more complex / enriched model.
% -------------------------------------------------------------------------

% 0) Base: self vs opp (no bias)  -- who explains behavior better?
[cmp_selfBase_vs_oppBase, sum_selfBase_vs_oppBase] = ...
    local_buildModelComparison(allResults.behaviorFits, ...
                               'OppOnly', 'SelfOnly');

% 1) Self: no bias vs BUY bias   (SelfOnly vs SelfOnlyBuyBias)
[cmp_noBias_vs_buyBias, sum_noBias_vs_buyBias] = ...
    local_buildModelComparison(allResults.behaviorFits, ...
                               'SelfOnly', 'SelfOnlyBuyBias');

% 1b) Opp: no bias vs BUY bias   (OppOnly vs OppOnlyBuyBias)
[cmp_opp_noBias_vs_buyBias, sum_opp_noBias_vs_buyBias] = ...
    local_buildModelComparison(allResults.behaviorFits, ...
                               'OppOnly', 'OppOnlyBuyBias');

% 2) Self BUY bias vs +Bubble sensitive   (SelfOnlyBuyBias vs SelfBubble)
[cmp_buyBias_vs_bubble, sum_buyBias_vs_bubble] = ...
    local_buildModelComparison(allResults.behaviorFits, ...
                               'SelfOnlyBuyBias', 'SelfBubble');

% 3) Self BUY bias vs +Opp Q term        (SelfOnlyBuyBias vs SelfOpp)
[cmp_buyBias_vs_opp, sum_buyBias_vs_opp] = ...
    local_buildModelComparison(allResults.behaviorFits, ...
                               'SelfOnlyBuyBias', 'SelfOpp');

% 4) Opp-only BUY bias vs Self+Opp coupled model
[cmp_oppOnlyBias_vs_selfOpp, sum_oppOnlyBias_vs_selfOpp] = ...
    local_buildModelComparison(allResults.behaviorFits, ...
                               'OppOnlyBuyBias', 'SelfOpp');

% ---------------- Store in allResults for later inspection ---------------
allResults.modelComparison.selfBase_vs_oppBase      = cmp_selfBase_vs_oppBase;
allResults.modelSummary.selfBase_vs_oppBase         = sum_selfBase_vs_oppBase;

allResults.modelComparison.noBias_vs_buyBias        = cmp_noBias_vs_buyBias;
allResults.modelSummary.noBias_vs_buyBias           = sum_noBias_vs_buyBias;

allResults.modelComparison.opp_noBias_vs_buyBias    = cmp_opp_noBias_vs_buyBias;
allResults.modelSummary.opp_noBias_vs_buyBias       = sum_opp_noBias_vs_buyBias;

allResults.modelComparison.buyBias_vs_bubble        = cmp_buyBias_vs_bubble;
allResults.modelSummary.buyBias_vs_bubble           = sum_buyBias_vs_bubble;

allResults.modelComparison.buyBias_vs_opp           = cmp_buyBias_vs_opp;
allResults.modelSummary.buyBias_vs_opp              = sum_buyBias_vs_opp;

allResults.modelComparison.oppOnlyBias_vs_selfOpp   = cmp_oppOnlyBias_vs_selfOpp;
allResults.modelSummary.oppOnlyBias_vs_selfOpp      = sum_oppOnlyBias_vs_selfOpp;

fprintf('  -> selfBase vs oppBase rows         : %d\n', height(cmp_selfBase_vs_oppBase));
fprintf('  -> self noBias vs BUY bias rows     : %d\n', height(cmp_noBias_vs_buyBias));
fprintf('  -> opp  noBias vs BUY bias rows     : %d\n', height(cmp_opp_noBias_vs_buyBias));
fprintf('  -> self BUY bias vs Bubble rows     : %d\n', height(cmp_buyBias_vs_bubble));
fprintf('  -> self BUY bias vs Opp rows        : %d\n', height(cmp_buyBias_vs_opp));
fprintf('  -> opp-only BUY bias vs Self+Opp rows: %d\n', height(cmp_oppOnlyBias_vs_selfOpp));


%% =========== PART C: Q & DeltaQ TIMECOURSES AND TERMINAL ADVANTAGE ======
% (C1) QTimecourse:
%       mean and SEM of Q(action, TrialInBlock) grouped by
%       [SetName, MonkeyRaw, Condition, Action, TrialInBlock].
%
% (C2) DeltaQTimecourse:
%       mean and SEM of RL_DeltaQ (decision advantage), grouped by
%       [SetName, MonkeyRaw, Condition, ChosenAction, TrialInBlock].
%
% (C3) terminalQAdvantage:
%       for a fixed TrialInBlock (default 15), compute
%         BuyAdvantage = Q_buy - max(Q_other actions)
%       per block and summarize by condition/monkey.

fprintf('\n[PART C] Building Q and DeltaQ timecourse and terminal-advantage tables...\n');

allResults.QTimecourse = local_buildQTimecourse(allResults, setNames, BUY_NON_BUY_ONLY);
allResults.DeltaQTimecourse = local_buildDeltaQTimecourse(allResults, setNames);
allResults.terminalQAdvantage = local_buildTerminalQAdvantage(allResults, setNames, BUY_NON_BUY_ONLY);

fprintf('  -> QTimecourse rows          : %d\n', height(allResults.QTimecourse));
fprintf('  -> DeltaQTimecourse rows     : %d\n', height(allResults.DeltaQTimecourse));
fprintf('  -> terminalQAdvantage rows   : %d\n', height(allResults.terminalQAdvantage));

fprintf('\n=== Advanced RL analysis complete. New fields added to allResults. ===\n');

%% Make plots:
% local_plotBehaviorParamComparisonBySet(allResults.behaviorFits, paramNamesPerModel);
local_plotBehaviorParamByMonkey(allResults.behaviorFits, paramNamesPerModel);
local_plotModelBIC_AllModelsByCondition(allResults.behaviorFits, MODEL_TO_RUN);
%%
local_plotQTimecourse(allResults.QTimecourse);
local_plotDeltaQTimecourse(allResults.DeltaQTimecourse);
local_plotTerminalQAdvantage(allResults.terminalQAdvantage);
% local_plotModelBIC_EvidenceByCondition(cmp_noBias_vs_buyBias, ...
%                                        'SelfOnly', 'SelfOnlyBuyBias');
% local_plotModelBIC_EvidenceByCondition(cmp_buyBias_vs_bubble, ...
%                                        'SelfOnlyBuyBias', 'SelfBubble');
% local_plotModelBIC_EvidenceByCondition(cmp_buyBias_vs_opp, ...
%                                        'SelfOnlyBuyBias', 'SelfOpp');

% % Optional bootstrap visualization (kept commented so it can be re-enabled)
% for s = 1:numel(setNames)
%     setName = setNames{s};
%     if isfield(allResults.behaviorBootstrap, setName)
%         visualizeBehaviorBootstrap(allResults.behaviorBootstrap.(setName), ...
%             paramNamesPerModel, setName, []);
%     end
% end

%% Save all figures to folder
outDirPlots = fullfile(OUTDIR,'RL_results');
if ~exist(outDirPlots,'dir')
    mkdir(outDirPlots);
end

figHandles = findall(0,'Type','figure');

for i = 1:numel(figHandles)
    h = figHandles(i);
    % Use figure name if it exists, else a generic name
    figName = get(h,'Name');
    if isempty(figName)
        figName = sprintf('Figure_%02d', i);
    end
    % Sanitize filename (remove illegal characters)
    figName = regexprep(figName,'[^\w\-]','_');

    % Save as .fig and .png
    saveas(h, fullfile(outDirPlots, [figName '.fig']));
    saveas(h, fullfile(outDirPlots, [figName '.png']));
end

%% ========================================================================
%                               LOCAL HELPERS
% ========================================================================

function out = local_mapChoicesToLabels(x)
    if iscell(x)
        x = x(:);
    end
    if isnumeric(x)
        out = strings(numel(x),1); out(:) = missing;
        out(x==1) = "BUY"; out(x==2) = "HOLD"; out(x==3) = "SELL";
    elseif isstring(x)
        out = upper(x(:));
    elseif iscellstr(x)
        out = upper(string(x(:)));
    else
        out = strings(numel(x),1); out(:) = missing;
    end
end

function Z = local_get_Z_fr_allNeurons(C, f, target_event, target_win, baseline_event, baseline_win)
    target_interval   = (target_win ./ 1000) ./ C.dt;
    baseline_interval = (baseline_win ./ 1000) ./ C.dt;
    eps_std = 1e-6;

    try
        target_anchor   = C.eventTables{f}.(target_event){1,1}(:);
        target_ts       = target_anchor + target_interval(1);
        target_te       = target_anchor + target_interval(2);

        baseline_anchor = C.eventTables{f}.(baseline_event){1,1}(:);
        baseline_ts     = baseline_anchor + baseline_interval(1);
        baseline_te     = baseline_anchor + baseline_interval(2);
    catch
        Z = [];
        return;
    end

    Ssel     = C.S{f};
    nNeurons = size(Ssel, 1);
    nTrials  = numel(target_ts);

    Z = nan(nNeurons, nTrials);

    baseline_values = nan(nNeurons, nTrials);
    durBase = max((baseline_win(2) - baseline_win(1))/1000/C.dt, eps_std);
    for t = 1:nTrials
        a = baseline_ts(t); b = baseline_te(t);
        if ~isnan(a) && ~isnan(b) && b >= a
            baseline_values(:,t) = sum(Ssel(:,a:b),2) ./ durBase;
        end
    end
    baseline_mean = mean(baseline_values, 2, 'omitnan');
    baseline_std  = std(baseline_values, 0, 2, 'omitnan');
    baseline_std  = max(baseline_std, eps_std);

    durTarget = max((target_win(2) - target_win(1))/1000/C.dt, eps_std);
    for t = 1:nTrials
        a = target_ts(t); b = target_te(t);
        if ~isnan(a) && ~isnan(b) && b >= a
            fr_raw = sum(Ssel(:,a:b),2) ./ durTarget;
            Z(:,t) = (fr_raw - baseline_mean) ./ baseline_std;
        end
    end
end

function [thetaHat, negLL, nParams] = local_fit_one_model(modelName, sessData, BUY_NON_BUY_ONLY, NUM_RESTARTS, optimOpts)
    [theta0, lb, ub, objFun] = local_get_model_spec(modelName, sessData, BUY_NON_BUY_ONLY);
    nParams = numel(theta0);
    bestLL = inf; bestTheta = theta0;

    for r = 1:max(1,NUM_RESTARTS)
        if all(isfinite(lb)) && all(isfinite(ub))
            x0 = lb + rand(1,nParams).*(ub - lb);
        else
            x0 = theta0(:).';
        end
        try
            if exist('fmincon','file') == 2 && all(isfinite(lb)) && all(isfinite(ub))
                [theta, fval] = fmincon(objFun, x0, [],[],[],[], lb, ub, [], optimOpts);
            else
                [theta, fval] = fminsearch(objFun, x0);
            end
        catch
            [theta, fval] = fminsearch(objFun, x0);
        end
        if fval < bestLL
            bestLL = fval; bestTheta = theta;
        end
    end

    thetaHat = bestTheta(:).';
    negLL    = bestLL;
end

function [theta0, lb, ub, objFun] = local_get_model_spec(modelName, sessData, BUY_NON_BUY_ONLY)
    switch modelName
        case 'SelfOnly'
            % theta = [alpha, beta, omega]
            theta0 = [0.2, 3.0, 0.5];
            lb     = [0.0, 0.0, 0.0];
            ub     = [1.0, 20.0, 1.0];
            objFun = @(theta) local_rlNegLogLik_SelfOnly(theta, sessData, BUY_NON_BUY_ONLY);

        case 'SelfOnlyBuyBias'
            % theta = [alpha, beta, omega, kappaBuy]
            theta0 = [0.2, 3.0, 0.5, 0.0];
            lb     = [0.0, 0.0, 0.0, -5.0];
            ub     = [1.0, 20.0, 1.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_SelfOnlyBuyBias(theta, sessData, BUY_NON_BUY_ONLY);

        case 'OppOnly'
            % theta = [alpha_opp, beta]
            theta0 = [0.2, 3.0];
            lb     = [0.0, 0.0];
            ub     = [1.0, 20.0];
            objFun = @(theta) local_rlNegLogLik_OppOnly(theta, sessData, BUY_NON_BUY_ONLY);

        case 'OppOnlyBuyBias'
            % theta = [alpha_opp, beta, kappaOppBuy]
            theta0 = [0.2, 3.0, 0.0];
            lb     = [0.0, 0.0, -5.0];
            ub     = [1.0, 20.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_OppOnlyBuyBias(theta, sessData, BUY_NON_BUY_ONLY);

        case 'SelfOpp'
            % theta = [alpha_self, alpha_opp, beta, omega_self, gamma, kappaBuy]
            % omega_self in [0,1] weights R2 (dividend) for the subject.
            % Opponent portfolio is z-scored within-session; gamma scales Q_opp in choice.
            theta0 = [0.2, 0.2, 3.0, 0.5, 0.5, 0.0];
            lb     = [0.0, 0.0, 0.0, 0.0, -5.0, -5.0];
            ub     = [1.0, 1.0, 20.0, 1.0,  5.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_SelfOpp(theta, sessData, BUY_NON_BUY_ONLY);

        case 'SelfBubble'
            % theta = [alpha, beta, omega_nonBub, omega_bub, kappa_nonBub, kappa_bub]
            theta0 = [0.2, 3.0, 0.5, 0.5, 0.0, 0.0];
            lb     = [0.0, 0.0, 0.0, 0.0, -5.0, -5.0];
            ub     = [1.0, 20.0, 1.0, 1.0,  5.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_SelfBubble(theta, sessData, BUY_NON_BUY_ONLY);

        otherwise
            error('Unknown modelName: %s', modelName);
    end
end

%%
function [negLL, latents] = local_rlNegLogLik_SelfOnlyBuyBias(theta, sessData, buyNonBuyOnly)
    alpha = max(min(theta(1),1),0);
    beta  = theta(2);
    omega = min(max(theta(3),0),1);   % clamp omega to [0,1]
    kappa = theta(4);

    if buyNonBuyOnly
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);

    choiceStr = string(sessData.subjectChoice);
    actIdx = zeros(height(sessData),1);
    for a = 1:nA
        actIdx(choiceStr == actionLabels{a}) = a;
    end
    valid = actIdx > 0 & ~isnan(actIdx);

    blockID = sessData.BlockID;
    uBlocks = unique(blockID);

    T = height(sessData);
    Qhist = nan(T, nA);
    value = nan(T, 1);
    delta = nan(T, 1);

    negLL = 0;
    tCounter = 0;

    for b = uBlocks.'
        idxBlock = find(blockID == b);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            tCounter = tCounter + 1;

            if ~valid(t)
                Qhist(tCounter,:) = Q;
                continue;
            end

            a = actIdx(t);
            Rc = sessData.choiceOutcomeReward(t);             % R1: choice-linked reward
            Rp = sessData.portfolioRewardFromDividend(t);     % R2: dividend-based reward
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = (1 - omega) * Rc + omega * Rp;               % convex combination in raw units

            Vt = Q(a);
            value(tCounter) = Vt;
            delta(tCounter) = r - Vt;

            logits = beta * Q;
            logits(1) = logits(1) + kappa;

            maxLogit = max(logits);
            p = exp(logits - maxLogit); p = p ./ sum(p);
            pChosen = max(p(a), realmin);

            negLL = negLL - log(pChosen);

            Q(a) = Q(a) + alpha * (r - Vt);
            Qhist(tCounter,:) = Q;
        end
    end

    if nargout > 1
        latents = struct('value',value,'delta',delta,'Q',Qhist);
    end
end

function [negLL, latents] = local_rlNegLogLik_SelfOnly(theta, sessData, buyNonBuyOnly)
    alpha = max(min(theta(1),1),0);
    beta  = theta(2);
    omega = min(max(theta(3),0),1);

    if buyNonBuyOnly
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);

    choiceStr = string(sessData.subjectChoice);
    actIdx = zeros(height(sessData),1);
    for a = 1:nA
        actIdx(choiceStr == actionLabels{a}) = a;
    end
    valid = actIdx > 0 & ~isnan(actIdx);

    blockID = sessData.BlockID;
    uBlocks = unique(blockID);

    T = height(sessData);
    Qhist = nan(T, nA);
    value = nan(T, 1);
    delta = nan(T, 1);

    negLL = 0;
    tCounter = 0;

    for b = uBlocks.'
        idxBlock = find(blockID == b);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            tCounter = tCounter + 1;

            if ~valid(t)
                Qhist(tCounter,:) = Q;
                continue;
            end

            a = actIdx(t);
            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = (1 - omega) * Rc + omega * Rp;

            Vt = Q(a);
            value(tCounter) = Vt;
            delta(tCounter) = r - Vt;

            logits = beta * Q;      % no kappa term
            maxLogit = max(logits);
            p = exp(logits - maxLogit); p = p ./ sum(p);
            pChosen = max(p(a), realmin);

            negLL = negLL - log(pChosen);

            Q(a) = Q(a) + alpha * (r - Vt);
            Qhist(tCounter,:) = Q;
        end
    end

    if nargout > 1
        latents = struct('value',value,'delta',delta,'Q',Qhist);
    end
end

function [negLL, latents] = local_rlNegLogLik_OppOnly(theta, sessData, buyNonBuyOnly)
    % OppOnly:
    %   - Track Q over opponent's actions using opponent portfolio reward
    %   - Subject choice policy depends only on Q_opp
    %   - theta = [alpha_opp, beta]

    alpha_opp = max(min(theta(1),1),0);
    beta      = theta(2);

    if buyNonBuyOnly
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);

    T = height(sessData);

    % Subject choices (target for likelihood)
    choiceStrSubj = string(sessData.subjectChoice);
    actIdxSubj = zeros(T,1);
    for a = 1:nA
        actIdxSubj(choiceStrSubj == actionLabels{a}) = a;
    end
    validSubj = actIdxSubj > 0 & ~isnan(actIdxSubj);

    % Opponent choices (drive Q updates)
    choiceStrOpp = string(sessData.opponentChoice);
    actIdxOpp = zeros(T,1);
    for a = 1:nA
        actIdxOpp(choiceStrOpp == actionLabels{a}) = a;
    end
    validOpp = actIdxOpp > 0 & ~isnan(actIdxOpp);

    % Reward: normalized opponent portfolio value
    Ropp_raw = sessData.opponentPortfolioValue;
    Ropp_raw(~isfinite(Ropp_raw)) = NaN;
    validR = ~isnan(Ropp_raw);
    if any(validR)
        muR    = mean(Ropp_raw(validR));
        sigmaR = std(Ropp_raw(validR));
        if sigmaR < 1e-6
            sigmaR = 1.0;
        end
        RoppNorm = (Ropp_raw - muR) ./ sigmaR;
    else
        RoppNorm = zeros(size(Ropp_raw));
    end

    blockID = sessData.BlockID;
    uBlocks = unique(blockID);

    Qhist = nan(T, nA);
    value = nan(T, 1);
    delta = nan(T, 1);

    negLL   = 0;
    tCounter = 0;

    for b = uBlocks.'
        idxBlock  = find(blockID == b);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            tCounter = tCounter + 1;

            % Likelihood for subject choice given current Q_opp
            if validSubj(t)
                aSubj = actIdxSubj(t);

                Vt = Q(aSubj);
                value(tCounter) = Vt;

                logits = beta * Q;
                maxLogit = max(logits);
                p = exp(logits - maxLogit); 
                p = p ./ sum(p);
                pChosen = max(p(aSubj), realmin);

                negLL = negLL - log(pChosen);
            end

            % RL update from opponent's action and portfolio reward
            if validOpp(t)
                aOpp = actIdxOpp(t);
                r    = RoppNorm(t);
                if isnan(r), r = 0; end

                Vt_opp = Q(aOpp);
                delta(tCounter) = r - Vt_opp;

                Q(aOpp) = Q(aOpp) + alpha_opp * (r - Vt_opp);
            end

            Qhist(tCounter,:) = Q;
        end
    end

    if nargout > 1
        latents = struct('value',value,'delta',delta,'Q',Qhist);
    end
end

function [negLL, latents] = local_rlNegLogLik_OppOnlyBuyBias(theta, sessData, buyNonBuyOnly)
    % OppOnlyBuyBias:
    %   - Same as OppOnly, but add kappaOppBuy bias on BUY logit
    %   - theta = [alpha_opp, beta, kappaOppBuy]

    alpha_opp   = max(min(theta(1),1),0);
    beta        = theta(2);
    kappaOppBuy = theta(3);

    if buyNonBuyOnly
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA     = numel(actionLabels);
    buyIdx = 1;  % BUY index

    T = height(sessData);

    % Subject choices (target for likelihood)
    choiceStrSubj = string(sessData.subjectChoice);
    actIdxSubj = zeros(T,1);
    for a = 1:nA
        actIdxSubj(choiceStrSubj == actionLabels{a}) = a;
    end
    validSubj = actIdxSubj > 0 & ~isnan(actIdxSubj);

    % Opponent choices (drive Q updates)
    choiceStrOpp = string(sessData.opponentChoice);
    actIdxOpp = zeros(T,1);
    for a = 1:nA
        actIdxOpp(choiceStrOpp == actionLabels{a}) = a;
    end
    validOpp = actIdxOpp > 0 & ~isnan(actIdxOpp);

    % Reward: normalized opponent portfolio value
    Ropp_raw = sessData.opponentPortfolioValue;
    Ropp_raw(~isfinite(Ropp_raw)) = NaN;
    validR = ~isnan(Ropp_raw);
    if any(validR)
        muR    = mean(Ropp_raw(validR));
        sigmaR = std(Ropp_raw(validR));
        if sigmaR < 1e-6
            sigmaR = 1.0;
        end
        RoppNorm = (Ropp_raw - muR) ./ sigmaR;
    else
        RoppNorm = zeros(size(Ropp_raw));
    end

    blockID = sessData.BlockID;
    uBlocks = unique(blockID);

    Qhist = nan(T, nA);
    value = nan(T, 1);
    delta = nan(T, 1);

    negLL   = 0;
    tCounter = 0;

    for b = uBlocks.'
        idxBlock  = find(blockID == b);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            tCounter = tCounter + 1;

            % Likelihood for subject choice given current Q_opp (+ BUY bias)
            if validSubj(t)
                aSubj = actIdxSubj(t);

                Vt = Q(aSubj);
                value(tCounter) = Vt;

                logits = beta * Q;
                logits(buyIdx) = logits(buyIdx) + kappaOppBuy;

                maxLogit = max(logits);
                p = exp(logits - maxLogit); 
                p = p ./ sum(p);
                pChosen = max(p(aSubj), realmin);

                negLL = negLL - log(pChosen);
            end

            % RL update from opponent's action and portfolio reward
            if validOpp(t)
                aOpp = actIdxOpp(t);
                r    = RoppNorm(t);
                if isnan(r), r = 0; end

                Vt_opp = Q(aOpp);
                delta(tCounter) = r - Vt_opp;

                Q(aOpp) = Q(aOpp) + alpha_opp * (r - Vt_opp);
            end

            Qhist(tCounter,:) = Q;
        end
    end

    if nargout > 1
        latents = struct('value',value,'delta',delta,'Q',Qhist);
    end
end

function [negLL, latents] = local_rlNegLogLik_SelfBubble(theta, sessData, buyNonBuyOnly)
    alpha          = max(min(theta(1),1),0);
    beta           = theta(2);
    omega_nonBub   = min(max(theta(3),0),1);  % weight on R2 in non-bubble markets
    omega_bub      = min(max(theta(4),0),1);  % weight on R2 in bubble markets
    kappa_nonBub   = theta(5);
    kappa_bub      = theta(6);

    if buyNonBuyOnly
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);

    choiceStr = string(sessData.subjectChoice);
    actIdx = zeros(height(sessData),1);
    for a = 1:nA
        actIdx(choiceStr == actionLabels{a}) = a;
    end
    valid = actIdx > 0 & ~isnan(actIdx);

    bubbleFlag = sessData.BubbleFlag;
    blockID = sessData.BlockID;
    uBlocks = unique(blockID);

    T = height(sessData);
    Qhist = nan(T, nA);
    value = nan(T, 1);
    delta = nan(T, 1);

    negLL = 0; tCounter = 0;

    for b = uBlocks.'
        idxBlock = find(blockID == b);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k); tCounter = tCounter + 1;

            if ~valid(t)
                Qhist(tCounter,:) = Q; continue;
            end

            a = actIdx(t);
            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            
            if bubbleFlag(t) == 1
                omega = omega_bub;  kappa = kappa_bub;
            else
                omega = omega_nonBub; kappa = kappa_nonBub;
            end
            
            r = (1 - omega) * Rc + omega * Rp;

            Vt = Q(a);
            value(tCounter) = Vt;
            delta(tCounter) = r - Vt;

            logits = beta * Q;
            logits(1) = logits(1) + kappa;

            maxLogit = max(logits);
            p = exp(logits - maxLogit); p = p ./ sum(p);
            pChosen = max(p(a), realmin);

            negLL = negLL - log(pChosen);

            Q(a) = Q(a) + alpha * (r - Vt);
            Qhist(tCounter,:) = Q;
        end
    end

    if nargout > 1
        latents = struct('value',value,'delta',delta,'Q',Qhist);
    end
end

function [negLL, latentsSelf, latentsOpp] = local_rlNegLogLik_SelfOpp(theta, sessData, buyNonBuyOnly)
    % theta = [alpha_self, alpha_opp, beta, omega_self, gamma, kappaBuy]
    alpha_self = max(min(theta(1),1),0);
    alpha_opp  = max(min(theta(2),1),0);
    beta       = theta(3);
    omega_self = min(max(theta(4),0),1);   % subject R2 weight
    gamma      = theta(5);                 % coupling from Q_opp into choice
    kappa_buy  = theta(6);

    if buyNonBuyOnly
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);

    choiceStrSelf = string(sessData.subjectChoice);
    actIdxSelf = zeros(height(sessData),1);
    for a = 1:nA
        actIdxSelf(choiceStrSelf == actionLabels{a}) = a;
    end
    validSelf = actIdxSelf > 0 & ~isnan(actIdxSelf);

    choiceStrOpp = string(sessData.opponentChoice);
    actIdxOpp = zeros(height(sessData),1);
    for a = 1:nA
        actIdxOpp(choiceStrOpp == actionLabels{a}) = a;
    end
    validOpp = actIdxOpp > 0 & ~isnan(actIdxOpp);

    blockID = sessData.BlockID;
    uBlocks = unique(blockID);

    T = height(sessData);
    Q_self_hist = nan(T, nA);
    Q_opp_hist  = nan(T, nA);
    valueSelf   = nan(T,1);
    deltaSelf   = nan(T,1);
    valueOpp    = nan(T,1);
    deltaOpp    = nan(T,1);

    % ---------- normalize opponentPortfolioValue within-session ----------
    Ropp_raw = sessData.opponentPortfolioValue;
    Ropp_raw(~isfinite(Ropp_raw)) = NaN;
    validR  = ~isnan(Ropp_raw);
    if any(validR)
        muR    = mean(Ropp_raw(validR));
        sigmaR = std(Ropp_raw(validR));
        if sigmaR < 1e-6
            sigmaR = 1.0;
        end
        RoppNorm = (Ropp_raw - muR) ./ sigmaR;
    else
        RoppNorm = zeros(size(Ropp_raw));
    end
    % --------------------------------------------------------------------

    negLL = 0; 
    tCounter = 0;

    for b = uBlocks.'
        idxBlock  = find(blockID == b);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q_self = zeros(nA,1);
        Q_opp  = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k); 
            tCounter = tCounter + 1;

            % ------------------ subject update + likelihood ------------------
            if validSelf(t)
                aSelf = actIdxSelf(t);

                Rc = sessData.choiceOutcomeReward(t);
                Rp = sessData.portfolioRewardFromDividend(t);
                if isnan(Rc), Rc = 0; end
                if isnan(Rp), Rp = 0; end

                rSelf = (1 - omega_self) * Rc + omega_self * Rp;

                VtSelf = Q_self(aSelf);
                valueSelf(tCounter) = VtSelf;
                deltaSelf(tCounter) = rSelf - VtSelf;

                V = Q_self + gamma * Q_opp;        % integrate opponent Q
                logits = beta * V;
                logits(1) = logits(1) + kappa_buy; % BUY bias

                maxLogit = max(logits);
                p = exp(logits - maxLogit); 
                p = p ./ sum(p);
                pChosen = max(p(aSelf), realmin);

                negLL = negLL - log(pChosen);

                Q_self(aSelf) = Q_self(aSelf) + alpha_self * (rSelf - VtSelf);
            end

            % ------------------ opponent update (normalized Ropp) -----------
            if validOpp(t)
                aOpp = actIdxOpp(t);
                if aOpp > 0
                    Rp_opp = RoppNorm(t);
                    if isnan(Rp_opp), Rp_opp = 0; end
                    r_opp = Rp_opp;

                    VtOpp = Q_opp(aOpp);
                    valueOpp(tCounter) = VtOpp;
                    deltaOpp(tCounter) = r_opp - VtOpp;

                    Q_opp(aOpp) = Q_opp(aOpp) + alpha_opp * (r_opp - VtOpp);
                end
            end

            Q_self_hist(tCounter,:) = Q_self;
            Q_opp_hist(tCounter,:)  = Q_opp;
        end
    end

    if nargout > 1
        latentsSelf = struct('value',valueSelf,'delta',deltaSelf,'Q',Q_self_hist);
        latentsOpp  = struct('value',valueOpp, 'delta',deltaOpp, 'Q',Q_opp_hist);
    end
end

%%
function behaviorBootstrap = local_rl_parametric_block_bootstrap( ...
    trialTable, behaviorFitSummary, BUY_NON_BUY_ONLY, modelsForBootstrap, N_BOOT, FIT_SCOPE)

    if nargin < 6 || isempty(FIT_SCOPE), FIT_SCOPE = 'per_monkey'; end
    if nargin < 5 || isempty(N_BOOT),   N_BOOT   = 200; end
    if nargin < 4 || isempty(modelsForBootstrap)
        modelsForBootstrap = {'SelfOnlyBuyBias','SelfBubble','SelfOpp'};
    end

    rng(0);

    behaviorBootstrap = table( ...
        string.empty(0,1), ...  % SetName
        string.empty(0,1), ...  % Event
        string.empty(0,1), ...  % Monkey
        string.empty(0,1), ...  % Condition
        string.empty(0,1), ...  % DateID
        string.empty(0,1), ...  % ModelName
        cell(0,1), ...          % ThetaHat
        cell(0,1), ...          % ThetaBoot
        double.empty(0,1), ...  % NumBoot
        'VariableNames', {'SetName','Event','Monkey','Condition','DateID', ...
                          'ModelName','ThetaHat','ThetaBoot','NumBoot'});

    if isempty(trialTable) || isempty(behaviorFitSummary)
        return;
    end

    if ~ismember('BlockID', trialTable.Properties.VariableNames)
        error('trialTable must contain BlockID.');
    end
    if ~ismember('TrialInBlock', trialTable.Properties.VariableNames)
        error('trialTable must contain TrialInBlock.');
    end

    switch lower(FIT_SCOPE)
        case 'per_monkey'
            groupVars = {'Monkey','Condition','Date'};
        case 'pooled'
            groupVars = {'Condition','Date'};
        otherwise
            error('Unknown FIT_SCOPE: %s', FIT_SCOPE);
    end

    [mcGroup, mcKey] = findgroups(trialTable(:, groupVars));
    nGroups = height(mcKey);

    for g = 1:nGroups
        if ismember('Monkey', groupVars)
            monkeyID    = char(mcKey.Monkey(g));
        else
            monkeyID    = 'ALL';
        end
        conditionID = char(mcKey.Condition(g));
        dateID   = mcKey.Date(g);

        rowMaskMC  = (mcGroup == g);
        sessDataMC = sortrows(trialTable(rowMaskMC,:), {'Date','BlockID','TrialInBlock'});

        for m = 1:numel(modelsForBootstrap)
            modelName = modelsForBootstrap{m};

            fitMask = strcmp(behaviorFitSummary.Monkey,   monkeyID) & ...
                      strcmp(behaviorFitSummary.Condition, conditionID) & ...
                      behaviorFitSummary.DateID == dateID & ...
                      strcmp(behaviorFitSummary.ModelName, modelName);
            if ~any(fitMask)
                continue;
            end

            fitRow   = behaviorFitSummary(find(fitMask,1), :);
            thetaHat = fitRow.Theta{1};
            P        = numel(thetaHat);

            thetaBoot = nan(N_BOOT, P);
    
            switch modelName
                case 'SelfOnly'
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_SelfOnly(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_SelfOnly(th, sd, BUY_NON_BUY_ONLY);
                    lb = [0.0, 0.0, 0]; 
                    ub = [1.0, 20.0, 1.0];
            
                case 'SelfOnlyBuyBias'
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_SelfOnlyBuyBias(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_SelfOnlyBuyBias(th, sd, BUY_NON_BUY_ONLY);
                    lb = [0.0, 0.0, 0, -5.0]; 
                    ub = [1.0, 20.0, 1.0, 5.0];            

                case 'OppOnly'
                    % theta = [alpha_opp, beta]
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_OppOnly(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_OppOnly(th, sd, BUY_NON_BUY_ONLY);
                    lb = [0.0, 0.0]; 
                    ub = [1.0, 20.0];

                case 'OppOnlyBuyBias'
                    % theta = [alpha_opp, beta, kappaOppBuy]
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_OppOnlyBuyBias(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_OppOnlyBuyBias(th, sd, BUY_NON_BUY_ONLY);
                    lb = [0.0, 0.0, -5.0]; 
                    ub = [1.0, 20.0,  5.0];

                case 'SelfBubble'
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_SelfBubble(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_SelfBubble(th, sd, BUY_NON_BUY_ONLY);
                    lb = [0.0, 0.0, 0.0, 0.0, -5.0, -5.0]; 
                    ub = [1.0, 20.0, 1.0, 1.0, 5.0, 5.0];
            
                case 'SelfOpp'
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_SelfOpp(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_SelfOpp(th, sd, BUY_NON_BUY_ONLY);
                    % theta = [alpha_self, alpha_opp, beta, omega_self, gamma, kappaBuy]
                    lb = [0.0, 0.0, 0.0, 0.0, -5.0, -5.0]; 
                    ub = [1.0, 1.0, 20.0, 1.0,  5.0,  5.0];
            
                otherwise
                    error('Unknown modelName: %s', modelName);
            end

            parfor b = 1:N_BOOT
                simChoices = simFun(thetaHat, sessDataMC);
                sessDataSim = sessDataMC;
                sessDataSim.subjectChoice = categorical(simChoices, categories(sessDataMC.subjectChoice));

                objFun = objFunMaker(sessDataSim);
                x0 = thetaHat(:).';

                try
                    if exist('fmincon','file') == 2
                        thetaFit = fmincon(objFun, x0, [],[],[],[], lb, ub, []);
                    else
                        thetaFit = fminsearch(objFun, x0);
                    end
                catch
                    thetaFit = fminsearch(objFun, x0);
                end
                thetaBoot(b,:) = thetaFit(:).';
            end

            newRow = table( ...
                string(fitRow.SetName), ...
                string(fitRow.Event), ...
                string(monkeyID), ...
                string(conditionID), ...
                string(dateID), ...
                string(modelName), ...
                {thetaHat}, ...
                {thetaBoot}, ...
                N_BOOT, ...
                'VariableNames', behaviorBootstrap.Properties.VariableNames);
            behaviorBootstrap = [behaviorBootstrap; newRow]; %#ok<AGROW>
        end
    end
end

function simChoices = local_simulateChoices_SelfOnly(theta, sessData, BUY_NON_BUY_ONLY)
    % theta = [alpha, beta, omega]
    alpha = max(min(theta(1),1),0);
    beta  = theta(2);
    omega = theta(3);

    if BUY_NON_BUY_ONLY
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);

    T = height(sessData);
    origChoiceStr = string(sessData.subjectChoice);
    valid = ~ismissing(origChoiceStr);

    simChoices = strings(T,1); 
    simChoices(:) = missing;

    blockID = sessData.BlockID; 
    uBlocks = unique(blockID,'stable');

    for bb = 1:numel(uBlocks)
        bID      = uBlocks(bb);
        idxBlock = find(blockID == bID);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            if ~valid(t), continue; end

            % Softmax over Q, no kappa term
            logits = beta * Q;
            maxLogit = max(logits);
            p = exp(logits - maxLogit); 
            p = p / sum(p);

            u = rand; 
            cum = cumsum(p);
            aSim = find(u <= cum, 1); 
            if isempty(aSim), aSim = 1; end
            simChoices(t) = string(actionLabels{aSim});

            % Reward and Q update
            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = (1 - omega) * Rc + omega * Rp;

            Vt = Q(aSim);
            Q(aSim) = Q(aSim) + alpha * (r - Vt);
        end
    end
end

function simChoices = local_simulateChoices_SelfOnlyBuyBias(theta, sessData, BUY_NON_BUY_ONLY)
    alpha = max(min(theta(1),1),0);
    beta  = theta(2);
    omega = theta(3);
    kappa = theta(4);

    if BUY_NON_BUY_ONLY
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);
    buyIdx = 1;

    T = height(sessData);
    origChoiceStr = string(sessData.subjectChoice);
    valid = ~ismissing(origChoiceStr);

    simChoices = strings(T,1); simChoices(:) = missing;
    blockID = sessData.BlockID; uBlocks = unique(blockID,'stable');

    for bb = 1:numel(uBlocks)
        bID = uBlocks(bb);
        idxBlock = find(blockID == bID);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);
        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            if ~valid(t), continue; end

            logits = beta * Q; logits(buyIdx) = logits(buyIdx) + kappa;
            maxLogit = max(logits);
            p = exp(logits - maxLogit); p = p / sum(p);
            u = rand; cum = cumsum(p);
            aSim = find(u <= cum, 1); if isempty(aSim), aSim = buyIdx; end
            simChoices(t) = string(actionLabels{aSim});

            % Reward and Q update
            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = (1 - omega) * Rc + omega * Rp;

            Vt = Q(aSim);
            Q(aSim) = Q(aSim) + alpha * (r - Vt);
        end
    end
end

function simChoices = local_simulateChoices_SelfBubble(theta, sessData, BUY_NON_BUY_ONLY)
    alpha = max(min(theta(1),1),0);
    beta            = theta(2);
    omega_nonBubble = theta(3);
    omega_bubble    = theta(4);
    kappa_nonBubble = theta(5);
    kappa_bubble    = theta(6);

    if BUY_NON_BUY_ONLY
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);
    buyIdx = 1;

    T = height(sessData);
    origChoiceStr = string(sessData.subjectChoice);
    valid = ~ismissing(origChoiceStr);

    simChoices = strings(T,1); simChoices(:) = missing;
    blockID = sessData.BlockID; uBlocks = unique(blockID,'stable');

    for bb = 1:numel(uBlocks)
        bID = uBlocks(bb);
        idxBlock = find(blockID == bID);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            if ~valid(t), continue; end

            isBubble = logical(sessData.BubbleFlag(t));
            if isBubble
                omega = omega_bubble; kappa = kappa_bubble;
            else
                omega = omega_nonBubble; kappa = kappa_nonBubble;
            end

            logits = beta * Q; logits(buyIdx) = logits(buyIdx) + kappa;
            maxLogit = max(logits);
            p = exp(logits - maxLogit); p = p / sum(p);
            u = rand; cum = cumsum(p);
            aSim = find(u <= cum, 1); if isempty(aSim), aSim = buyIdx; end

            simChoices(t) = string(actionLabels{aSim});

            % Reward and Q update
            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = (1 - omega) * Rc + omega * Rp;

            Vt = Q(aSim);
            Q(aSim) = Q(aSim) + alpha * (r - Vt);
        end
    end
end

function simChoices = local_simulateChoices_SelfOpp(theta, sessData, BUY_NON_BUY_ONLY)
    % theta = [alpha_self, alpha_opp, beta, omega_self, gamma, kappaBuy]
    alpha_self = max(min(theta(1),1),0);
    alpha_opp  = max(min(theta(2),1),0);
    beta       = theta(3);
    omega_self = min(max(theta(4),0),1);
    gamma      = theta(5);
    kappaBuy   = theta(6);

    if BUY_NON_BUY_ONLY
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA     = numel(actionLabels);
    buyIdx = 1;

    T = height(sessData);
    origChoiceStr = string(sessData.subjectChoice);
    valid = ~ismissing(origChoiceStr);

    oppChoiceStr = string(sessData.opponentChoice);
    oppValid     = ~ismissing(oppChoiceStr);

    simChoices = strings(T,1); 
    simChoices(:) = missing;
    blockID = sessData.BlockID; 
    uBlocks = unique(blockID,'stable');

    % ---------- normalize opponentPortfolioValue within-session ----------
    Ropp_raw = sessData.opponentPortfolioValue;
    Ropp_raw(~isfinite(Ropp_raw)) = NaN;
    validR  = ~isnan(Ropp_raw);
    if any(validR)
        muR    = mean(Ropp_raw(validR));
        sigmaR = std(Ropp_raw(validR));
        if sigmaR < 1e-6
            sigmaR = 1.0;
        end
        RoppNorm = (Ropp_raw - muR) ./ sigmaR;
    else
        RoppNorm = zeros(size(Ropp_raw));
    end
    % --------------------------------------------------------------------

    for bb = 1:numel(uBlocks)
        bID      = uBlocks(bb);
        idxBlock = find(blockID == bID);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q_self = zeros(nA,1);
        Q_opp  = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            if ~valid(t), continue; end

            % choice policy using current Q_self and Q_opp
            V = Q_self + gamma * Q_opp;
            logits = beta * V;
            logits(buyIdx) = logits(buyIdx) + kappaBuy;

            maxLogit = max(logits);
            p = exp(logits - maxLogit); 
            p = p / sum(p);

            u = rand; 
            cum = cumsum(p);
            aSim = find(u <= cum, 1); 
            if isempty(aSim), aSim = buyIdx; end
            simChoices(t) = string(actionLabels{aSim});

            % subject reward
            Rc = sessData.choiceOutcomeReward(t);
            Rp_self = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp_self), Rp_self = 0; end
            r_self = (1 - omega_self) * Rc + omega_self * Rp_self;

            Vt_self = Q_self(aSim);
            Q_self(aSim) = Q_self(aSim) + alpha_self * (r_self - Vt_self);

            % opponent learning on normalized portfolio
            if oppValid(t)
                oppStr = oppChoiceStr(t);
                aOpp = 0;
                for a = 1:nA
                    if oppStr == actionLabels{a}
                        aOpp = a; 
                        break; 
                    end
                end
                if aOpp > 0
                    Rp_opp = RoppNorm(t);
                    if isnan(Rp_opp), Rp_opp = 0; end
                    r_opp = Rp_opp;

                    Vt_opp = Q_opp(aOpp);
                    Q_opp(aOpp) = Q_opp(aOpp) + alpha_opp * (r_opp - Vt_opp);
                end
            end
        end
    end
end

function simChoices = local_simulateChoices_OppOnly(theta, sessData, BUY_NON_BUY_ONLY)
    % theta = [alpha_opp, beta]
    alpha_opp = max(min(theta(1),1),0);
    beta      = theta(2);

    if BUY_NON_BUY_ONLY
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA = numel(actionLabels);

    T = height(sessData);
    simChoices = strings(T,1); 
    simChoices(:) = missing;

    % Opponent choices
    choiceStrOpp = string(sessData.opponentChoice);
    actIdxOpp = zeros(T,1);
    for a = 1:nA
        actIdxOpp(choiceStrOpp == actionLabels{a}) = a;
    end
    oppValid = actIdxOpp > 0 & ~isnan(actIdxOpp);

    % Reward: normalized opponent portfolio value
    Ropp_raw = sessData.opponentPortfolioValue;
    Ropp_raw(~isfinite(Ropp_raw)) = NaN;
    validR = ~isnan(Ropp_raw);
    if any(validR)
        muR    = mean(Ropp_raw(validR));
        sigmaR = std(Ropp_raw(validR));
        if sigmaR < 1e-6
            sigmaR = 1.0;
        end
        RoppNorm = (Ropp_raw - muR) ./ sigmaR;
    else
        RoppNorm = zeros(size(Ropp_raw));
    end

    blockID = sessData.BlockID;
    uBlocks = unique(blockID,'stable');

    for bb = 1:numel(uBlocks)
        bID      = uBlocks(bb);
        idxBlock = find(blockID == bID);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);

            % Subject choice from Q_opp
            logits = beta * Q;
            maxLogit = max(logits);
            p = exp(logits - maxLogit); 
            p = p ./ sum(p);

            u = rand;
            cum = cumsum(p);
            aSim = find(u <= cum, 1);
            if isempty(aSim), aSim = 1; end
            simChoices(t) = string(actionLabels{aSim});

            % RL update from opponent's action and portfolio reward
            if oppValid(t)
                aOpp = actIdxOpp(t);
                r    = RoppNorm(t);
                if isnan(r), r = 0; end

                Vt_opp = Q(aOpp);
                Q(aOpp) = Q(aOpp) + alpha_opp * (r - Vt_opp);
            end
        end
    end
end

function simChoices = local_simulateChoices_OppOnlyBuyBias(theta, sessData, BUY_NON_BUY_ONLY)
    % theta = [alpha_opp, beta, kappaOppBuy]
    alpha_opp   = max(min(theta(1),1),0);
    beta        = theta(2);
    kappaOppBuy = theta(3);

    if BUY_NON_BUY_ONLY
        actionLabels = {'BUY','Non-BUY'};
    else
        actionLabels = {'BUY','HOLD','SELL'};
    end
    nA     = numel(actionLabels);
    buyIdx = 1;

    T = height(sessData);
    simChoices = strings(T,1); 
    simChoices(:) = missing;

    % Opponent choices
    choiceStrOpp = string(sessData.opponentChoice);
    actIdxOpp = zeros(T,1);
    for a = 1:nA
        actIdxOpp(choiceStrOpp == actionLabels{a}) = a;
    end
    oppValid = actIdxOpp > 0 & ~isnan(actIdxOpp);

    % Reward: normalized opponent portfolio value
    Ropp_raw = sessData.opponentPortfolioValue;
    Ropp_raw(~isfinite(Ropp_raw)) = NaN;
    validR = ~isnan(Ropp_raw);
    if any(validR)
        muR    = mean(Ropp_raw(validR));
        sigmaR = std(Ropp_raw(validR));
        if sigmaR < 1e-6
            sigmaR = 1.0;
        end
        RoppNorm = (Ropp_raw - muR) ./ sigmaR;
    else
        RoppNorm = zeros(size(Ropp_raw));
    end

    blockID = sessData.BlockID;
    uBlocks = unique(blockID,'stable');

    for bb = 1:numel(uBlocks)
        bID      = uBlocks(bb);
        idxBlock = find(blockID == bID);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);

            % Subject choice from Q_opp with BUY bias
            logits = beta * Q;
            logits(buyIdx) = logits(buyIdx) + kappaOppBuy;

            maxLogit = max(logits);
            p = exp(logits - maxLogit); 
            p = p ./ sum(p);

            u = rand;
            cum = cumsum(p);
            aSim = find(u <= cum, 1);
            if isempty(aSim), aSim = buyIdx; end
            simChoices(t) = string(actionLabels{aSim});

            % RL update from opponent's action and portfolio reward
            if oppValid(t)
                aOpp = actIdxOpp(t);
                r    = RoppNorm(t);
                if isnan(r), r = 0; end

                Vt_opp = Q(aOpp);
                Q(aOpp) = Q(aOpp) + alpha_opp * (r - Vt_opp);
            end
        end
    end
end

%%
function paramTable = local_behaviorParamTable(fitTable, paramNamesPerModel)
    if isempty(fitTable)
        paramTable = table();
        return;
    end

    rows = table();

    for i = 1:height(fitTable)
        modelName = char(fitTable.ModelName(i));
        theta     = fitTable.Theta{i};
        P         = numel(theta);

        if isfield(paramNamesPerModel, modelName)
            names = paramNamesPerModel.(modelName);
            if numel(names) ~= P
                names = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
            end
        else
            names = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
        end

        for p = 1:P
            newRow = table( ...
                fitTable.SetName(i), ...
                fitTable.Event(i), ...
                fitTable.Monkey(i), ...
                fitTable.Condition(i), ...
                fitTable.DateID(i), ...
                fitTable.ModelName(i), ...
                uint32(p), ...
                string(names{p}), ...
                theta(p), ...
                'VariableNames', {'SetName','Event','Monkey','Condition','DateID', ...
                                  'ModelName','ParamIndex','ParamName','Value'});
            rows = [rows; newRow]; %#ok<AGROW>
        end
    end

    paramTable = rows;
end

function local_plotBehaviorParamComparisonBySet(fitTable, paramNamesPerModel)
% local_plotBehaviorParamComparisonBySet
% -------------------------------------------------------------------------
% New version: mean ± STD bars per condition instead of box+scatter.
%
% For each Set / Model / Param:
%   - X-axis: Condition (ordered by keywords)
%   - Y-axis: mean(parameter value across sessions)
%   - Error bars: ± STD (not SEM)
%   - X tick labels: "<Condition>\n(n=#sessions)"
%   - Horizontal dashed line at y = 0
%   - Pairwise significance (same logic as before) drawn on top.
%
% Also supports OT vs Saline paired plots (same layout as before).
% -------------------------------------------------------------------------

    if isempty(fitTable), return; end

    sets   = unique(fitTable.SetName, 'stable');
    events = unique(fitTable.Event,   'stable');
    if isempty(events)
        return;
    end
    eventUse = events(1);  % single-event pipeline; adjust if multi-event

    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    % ===== Precompute global y-limits per (Model, Param) across ALL sets =====
    modelListGlobal = unique(fitTable.ModelName, 'stable');
    yLimsByModel = struct();

    for miGlobal = 1:numel(modelListGlobal)
        modelNameGlobal = modelListGlobal(miGlobal);
        subAll = fitTable(fitTable.ModelName == modelNameGlobal & ...
                          fitTable.Event     == eventUse, :);
        if isempty(subAll)
            continue;
        end

        thetaExample = subAll.Theta{1};
        P = numel(thetaExample);

        yMin = inf(1,P);
        yMax = -inf(1,P);

        for r = 1:height(subAll)
            theta = subAll.Theta{r};
            if numel(theta) ~= P
                continue;
            end
            thetaRow = theta(:).';
            yMin = min(yMin, thetaRow);
            yMax = max(yMax, thetaRow);
        end

        fld = matlab.lang.makeValidName(char(modelNameGlobal));
        % Slight padding on lower bound; upper bound scaled a bit
        yLimsByModel.(fld) = [yMin - max(0.1, abs(yMin)); yMax * 1.5];  % 2 x P
    end

    % ===== Detect OT and Saline sets (by name substring) =====
    otIdx  = find(contains(sets, "OT",      'IgnoreCase', true), 1, 'first');
    salIdx = find(contains(sets, "Saline",  'IgnoreCase', true), 1, 'first');
    pairPossible = ~isempty(otIdx) && ~isempty(salIdx);

    % =====================================================================
    % If we have both OT and Saline, make paired figures (OT vs Saline).
    % Otherwise, fall back to per-set plotting.
    % =====================================================================

    if pairPossible
        setOT   = sets(otIdx);
        setSal  = sets(salIdx);

        for miGlobal = 1:numel(modelListGlobal)
            modelNameGlobal = modelListGlobal(miGlobal);

            subOT  = fitTable(fitTable.SetName == setOT  & ...
                              fitTable.Event   == eventUse & ...
                              fitTable.ModelName == modelNameGlobal, :);
            subSal = fitTable(fitTable.SetName == setSal & ...
                              fitTable.Event   == eventUse & ...
                              fitTable.ModelName == modelNameGlobal, :);

            if isempty(subOT) && isempty(subSal)
                continue;
            end

            thetaExample = [];
            if ~isempty(subOT)
                thetaExample = subOT.Theta{1};
            elseif ~isempty(subSal)
                thetaExample = subSal.Theta{1};
            end
            if isempty(thetaExample)
                continue;
            end
            P = numel(thetaExample);

            modelNameChar = char(modelNameGlobal);
            if isfield(paramNamesPerModel, modelNameChar)
                paramNames = paramNamesPerModel.(modelNameChar);
                if numel(paramNames) ~= P
                    paramNames = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
                end
            else
                paramNames = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
            end

            fld = matlab.lang.makeValidName(modelNameChar);
            useGlobalY = isfield(yLimsByModel, fld);
            if useGlobalY
                yMat = yLimsByModel.(fld);  % 2 x P
            else
                yMat = nan(2,P);
            end

            figName = sprintf('RL Params (OT vs Saline) - Model %s', modelNameChar);
            f = figure('Name', figName, 'Color','w');
            try
                t = tiledlayout(f, P, 2, 'TileSpacing','compact', 'Padding','compact');
            catch
                t = [];
            end

            for p = 1:P
                % ---------- Column 1: OT ----------
                subT = subOT;
                if ~isempty(subT)
                    condLevelsRaw = unique(string(subT.Condition), 'stable');

                    % Order conditions by keyword (works for "AI" and "OT AI", etc.)
                    ordered = strings(0,1);
                    tmpLevels = condLevelsRaw;
                    for kk = 1:numel(baseOrderKeywords)
                        key = baseOrderKeywords(kk);
                        mask = contains(tmpLevels, key, 'IgnoreCase', true);
                        if any(mask)
                            ordered = [ordered; tmpLevels(mask)]; %#ok<AGROW>
                            tmpLevels(mask) = "";
                        end
                    end
                    leftover     = tmpLevels(tmpLevels ~= "");
                    orderedConds = [ordered; leftover];

                    condCatsAll = categorical(string(subT.Condition), orderedConds, 'Ordinal', true);
                    condIdxAll  = double(condCatsAll);
                    valsAll     = cellfun(@(theta) theta(p), subT.Theta);
                    K           = numel(orderedConds);

                    if ~isempty(t)
                        ax1 = nexttile(t, (p-1)*2 + 1);
                    else
                        subplot(P,2,(p-1)*2+1);
                        ax1 = gca;
                    end
                    hold(ax1, 'on');

                    % Mean, STD, N per condition
                    meanVals = nan(1,K);
                    stdVals  = nan(1,K);
                    nVals    = zeros(1,K);
                    for k = 1:K
                        maskK = (condIdxAll == k);
                        thisVals = valsAll(maskK);
                        nVals(k) = numel(thisVals);
                        if nVals(k) > 0
                            meanVals(k) = mean(thisVals, 'omitnan');
                            stdVals(k)  = std(thisVals, 0, 'omitnan');
                        end
                    end

                    bar(ax1, 1:K, meanVals, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
                    errorbar(ax1, 1:K, meanVals, stdVals, ...
                             'k', 'LineStyle','none', 'LineWidth',1);

                    yline(ax1, 0, 'k--', 'LineWidth',1);

                    % X tick labels: Condition + (n=#) on second line
                    condStr = cellstr(orderedConds);
                    tickLabels = cell(K,1);
                    for k = 1:K
                        tickLabels{k} = sprintf('%s\\newline(n=%d)', condStr{k}, nVals(k));
                    end
                    set(ax1, 'XTick', 1:K, ...
                             'XTickLabel', tickLabels, ...
                             'TickLabelInterpreter', 'tex');

                    if p == P
                        xlabel(ax1, 'Condition');
                    else
                        xlabel(ax1, '');
                    end
                    ylabel(ax1, paramNames{p}, 'Interpreter','none');
                    title(ax1, sprintf('%s', string(setOT)), 'Interpreter','none');

                    % Apply global y-limits then significance bars
                    if useGlobalY
                        yLo = yMat(1,p); yHi = yMat(2,p);
                        if isfinite(yLo) && isfinite(yHi) && (yHi > yLo)
                            ylim(ax1, [yLo yHi]);
                        end
                    end

                    if K >= 2
                        [sigPairs, pvals] = local_computeSig(valsAll, condIdxAll, K);
                        local_plotSignificanceLines(ax1, sigPairs, pvals);
                    end

                    hold(ax1, 'off');
                else
                    if ~isempty(t)
                        ax1 = nexttile(t, (p-1)*2 + 1);
                    else
                        subplot(P,2,(p-1)*2+1);
                        ax1 = gca;
                    end
                    axis(ax1, 'off');
                    title(ax1, sprintf('%s (no data)', string(setOT)), 'Interpreter','none');
                end

                % ---------- Column 2: Saline ----------
                subT = subSal;
                if ~isempty(subT)
                    condLevelsRaw = unique(string(subT.Condition), 'stable');

                    ordered = strings(0,1);
                    tmpLevels = condLevelsRaw;
                    for kk = 1:numel(baseOrderKeywords)
                        key = baseOrderKeywords(kk);
                        mask = contains(tmpLevels, key, 'IgnoreCase', true);
                        if any(mask)
                            ordered = [ordered; tmpLevels(mask)]; %#ok<AGROW>
                            tmpLevels(mask) = "";
                        end
                    end
                    leftover     = tmpLevels(tmpLevels ~= "");
                    orderedConds = [ordered; leftover];

                    condCatsAll = categorical(string(subT.Condition), orderedConds, 'Ordinal', true);
                    condIdxAll  = double(condCatsAll);
                    valsAll     = cellfun(@(theta) theta(p), subT.Theta);
                    K           = numel(orderedConds);

                    if ~isempty(t)
                        ax2 = nexttile(t, (p-1)*2 + 2);
                    else
                        subplot(P,2,(p-1)*2+2);
                        ax2 = gca;
                    end
                    hold(ax2, 'on');

                    meanVals = nan(1,K);
                    stdVals  = nan(1,K);
                    nVals    = zeros(1,K);
                    for k = 1:K
                        maskK = (condIdxAll == k);
                        thisVals = valsAll(maskK);
                        nVals(k) = numel(thisVals);
                        if nVals(k) > 0
                            meanVals(k) = mean(thisVals, 'omitnan');
                            stdVals(k)  = std(thisVals, 0, 'omitnan');
                        end
                    end

                    bar(ax2, 1:K, meanVals, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
                    errorbar(ax2, 1:K, meanVals, stdVals, ...
                             'k', 'LineStyle','none', 'LineWidth',1);

                    yline(ax2, 0, 'k--', 'LineWidth',1);

                    condStr = cellstr(orderedConds);
                    tickLabels = cell(K,1);
                    for k = 1:K
                        tickLabels{k} = sprintf('%s\\newline(n=%d)', condStr{k}, nVals(k));
                    end
                    set(ax2, 'XTick', 1:K, ...
                             'XTickLabel', tickLabels, ...
                             'TickLabelInterpreter', 'tex');

                    if p == P
                        xlabel(ax2, 'Condition');
                    else
                        xlabel(ax2, '');
                    end
                    ylabel(ax2, '');
                    title(ax2, sprintf('%s', string(setSal)), 'Interpreter','none');

                    if useGlobalY
                        yLo = yMat(1,p); yHi = yMat(2,p);
                        if isfinite(yLo) && isfinite(yHi) && (yHi > yLo)
                            ylim(ax2, [yLo yHi]);
                        end
                    end

                    if K >= 2
                        [sigPairs, pvals] = local_computeSig(valsAll, condIdxAll, K);
                        local_plotSignificanceLines(ax2, sigPairs, pvals);
                    end

                    hold(ax2, 'off');
                else
                    if ~isempty(t)
                        ax2 = nexttile(t, (p-1)*2 + 2);
                    else
                        subplot(P,2,(p-1)*2+2);
                        ax2 = gca;
                    end
                    axis(ax2, 'off');
                    title(ax2, sprintf('%s (no data)', string(setSal)), 'Interpreter','none');
                end
            end

            if ~isempty(t)
                sgtitle(t, sprintf('OT vs Saline | Model: %s', modelNameChar), ...
                        'Interpreter','none');
            else
                suptitle(sprintf('OT vs Saline | Model: %s', modelNameChar));
            end
        end

        % return;  % do not fall back to per-set plotting when OT/Saline pair exists
    end

    % =====================================================================
    % Fallback: per-set plotting (when OT/Saline pair not present)
    % =====================================================================

    for si = 1:numel(sets)
        setName = sets(si);
        if ~strcmp(setName, "baseline")
            continue;
        end
        subSetSet = fitTable(fitTable.SetName == setName & fitTable.Event == eventUse, :);
        if isempty(subSetSet), continue; end

        models = unique(subSetSet.ModelName, 'stable');
        for mi = 1:numel(models)
            modelName = models(mi);
            subT = subSetSet(subSetSet.ModelName == modelName, :);
            if isempty(subT), continue; end

            thetaExample = subT.Theta{1};
            P = numel(thetaExample);

            modelNameChar = char(modelName);
            if isfield(paramNamesPerModel, modelNameChar)
                paramNames = paramNamesPerModel.(modelNameChar);
                if numel(paramNames) ~= P
                    paramNames = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
                end
            else
                paramNames = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
            end

            % condition ordering by keywords
            condLevelsRaw = unique(string(subT.Condition), 'stable');
            ordered = strings(0,1);
            tmpLevels = condLevelsRaw;
            for kk = 1:numel(baseOrderKeywords)
                key = baseOrderKeywords(kk);
                mask = contains(tmpLevels, key, 'IgnoreCase', true);
                if any(mask)
                    ordered = [ordered; tmpLevels(mask)]; %#ok<AGROW>
                    tmpLevels(mask) = "";
                end
            end
            leftover     = tmpLevels(tmpLevels ~= "");
            orderedConds = [ordered; leftover];

            condCatsAll = categorical(string(subT.Condition), orderedConds, 'Ordinal', true);
            condIdxAll  = double(condCatsAll);
            K           = numel(orderedConds);

            fld = matlab.lang.makeValidName(modelNameChar);
            useGlobalY = isfield(yLimsByModel, fld);
            if useGlobalY
                yMat = yLimsByModel.(fld);
            else
                yMat = nan(2,P);
            end

            figure('Name', sprintf('RL Parameters by Condition - %s - %s', setName, modelNameChar), ...
                   'Color','w');
            nCols = min(3,P);
            nRows = ceil(P/nCols);

            for p = 1:P
                ax = subplot(nRows,nCols,p); hold(ax, 'on');

                valsAll = cellfun(@(theta) theta(p), subT.Theta);

                % Mean, STD, N per condition
                meanVals = nan(1,K);
                stdVals  = nan(1,K);
                nVals    = zeros(1,K);
                for k = 1:K
                    maskK = (condIdxAll == k);
                    thisVals = valsAll(maskK);
                    nVals(k) = numel(thisVals);
                    if nVals(k) > 0
                        meanVals(k) = mean(thisVals, 'omitnan');
                        stdVals(k)  = std(thisVals, 0, 'omitnan');
                    end
                end

                bar(ax, 1:K, meanVals, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
                errorbar(ax, 1:K, meanVals, stdVals, ...
                         'k', 'LineStyle','none', 'LineWidth',1);

                yline(ax, 0, 'k--', 'LineWidth',1);

                condStr = cellstr(orderedConds);
                tickLabels = cell(K,1);
                for k = 1:K
                    tickLabels{k} = sprintf('%s\\newline(n=%d)', condStr{k}, nVals(k));
                end
                set(ax, 'XTick', 1:K, ...
                        'XTickLabel', tickLabels, ...
                        'TickLabelInterpreter', 'tex');

                title(ax, paramNames{p}, 'Interpreter','none');
                xlabel(ax, 'Condition');
                ylabel(ax, 'Param value');

                if useGlobalY
                    yLo = yMat(1,p); yHi = yMat(2,p);
                    if isfinite(yLo) && isfinite(yHi) && (yHi > yLo)
                        ylim(ax, [yLo yHi]);
                    end
                end

                if K >= 2
                    [sigPairs, pvals] = local_computeSig(valsAll, condIdxAll, K);
                    local_plotSignificanceLines(ax, sigPairs, pvals);
                end

                hold(ax, 'off');
            end

            sgtitle(sprintf('Set: %s | Model: %s', setName, modelNameChar), ...
                    'Interpreter','none');
        end
    end
end

function [sigPairs, pvals] = local_computeSig(vals, condIdx, K)
    sigPairs = [];
    pvals    = [];

    if K < 2
        return;
    end

    if K == 2
        x1 = vals(condIdx == 1);
        x2 = vals(condIdx == 2);
        if numel(x1) > 1 && numel(x2) > 1
            [~,p] = ttest2(x1, x2);
            if ~isnan(p) && p < 0.05
                sigPairs = [1 2];
                pvals    = p;
            end
        end
    else
        [~,~,stats] = anova1(vals, condIdx, 'off');
        c = multcompare(stats, 'Display','off');
        sigMask = c(:,6) < 0.05;
        sigPairs = c(sigMask,1:2);
        pvals    = c(sigMask,6);
    end
end

function local_plotSignificanceLines(ax, sigPairs, pvals)
    if isempty(sigPairs)
        return;
    end

    axes(ax);
    hold(ax,'on');

    yLim   = get(ax,'YLim');
    yRange = diff(yLim);
    yBase  = yLim(2);
    step   = 0.06 * yRange;
    tick   = 0.015 * yRange;

    nPairs = size(sigPairs,1);

    for i = 1:nPairs
        x1 = sigPairs(i,1);
        x2 = sigPairs(i,2);
        if x1 == x2, continue; end

        y = yBase + i * step;

        plot(ax, [x1 x2], [y y], 'k-', 'LineWidth', 1);
        plot(ax, [x1 x1], [y-tick y], 'k-', 'LineWidth', 1);
        plot(ax, [x2 x2], [y-tick y], 'k-', 'LineWidth', 1);

        starText = local_pval2star(pvals(i));
        text(mean([x1 x2]), y+tick*0.5, starText, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Color','r','FontSize',12,'FontWeight','bold');
    end

    set(ax,'YLim',[yLim(1), yBase + (nPairs+2)*step]);
end

function s = local_pval2star(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end

function local_plotTerminalQDistributions(terminalQTable, cond_sets, setNames)
    if isempty(terminalQTable)
        return;
    end

    useBoxchart = exist('boxchart','file') == 2;

    terminalQTable.SetName   = string(terminalQTable.SetName);
    terminalQTable.Event     = string(terminalQTable.Event);
    terminalQTable.MonkeyRaw = string(terminalQTable.MonkeyRaw);
    terminalQTable.Condition = string(terminalQTable.Condition);
    terminalQTable.Action    = string(terminalQTable.Action);

    % ---------- GLOBAL Y-LIMITS ACROSS ALL FIGURES ----------
    qVals = double(terminalQTable.Q_self);
    yMinData = min(qVals, [], 'omitnan');
    yMaxData = max(qVals, [], 'omitnan');

    if ~isfinite(yMinData) || ~isfinite(yMaxData)
        return;
    end

    if yMinData == yMaxData
        yMinData = yMinData - 0.5;
        yMaxData = yMaxData + 0.5;
    end

    yRange   = yMaxData - yMinData;
    yMinPlot = yMinData - 0.05 * yRange;
    yMaxPlot = yMaxData + 0.20 * yRange;  % extra headroom for sig bars

    for s = 1:numel(setNames)
        setName = setNames{s};
        subSet = terminalQTable(terminalQTable.SetName == string(setName), :);
        if isempty(subSet)
            continue;
        end

        monkeys = unique(subSet.MonkeyRaw, 'stable');
        for mi = 1:numel(monkeys)
            mID = monkeys(mi);
            subM = subSet(subSet.MonkeyRaw == mID, :);
            if isempty(subM)
                continue;
            end

            % Ordered actions
            uActions = unique(subM.Action, 'stable');
            uActions = intersect(["BUY","HOLD","SELL","Non-BUY"], uActions, 'stable');
            if isempty(uActions)
                uActions = unique(subM.Action, 'stable');
            end
            nActions = numel(uActions);

            % Ordered conditions using cond_sets if possible
            allConds = unique(subM.Condition, 'stable');
            if s <= numel(cond_sets) && ~isempty(cond_sets{s})
                desired = string(cond_sets{s});
                allCondsNorm = regexprep(allConds, '\s+', '');
                desiredNorm  = regexprep(desired, '\s+', '');
                present = desired(ismember(desiredNorm, allCondsNorm));
                leftover = setdiff(allConds, present, 'stable');
                uConds = [present; leftover];
            else
                uConds = allConds;
            end
            nConds = numel(uConds);

            actionCats = categorical(subM.Action, uActions);
            condCats   = categorical(subM.Condition, uConds);
            actIdx  = double(actionCats);   % 1..nActions
            condIdx = double(condCats);     % 1..nConds

            % Layout: condition groups × actions within each condition
            x = (condIdx - 1) * (nActions + 1) + actIdx;

            figName = sprintf('Terminal Q (Trial 15) | %s | Monkey %s', setName, mID);
            figure('Name', figName, 'Color','w'); hold on;

            if useBoxchart
                for ic = 1:nConds
                    for ia = 1:nActions
                        mask = condIdx == ic & actIdx == ia;
                        if ~any(mask)
                            continue;
                        end
                        thisX = (ic - 1) * (nActions + 1) + ia;
                        boxchart(repmat(thisX, sum(mask), 1), subM.Q_self(mask));
                    end
                end
            else
                for ic = 1:nConds
                    for ia = 1:nActions
                        mask = condIdx == ic & actIdx == ia;
                        if ~any(mask)
                            continue;
                        end
                        thisX = (ic - 1) * (nActions + 1) + ia;
                        boxplot(subM.Q_self(mask), repmat(thisX, sum(mask), 1), ...
                                'Positions', thisX, 'Widths', 0.7, 'Colors','k');
                    end
                end
            end

            % X-ticks and labels: Condition-Action, e.g. AI-BUY, AI-HOLD, ...
            xticks   = [];
            xtlabels = strings(0,1);
            for ic = 1:nConds
                for ia = 1:nActions
                    xticks(end+1)   = (ic - 1) * (nActions + 1) + ia; %#ok<AGROW>
                    xtlabels(end+1) = sprintf('%s-%s', uConds(ic), uActions(ia)); %#ok<AGROW>
                end
            end
            set(gca, 'XTick', xticks, 'XTickLabel', xtlabels, 'XTickLabelRotation', 45);
            ylabel('Terminal Q_{self}');
            title(sprintf('Set: %s | Monkey: %s (Trial 15)', setName, mID), 'Interpreter','none');

            % ---------- SIGNIFICANCE: BETWEEN ACTIONS WITHIN EACH CONDITION ----------
            for ic = 1:nConds
                maskCond = condIdx == ic;
                if sum(maskCond) < 2
                    continue;
                end

                valsCond = subM.Q_self(maskCond);

                actionCatsCond   = categorical(subM.Action(maskCond), uActions);
                actionLevelsCond = categories(actionCatsCond);
                Kact             = numel(actionLevelsCond);
                actionIdxCond    = double(actionCatsCond); % 1..Kact

                % x-positions for these actions in this condition (match layout)
                xPosAct = nan(Kact,1);
                for k = 1:Kact
                    actName   = actionLevelsCond{k};
                    globalIdx = find(uActions == actName, 1, 'first');
                    if isempty(globalIdx)
                        continue;
                    end
                    xPosAct(k) = (ic - 1) * (nActions + 1) + globalIdx;
                end

                [sigPairs, pvals] = local_computeSig(valsCond, actionIdxCond, Kact);
                local_plotSignificanceLines_custom(gca, sigPairs, pvals, xPosAct, yMinPlot, yMaxPlot);
            end

            % ---------- APPLY GLOBAL Y-LIMITS ----------
            ylim([yMinPlot, yMaxPlot]);

            hold off;
        end
    end
end

function local_plotSignificanceLines_custom(ax, sigPairs, pvals, xPos, yMinPlot, yMaxPlot)
    if isempty(sigPairs) || isempty(xPos)
        return;
    end

    axes(ax);
    hold(ax,'on');

    yRange = yMaxPlot - yMinPlot;
    if yRange <= 0
        yRange = 1;
    end

    % Reserve top ~15% of axis for significance bars
    yBase = yMaxPlot - 0.15 * yRange;
    step  = 0.06 * yRange;
    tick  = 0.015 * yRange;

    nPairs = size(sigPairs,1);

    for i = 1:nPairs
        i1 = sigPairs(i,1);
        i2 = sigPairs(i,2);
        if i1 == i2 || i1 < 1 || i2 < 1 || i1 > numel(xPos) || i2 > numel(xPos)
            continue;
        end

        x1 = xPos(i1);
        x2 = xPos(i2);
        if isnan(x1) || isnan(x2)
            continue;
        end

        y = yBase + (i-1) * step;
        if y + tick > yMaxPlot
            y = yMaxPlot - 2*tick;
        end

        plot(ax, [x1 x2], [y y], 'k-', 'LineWidth', 1);
        plot(ax, [x1 x1], [y-tick y], 'k-', 'LineWidth', 1);
        plot(ax, [x2 x2], [y-tick y], 'k-', 'LineWidth', 1);

        starText = local_pval2star(pvals(i));
        text(mean([x1 x2]), y + tick, starText, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');
    end
end

function terminalQTable = local_buildTerminalQTable(allResults, setNames, canonicalEvt, BUY_NON_BUY_ONLY)
    terminalQTable = table( ...
        string.empty(0,1), ... % SetName
        string.empty(0,1), ... % Event
        string.empty(0,1), ... % MonkeyRaw
        string.empty(0,1), ... % Monkey
        string.empty(0,1), ... % Condition
        double.empty(0,1), ... % Session
        double.empty(0,1), ... % BlockID
        double.empty(0,1), ... % TrialInBlock
        string.empty(0,1), ... % Action
        double.empty(0,1), ... % Q_self
        'VariableNames', {'SetName','Event','MonkeyRaw','Monkey','Condition', ...
                          'Session','BlockID','TrialInBlock','Action','Q_self'});

    for s = 1:numel(setNames)
        setName = setNames{s};

        if ~isfield(allResults, setName) || ...
           ~isfield(allResults.(setName), 'behaviorTrialTableWithRL')
            continue;
        end

        trialTable = allResults.(setName).behaviorTrialTableWithRL;
        if isempty(trialTable) || ...
           ~ismember('BlockID', trialTable.Properties.VariableNames) || ...
           ~ismember('TrialInBlock', trialTable.Properties.VariableNames)
            continue;
        end

        idxTerm = find(trialTable.TrialInBlock == 15);
        if isempty(idxTerm)
            continue;
        end

        for ii = 1:numel(idxTerm)
            r   = idxTerm(ii);
            row = trialTable(r,:);

            if BUY_NON_BUY_ONLY
                actions = {'BUY','Non-BUY'};
                qVals   = [row.RL_Q_self_buy, row.RL_Q_self_hold];
            else
                actions = {'BUY','HOLD','SELL'};
                qVals   = [row.RL_Q_self_buy, row.RL_Q_self_hold, row.RL_Q_self_sell];
            end

            for a = 1:numel(actions)
                if isnan(qVals(a))
                    continue;
                end

                newRow = table( ...
                    string(setName), ...
                    string(canonicalEvt), ...
                    string(row.MonkeyRaw), ...
                    string(row.Monkey), ...
                    string(row.Condition), ...
                    double(row.Session), ...
                    double(row.BlockID), ...
                    double(row.TrialInBlock), ...
                    string(actions{a}), ...
                    double(qVals(a)), ...
                    'VariableNames', terminalQTable.Properties.VariableNames);

                terminalQTable = [terminalQTable; newRow]; %#ok<AGROW>
            end
        end
    end
end

function allResults = local_addDeltaQ(allResults, setNames, TARGET_EVTS, BUY_NON_BUY_ONLY)
% local_addDeltaQ
% -------------------------------------------------------------------------
% For each setName:
%   - Read allResults.(setName).behaviorTrialTableWithRL
%   - Using RL_Q_self_* fields and subjectChoice, compute:
%       * RL_Q_self_chosen   : Q-value of chosen option at decision time
%       * RL_Q_self_otherMax : maximum Q among unchosen options
%       * RL_DeltaQ          : chosenQ - otherMaxQ
%   - Attach these three columns to:
%       * behaviorTrialTableWithRL
%       * each event-specific trialTableWithRL (checking row alignment)
%
% This is done without touching any RL fitting code.

    if nargin < 4
        error('local_addDeltaQ requires allResults, setNames, TARGET_EVTS, BUY_NON_BUY_ONLY.');
    end

    for s = 1:numel(setNames)
        setName = setNames{s};
        if ~isfield(allResults, setName)
            continue;
        end

        if ~isfield(allResults.(setName), 'behaviorTrialTableWithRL')
            % No RL latents stored for this set; skip.
            continue;
        end

        T = allResults.(setName).behaviorTrialTableWithRL;
        nTrials = height(T);
        if nTrials == 0
            continue;
        end

        % Preallocate
        Q_chosen   = nan(nTrials,1);
        Q_otherMax = nan(nTrials,1);
        DeltaQ     = nan(nTrials,1);

        % For convenience
        choiceStr = string(T.subjectChoice);

        if BUY_NON_BUY_ONLY
            % Binary case: BUY vs Non-BUY.
            % Convention in existing code:
            %   RL_Q_self_buy  : Q(BUY)
            %   RL_Q_self_hold : Q(Non-BUY) when BUY_NON_BUY_ONLY == 1
            Q_buy    = T.RL_Q_self_buy;
            Q_nonBuy = T.RL_Q_self_hold;

            for t = 1:nTrials
                ch = choiceStr(t);
                if ch == "BUY"
                    qc = Q_buy(t);
                    qo = Q_nonBuy(t);
                elseif ch == "Non-BUY"
                    qc = Q_nonBuy(t);
                    qo = Q_buy(t);
                else
                    qc = NaN;
                    qo = NaN;
                end
                Q_chosen(t)   = qc;
                Q_otherMax(t) = qo;
                DeltaQ(t)     = qc - qo;
            end
        else
            % Three-action case: BUY, HOLD, SELL
            Q_buy  = T.RL_Q_self_buy;
            Q_hold = T.RL_Q_self_hold;
            Q_sell = T.RL_Q_self_sell;

            for t = 1:nTrials
                ch = choiceStr(t);

                if ch == "BUY"
                    qc = Q_buy(t);
                    qo = max([Q_hold(t), Q_sell(t)]);
                elseif ch == "HOLD"
                    qc = Q_hold(t);
                    qo = max([Q_buy(t), Q_sell(t)]);
                elseif ch == "SELL"
                    qc = Q_sell(t);
                    qo = max([Q_buy(t), Q_hold(t)]);
                else
                    qc = NaN;
                    qo = NaN;
                end

                Q_chosen(t)   = qc;
                Q_otherMax(t) = qo;
                DeltaQ(t)     = qc - qo;
            end
        end

        % Attach to canonical behavior table
        T.RL_Q_self_chosen   = Q_chosen;
        T.RL_Q_self_otherMax = Q_otherMax;
        T.RL_DeltaQ          = DeltaQ;

        allResults.(setName).behaviorTrialTableWithRL = T;

        % Propagate to all event-specific trialTableWithRL tables for this set
        for pair = 1:numel(TARGET_EVTS)
            evt = TARGET_EVTS{pair};
            if ~isfield(allResults.(setName), evt)
                continue;
            end
            if ~isfield(allResults.(setName).(evt), 'trialTableWithRL')
                continue;
            end

            Tevt = allResults.(setName).(evt).trialTableWithRL;
            if height(Tevt) ~= nTrials
                warning(['DeltaQ propagation skipped for set "%s", event "%s": ', ...
                         'row mismatch between behaviorTrialTableWithRL (%d) and trialTableWithRL (%d).'], ...
                         setName, string(evt), nTrials, height(Tevt));
                continue;
            end

            Tevt.RL_Q_self_chosen   = Q_chosen;
            Tevt.RL_Q_self_otherMax = Q_otherMax;
            Tevt.RL_DeltaQ          = DeltaQ;

            allResults.(setName).(evt).trialTableWithRL = Tevt;
        end
    end
end

function [kappaComparisonTable, kappaSummaryByCond] = local_buildKappaComparison(behaviorFitSummary, behaviorParamTable)
% local_buildKappaComparison
% -------------------------------------------------------------------------
% Construct a session-level table that quantifies the benefit of including
% kappa (BUY-bias) by comparing:
%   SelfOnly  (no kappa)
%   vs
%   SelfOnlyBuyBias        (with kappa)
%
% For each (SetName, Monkey, Condition, Session) with both models present:
%   - LL0   = -NegLogLik(SelfOnly)
%   - LL1   = -NegLogLik(SelfOnlyBuyBias)
%   - dLL   = LL1 - LL0
%   - LRT   = 2 * dLL
%   - df    = NumParams(SelfOnlyBuyBias) - NumParams(SelfOnly)
%   - p_LRT = 1 - chi2cdf(LRT, df)
%   - dAIC  = AIC_noBias   - AIC_withKappa
%   - dBIC  = BIC_noBias   - BIC_withKappa
% and we also attach kappaHat (from behaviorParamTable).
%
% A positive dBIC means BIC is lower for SelfOnlyBuyBias (with kappa), i.e.
%   favoring inclusion of kappa.
%
% We also build a summary table aggregated by (SetName, Monkey, Condition):
%   - mean dBIC, mean dAIC, mean LRT, mean p_LRT
%   - fraction of sessions with dBIC > 0 (BIC favors kappa)
%   - N_sessions, mean kappaHat, etc.

    if isempty(behaviorFitSummary)
        kappaComparisonTable = table();
        kappaSummaryByCond   = table();
        return;
    end

    % Focus on canonical event label (just for bookkeeping)
    events = unique(behaviorFitSummary.Event, 'stable');
    if ~isempty(events)
        canonicalEvt = events(1);
    else
        canonicalEvt = "";
    end

    % Filter to relevant models
    mask0 = behaviorFitSummary.ModelName == "SelfOnly";
    mask1 = behaviorFitSummary.ModelName == "SelfOnlyBuyBias";

    fits0 = behaviorFitSummary(mask0,:);
    fits1 = behaviorFitSummary(mask1,:);

    if isempty(fits0) || isempty(fits1)
        warning('No fits found for SelfOnly or SelfOnlyBuyBias. kappaComparisonTable will be empty.');
        kappaComparisonTable = table();
        kappaSummaryByCond   = table();
        return;
    end

    rows = [];

    for i = 1:height(fits0)
        row0 = fits0(i,:);

        % Match SelfOnlyBuyBias fit for the same SetName/Monkey/Condition/Session
        maskMatch = (fits1.SetName   == row0.SetName)   & ...
                    (fits1.Monkey    == row0.Monkey)    & ...
                    (fits1.Condition == row0.Condition) & ...
                    (fits1.DateID   == row0.DateID);

        idx = find(maskMatch, 1, 'first');
        if isempty(idx)
            continue;
        end
        row1 = fits1(idx,:);

        % Make sure number of trials matches
        if row0.NumTrials ~= row1.NumTrials
            warning('Mismatched NumTrials for SelfOnly vs SelfOnlyBuyBias in Set=%s Monkey=%s Cond=%s Sess=%d. Skipping.', ...
                    char(row0.SetName), char(row0.Monkey), char(row0.Condition), row0.DateID);
            continue;
        end

        % Compute LL, LRT, AIC / BIC differences
        LL0 = -row0.NegLogLik;
        LL1 = -row1.NegLogLik;

        dLL  = LL1 - LL0;
        LRT  = 2*dLL;
        df   = row1.NumParams - row0.NumParams;

        if df <= 0
            p_LRT = NaN;
        else
            % chi2cdf requires Statistics Toolbox; if missing, this call may error.
            % We catch and set NaN if that happens.
            try
                p_LRT = 1 - chi2cdf(LRT, df);
            catch
                warning('chi2cdf not available. Setting p_LRT to NaN.');
                p_LRT = NaN;
            end
        end

        dAIC = row0.AIC - row1.AIC;
        dBIC = row0.BIC - row1.BIC;

        % BIC-based Bayes factor (approximation):
        %   DeltaBIC = BIC(noBias) - BIC(withKappa)
        %   log BF (kappa vs noBias) ~= 0.5 * DeltaBIC
        logBF_BIC = 0.5 * dBIC;
        BF_BIC    = exp(logBF_BIC);
        
        % Evidence strength based on |DeltaBIC|
        absDBIC = abs(dBIC);
        if absDBIC < 2
            evidenceStrength = "none";        % negligible
        elseif absDBIC < 6
            evidenceStrength = "positive";    % positive/substantial
        elseif absDBIC < 10
            evidenceStrength = "strong";      % strong
        else
            evidenceStrength = "very_strong"; % very strong
        end
        
        % Direction: which model is favored?
        if dBIC > 0
            evidenceDirection = "favor_kappa";     % BIC lower for SelfOnlyBuyBias (with kappa)
        elseif dBIC < 0
            evidenceDirection = "favor_noKappa";   % BIC lower for SelfOnly
        else
            evidenceDirection = "tie";
        end
        
        evidenceLabel = evidenceStrength + "_" + evidenceDirection;

        % Positive dAIC/dBIC => SelfOnlyBuyBias (with kappa) is better (lower IC).
        % Extract kappaHat from behaviorParamTable (may be missing).
        kappaHat = NaN;
        maskParam = (behaviorParamTable.SetName   == row0.SetName)   & ...
                    (behaviorParamTable.Monkey    == row0.Monkey)    & ...
                    (behaviorParamTable.Condition == row0.Condition) & ...
                    (behaviorParamTable.DateID   == row0.DateID)   & ...
                    (behaviorParamTable.ModelName == "SelfOnlyBuyBias")     & ...
                    (behaviorParamTable.ParamName == "kappaBuy");
        idxParam = find(maskParam, 1, 'first');
        if ~isempty(idxParam)
            kappaHat = behaviorParamTable.Value(idxParam);
        end

        newRow = table( ...
            row0.SetName, ...
            canonicalEvt, ...
            row0.Monkey, ...
            row0.Condition, ...
            row0.DateID, ...
            row0.NumTrials, ...
            LL0, LL1, dLL, LRT, df, p_LRT, ...
            row0.AIC, row1.AIC, dAIC, ...
            row0.BIC, row1.BIC, dBIC, ...
            logBF_BIC, BF_BIC, ...
            evidenceStrength, evidenceDirection, evidenceLabel, ...
            kappaHat, ...
            'VariableNames', {'SetName','Event','Monkey','Condition','DateID', ...
                              'NumTrials','LL_noBias','LL_withKappa','DeltaLL','LRT','df','p_LRT', ...
                              'AIC_noBias','AIC_withKappa','DeltaAIC', ...
                              'BIC_noBias','BIC_withKappa','DeltaBIC', ...
                              'logBF_BIC','BF_BIC', ...
                              'EvidenceStrength','EvidenceDirection','EvidenceLabel', ...
                              'kappaHat'});

        rows = [rows; newRow]; %#ok<AGROW>
    end

    kappaComparisonTable = rows;

    % ------------------------- Summaries by condition ---------------------
    if isempty(kappaComparisonTable)
        kappaSummaryByCond = table();
        return;
    end

    T = kappaComparisonTable;

    [G, key] = findgroups(T(:,{'SetName','Monkey','Condition'}));

    N_sessions     = splitapply(@numel,        T.DateID,   G);
    mean_dBIC      = splitapply(@(x) mean(x,'omitnan'), T.DeltaBIC, G);
    mean_dAIC      = splitapply(@(x) mean(x,'omitnan'), T.DeltaAIC, G);
    mean_LRT       = splitapply(@(x) mean(x,'omitnan'), T.LRT,      G);
    mean_p_LRT     = splitapply(@(x) mean(x,'omitnan'), T.p_LRT,    G);
    mean_kappaHat  = splitapply(@(x) mean(x,'omitnan'), T.kappaHat, G);
    frac_BIC_better = splitapply(@(x) mean(x > 0), T.DeltaBIC, G);  % fraction where BIC favors kappa
    frac_AIC_better = splitapply(@(x) mean(x > 0), T.DeltaAIC, G);  % fraction where AIC favors kappa

    % Bayes-factor based summaries
    mean_logBF_BIC = splitapply(@(x) mean(x,'omitnan'), T.logBF_BIC, G);
    frac_BF_gt1    = splitapply(@(x) mean(x > 0),       T.logBF_BIC, G);   % BF>1 for kappa
    frac_strongBIC = splitapply(@(d) mean(d > 6),       T.DeltaBIC,  G);   % strong+ evidence
    frac_veryStrongBIC = splitapply(@(d) mean(d > 10),  T.DeltaBIC,  G);

    kappaSummaryByCond = key;
    kappaSummaryByCond.N_sessions          = N_sessions;
    kappaSummaryByCond.meanDeltaBIC        = mean_dBIC;
    kappaSummaryByCond.meanDeltaAIC        = mean_dAIC;
    kappaSummaryByCond.meanLRT             = mean_LRT;
    kappaSummaryByCond.mean_p_LRT          = mean_p_LRT;
    kappaSummaryByCond.meanKappaHat        = mean_kappaHat;
    kappaSummaryByCond.fracBIC_favorsKappa = frac_BIC_better;
    kappaSummaryByCond.fracAIC_favorsKappa = frac_AIC_better;
    
    % New BF-based summary columns
    kappaSummaryByCond.meanLogBF_BIC           = mean_logBF_BIC;
    kappaSummaryByCond.fracBF_gt1_forKappa     = frac_BF_gt1;
    kappaSummaryByCond.fracStrongEvidenceBIC   = frac_strongBIC;
    kappaSummaryByCond.fracVeryStrongEvidenceBIC = frac_veryStrongBIC;
end

function Qtc = local_buildQTimecourse(allResults, setNames, BUY_NON_BUY_ONLY)
% local_buildQTimecourse
% -------------------------------------------------------------------------
% Build a long-form table summarizing how Q-values evolve over TrialInBlock.
% For each setName:
%   - We use allResults.(setName).behaviorTrialTableWithRL as the canonical
%     per-trial RL latent table (already aligned to canonicalEvt).
%   - For each action (BUY/HOLD/SELL or BUY/Non-BUY, depending on
%     BUY_NON_BUY_ONLY), we create rows:
%       [SetName, MonkeyRaw, Condition, Action, TrialInBlock, Q, BlockID]
%   - We then group by
%       (SetName, MonkeyRaw, Condition, Action, TrialInBlock)
%     and compute:
%       meanQ, semQ, NumBlocks (unique BlockID count)
%
% The output Qtc is a table with columns:
%   SetName, MonkeyRaw, Condition, Action, TrialInBlock, meanQ, semQ, NumBlocks

    rows = table( ...
        string.empty(0,1), ... % SetName
        string.empty(0,1), ... % MonkeyRaw
        string.empty(0,1), ... % Condition
        string.empty(0,1), ... % Action
        double.empty(0,1), ... % TrialInBlock
        double.empty(0,1), ... % Q
        double.empty(0,1), ... % BlockID
        'VariableNames', {'SetName','MonkeyRaw','Condition','Action','TrialInBlock','Q','BlockID'});

    for s = 1:numel(setNames)
        setName = setNames{s};
        if ~isfield(allResults, setName)
            continue;
        end

        if ~isfield(allResults.(setName), 'behaviorTrialTableWithRL')
            continue;
        end

        T = allResults.(setName).behaviorTrialTableWithRL;
        if isempty(T)
            continue;
        end

        if ~ismember('BlockID', T.Properties.VariableNames)
            warning('BlockID missing in behaviorTrialTableWithRL for set "%s". Skipping.', setName);
            continue;
        end

        if ~ismember('MonkeyRaw', T.Properties.VariableNames)
            % Fall back to 'Monkey' if MonkeyRaw is absent
            monkeyRawStr = string(T.Monkey);
        else
            monkeyRawStr = string(T.MonkeyRaw);
        end

        nTrials = height(T);

        if BUY_NON_BUY_ONLY
            actions = {'BUY','Non-BUY'};
            Qmat    = [T.RL_Q_self_buy, T.RL_Q_self_hold];
        else
            actions = {'BUY','HOLD','SELL'};
            Qmat    = [T.RL_Q_self_buy, T.RL_Q_self_hold, T.RL_Q_self_sell];
        end

        for a = 1:numel(actions)
            actionName = actions{a};
            qVals      = Qmat(:,a);

            newRows = table( ...
                repmat(string(setName), nTrials, 1), ...
                monkeyRawStr, ...
                string(T.Condition), ...
                repmat(string(actionName), nTrials, 1), ...
                double(T.TrialInBlock), ...
                double(qVals), ...
                double(T.BlockID), ...
                'VariableNames', rows.Properties.VariableNames);

            rows = [rows; newRows]; %#ok<AGROW>
        end
    end

    if height(rows) == 0
        Qtc = rows([],:);
        return;
    end

    [G, key] = findgroups(rows(:,{'SetName','MonkeyRaw','Condition','Action','TrialInBlock'}));

    meanQ    = splitapply(@(x) mean(x, 'omitnan'), rows.Q, G);
    semQ     = splitapply(@(x) std(x, 'omitnan')./sqrt(sum(~isnan(x))), rows.Q, G);
    nBlocks  = splitapply(@(b) numel(unique(b(~isnan(b)))), rows.BlockID, G);

    Qtc = key;
    Qtc.meanQ    = meanQ;
    Qtc.semQ     = semQ;
    Qtc.NumBlocks = nBlocks;
end

function DQtc = local_buildDeltaQTimecourse(allResults, setNames)
% local_buildDeltaQTimecourse
% -------------------------------------------------------------------------
% Build a long-form table summarizing the time-course of RL_DeltaQ
% (decision advantage) as a function of TrialInBlock.
%
% For each setName, we read:
%   allResults.(setName).behaviorTrialTableWithRL
% which should contain:
%   RL_DeltaQ, subjectChoice, MonkeyRaw, Condition, TrialInBlock, BlockID
%
% We construct rows:
%   [SetName, MonkeyRaw, Condition, ChosenAction, TrialInBlock, DeltaQ, BlockID]
%
% Then group by:
%   (SetName, MonkeyRaw, Condition, ChosenAction, TrialInBlock)
% and compute:
%   meanDeltaQ, semDeltaQ, NumBlocks (number of unique BlockID per point)
%
% Output table DQtc columns:
%   SetName, MonkeyRaw, Condition, ChosenAction, TrialInBlock,
%   meanDeltaQ, semDeltaQ, NumBlocks

    rows = table( ...
        string.empty(0,1), ... % SetName
        string.empty(0,1), ... % MonkeyRaw
        string.empty(0,1), ... % Condition
        string.empty(0,1), ... % ChosenAction
        double.empty(0,1), ... % TrialInBlock
        double.empty(0,1), ... % DeltaQ
        double.empty(0,1), ... % BlockID
        'VariableNames', {'SetName','MonkeyRaw','Condition','ChosenAction','TrialInBlock','DeltaQ','BlockID'});

    for s = 1:numel(setNames)
        setName = setNames{s};
        if ~isfield(allResults, setName)
            continue;
        end

        if ~isfield(allResults.(setName), 'behaviorTrialTableWithRL')
            continue;
        end

        T = allResults.(setName).behaviorTrialTableWithRL;
        if isempty(T)
            continue;
        end

        if ~ismember('RL_DeltaQ', T.Properties.VariableNames)
            warning('RL_DeltaQ missing in behaviorTrialTableWithRL for set "%s". Skipping.', setName);
            continue;
        end

        if ~ismember('BlockID', T.Properties.VariableNames)
            warning('BlockID missing in behaviorTrialTableWithRL for set "%s". Skipping.', setName);
            continue;
        end

        if ~ismember('MonkeyRaw', T.Properties.VariableNames)
            monkeyRawStr = string(T.Monkey);
        else
            monkeyRawStr = string(T.MonkeyRaw);
        end

        nTrials = height(T);

        newRows = table( ...
            repmat(string(setName), nTrials, 1), ...
            monkeyRawStr, ...
            string(T.Condition), ...
            string(T.subjectChoice), ...
            double(T.TrialInBlock), ...
            double(T.RL_DeltaQ), ...
            double(T.BlockID), ...
            'VariableNames', rows.Properties.VariableNames);

        rows = [rows; newRows]; %#ok<AGROW>
    end

    if height(rows) == 0
        DQtc = rows([],:);
        return;
    end

    [G, key] = findgroups(rows(:,{'SetName','MonkeyRaw','Condition','ChosenAction','TrialInBlock'}));

    meanDQ = splitapply(@(x) mean(x, 'omitnan'), rows.DeltaQ, G);
    semDQ  = splitapply(@(x) std(x, 'omitnan')./sqrt(sum(~isnan(x))), rows.DeltaQ, G);
    nBlocks = splitapply(@(b) numel(unique(b(~isnan(b)))), rows.BlockID, G);

    DQtc = key;
    DQtc.meanDeltaQ = meanDQ;
    DQtc.semDeltaQ  = semDQ;
    DQtc.NumBlocks  = nBlocks;
end

function Tadv = local_buildTerminalQAdvantage(allResults, setNames, BUY_NON_BUY_ONLY)
% local_buildTerminalQAdvantage
% -------------------------------------------------------------------------
% Construct a trial-level table capturing the terminal Buy advantage:
%
%   BuyAdvantage = Q_buy - maxQ_other
%
% for blocks at a specified TrialInBlock (default 15).
%
% We use allResults.(setName).behaviorTrialTableWithRL, which contains
% RL_Q_self_buy/hold/sell (for the chosen RL model used for neural latents).
%
% For each setName:
%   - find rows where TrialInBlock == terminalTrial (default 15)
%   - compute Q_buy, Q_otherMax, BuyAdvantage
%   - store:
%       SetName, MonkeyRaw, Monkey, Condition, Session, BlockID,
%       TrialInBlock, Q_buy, Q_otherMax, BuyAdvantage

    terminalTrial = 15;

    Tadv = table( ...
        string.empty(0,1), ... % SetName
        string.empty(0,1), ... % MonkeyRaw
        string.empty(0,1), ... % Monkey
        string.empty(0,1), ... % Condition
        double.empty(0,1), ... % Session
        double.empty(0,1), ... % BlockID
        double.empty(0,1), ... % TrialInBlock
        double.empty(0,1), ... % Q_buy
        double.empty(0,1), ... % Q_otherMax
        double.empty(0,1), ... % BuyAdvantage
        'VariableNames', {'SetName','MonkeyRaw','Monkey','Condition', ...
                          'Session','BlockID','TrialInBlock', ...
                          'Q_buy','Q_otherMax','BuyAdvantage'});

    for s = 1:numel(setNames)
        setName = setNames{s};
        if ~isfield(allResults, setName)
            continue;
        end

        if ~isfield(allResults.(setName), 'behaviorTrialTableWithRL')
            continue;
        end

        T = allResults.(setName).behaviorTrialTableWithRL;
        if isempty(T)
            continue;
        end

        idxTerm = find(T.TrialInBlock == terminalTrial);
        if isempty(idxTerm)
            continue;
        end

        if ~ismember('MonkeyRaw', T.Properties.VariableNames)
            monkeyRawStr = string(T.Monkey);
        else
            monkeyRawStr = string(T.MonkeyRaw);
        end

        for j = 1:numel(idxTerm)
            r = idxTerm(j);
            Qb  = T.RL_Q_self_buy(r);

            if BUY_NON_BUY_ONLY
                % Binary case: other is Non-BUY (stored in hold column)
                Q_other = T.RL_Q_self_hold(r);
            else
                Qh  = T.RL_Q_self_hold(r);
                Qs  = T.RL_Q_self_sell(r);
                Q_other = max([Qh, Qs]);
            end

            adv = Qb - Q_other;

            newRow = table( ...
                string(setName), ...
                monkeyRawStr(r), ...
                string(T.Monkey(r)), ...
                string(T.Condition(r)), ...
                double(T.Session(r)), ...
                double(T.BlockID(r)), ...
                double(T.TrialInBlock(r)), ...
                double(Qb), ...
                double(Q_other), ...
                double(adv), ...
                'VariableNames', Tadv.Properties.VariableNames);

            Tadv = [Tadv; newRow]; %#ok<AGROW>
        end
    end
end

function local_plotKappaBIC_EvidenceByCondition(kappaComparisonTable)
% local_plotKappaBIC_EvidenceByCondition
% -------------------------------------------------------------------------
% Visualize BIC-based kappa evidence across conditions, per (SetName, Monkey).
%
% Uses:
%   DeltaBIC  = BIC(noBias) - BIC(withKappa)
%   logBF_BIC = 0.5 * DeltaBIC
%   BF_BIC    = exp(logBF_BIC)
%
% For each (SetName, Monkey) combination, creates a figure with:
%   - Top panel:  mean DeltaBIC by Condition (bar + SEM)
%   - Bottom:     mean log10(BF_BIC) by Condition (bar + SEM)
%
% Positive DeltaBIC / log10(BF_BIC) => evidence favoring kappa model.

    if isempty(kappaComparisonTable)
        fprintf('  [kappaBFPlot] kappaComparisonTable is empty; no plots generated.\n');
        return;
    end

    % Required columns
    requiredVars = {'DeltaBIC','logBF_BIC','BF_BIC','Condition','SetName','Monkey'};
    for i = 1:numel(requiredVars)
        if ~ismember(requiredVars{i}, kappaComparisonTable.Properties.VariableNames)
            error('kappaComparisonTable must contain column "%s".', requiredVars{i});
        end
    end

    sets = unique(kappaComparisonTable.SetName, 'stable');
    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    for si = 1:numel(sets)
        setName = sets(si);
        subSet  = kappaComparisonTable(kappaComparisonTable.SetName == setName, :);
        if isempty(subSet)
            continue;
        end

        monkeys = unique(subSet.Monkey, 'stable');
        for mi = 1:numel(monkeys)
            mID  = monkeys(mi);
            subM = subSet(subSet.Monkey == mID, :);
            if isempty(subM)
                continue;
            end

            condLevelsRaw = unique(string(subM.Condition), 'stable');

            % Order conditions by keywords (AI/Replay/Decoy/Live first)
            ordered = strings(0,1);
            tmp     = condLevelsRaw;
            for kk = 1:numel(baseOrderKeywords)
                key  = baseOrderKeywords(kk);
                mask = contains(tmp, key, 'IgnoreCase', true);
                if any(mask)
                    ordered = [ordered; tmp(mask)]; %#ok<AGROW>
                    tmp(mask) = "";
                end
            end
            leftover     = tmp(tmp ~= "");
            orderedConds = [ordered; leftover];
            nConds       = numel(orderedConds);

            if nConds == 0
                continue;
            end

            % Aggregate per condition
            meanDeltaBIC   = nan(nConds,1);
            semDeltaBIC    = nan(nConds,1);
            meanLog10BF    = nan(nConds,1);
            semLog10BF     = nan(nConds,1);

            for ci = 1:nConds
                cond = orderedConds(ci);
                mask = string(subM.Condition) == cond;
                if ~any(mask)
                    continue;
                end

                dBIC   = subM.DeltaBIC(mask);
                logBF  = subM.logBF_BIC(mask);
                log10BF = logBF ./ log(10);  % convert natural log to log10

                meanDeltaBIC(ci) = mean(dBIC,   'omitnan');
                semDeltaBIC(ci)  = std(dBIC,    'omitnan') ./ sqrt(sum(~isnan(dBIC)));

                meanLog10BF(ci)  = mean(log10BF,'omitnan');
                semLog10BF(ci)   = std(log10BF, 'omitnan') ./ sqrt(sum(~isnan(log10BF)));
            end

            x = 1:nConds;

            figName = sprintf('BIC/BF Evidence (kappa) | Set=%s | Monkey=%s', ...
                               char(setName), char(mID));
            figure('Name', figName, 'Color','w');

            % ----------------- Top: mean DeltaBIC by condition -----------------
            subplot(2,1,1); hold on;
            bar(x, meanDeltaBIC);
            errorbar(x, meanDeltaBIC, semDeltaBIC, 'k.', 'LineStyle','none', 'LineWidth',1);
            yline(0,'k--','LineWidth',1);

            set(gca, 'XTick', x, 'XTickLabel', orderedConds, 'XTickLabelRotation', 45);
            ylabel('\DeltaBIC (noBias - with \kappa)');
            title(sprintf('BIC difference by condition\nSet: %s | Monkey: %s', ...
                  char(setName), char(mID)), 'Interpreter','none');
            hold off;

            % ----------------- Bottom: mean log10(BF_BIC) by condition --------
            subplot(2,1,2); hold on;
            bar(x, meanLog10BF);
            errorbar(x, meanLog10BF, semLog10BF, 'k.', 'LineStyle','none', 'LineWidth',1);
            yline(0,'k--','LineWidth',1);

            set(gca, 'XTick', x, 'XTickLabel', orderedConds, 'XTickLabelRotation', 45);
            ylabel('log_{10} BF_{BIC} (kappa vs noBias)');
            title('Bayes-factor evidence (BIC-based)', 'Interpreter','none');
            hold off;
        end
    end
end

function local_plotQTimecourse(Qtc)
% local_plotQTimecourse
% -------------------------------------------------------------------------
% Simple diagnostic plots of mean Q(action, TrialInBlock) per condition.
% For each (SetName, MonkeyRaw) pair, we make a figure with one pair of
% subplots per condition:
%   - Top:   meanQ vs TrialInBlock for each action (with SEM)
%   - Bottom:numBlocks vs TrialInBlock as a bar chart
%
% All top panels share the same y-limits; all bottom panels share the same
% y-limits.

    if isempty(Qtc)
        fprintf('  [QTimecourse] QTimecourse table is empty; no plots generated.\n');
        return;
    end

    if ~ismember('NumBlocks', Qtc.Properties.VariableNames)
        error('Qtc must contain a NumBlocks column. Build Qtc with block counts first.');
    end

    sets = unique(Qtc.SetName, 'stable');
    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    for si = 1:numel(sets)
        setName = sets(si);
        subSet  = Qtc(Qtc.SetName == setName, :);

        if isempty(subSet)
            continue;
        end

        monkeys = unique(subSet.MonkeyRaw, 'stable');

        for mi = 1:numel(monkeys)
            mID  = monkeys(mi);
            subM = subSet(subSet.MonkeyRaw == mID, :);
            if isempty(subM)
                continue;
            end

            % Ordered conditions
            condLevelsRaw = unique(string(subM.Condition), 'stable');
            ordered = strings(0,1);
            tmp     = condLevelsRaw;
            for kk = 1:numel(baseOrderKeywords)
                key  = baseOrderKeywords(kk);
                mask = contains(tmp, key, 'IgnoreCase', true);
                if any(mask)
                    ordered = [ordered; tmp(mask)]; %#ok<AGROW>
                    tmp(mask) = "";
                end
            end
            leftover     = tmp(tmp ~= "");
            orderedConds = [ordered; leftover];

            nConds = numel(orderedConds);
            nCols  = 4;
            nRows  = ceil(nConds / nCols);
            bigRows = 2 * nRows;   % each condition gets 2 stacked rows (Q + blocks)

            % Global y-limits for Q across all conditions/actions in this figure
            qLow  = subM.meanQ  - subM.semQ;
            qHigh = subM.meanQ  + subM.semQ;
            yMin = min(qLow,  [], 'omitnan');
            yMax = max(qHigh, [], 'omitnan');
            if ~isfinite(yMin) || ~isfinite(yMax)
                yMin = -1; yMax = 1;
            end
            if yMin == yMax
                yMin = yMin - 0.5;
                yMax = yMax + 0.5;
            end
            pad  = 0.05 * (yMax - yMin);
            yMin = yMin - pad;
            yMax = yMax + pad;

            % Global y-limits for NumBlocks across all conditions in this figure
            blockMax = max(subM.NumBlocks, [], 'omitnan');
            if ~isfinite(blockMax) || blockMax <= 0
                blockMax = 1;
            end
            blockYMax = blockMax * 1.1;

            figure('Name', sprintf('Q Timecourse | Set=%s | Monkey=%s', ...
                                   char(setName), char(mID)), ...
                   'Color','w');

            for ci = 1:nConds
                cond = orderedConds(ci);
                subC = subM(string(subM.Condition) == cond, :);
                if isempty(subC)
                    continue;
                end

                % Row/col position within the condition grid
                r0 = ceil(ci / nCols);
                c0 = mod(ci-1, nCols) + 1;

                % Linear indices for top (Q) and bottom (NumBlocks) plots
                topIdx    = (2*(r0-1))   * nCols + c0;   % row = 2*r0-1
                bottomIdx = (2*(r0-1)+1) * nCols + c0;   % row = 2*r0

                % ---------- Top panel: Q timecourse ----------
                subplot(bigRows, nCols, topIdx); hold on;

                actions = unique(subC.Action, 'stable');
                for ai = 1:numel(actions)
                    act  = actions(ai);
                    subA = subC(subC.Action == act, :);

                    [trialsSorted, idxSort] = sort(subA.TrialInBlock);
                    meanQ = subA.meanQ(idxSort);
                    semQ  = subA.semQ(idxSort);

                    plot(trialsSorted, meanQ, '-', 'LineWidth', 1.5);
                    errorbar(trialsSorted, meanQ, semQ, ...
                        'LineStyle','none','Color','k','HandleVisibility','off');
                end

                xlabel('Trial');
                ylabel('Q_{self}');
                title(sprintf('%s', char(cond)), 'Interpreter','none');
                legend(cellstr(actions), 'Location','northwest', 'Box','off');
                ylim([yMin, yMax]);
                hold off;

                % ---------- Bottom panel: NumBlocks per trial ----------
                subplot(bigRows, nCols, bottomIdx); hold on;

                % Aggregate NumBlocks over actions per TrialInBlock
                [Gtr, keyTr] = findgroups(subC.TrialInBlock);
                blocksPerTrial = splitapply(@(x) max(x, [], 'omitnan'), subC.NumBlocks, Gtr);

                [trialsSorted2, idx2] = sort(keyTr);
                blocksSorted = blocksPerTrial(idx2);

                bar(trialsSorted2, blocksSorted, EdgeColor="none");
                xlabel('Trial');
                ylabel('# Blocks');
                ylim([0, blockYMax]);
                title('Blocks per trial', 'Interpreter','none');

                hold off;
            end

            sgtitle(sprintf('Q timecourse | Set: %s | Monkey: %s', ...
                    char(setName), char(mID)), 'Interpreter','none');
        end
    end
end

function local_plotDeltaQTimecourse(DQtc)
% local_plotDeltaQTimecourse
% -------------------------------------------------------------------------
% Diagnostic plots of mean RL_DeltaQ over TrialInBlock, broken down
% by chosen action.
%
% For each (SetName, MonkeyRaw) pair, we create a figure with one pair of
% subplots per condition:
%   - Top:   meanDeltaQ vs TrialInBlock for each ChosenAction (with SEM)
%   - Bottom:numBlocks vs TrialInBlock as a bar chart
%
% All top panels share the same y-limits; all bottom panels share the same
% y-limits.

    if isempty(DQtc)
        fprintf('  [DeltaQTimecourse] DeltaQTimecourse table is empty; no plots.\n');
        return;
    end

    if ~ismember('NumBlocks', DQtc.Properties.VariableNames)
        error('DQtc must contain a NumBlocks column. Build DQtc with block counts first.');
    end

    sets = unique(DQtc.SetName, 'stable');
    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    for si = 1:numel(sets)
        setName = sets(si);
        subSet  = DQtc(DQtc.SetName == setName, :);
        if isempty(subSet)
            continue;
        end

        monkeys = unique(subSet.MonkeyRaw, 'stable');

        for mi = 1:numel(monkeys)
            mID  = monkeys(mi);
            subM = subSet(subSet.MonkeyRaw == mID, :);
            if isempty(subM)
                continue;
            end

            % Ordered conditions
            condLevelsRaw = unique(string(subM.Condition), 'stable');
            ordered = strings(0,1);
            tmp     = condLevelsRaw;
            for kk = 1:numel(baseOrderKeywords)
                key  = baseOrderKeywords(kk);
                mask = contains(tmp, key, 'IgnoreCase', true);
                if any(mask)
                    ordered = [ordered; tmp(mask)]; %#ok<AGROW>
                    tmp(mask) = "";
                end
            end
            leftover     = tmp(tmp ~= "");
            orderedConds = [ordered; leftover];

            nConds = numel(orderedConds);
            nCols  = 4;
            nRows  = ceil(nConds / nCols);
            bigRows = 2 * nRows;   % each condition gets 2 stacked rows (DeltaQ + blocks)

            % Global y-limits for DeltaQ across all conditions/actions
            dqLow  = subM.meanDeltaQ  - subM.semDeltaQ;
            dqHigh = subM.meanDeltaQ  + subM.semDeltaQ;
            yMin = min(dqLow,  [], 'omitnan');
            yMax = max(dqHigh, [], 'omitnan');
            if ~isfinite(yMin) || ~isfinite(yMax)
                yMin = -1; yMax = 1;
            end
            if yMin == yMax
                yMin = yMin - 0.5;
                yMax = yMax + 0.5;
            end
            pad  = 0.05 * (yMax - yMin);
            yMin = yMin - pad;
            yMax = yMax + pad;

            % Global y-limits for NumBlocks across all conditions
            blockMax = max(subM.NumBlocks, [], 'omitnan');
            if ~isfinite(blockMax) || blockMax <= 0
                blockMax = 1;
            end
            blockYMax = blockMax * 1.1;

            figure('Name', sprintf('DeltaQ Timecourse | Set=%s | Monkey=%s', ...
                                   char(setName), char(mID)), ...
                   'Color','w');

            for ci = 1:nConds
                cond = orderedConds(ci);
                subC = subM(string(subM.Condition) == cond, :);
                if isempty(subC)
                    continue;
                end

                % Row/col position within the condition grid
                r0 = ceil(ci / nCols);
                c0 = mod(ci-1, nCols) + 1;

                % Linear indices for top (DeltaQ) and bottom (NumBlocks) plots
                topIdx    = (2*(r0-1))   * nCols + c0;   % row = 2*r0-1
                bottomIdx = (2*(r0-1)+1) * nCols + c0;   % row = 2*r0

                % ---------- Top panel: DeltaQ timecourse ----------
                subplot(bigRows, nCols, topIdx); hold on;

                actions = unique(subC.ChosenAction, 'stable');
                for ai = 1:numel(actions)
                    act  = actions(ai);
                    subA = subC(subC.ChosenAction == act, :);

                    [trialsSorted, idxSort] = sort(subA.TrialInBlock);
                    meanDQ = subA.meanDeltaQ(idxSort);
                    semDQ  = subA.semDeltaQ(idxSort);

                    plot(trialsSorted, meanDQ, '-', 'LineWidth', 1.5);
                    errorbar(trialsSorted, meanDQ, semDQ, ...
                        'LineStyle','none','Color','k','HandleVisibility','off');
                end

                xlabel('Trial');
                ylabel('\DeltaQ (chosen - best other)');
                title(sprintf('%s', char(cond)), 'Interpreter','none');
                legend(cellstr(actions), 'Location','northwest', 'Box','off');
                ylim([yMin, yMax]);
                hold off;

                % ---------- Bottom panel: NumBlocks per trial ----------
                subplot(bigRows, nCols, bottomIdx); hold on;

                % Aggregate NumBlocks over actions per TrialInBlock
                [Gtr, keyTr] = findgroups(subC.TrialInBlock);
                blocksPerTrial = splitapply(@(x) max(x, [], 'omitnan'), subC.NumBlocks, Gtr);

                [trialsSorted2, idx2] = sort(keyTr);
                blocksSorted = blocksPerTrial(idx2);

                bar(trialsSorted2, blocksSorted, EdgeColor="none");
                xlabel('Trial');
                ylabel('# Blocks');
                ylim([0, blockYMax]);
                title('Blocks per trial', 'Interpreter','none');

                hold off;
            end

            sgtitle(sprintf('DeltaQ timecourse | Set: %s | Monkey: %s', ...
                    char(setName), char(mID)), 'Interpreter','none');
        end
    end
end

function local_plotTerminalQAdvantage(Tadv)
% local_plotTerminalQAdvantage
% -------------------------------------------------------------------------
% Visualization of terminal BuyAdvantage (Q_buy - max Q_other) across
% conditions for each (SetName, MonkeyRaw).
%
% For each SetName and MonkeyRaw:
%   - Create a figure showing boxplot/boxchart of BuyAdvantage grouped
%     by Condition.

    if isempty(Tadv)
        fprintf('  [terminalQAdvantage] table is empty; no plots.\n');
        return;
    end

    useBoxchart = exist('boxchart','file') == 2;
    sets = unique(Tadv.SetName, 'stable');
    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    for si = 1:numel(sets)
        setName = sets(si);
        subSet  = Tadv(Tadv.SetName == setName, :);
        if isempty(subSet)
            continue;
        end

        monkeys = unique(subSet.MonkeyRaw, 'stable');
        for mi = 1:numel(monkeys)
            mID  = monkeys(mi);
            subM = subSet(subSet.MonkeyRaw == mID, :);
            if isempty(subM)
                continue;
            end

            condLevelsRaw = unique(string(subM.Condition), 'stable');
            ordered = strings(0,1);
            tmp     = condLevelsRaw;
            for kk = 1:numel(baseOrderKeywords)
                key  = baseOrderKeywords(kk);
                mask = contains(tmp, key, 'IgnoreCase', true);
                if any(mask)
                    ordered = [ordered; tmp(mask)]; %#ok<AGROW>
                    tmp(mask) = "";
                end
            end
            leftover     = tmp(tmp ~= "");
            orderedConds = [ordered; leftover];

            condCats = categorical(string(subM.Condition), orderedConds, 'Ordinal', true);

            figure('Name', sprintf('Terminal BuyAdvantage | Set=%s | Monkey=%s', ...
                                   char(setName), char(mID)), ...
                   'Color','w');
            hold on;

            if useBoxchart
                boxchart(condCats, subM.BuyAdvantage);
            else
                boxplot(subM.BuyAdvantage, condCats);
            end

            yline(0,'k--','LineWidth',1);
            xlabel('Condition');
            ylabel('BuyAdvantage = Q_{buy} - max(Q_{other})');
            title(sprintf('Terminal Buy advantage (@ trial 15)\nSet: %s | Monkey: %s', ...
                  char(setName), char(mID)), 'Interpreter','none');

            hold off;
        end
    end
end

function [modelComparisonTable, modelSummaryByCond] = local_buildModelComparison(behaviorFitSummary, modelNameA, modelNameB)
% General pairwise model comparison:
%   modelNameA = "reference" (simpler) model
%   modelNameB = "alternative" (more complex) model
%
% For each (SetName, Monkey, Condition, DateID) where both models are fit:
%   DeltaBIC  = BIC_A - BIC_B
%   logBF_BIC = 0.5 * DeltaBIC  (log BF in favor of model B vs A)
%   BF_BIC    = exp(logBF_BIC)
%
% Positive DeltaBIC / logBF_BIC => evidence favoring model B.

    if isempty(behaviorFitSummary)
        modelComparisonTable = table();
        modelSummaryByCond   = table();
        return;
    end

    % Use canonical event tag just for bookkeeping
    events = unique(behaviorFitSummary.Event, 'stable');
    if ~isempty(events)
        canonicalEvt = events(1);
    else
        canonicalEvt = "";
    end

    maskA = behaviorFitSummary.ModelName == string(modelNameA);
    maskB = behaviorFitSummary.ModelName == string(modelNameB);

    fitsA = behaviorFitSummary(maskA,:);
    fitsB = behaviorFitSummary(maskB,:);

    if isempty(fitsA) || isempty(fitsB)
        warning('No fits found for %s or %s. modelComparisonTable will be empty.', ...
                char(modelNameA), char(modelNameB));
        modelComparisonTable = table();
        modelSummaryByCond   = table();
        return;
    end

    rows = [];

    for i = 1:height(fitsA)
        rowA = fitsA(i,:);

        % Match model B fit for same Set/Monkey/Condition/Session
        maskMatch = (fitsB.SetName   == rowA.SetName)   & ...
                    (fitsB.Monkey    == rowA.Monkey)    & ...
                    (fitsB.Condition == rowA.Condition) & ...
                    (fitsB.DateID    == rowA.DateID);

        idx = find(maskMatch, 1, 'first');
        if isempty(idx)
            continue;
        end
        rowB = fitsB(idx,:);

        if rowA.NumTrials ~= rowB.NumTrials
            warning('NumTrials mismatch for %s vs %s in Set=%s Monkey=%s Cond=%s Sess=%d. Skipping.', ...
                    char(modelNameA), char(modelNameB), ...
                    char(rowA.SetName), char(rowA.Monkey), char(rowA.Condition), rowA.DateID);
            continue;
        end

        % Likelihoods and IC differences
        LL_A = -rowA.NegLogLik;
        LL_B = -rowB.NegLogLik;

        dLL = LL_B - LL_A;               % improvement when going A -> B
        LRT = 2 * dLL;
        df  = rowB.NumParams - rowA.NumParams;

        if df <= 0
            p_LRT = NaN;
        else
            try
                p_LRT = 1 - chi2cdf(LRT, df);
            catch
                warning('chi2cdf not available. Setting p_LRT to NaN.');
                p_LRT = NaN;
            end
        end

        dAIC = rowA.AIC - rowB.AIC;      % >0 => AIC favors B
        dBIC = rowA.BIC - rowB.BIC;      % >0 => BIC favors B

        logBF_BIC = 0.5 * dBIC;          % log BF (B vs A)
        BF_BIC    = exp(logBF_BIC);

        % Evidence strength based on |DeltaBIC|
        absDBIC = abs(dBIC);
        if absDBIC < 2
            evidenceStrength = "none";
        elseif absDBIC < 6
            evidenceStrength = "positive";
        elseif absDBIC < 10
            evidenceStrength = "strong";
        else
            evidenceStrength = "very_strong";
        end

        if dBIC > 0
            evidenceDirection = "favor_modelB";
        elseif dBIC < 0
            evidenceDirection = "favor_modelA";
        else
            evidenceDirection = "tie";
        end

        evidenceLabel = evidenceStrength + "_" + evidenceDirection;

        newRow = table( ...
            rowA.SetName, ...
            canonicalEvt, ...
            rowA.Monkey, ...
            rowA.Condition, ...
            rowA.DateID, ...
            rowA.NumTrials, ...
            string(modelNameA), ...
            string(modelNameB), ...
            LL_A, LL_B, dLL, LRT, df, p_LRT, ...
            rowA.AIC, rowB.AIC, dAIC, ...
            rowA.BIC, rowB.BIC, dBIC, ...
            logBF_BIC, BF_BIC, ...
            evidenceStrength, evidenceDirection, evidenceLabel, ...
            'VariableNames', {'SetName','Event','Monkey','Condition','DateID', ...
                              'NumTrials','ModelA','ModelB', ...
                              'LL_A','LL_B','DeltaLL','LRT','df','p_LRT', ...
                              'AIC_A','AIC_B','DeltaAIC', ...
                              'BIC_A','BIC_B','DeltaBIC', ...
                              'logBF_BIC','BF_BIC', ...
                              'EvidenceStrength','EvidenceDirection','EvidenceLabel'});

        rows = [rows; newRow]; %#ok<AGROW>
    end

    modelComparisonTable = rows;

    % ---------- Summaries by condition (like kappaSummaryByCond) ----------
    if isempty(modelComparisonTable)
        modelSummaryByCond = table();
        return;
    end

    T = modelComparisonTable;
    [G, key] = findgroups(T(:,{'SetName','Monkey','Condition'}));

    N_sessions      = splitapply(@numel, T.DateID, G);
    mean_dBIC       = splitapply(@(x) mean(x,'omitnan'), T.DeltaBIC,  G);
    mean_dAIC       = splitapply(@(x) mean(x,'omitnan'), T.DeltaAIC,  G);
    mean_LRT        = splitapply(@(x) mean(x,'omitnan'), T.LRT,       G);
    mean_p_LRT      = splitapply(@(x) mean(x,'omitnan'), T.p_LRT,     G);
    frac_BIC_favorB = splitapply(@(x) mean(x > 0),       T.DeltaBIC,  G); % BIC favors B
    frac_AIC_favorB = splitapply(@(x) mean(x > 0),       T.DeltaAIC,  G); % AIC favors B

    mean_logBF_BIC  = splitapply(@(x) mean(x,'omitnan'), T.logBF_BIC, G);
    frac_BF_gt1_B   = splitapply(@(x) mean(x > 0),       T.logBF_BIC, G); % BF>1 for B
    frac_strongBIC  = splitapply(@(d) mean(d > 6),       T.DeltaBIC,  G);
    frac_veryStrongBIC = splitapply(@(d) mean(d > 10),   T.DeltaBIC,  G);

    modelSummaryByCond = key;
    modelSummaryByCond.ModelA = repmat(string(modelNameA), height(key), 1);
    modelSummaryByCond.ModelB = repmat(string(modelNameB), height(key), 1);

    modelSummaryByCond.N_sessions           = N_sessions;
    modelSummaryByCond.meanDeltaBIC         = mean_dBIC;
    modelSummaryByCond.meanDeltaAIC         = mean_dAIC;
    modelSummaryByCond.meanLRT              = mean_LRT;
    modelSummaryByCond.mean_p_LRT           = mean_p_LRT;
    modelSummaryByCond.fracBIC_favorsModelB = frac_BIC_favorB;
    modelSummaryByCond.fracAIC_favorsModelB = frac_AIC_favorB;

    modelSummaryByCond.meanLogBF_BIC            = mean_logBF_BIC;
    modelSummaryByCond.fracBF_gt1_forModelB     = frac_BF_gt1_B;
    modelSummaryByCond.fracStrongEvidenceBIC    = frac_strongBIC;
    modelSummaryByCond.fracVeryStrongEvidenceBIC= frac_veryStrongBIC;
end

function local_plotModelBIC_EvidenceByCondition(modelComparisonTable, modelNameA, modelNameB)
% General version of local_plotKappaBIC_EvidenceByCondition.
% Uses DeltaBIC and logBF_BIC from modelComparisonTable.
% Annotates, per condition, whether the improvement is significant (ΔBIC ≠ 0).

    if isempty(modelComparisonTable)
        fprintf('  [modelBFPlot] modelComparisonTable is empty; no plots.\n');
        return;
    end

    requiredVars = {'DeltaBIC','logBF_BIC','Condition','SetName','Monkey'};
    for i = 1:numel(requiredVars)
        if ~ismember(requiredVars{i}, modelComparisonTable.Properties.VariableNames)
            error('modelComparisonTable must contain column "%s".', requiredVars{i});
        end
    end

    sets = unique(modelComparisonTable.SetName, 'stable');
    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    for si = 1:numel(sets)
        setName = sets(si);
        subSet  = modelComparisonTable(modelComparisonTable.SetName == setName, :);
        if isempty(subSet), continue; end

        monkeys = unique(subSet.Monkey, 'stable');
        for mi = 1:numel(monkeys)
            mID  = monkeys(mi);
            subM = subSet(subSet.Monkey == mID, :);
            if isempty(subM), continue; end

            condLevelsRaw = unique(string(subM.Condition), 'stable');

            ordered = strings(0,1);
            tmp     = condLevelsRaw;
            for kk = 1:numel(baseOrderKeywords)
                key  = baseOrderKeywords(kk);
                mask = contains(tmp, key, 'IgnoreCase', true);
                if any(mask)
                    ordered = [ordered; tmp(mask)]; %#ok<AGROW>
                    tmp(mask) = "";
                end
            end
            leftover     = tmp(tmp ~= "");
            orderedConds = [ordered; leftover];
            nConds       = numel(orderedConds);
            if nConds == 0, continue; end

            % Session-level values
            condCatsAll     = categorical(string(subM.Condition), orderedConds, 'Ordinal', true);
            condIdxAll      = double(condCatsAll);
            valsAll_dBIC    = subM.DeltaBIC;
            valsAll_log10BF = subM.logBF_BIC ./ log(10);

            meanDeltaBIC = nan(nConds,1);
            semDeltaBIC  = nan(nConds,1);
            meanLog10BF  = nan(nConds,1);
            semLog10BF   = nan(nConds,1);

            for ci = 1:nConds
                cond = orderedConds(ci);
                mask = string(subM.Condition) == cond;
                if ~any(mask), continue; end

                dBIC    = subM.DeltaBIC(mask);
                logBF   = subM.logBF_BIC(mask);
                log10BF = logBF ./ log(10);

                meanDeltaBIC(ci) = mean(dBIC,   'omitnan');
                semDeltaBIC(ci)  = std(dBIC,    'omitnan') ./ sqrt(sum(~isnan(dBIC)));
                meanLog10BF(ci)  = mean(log10BF,'omitnan');
                semLog10BF(ci)   = std(log10BF, 'omitnan') ./ sqrt(sum(~isnan(log10BF)));
            end

            x = 1:nConds;
            figName = sprintf('BIC/BF Evidence | %s vs %s | Set=%s | Monkey=%s', ...
                              char(modelNameA), char(modelNameB), ...
                              char(setName), char(mID));
            figure('Name', figName, 'Color','w');

            % ---------- Top: DeltaBIC (A - B) ----------
            ax1 = subplot(2,1,1); hold(ax1,'on');
            bar(ax1, x, meanDeltaBIC);
            errorbar(ax1, x, meanDeltaBIC, semDeltaBIC, 'k.', ...
                     'LineStyle','none','LineWidth',1);
            yline(ax1, 0,'k--','LineWidth',1);
            set(ax1, 'XTick', x, 'XTickLabel', orderedConds, 'XTickLabelRotation', 45);
            ylabel(ax1, sprintf('\\DeltaBIC (%s - %s)', char(modelNameA), char(modelNameB)));
            title(ax1, sprintf('BIC difference by condition\nSet: %s | Monkey: %s', ...
                  char(setName), char(mID)), 'Interpreter','none');

            % --- Per-condition significance vs 0 on ΔBIC ---
            yl = ylim(ax1);
            yRange = yl(2) - yl(1);
            for ci = 1:nConds
                valsCond = valsAll_dBIC(condIdxAll == ci);
                valsCond = valsCond(~isnan(valsCond));
                if numel(valsCond) < 2
                    continue;
                end
                [~,p] = ttest(valsCond, 0, 'Tail','right');       % one-sided against 0. H1 = mean(DeltaBIC) > 0  (evidence favoring model B)
                if isnan(p) || p >= 0.05
                    continue;
                end

                if p < 0.001
                    star = '***';
                elseif p < 0.01
                    star = '**';
                else
                    star = '*';
                end

                yBarTop = meanDeltaBIC(ci) + semDeltaBIC(ci);
                if isnan(yBarTop)
                    yBarTop = meanDeltaBIC(ci);
                end
                yStar = yBarTop + 0.05*yRange;
                if yStar > yl(2)
                    yl(2) = yStar + 0.02*yRange;
                end
                text(ax1, ci, yStar, star, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',12,'FontWeight','bold');
            end
            ylim(ax1, yl);
            hold(ax1,'off');

            % ---------- Bottom: log10 BF (B vs A) ----------
            ax2 = subplot(2,1,2); hold(ax2,'on');
            bar(ax2, x, meanLog10BF);
            errorbar(ax2, x, meanLog10BF, semLog10BF, 'k.', ...
                     'LineStyle','none','LineWidth',1);
            yline(ax2, 0,'k--','LineWidth',1);
            set(ax2, 'XTick', x, 'XTickLabel', orderedConds, 'XTickLabelRotation', 45);
            ylabel(ax2, sprintf('log_{10} BF_{BIC} (%s vs %s)', char(modelNameB), char(modelNameA)));
            title(ax2, 'Bayes-factor evidence (BIC-based)', 'Interpreter','none');
            hold(ax2,'off');
        end
    end
end

function local_plotModelBIC_AllModelsByCondition(behaviorFitSummary, MODEL_TO_RUN)
% local_plotModelBIC_AllModelsByCondition
% -------------------------------------------------------------------------
% Plots mean BIC (± SEM) of multiple RL models by condition, separately for
% each (SetName, Monkey) pair.
%
% INPUTS
%   behaviorFitSummary : table allResults.behaviorFits
%                        (requires columns: SetName, Event, Monkey,
%                         Condition, ModelName, BIC)
%   MODEL_TO_RUN          : (optional) cellstr/string array of model names to
%                        include, in desired plotting order.
%                        If omitted, all models present are used.
%
% For each SetName × Monkey:
%   - X-axis: Condition (AI / Replay / Decoy / Live ordered if present)
%   - Bars:   mean BIC per model (grouped bars)
%   - Error bars: ± SEM across sessions
%   - Lower BIC = better model.
% -------------------------------------------------------------------------

    if isempty(behaviorFitSummary)
        fprintf('[local_plotModelBIC_AllModelsByCondition] behaviorFitSummary is empty.\n');
        return;
    end

    requiredVars = {'SetName','Event','Monkey','Condition','ModelName','BIC'};
    for i = 1:numel(requiredVars)
        if ~ismember(requiredVars{i}, behaviorFitSummary.Properties.VariableNames)
            error('behaviorFitSummary must contain column "%s".', requiredVars{i});
        end
    end

    % Use canonical event (first one) for plotting, as in rest of pipeline
    events = unique(behaviorFitSummary.Event, 'stable');
    if isempty(events)
        fprintf('[local_plotModelBIC_AllModelsByCondition] No Event labels found.\n');
        return;
    end
    canonicalEvt = events(1);

    % Global model list (optional override)
    if nargin < 2 || isempty(MODEL_TO_RUN)
        MODEL_TO_RUN = unique(behaviorFitSummary.ModelName, 'stable');
    end
    MODEL_TO_RUN = string(MODEL_TO_RUN(:));  % ensure string column vector
    nModels   = numel(MODEL_TO_RUN);

    % Condition ordering keywords (like other plotting helpers)
    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    sets = unique(behaviorFitSummary.SetName, 'stable');

    for si = 1:numel(sets)
        setName = sets(si);
        subSet  = behaviorFitSummary(behaviorFitSummary.SetName == setName & ...
                                     behaviorFitSummary.Event   == canonicalEvt, :);
        if isempty(subSet)
            continue;
        end

        monkeys = unique(subSet.Monkey, 'stable');

        for mi = 1:numel(monkeys)
            mID  = monkeys(mi);
            subM = subSet(subSet.Monkey == mID, :);
            if isempty(subM)
                continue;
            end

            % ----- Condition ordering -----
            condLevelsRaw = unique(string(subM.Condition), 'stable');
            ordered = strings(0,1);
            tmp     = condLevelsRaw;
            for kk = 1:numel(baseOrderKeywords)
                key  = baseOrderKeywords(kk);
                mask = contains(tmp, key, 'IgnoreCase', true);
                if any(mask)
                    ordered = [ordered; tmp(mask)]; %#ok<AGROW>
                    tmp(mask) = "";
                end
            end
            leftover     = tmp(tmp ~= "");
            orderedConds = [ordered; leftover];
            nConds       = numel(orderedConds);
            if nConds == 0
                continue;
            end

            % ----- Aggregate BIC by (Condition, ModelName) -----
            condStrAll  = string(subM.Condition);
            modelStrAll = string(subM.ModelName);

            [G, key] = findgroups(table(condStrAll, modelStrAll, ...
                                 'VariableNames', {'Condition','ModelName'}));

            meanBIC = splitapply(@(x) mean(x,'omitnan'), subM.BIC, G);
            semBIC  = splitapply(@(x) std(x,'omitnan')./sqrt(sum(~isnan(x))), ...
                                  subM.BIC, G);

            % ----- Fill condition × model matrices -----
            BIC_mean_mat = nan(nConds, nModels);
            BIC_sem_mat  = nan(nConds, nModels);

            for r = 1:height(key)
                condR  = string(key.Condition(r));
                modelR = string(key.ModelName(r));

                ci = find(orderedConds == condR, 1, 'first');
                mi2 = find(MODEL_TO_RUN   == modelR, 1, 'first');
                if isempty(ci) || isempty(mi2)
                    continue;
                end

                BIC_mean_mat(ci, mi2) = meanBIC(r);
                BIC_sem_mat(ci, mi2)  = semBIC(r);
            end

            % If everything is NaN, skip
            if all(isnan(BIC_mean_mat), 'all')
                continue;
            end

            % ----- Plot grouped bars with error bars -----
            figure('Name', sprintf('BIC by Condition | Set=%s | Monkey=%s', ...
                                   char(setName), char(mID)), ...
                   'Color','w');

            ax = axes; hold(ax, 'on');

            % grouped bar: each row = condition, each column = model
            bh = bar(ax, 1:nConds, BIC_mean_mat, 'grouped');
            % bar returns an array of bar objects, one per model

            % Error bars
            for m = 1:nModels
                thisBar = bh(m);
                % XEndPoints property available in recent MATLAB; fallback if absent
                if isprop(thisBar, 'XEndPoints')
                    xCenters = thisBar.XEndPoints;
                else
                    % approximate centers from XData
                    xCenters = thisBar.XData + thisBar.XOffset;
                end
                yVals = BIC_mean_mat(:,m);
                yErr  = BIC_sem_mat(:,m);

                valid = ~isnan(yVals);
                if any(valid)
                    errorbar(ax, xCenters(valid), yVals(valid), yErr(valid), ...
                             'k', 'LineStyle', 'none', 'LineWidth', 1);
                end
            end

            % Axis labels / legend
            set(ax, 'XTick', 1:nConds, ...
                    'XTickLabel', orderedConds, ...
                    'XTickLabelRotation', 45);
            ylabel(ax, 'BIC');
            xlabel(ax, 'Condition');
            title(ax, sprintf('BIC by Condition\nSet: %s | Monkey: %s', ...
                              char(setName), char(mID)), ...
                  'Interpreter','none');

            legend(ax, cellstr(MODEL_TO_RUN), 'Location','northoutside', ...
                   'Orientation','horizontal', 'Box','off');

            % Make sure y-axis has a bit of headroom
            yAll = BIC_mean_mat(:);
            yAll = yAll(isfinite(yAll));
            if ~isempty(yAll)
                yMin = min(yAll);
                yMax = max(yAll);
                if yMin == yMax
                    yMin = yMin - 1;
                    yMax = yMax + 1;
                end
                pad  = 0.05 * (yMax - yMin);
                ylim(ax, [yMin - pad, yMax + pad]);
            end

            box(ax, 'off');
            hold(ax, 'off');
        end
    end
end

function local_plotBehaviorParamByMonkey(fitTable, paramNamesPerModel)

    if isempty(fitTable), return; end

    sets   = unique(fitTable.SetName, 'stable');
    events = unique(fitTable.Event,   'stable');
    if isempty(events), return; end
    eventUse = events(1);  % canonical event

    baseOrderKeywords = ["AI","Replay","Decoy","Live"];

    for si = 1:numel(sets)
        setName = sets(si);
        subSetSet = fitTable(fitTable.SetName == setName & ...
                             fitTable.Event   == eventUse, :);
        if isempty(subSetSet), continue; end

        monkeys = unique(subSetSet.Monkey, 'stable');

        for mi = 1:numel(monkeys)
            monkeyID = monkeys(mi);
            subM = subSetSet(subSetSet.Monkey == monkeyID, :);
            if isempty(subM), continue; end

            models = unique(subM.ModelName, 'stable');

            for mj = 1:numel(models)
                modelName = models(mj);
                subT = subM(subM.ModelName == modelName, :);
                if isempty(subT), continue; end

                thetaExample = subT.Theta{1};
                P = numel(thetaExample);

                modelNameChar = char(modelName);
                if isfield(paramNamesPerModel, modelNameChar)
                    paramNames = paramNamesPerModel.(modelNameChar);
                    if numel(paramNames) ~= P
                        paramNames = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
                    end
                else
                    paramNames = arrayfun(@(k)sprintf('param%d',k), 1:P, 'UniformOutput', false);
                end

                condLevelsRaw = unique(string(subT.Condition), 'stable');
                ordered = strings(0,1);
                tmp     = condLevelsRaw;
                for kk = 1:numel(baseOrderKeywords)
                    key  = baseOrderKeywords(kk);
                    mask = contains(tmp, key, 'IgnoreCase', true);
                    if any(mask)
                        ordered = [ordered; tmp(mask)]; %#ok<AGROW>
                        tmp(mask) = "";
                    end
                end
                leftover     = tmp(tmp ~= "");
                orderedConds = [ordered; leftover];

                condCatsAll = categorical(string(subT.Condition), orderedConds, 'Ordinal', true);
                condIdxAll  = double(condCatsAll);
                K           = numel(orderedConds);
                if K == 0, continue; end

                f = figure('Name', sprintf('RL Params by Condition - %s - Monkey %s - Model %s', ...
                                           string(setName), string(monkeyID), modelNameChar), ...
                           'Color','w');
                nCols = min(3,P);
                nRows = ceil(P / nCols);

                for p = 1:P
                    ax = subplot(nRows,nCols,p); hold(ax, 'on');

                    valsAll = cellfun(@(theta) theta(p), subT.Theta);

                    meanVals = nan(1,K);
                    stdVals  = nan(1,K);
                    nVals    = zeros(1,K);

                    for k = 1:K
                        maskK   = (condIdxAll == k);
                        these   = valsAll(maskK);
                        nVals(k)= numel(these);
                        if nVals(k) > 0
                            meanVals(k) = mean(these, 'omitnan');
                            stdVals(k)  = std(these, 0, 'omitnan');
                        end
                    end

                    bar(ax, 1:K, meanVals, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor','k');
                    errorbar(ax, 1:K, meanVals, stdVals, ...
                             'k', 'LineStyle','none', 'LineWidth',1);

                    yline(ax, 0, 'k--', 'LineWidth',1);

                    condStr    = cellstr(orderedConds);
                    tickLabels = cell(K,1);
                    for k = 1:K
                        tickLabels{k} = sprintf('%s\\newline(n=%d)', condStr{k}, nVals(k));
                    end
                    set(ax, 'XTick', 1:K, ...
                            'XTickLabel', tickLabels, ...
                            'TickLabelInterpreter', 'tex');

                    xlabel(ax, 'Condition');
                    ylabel(ax, 'Param value');
                    title(ax, paramNames{p}, 'Interpreter','none');

                    if K >= 2
                        [sigPairs, pvals] = local_computeSig(valsAll, condIdxAll, K);
                        local_plotSignificanceLines(ax, sigPairs, pvals);
                    end

                    hold(ax, 'off');
                end

                if exist('sgtitle','file') == 2
                    sgtitle(f, sprintf('Set: %s | Monkey: %s | Model: %s', ...
                        string(setName), string(monkeyID), modelNameChar), 'Interpreter','none');
                end
            end
        end
    end
end
