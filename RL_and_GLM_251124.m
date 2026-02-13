%% Buy Bias Analysis (Cleaned v3): RL-based Behavior & Neural Modeling
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
%   - TARGET_EVTS is ONLY used where neural event timing matters
%     (firing-rate windows, FR extraction, GLMs, clustering).
%   - RL fitting and RL latents are now event-agnostic:
%       * A single "canonical" event (TARGET_EVTS{1}) is used to build the
%         behavior table for RL.
%       * RL parameters are fit per (Set × Monkey × Condition × Session),
%         independent of which neural event is being analyzed.
%       * Once RL latents are computed for this canonical behavior table,
%         they are copied into every event-specific trial table for GLM.
%   - A simpler RL model 'SelfOnlyNoBias' (no BUY-bias kappa) is added
%     alongside 'SelfOnly', 'SelfBubble', and 'SelfOpp'.
%   - Progress messages (fprintf) are sprinkled throughout long loops
%     so you can monitor run progress in the MATLAB console.
% =========================================================================

clear; clc; close all;
warning('off');  % suppress noisy warnings (fitlm rank deficiency, etc.)

%% ======================= PATHS & CONDITION SETS ==========================
OUTDIR = 'C:\Users\plattlab\Tim\Stock_market_experiment\Tim\condition_packs_251118';

% Condition sets to analyze (each cell-array is a "set")
cond_sets = {
    {'AI','Replay','Decoy','Live'}, ...
    % {'OT AI','OT Replay','OT Decoy','OT Live'}, ...
    % {'Saline AI','Saline Replay','Saline Decoy','Saline Live'}
    };

setNames = {
    'baseline', ...
    % 'OT', ...
    % 'Saline'
    };

% Target events and windows for neural analysis (ms).
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

% RL model used to generate latents for neural GLM:
%   - Fixed model name:
%       'SelfOnly' | 'SelfOnlyNoBias' | 'SelfBubble' | 'SelfOpp'
%   - Automatic: 'best_by_AIC' | 'best_by_BIC' (per group, per set)
RL_MODEL_FOR_NEURAL = 'SelfOnlyNoBias';

% Number of neuron clusters for coefficient-based clustering.
% Interpreted as K_max (max number of clusters to test per group).
N_NEURON_CLUSTERS = 5;

% Bootstrap settings
N_BOOT = 100;
modelsForBootstrap = {'SelfOnly','SelfBubble','SelfOpp'};  % can add 'SelfOnlyNoBias' if desired

% Optimization controls
NUM_RESTARTS = 3;  % restarts per fit group per model

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

    fprintf('\n[PART 0] Event %s (%d/%d)\n', string(TARGET_EVT), pair, numel(TARGET_EVTS));

    for s = 1:numel(cond_sets)
        conds   = cond_sets{s};
        setName = setNames{s};

        fprintf('  Set: %s | Event: %s\n', setName, string(TARGET_EVT));

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

                % Portfolio value at analysis event (pre or post)
                if string(TARGET_EVT) == "m1rew2on"
                    portfolioValueAtEvent = double(etab.postPortfolio{1}(:));
                else
                    portfolioValueAtEvent = double(etab.prePortfolio{1}(:));
                end

                % Opponent portfolio (per trial)
                opponentPortfolioValue = etab.postPortfolioOpp{1}(:);

                % Reward components
                choiceOutcomeReward = etab.fb1Off{1}(:) - etab.fb1On{1}(:);
                assert(sum(choiceOutcomeReward < 0) == 0);

                % Dividend-based portfolio reward
                % portfolioRewardFromDividend = etab.postPortfolio{1}(:) .* etab.divPerShare{1}(:);
                portfolioRewardFromDividend = etab.tfb2Off{1}(:) - etab.tfb2On{1}(:);
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
        [blockID, ~] = findgroups(trialTable(:,{'Monkey','Condition','Session','Market'}));
        trialTable.BlockID = blockID;

        trialTable.TrialInBlock = nan(height(trialTable),1);
        uBlocks = unique(blockID);
        for b = uBlocks.'
            idxBlock = find(blockID == b);
            [~, sortIdx] = sort(trialTable.Trial(idxBlock));
            trialTable.TrialInBlock(idxBlock(sortIdx)) = (1:numel(idxBlock)).';
        end

        % Bubble flag (1 if Bubble name contains 'bubble')
        bubbleCats   = categories(trialTable.Bubble);
        isBubbleCat  = contains(bubbleCats, 'bubble', 'IgnoreCase', true);
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

        fprintf('    -> trialTableRaw: %d trials, %d neurons\n', ...
                height(trialTable), numel(neurons));
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
    double.empty(0,1), ...  % Session
    string.empty(0,1), ...  % ModelName
    double.empty(0,1), ...  % NumTrials
    double.empty(0,1), ...  % NumParams
    double.empty(0,1), ...  % NegLogLik
    double.empty(0,1), ...  % AIC
    double.empty(0,1), ...  % BIC
    cell(0,1), ...          % Theta
    'VariableNames', {'SetName','Event','Monkey','Condition','Session','ModelName', ...
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
            groupVars = {'Monkey','Condition','Session'};
        case 'pooled'
            groupVars = {'Condition','Session'};
        otherwise
            error('Unknown FIT_SCOPE: %s', FIT_SCOPE);
    end

    [mcGroup, mcKey] = findgroups(trialTable(:,groupVars));
    nGroupsMC = max(mcGroup);
    fprintf('    -> %d behavior groups for RL\n', nGroupsMC);

    groupTables = cell(nGroupsMC,1);
    modelList   = {'SelfOnly','SelfOnlyNoBias','SelfBubble','SelfOpp'};

    parfor g = 1:nGroupsMC
        idxMC   = (mcGroup == g);
        sessData = trialTable(idxMC,:);
        sessData = sortrows(sessData, {'Session','BlockID','TrialInBlock'});

        if ismember('Monkey', groupVars)
            monkeyID = char(mcKey.Monkey(g));
        else
            monkeyID = 'ALL';
        end
        conditionID = char(mcKey.Condition(g));
        sessionID   = mcKey.Session(g);

        localTbl = behaviorFitSummary([],:);  % empty template

        if height(sessData) < 10
            groupTables{g} = localTbl;
            continue;
        end

        for m = 1:numel(modelList)
            modelName = modelList{m};
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
                sessionID, ...
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
paramNamesPerModel.SelfOnlyNoBias = {'alpha','beta','omega'};
paramNamesPerModel.SelfOnly       = {'alpha','beta','omega','kappaBuy'};
paramNamesPerModel.SelfBubble     = {'alpha','beta','omega_nonBubble','omega_bubble','kappa_nonBubble','kappa_bubble'};
paramNamesPerModel.SelfOpp        = {'alpha_self','alpha_opp','beta','omega_self','omega_opp','gamma','kappaBuy'};

allResults.behaviorParamTable = local_behaviorParamTable(allResults.behaviorFits, paramNamesPerModel);
fprintf('  -> behaviorParamTable rows: %d\n', height(allResults.behaviorParamTable));
%%
% local_plotBehaviorParamComparisonBySet(allResults.behaviorFits, paramNamesPerModel);

% % Optional bootstrap visualization (kept commented so it can be re-enabled)
% for s = 1:numel(setNames)
%     setName = setNames{s};
%     if isfield(allResults.behaviorBootstrap, setName)
%         visualizeBehaviorBootstrap(allResults.behaviorBootstrap.(setName), ...
%             paramNamesPerModel, setName, []);
%     end
% end

%% =========================================================================
% PART 1c-2: RECOMPUTE RL LATENTS FOR NEURAL ANALYSIS
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
            groupVars = {'Monkey','Condition','Session'};
        case 'pooled'
            groupVars = {'Condition','Session'};
        otherwise
            error('Unknown FIT_SCOPE: %s', FIT_SCOPE);
    end

    [mcGroup, mcKey] = findgroups(trialTable(:,groupVars));
    nGroupsMC = max(mcGroup);
    fprintf('    -> %d groups for RL latents\n', nGroupsMC);

    behaviorFits_local        = allResults.behaviorFits;
    RL_MODEL_FOR_NEURAL_local = RL_MODEL_FOR_NEURAL;

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
        sessData = sortrows(sessData, {'Session','BlockID','TrialInBlock'});
        if isempty(sessData)
            continue;
        end

        if ismember('Monkey', groupVars)
            monkeyID = char(mcKey.Monkey(g));
        else
            monkeyID = 'ALL';
        end
        conditionID = char(mcKey.Condition(g));
        sessionID   = mcKey.Session(g);

        modelNameForGroup = RL_MODEL_FOR_NEURAL_local;

        baseMask = strcmp(behaviorFits_local.SetName,   setName) & ...
                   strcmp(behaviorFits_local.Monkey,    monkeyID) & ...
                   strcmp(behaviorFits_local.Condition, conditionID) & ...
                   behaviorFits_local.Session == sessionID;

        if strcmpi(RL_MODEL_FOR_NEURAL_local,'best_by_AIC') || strcmpi(RL_MODEL_FOR_NEURAL_local,'best_by_BIC')
            fitsMC = behaviorFits_local(baseMask,:);
            if isempty(fitsMC)
                continue;
            end
            if strcmpi(RL_MODEL_FOR_NEURAL_local,'best_by_AIC')
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
            case 'SelfOnly'
                [~, latentsSelf] = local_rlNegLogLik_SelfOnly(thetaHat, sessData, BUY_NON_BUY_ONLY);
            case 'SelfOnlyNoBias'
                [~, latentsSelf] = local_rlNegLogLik_SelfOnlyNoBias(thetaHat, sessData, BUY_NON_BUY_ONLY);
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
%%
local_plotTerminalQDistributions(allResults.terminalQTable, cond_sets, setNames);

%% =========================================================================
% PART 2: NEURAL GLM (PER NEURON)
% =========================================================================
fprintf('\n=== PART 2: Neuron-wise GLM analysis ===\n');

neuronGLMResults = struct('SetName',{},'Event',{},'NeuronIndex',{}, ...
                          'Monkey',{},'Condition',{},'Session',{}, ...
                          'CoefNames',{},'CoefValues',{},'ModelFormula',{});

for pair = 1:numel(TARGET_EVTS)
    TARGET_EVT = TARGET_EVTS{pair};

    for s = 1:numel(cond_sets)
        setName = setNames{s};
        if ~isfield(allResults, setName) || ~isfield(allResults.(setName), TARGET_EVT) ...
                || ~isfield(allResults.(setName).(TARGET_EVT), 'trialTableWithRL')
            continue;
        end

        trialTable = allResults.(setName).(TARGET_EVT).trialTableWithRL;
        neurons    = allResults.(setName).(TARGET_EVT).neurons;

        if isempty(neurons)
            fprintf('  [PART 2] No neurons for set %s, event %s.\n', setName, TARGET_EVT);
            continue;
        end

        if BUY_NON_BUY_ONLY
            glmFormula = 'Z_fr ~ 1 + RL_predictionError + isBuy + portfolioValueAtEvent + opponentPortfolioValue_lag1 + Condition';
        else
            glmFormula = 'Z_fr ~ 1 + RL_predictionError + isBuy + isSell + portfolioValueAtEvent + opponentPortfolioValue_lag1 + Condition';
        end

        nNeurons = numel(neurons);
        tmpResults(nNeurons) = struct( ...
            'SetName',[], 'Event',[], 'NeuronIndex',[], ...
            'Monkey',[], 'Condition',[], 'Session',[], ...
            'CoefNames',[], 'CoefValues',[], 'ModelFormula',[]);

        fprintf('  [PART 2] Set %s | Event %s | %d neurons\n', setName, TARGET_EVT, nNeurons);

        parfor n = 1:nNeurons
            rows = neurons(n).trialRowIdx;
            if isempty(rows)
                continue;
            end

            glmData = table();
            glmData.Z_fr                        = neurons(n).FR(:);
            glmData.RL_predictionError          = trialTable.RL_predictionError(rows);
            glmData.isBuy                       = trialTable.isBuy(rows);
            if ~BUY_NON_BUY_ONLY
                glmData.isSell                  = trialTable.isSell(rows);
            end
            glmData.portfolioValueAtEvent       = trialTable.portfolioValueAtEvent(rows);
            glmData.opponentPortfolioValue_lag1 = trialTable.opponentPortfolioValue_lag1(rows);
            glmData.Condition                   = trialTable.Condition(rows);

            validRows = ~isnan(glmData.Z_fr);
            if sum(validRows) < 10
                continue;
            end
            glmData = glmData(validRows,:);

            try
                mdl = fitlm(glmData, glmFormula);
            catch
                continue;
            end

            coefTbl   = mdl.Coefficients;
            coefNames = coefTbl.Properties.RowNames;
            coefVals  = coefTbl.Estimate;

            r = struct();
            r.SetName      = setName;
            r.Event        = TARGET_EVT;
            r.NeuronIndex  = n;
            r.Monkey       = neurons(n).Monkey;
            r.Condition    = neurons(n).Condition;
            r.Session      = neurons(n).Session;
            r.CoefNames    = coefNames;
            r.CoefValues   = coefVals;
            r.ModelFormula = glmFormula;

            tmpResults(n) = r;
        end

        maskNonEmpty = ~arrayfun(@(x) isempty(x.SetName), tmpResults);
        neuronGLMResults = [neuronGLMResults, tmpResults(maskNonEmpty)]; %#ok<AGROW>
    end
end

allResults.neuronGLMResults = neuronGLMResults;
fprintf('  -> neuronGLMResults entries: %d\n', numel(neuronGLMResults));

%% =========================================================================
% PART 2a: OVERALL COEFFICIENTS PER MONKEY × CONDITION (OPTIONAL)
% =========================================================================
if ~isempty(neuronGLMResults)
    for s = 1:numel(setNames)
        setFilter = setNames{s};
        for pair = 1:numel(TARGET_EVTS)
            eventFilter = TARGET_EVTS{pair};

            mask = strcmp({neuronGLMResults.SetName}, setFilter) & ...
                   strcmp(cellstr({neuronGLMResults.Event}), cellstr(eventFilter));
            sub = neuronGLMResults(mask);
            if isempty(sub), continue; end

            coefNames = sub(1).CoefNames;
            nCoefs    = numel(coefNames);

            monkeys    = string({sub.Monkey})';
            conditions = string({sub.Condition})';
            groupsTbl  = table(monkeys, conditions, ...
                               'VariableNames',{'Monkey','Condition'});
            [grpID, grpKeyTbl] = findgroups(groupsTbl);

            % Example: only show one coefficient index (2) to limit plots;
            % this can be generalized back to 1:nCoefs if desired.
            for p = 2
                vals = cellfun(@(v)v(p), {sub.CoefValues})';

                figure('Color','w');
                try
                    boxchart(grpID, vals);
                catch
                    boxplot(vals, grpID);
                end

                ax = gca;
                ax.XTick = 1:height(grpKeyTbl);
                ax.XTickLabel = strcat(grpKeyTbl.Monkey,"_",grpKeyTbl.Condition);
                ax.XTickLabelRotation = 45;
                ax.TickLabelInterpreter = 'none';

                ylabel(coefNames{p}, 'Interpreter','none');
                title(sprintf('Coeff %s | Set=%s | Event=%s', ...
                      coefNames{p}, setFilter, eventFilter), ...
                      'Interpreter','none');
            end
        end
    end
end

%% =========================================================================
% PART 2b: NEURON COEFFICIENT CLUSTERING (BY MONKEY × CONDITION)
% =========================================================================
if ~isempty(allResults) && isfield(allResults,'neuronGLMResults')
    neuronGLMResults = allResults.neuronGLMResults;
else
    error('Expected allResults.neuronGLMResults to exist from PART 2.');
end

K_max = N_NEURON_CLUSTERS;

fprintf('\n=== PART 2b: Clustering neurons by condition-invariant coefficients ===\n');
fprintf('    (Per Set × Event × Monkey × Condition, K in [1, %d])\n', K_max);

allResults.neuronClusters = local_clusterNeuronCoefficientsByGroup_optimal( ...
    neuronGLMResults, K_max);

fprintf('=== PART 2b: Plotting neuron clusters (global y-limits across groups) ===\n');
%%
% local_plotNeuronClusters_condInvariant_optimal(allResults.neuronClusters);

fprintf('\n=== Analysis complete. Results are stored in allResults. ===\n');

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
            theta0 = [0.2, 3.0, 1.0, 0.0];          % [alpha, beta, omega, kappaBuy]
            lb     = [0.0, 0.0, 0, -5.0];
            ub     = [1.0, 20.0,  5.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_SelfOnly(theta, sessData, BUY_NON_BUY_ONLY);

        case 'SelfOnlyNoBias'
            theta0 = [0.2, 3.0, 1.0];              % [alpha, beta, omega], no kappa
            lb     = [0.0, 0.0, 0];
            ub     = [1.0, 20.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_SelfOnlyNoBias(theta, sessData, BUY_NON_BUY_ONLY);

        case 'SelfBubble'
            theta0 = [0.2, 3.0, 1.0, 1.0, 0.0, 0.0];
            lb     = [0.0, 0.0, 0, 0, -5.0, -5.0];
            ub     = [1.0, 20.0,  5.0,  5.0,  5.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_SelfBubble(theta, sessData, BUY_NON_BUY_ONLY);

        case 'SelfOpp'
            theta0 = [0.2, 0.2, 3.0, 1.0, 1.0, 0.5, 0.0];
            lb     = [0.0, 0.0, 0.0, 0, 0, -5.0, -5.0];
            ub     = [1.0, 1.0, 20.0,  5.0,  5.0,  5.0,  5.0];
            objFun = @(theta) local_rlNegLogLik_SelfOpp(theta, sessData, BUY_NON_BUY_ONLY);

        otherwise
            error('Unknown modelName: %s', modelName);
    end
end

function [negLL, latents] = local_rlNegLogLik_SelfOnly(theta, sessData, buyNonBuyOnly)
    alpha = max(min(theta(1),1),0);
    beta  = theta(2);
    omega = theta(3);
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
            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = Rc + omega * Rp;

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

function [negLL, latents] = local_rlNegLogLik_SelfOnlyNoBias(theta, sessData, buyNonBuyOnly)
    alpha = max(min(theta(1),1),0);
    beta  = theta(2);
    omega = theta(3);

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
            r = Rc + omega * Rp;

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

function [negLL, latents] = local_rlNegLogLik_SelfBubble(theta, sessData, buyNonBuyOnly)
    alpha          = max(min(theta(1),1),0);
    beta           = theta(2);
    omega_nonBub   = theta(3);
    omega_bub      = theta(4);
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

            r = Rc + omega * Rp;

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
    alpha_self = max(min(theta(1),1),0);
    alpha_opp  = max(min(theta(2),1),0);
    beta       = theta(3);
    omega_self = theta(4);
    omega_opp  = theta(5);
    gamma      = theta(6);
    kappa_buy  = theta(7);

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

    negLL = 0; tCounter = 0;

    for b = uBlocks.'
        idxBlock = find(blockID == b);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q_self = zeros(nA,1);
        Q_opp  = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k); tCounter = tCounter + 1;

            if validSelf(t)
                aSelf = actIdxSelf(t);

                Rc = sessData.choiceOutcomeReward(t);
                Rp = sessData.portfolioRewardFromDividend(t);
                if isnan(Rc), Rc = 0; end
                if isnan(Rp), Rp = 0; end

                rSelf = Rc + omega_self * Rp;

                VtSelf = Q_self(aSelf);
                valueSelf(tCounter) = VtSelf;
                deltaSelf(tCounter) = rSelf - VtSelf;

                V = Q_self + gamma * Q_opp;
                logits = beta * V;
                logits(1) = logits(1) + kappa_buy;

                maxLogit = max(logits);
                p = exp(logits - maxLogit); p = p ./ sum(p);
                pChosen = max(p(aSelf), realmin);
                negLL = negLL - log(pChosen);

                Q_self(aSelf) = Q_self(aSelf) + alpha_self * (rSelf - VtSelf);
            end

            if validOpp(t)
                aOpp = actIdxOpp(t);
                Ropp = sessData.opponentPortfolioValue(t);
                if isnan(Ropp), Ropp = 0; end
                rOpp = omega_opp * Ropp;

                VtOpp = Q_opp(aOpp);
                valueOpp(tCounter) = VtOpp;
                deltaOpp(tCounter) = rOpp - VtOpp;

                Q_opp(aOpp) = Q_opp(aOpp) + alpha_opp * (rOpp - VtOpp);
            end

            Q_self_hist(tCounter,:) = Q_self;
            Q_opp_hist(tCounter,:)  = Q_opp;
        end
    end

    if nargout > 1
        latentsSelf = struct('value',valueSelf,'delta',deltaSelf,'Q',Q_self_hist);
        latentsOpp  = struct('value',valueOpp,'delta',deltaOpp,'Q',Q_opp_hist);
    end
end

function behaviorBootstrap = local_rl_parametric_block_bootstrap( ...
    trialTable, behaviorFitSummary, BUY_NON_BUY_ONLY, modelsForBootstrap, N_BOOT, FIT_SCOPE)

    if nargin < 6 || isempty(FIT_SCOPE), FIT_SCOPE = 'per_monkey'; end
    if nargin < 5 || isempty(N_BOOT),   N_BOOT   = 200; end
    if nargin < 4 || isempty(modelsForBootstrap)
        modelsForBootstrap = {'SelfOnly','SelfBubble','SelfOpp'};
    end

    rng(0);

    behaviorBootstrap = table( ...
        string.empty(0,1), ...  % SetName
        string.empty(0,1), ...  % Event
        string.empty(0,1), ...  % Monkey
        string.empty(0,1), ...  % Condition
        double.empty(0,1), ...  % Session
        string.empty(0,1), ...  % ModelName
        cell(0,1), ...          % ThetaHat
        cell(0,1), ...          % ThetaBoot
        double.empty(0,1), ...  % NumBoot
        'VariableNames', {'SetName','Event','Monkey','Condition','Session', ...
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
            groupVars = {'Monkey','Condition','Session'};
        case 'pooled'
            groupVars = {'Condition','Session'};
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
        sessionID   = mcKey.Session(g);

        rowMaskMC  = (mcGroup == g);
        sessDataMC = sortrows(trialTable(rowMaskMC,:), {'Session','BlockID','TrialInBlock'});

        for m = 1:numel(modelsForBootstrap)
            modelName = modelsForBootstrap{m};

            fitMask = strcmp(behaviorFitSummary.Monkey,   monkeyID) & ...
                      strcmp(behaviorFitSummary.Condition, conditionID) & ...
                      behaviorFitSummary.Session == sessionID & ...
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
                    lb = [0.0, 0.0, 0, -5.0]; ub = [1.0, 20.0, 5.0, 5.0];
                case 'SelfBubble'
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_SelfBubble(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_SelfBubble(th, sd, BUY_NON_BUY_ONLY);
                    lb = [0.0, 0.0, -5.0, -5.0, -5.0, -5.0]; ub = [1.0, 20.0, 5.0, 5.0, 5.0, 5.0];
                case 'SelfOpp'
                    objFunMaker = @(sd) @(th) local_rlNegLogLik_SelfOpp(th, sd, BUY_NON_BUY_ONLY);
                    simFun      = @(th,sd) local_simulateChoices_SelfOpp(th, sd, BUY_NON_BUY_ONLY);
                    lb = [0.0, 0.0, 0.0, -5.0, -5.0, -5.0, -5.0]; ub = [1.0, 1.0, 20.0, 5.0, 5.0, 5.0, 5.0];
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
                sessionID, ...
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

            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = Rc + omega * Rp;

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

            Rc = sessData.choiceOutcomeReward(t);
            Rp = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp), Rp = 0; end
            r = Rc + omega * Rp;

            Vt = Q(aSim);
            Q(aSim) = Q(aSim) + alpha * (r - Vt);
        end
    end
end

function simChoices = local_simulateChoices_SelfOpp(theta, sessData, BUY_NON_BUY_ONLY)
    alpha_self = max(min(theta(1),1),0);
    alpha_opp  = max(min(theta(2),1),0);
    beta       = theta(3);
    omega_self = theta(4);
    omega_opp  = theta(5);
    gamma      = theta(6);
    kappaBuy   = theta(7);

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

    oppChoiceStr = string(sessData.opponentChoice);
    oppValid = ~ismissing(oppChoiceStr);

    simChoices = strings(T,1); simChoices(:) = missing;
    blockID = sessData.BlockID; uBlocks = unique(blockID,'stable');

    for bb = 1:numel(uBlocks)
        bID = uBlocks(bb);
        idxBlock = find(blockID == bID);
        blockData = sessData(idxBlock,:);
        [~, sortIdx] = sort(blockData.TrialInBlock);
        ordIdx = idxBlock(sortIdx);

        Q_self = zeros(nA,1);
        Q_opp  = zeros(nA,1);

        for k = 1:numel(ordIdx)
            t = ordIdx(k);
            if ~valid(t), continue; end

            V = Q_self + gamma * Q_opp;
            logits = beta * V; logits(buyIdx) = logits(buyIdx) + kappaBuy;
            maxLogit = max(logits);
            p = exp(logits - maxLogit); p = p / sum(p);
            u = rand; cum = cumsum(p);
            aSim = find(u <= cum, 1); if isempty(aSim), aSim = buyIdx; end
            simChoices(t) = string(actionLabels{aSim});

            Rc = sessData.choiceOutcomeReward(t);
            Rp_self = sessData.portfolioRewardFromDividend(t);
            if isnan(Rc), Rc = 0; end
            if isnan(Rp_self), Rp_self = 0; end
            r_self = Rc + omega_self * Rp_self;

            Vt_self = Q_self(aSim);
            Q_self(aSim) = Q_self(aSim) + alpha_self * (r_self - Vt_self);

            if oppValid(t)
                oppStr = oppChoiceStr(t);
                aOpp = 0;
                for a = 1:nA
                    if oppStr == actionLabels{a}, aOpp = a; break; end
                end
                if aOpp > 0
                    Rp_opp = sessData.opponentPortfolioValue(t);
                    if isnan(Rp_opp), Rp_opp = 0; end
                    r_opp = omega_opp * Rp_opp;

                    Vt_opp = Q_opp(aOpp);
                    Q_opp(aOpp) = Q_opp(aOpp) + alpha_opp * (r_opp - Vt_opp);
                end
            end
        end
    end
end

function paramTable = local_behaviorParamTable(fitTable, paramNamesPerModel)
    if isempty(fitTable)
        paramTable = table();
        return;
    end

    rows = [];

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
                fitTable.Session(i), ...
                fitTable.ModelName(i), ...
                uint32(p), ...
                string(names{p}), ...
                theta(p), ...
                'VariableNames', {'SetName','Event','Monkey','Condition','Session', ...
                                  'ModelName','ParamIndex','ParamName','Value'});
            rows = [rows; newRow]; %#ok<AGROW>
        end
    end

    paramTable = rows;
end

function local_plotBehaviorParamComparisonBySet(fitTable, paramNamesPerModel)
    if isempty(fitTable), return; end

    sets   = unique(fitTable.SetName, 'stable');
    events = unique(fitTable.Event,   'stable');
    if isempty(events)
        return;
    end
    eventUse = events(1);  % single-event pipeline; adjust if multi-event

    useBoxchart = exist('boxchart','file') == 2;
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
        yLimsByModel.(fld) = [yMin; yMax];  % 2 x P
    end

    % ===== Detect OT and Saline sets (by name substring) =====
    otIdx  = find(contains(sets, "OT",      'IgnoreCase', true), 1, 'first');
    salIdx = find(contains(sets, "Saline",  'IgnoreCase', true), 1, 'first');
    pairPossible = ~isempty(otIdx) && ~isempty(salIdx);

    % =====================================================================
    % If we have both OT and Saline, make paired figures (OT vs Saline).
    % Otherwise, fall back to per-set plotting (original-style).
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

            figName = sprintf('RL Params (OT vs Saline) - Model %s - Event %s', modelNameChar, eventUse);
            f = figure('Name', figName, 'Color','w');
            try
                t = tiledlayout(f, P, 2, 'TileSpacing','compact', 'Padding','compact');
            catch
                % older MATLAB without tiledlayout
                t = [];
            end

            for p = 1:P
                % ---------- Column 1: OT ----------
                subT = subOT;
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
                    leftover = tmpLevels(tmpLevels ~= "");
                    orderedConds = [ordered; leftover];

                    condCats = categorical(string(subT.Condition), orderedConds, 'Ordinal', true);
                    condIdx  = double(condCats);
                    vals     = cellfun(@(theta) theta(p), subT.Theta);
                    K        = numel(orderedConds);

                    if ~isempty(t)
                        ax1 = nexttile(t, (p-1)*2 + 1);
                    else
                        subplot(P,2,(p-1)*2+1);
                        ax1 = gca;
                    end
                    hold(ax1, 'on');

                    if useBoxchart
                        boxchart(ax1, condCats, vals);
                    else
                        boxplot(ax1, vals, condCats);
                    end

                    jitter = (rand(size(vals)) - 0.5) * 0.3;
                    scatter(ax1, condIdx + jitter, vals, 20, 'filled', ...
                            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);

                    yline(ax1, 0, 'k--');
                    if p == P
                        xlabel(ax1, 'Condition');
                    else
                        xlabel(ax1, '');
                    end
                    ylabel(ax1, paramNames{p}, 'Interpreter','none');
                    title(ax1, sprintf('%s', string(setOT)), 'Interpreter','none');

                    if K >= 2
                        [sigPairs, pvals] = local_computeSig(vals, condIdx, K);
                        local_plotSignificanceLines(ax1, sigPairs, pvals);
                    end

                    if useGlobalY
                        yLo = yMat(1,p); yHi = yMat(2,p);
                        if isfinite(yLo) && isfinite(yHi) && (yHi > yLo)
                            ylim(ax1, [yLo yHi]);
                        end
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
                    leftover = tmpLevels(tmpLevels ~= "");
                    orderedConds = [ordered; leftover];

                    condCats = categorical(string(subT.Condition), orderedConds, 'Ordinal', true);
                    condIdx  = double(condCats);
                    vals     = cellfun(@(theta) theta(p), subT.Theta);
                    K        = numel(orderedConds);

                    if ~isempty(t)
                        ax2 = nexttile(t, (p-1)*2 + 2);
                    else
                        subplot(P,2,(p-1)*2+2);
                        ax2 = gca;
                    end
                    hold(ax2, 'on');

                    if useBoxchart
                        boxchart(ax2, condCats, vals);
                    else
                        boxplot(ax2, vals, condCats);
                    end

                    jitter = (rand(size(vals)) - 0.5) * 0.3;
                    scatter(ax2, condIdx + jitter, vals, 20, 'filled', ...
                            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);

                    yline(ax2, 0, 'k--');
                    if p == P
                        xlabel(ax2, 'Condition');
                    else
                        xlabel(ax2, '');
                    end
                    ylabel(ax2, '');
                    title(ax2, sprintf('%s', string(setSal)), 'Interpreter','none');

                    if K >= 2
                        [sigPairs, pvals] = local_computeSig(vals, condIdx, K);
                        local_plotSignificanceLines(ax2, sigPairs, pvals);
                    end

                    if useGlobalY
                        yLo = yMat(1,p); yHi = yMat(2,p);
                        if isfinite(yLo) && isfinite(yHi) && (yHi > yLo)
                            ylim(ax2, [yLo yHi]);
                        end
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
                sgtitle(t, sprintf('OT vs Saline | Model: %s | Event: %s', modelNameChar, eventUse), ...
                        'Interpreter','none');
            else
                suptitle(sprintf('OT vs Saline | Model: %s | Event: %s', modelNameChar, eventUse));
            end
        end

        return;  % do not fall back to per-set plotting when OT/Saline pair exists
    end

    % =====================================================================
    % Fallback: original per-set plotting (when OT/Saline pair not present)
    % =====================================================================

    for si = 1:numel(sets)
        setName = sets(si);
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
            leftover = tmpLevels(tmpLevels ~= "");
            orderedConds = [ordered; leftover];

            condCatsAll = categorical(string(subT.Condition), orderedConds, 'Ordinal', true);
            condIdxAll  = double(condCatsAll);
            K = numel(orderedConds);

            fld = matlab.lang.makeValidName(modelNameChar);
            useGlobalY = isfield(yLimsByModel, fld);
            if useGlobalY
                yMat = yLimsByModel.(fld);
            else
                yMat = nan(2,P);
            end

            figure('Name', sprintf('RL Parameters by Condition - %s - %s', setName, modelNameChar), 'Color','w');
            nCols = min(3,P);
            nRows = ceil(P/nCols);

            for p = 1:P
                subplot(nRows,nCols,p); hold on;

                vals = cellfun(@(theta) theta(p), subT.Theta);
                condCats = condCatsAll;
                condIdx  = condIdxAll;

                if useBoxchart
                    boxchart(condCats, vals);
                else
                    boxplot(vals, condCats);
                end

                jitter = (rand(size(vals)) - 0.5) * 0.3;
                scatter(condIdx + jitter, vals, 20, 'filled', ...
                        'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);

                yline(0,'k--');
                title(paramNames{p}, 'Interpreter','none');
                xlabel('Condition');
                ylabel('Param value');

                if K >= 2
                    [sigPairs, pvals] = local_computeSig(vals, condIdx, K);
                    local_plotSignificanceLines(gca, sigPairs, pvals);
                end

                if useGlobalY
                    yLo = yMat(1,p); yHi = yMat(2,p);
                    if isfinite(yLo) && isfinite(yHi) && (yHi > yLo)
                        ylim([yLo yHi]);
                    end
                end

                hold off;
            end

            sgtitle(sprintf('Set: %s | Model: %s | Event: %s', setName, modelNameChar, eventUse), 'Interpreter','none');
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
        text(mean([x1 x2]), y + tick, starText, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');
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
                desired = regexprep(desired, '\s+', '');
                present = desired(ismember(desired, allConds));
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

function clusterStruct = local_clusterNeuronCoefficientsByGroup_optimal(neuronGLMResults, K_max)
    clusterStruct = struct('groupKey',{},'coefLabels',{},'coefMatrix',{}, ...
                           'coefMatrixZ',{},'clusterIdx',{},'bestK',{}, ...
                           'K_opt',{},'neuronIndices',{},'outlierNeuronIdx',{}, ...
                           'qualityMetrics',{});

    if isempty(neuronGLMResults) || K_max < 2
        return;
    end

    MIN_NEURONS_PER_GROUP   = 10;
    MIN_NEURONS_PER_CLUSTER = 5;
    Z_OUTLIER_THRESH        = 3.5;

    allCoefLabels = neuronGLMResults(1).CoefNames(:);
    nCoefs        = numel(allCoefLabels);

    isIntercept = strcmp(allCoefLabels, '(Intercept)');
    isCondition = startsWith(allCoefLabels, 'Condition');
    useMask     = ~(isIntercept | isCondition);

    featLabels  = allCoefLabels(useMask);
    nFeat       = numel(featLabels);
    if nFeat == 0
        fprintf('No condition-invariant coefficients found; skipping clustering.\n');
        return;
    end

    nNeurons = numel(neuronGLMResults);
    SetName   = strings(nNeurons,1);
    Event     = strings(nNeurons,1);
    Monkey    = strings(nNeurons,1);
    Condition = strings(nNeurons,1);

    for i = 1:nNeurons
        SetName(i)   = string(neuronGLMResults(i).SetName);
        Event(i)     = string(neuronGLMResults(i).Event);
        Monkey(i)    = string(neuronGLMResults(i).Monkey);
        Condition(i) = string(neuronGLMResults(i).Condition);
    end

    T = table(SetName, Event, Monkey, Condition);
    [grpID, grpKeyTbl] = findgroups(T);
    nGroups = height(grpKeyTbl);

    gCounter = 0;

    for g = 1:nGroups
        idxNeurons = find(grpID == g);
        nInGroup   = numel(idxNeurons);

        if nInGroup < MIN_NEURONS_PER_GROUP
            fprintf('  [CLUST] Skip group %s | %s | %s | %s: n=%d < %d\n', ...
                grpKeyTbl.SetName(g), grpKeyTbl.Event(g), ...
                grpKeyTbl.Monkey(g), grpKeyTbl.Condition(g), ...
                nInGroup, MIN_NEURONS_PER_GROUP);
            continue;
        end

        coefMatrix = nan(nInGroup, nFeat);
        for j = 1:nInGroup
            vals = neuronGLMResults(idxNeurons(j)).CoefValues;
            vals = vals(:).';
            if numel(vals) == nCoefs
                coefMatrix(j,:) = vals(useMask);
            end
        end

        validCoef = all(~isnan(coefMatrix),2);
        if sum(validCoef) < MIN_NEURONS_PER_GROUP
            fprintf('  [CLUST] Skip group %s | %s | %s | %s: NaNs, valid=%d < %d\n', ...
                grpKeyTbl.SetName(g), grpKeyTbl.Event(g), ...
                grpKeyTbl.Monkey(g), grpKeyTbl.Condition(g), ...
                sum(validCoef), MIN_NEURONS_PER_GROUP);
            continue;
        end

        coefMatrix = coefMatrix(validCoef,:);
        idxNeurons = idxNeurons(validCoef);
        nInGroup   = numel(idxNeurons);

        coefZ = zscore(coefMatrix, 0, 1);

        outlierMask = any(abs(coefZ) > Z_OUTLIER_THRESH, 2);

        nOutliers = sum(outlierMask);
        nKeep     = nInGroup - nOutliers;

        if nKeep < MIN_NEURONS_PER_GROUP
            fprintf('  [CLUST] Skip group %s | %s | %s | %s: keep=%d < %d after outliers\n', ...
                grpKeyTbl.SetName(g), grpKeyTbl.Event(g), ...
                grpKeyTbl.Monkey(g), grpKeyTbl.Condition(g), ...
                nKeep, MIN_NEURONS_PER_GROUP);
            continue;
        end

        coefKeep      = coefMatrix(~outlierMask,:);
        coefZKeep     = coefZ(~outlierMask,:);
        neuronIdxKeep = idxNeurons(~outlierMask);
        outlierNeuronIdx = idxNeurons(outlierMask);

        K_upper = min(K_max, floor(nKeep / MIN_NEURONS_PER_CLUSTER));
        if K_upper < 2
            fprintf('  [CLUST] Skip group %s | %s | %s | %s: not enough neurons for clustering (keep=%d).\n', ...
                grpKeyTbl.SetName(g), grpKeyTbl.Event(g), ...
                grpKeyTbl.Monkey(g), grpKeyTbl.Condition(g), ...
                nKeep);
            continue;
        end

        K_candidates = 2:K_upper;
        meanSilByK   = nan(size(K_candidates));
        bestIdx      = [];
        bestK        = NaN;
        bestSil      = -Inf;

        rng(0);

        for ki = 1:numel(K_candidates)
            K = K_candidates(ki);

            try
                idxK = kmeans(coefZKeep, K, ...
                              'Replicates', 10, ...
                              'MaxIter',    1000, ...
                              'Display',    'off');
            catch
                continue;
            end

            clusterSizes = accumarray(idxK, 1, [K 1]);
            if any(clusterSizes < MIN_NEURONS_PER_CLUSTER)
                continue;
            end

            try
                s = silhouette(coefZKeep, idxK);
                mSil = mean(s);
            catch
                mSil = NaN;
            end

            meanSilByK(ki) = mSil;

            if ~isnan(mSil) && mSil > bestSil
                bestSil = mSil;
                bestK   = K;
                bestIdx = idxK;
            end
        end

        if isnan(bestK) || isempty(bestIdx)
            fprintf('  [CLUST] Skip group %s | %s | %s | %s: no valid K.\n', ...
                grpKeyTbl.SetName(g), grpKeyTbl.Event(g), ...
                grpKeyTbl.Monkey(g), grpKeyTbl.Condition(g));
            continue;
        end

        finalClusterSizes = accumarray(bestIdx, 1, [bestK 1]);

        gCounter = gCounter + 1;

        clusterStruct(gCounter).groupKey = struct( ...
            'SetName',   grpKeyTbl.SetName(g), ...
            'Event',     grpKeyTbl.Event(g), ...
            'Monkey',    grpKeyTbl.Monkey(g), ...
            'Condition', grpKeyTbl.Condition(g));

        clusterStruct(gCounter).coefLabels       = featLabels;
        clusterStruct(gCounter).coefMatrix       = coefKeep;
        clusterStruct(gCounter).coefMatrixZ      = coefZKeep;
        clusterStruct(gCounter).clusterIdx       = bestIdx;
        clusterStruct(gCounter).bestK            = bestK;
        clusterStruct(gCounter).K_opt            = bestK;
        clusterStruct(gCounter).neuronIndices    = neuronIdxKeep;
        clusterStruct(gCounter).outlierNeuronIdx = outlierNeuronIdx;

        qMetrics = struct();
        qMetrics.K_candidates  = K_candidates;
        qMetrics.meanSilByK    = meanSilByK;
        qMetrics.K_opt         = bestK;
        qMetrics.meanSilOpt    = bestSil;
        qMetrics.clusterSizes  = finalClusterSizes;
        qMetrics.nTotal        = nInGroup;
        qMetrics.nKeep         = nKeep;
        qMetrics.nOutliers     = nOutliers;

        clusterStruct(gCounter).qualityMetrics = qMetrics;

        clusterSizesStr = sprintf('%d ', finalClusterSizes);
        fprintf('  [CLUST] Group %s | %s | %s | %s: total=%d, keep=%d, outliers=%d, K_opt=%d, meanSil=%.3f, sizes=[%s]\n', ...
            grpKeyTbl.SetName(g), grpKeyTbl.Event(g), ...
            grpKeyTbl.Monkey(g), grpKeyTbl.Condition(g), ...
            qMetrics.nTotal, qMetrics.nKeep, qMetrics.nOutliers, ...
            qMetrics.K_opt, qMetrics.meanSilOpt, strtrim(clusterSizesStr));
    end
end

function local_plotNeuronClusters_condInvariant_optimal(clusterStruct)
    if isempty(clusterStruct)
        return;
    end

    yMin = inf;
    yMax = -inf;

    for g = 1:numel(clusterStruct)
        Z   = clusterStruct(g).coefMatrixZ;
        idx = clusterStruct(g).clusterIdx;

        if isempty(Z) || isempty(idx)
            continue;
        end

        K_here = max(idx);
        for k = 1:K_here
            inCluster = (idx == k);
            if ~any(inCluster)
                continue;
            end
            m  = mean(Z(inCluster,:), 1);
            se = std(Z(inCluster,:), 0, 1) ./ sqrt(sum(inCluster));
            yMin = min(yMin, min(m - se));
            yMax = max(yMax, max(m + se));
        end
    end

    if ~isfinite(yMin) || ~isfinite(yMax)
        warning('No valid coefficients found for plotting clusters.');
        return;
    end

    margin = 0.1 * max(abs([yMin yMax]));
    yMin = yMin - margin;
    yMax = yMax + margin;

    for g = 1:numel(clusterStruct)
        coefLabels = clusterStruct(g).coefLabels;
        Z          = clusterStruct(g).coefMatrixZ;
        idx        = clusterStruct(g).clusterIdx;

        if isempty(Z) || isempty(idx)
            continue;
        end

        nFeat  = numel(coefLabels);
        K_here = max(idx);

        setName   = string(clusterStruct(g).groupKey.SetName);
        eventName = string(clusterStruct(g).groupKey.Event);
        monkey    = string(clusterStruct(g).groupKey.Monkey);
        cond      = string(clusterStruct(g).groupKey.Condition);
        K_opt     = clusterStruct(g).K_opt;

        if monkey ~= "1" || eventName ~= "m1rew2on"
            continue;
        end        

        figName = sprintf('Neuron clusters | %s | %s | %s | %s', ...
            setName, eventName, monkey, cond);

        figure('Name', figName, 'Color','w');

        for k = 1:K_here
            subplot(K_here,1,k); hold on;

            inCluster = (idx == k);
            if ~any(inCluster)
                continue;
            end

            m  = mean(Z(inCluster,:), 1);
            se = std(Z(inCluster,:), 0, 1) ./ sqrt(sum(inCluster));

            x = 1:nFeat;
            bar(x, m);
            errorbar(x, m, se, 'k', 'LineStyle','none', 'LineWidth', 1.0);

            set(gca, 'XTick', x, ...
                     'XTickLabel', coefLabels, ...
                     'XTickLabelRotation', 45);
            ylabel(sprintf('Cluster %d (n = %d)', k, sum(inCluster)));

            xline(0,'k-');
            yline(0,'k--');

            ylim([yMin yMax]);
            box on; grid on;
        end

        sgtitle(sprintf('%s | %s\n%s | %s\n(K_{opt} = %d)', ...
            setName, eventName, monkey, cond, K_opt));
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
