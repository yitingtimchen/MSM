% File: d20251123_fit_glme_with_neural_features.m

%% Buy Bias Analysis: Live vs Other Contexts
% This script runs a GLME on neural firing rates and behavior, comparing
% Live vs other task contexts. It has been updated to:
%   (1) Preserve condition names like "AI", "OT_AI", "Saline_AI", etc.
%       by normalizing struct field names with underscores instead of
%       collapsing spaces (no more "OTAI", "OTReplay", ...).
%   (2) Shift history-dependent behavioral variables (previous choice,
%       opponent choice, previous reward, opponent portfolio) within
%       each 15-trial block rather than across the full 90-trial session.
%       This yields 6 NaNs per session for those "previous" variables
%       (one at the start of each 15-trial block).

clear; clc; close all;

% ------------------ PATHS & GLOBAL OPTIONS ------------------

% Path to condition packs
% OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs';
% OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only';
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_251118';

% If set to 1, collapse all non-Live conditions into a single "not_Live" level
LIVE_NON_LIVE_ONLY = 0;

% ------------------ NEURAL ANALYSIS SETTINGS ------------------

% Target event(s) and windows for neural analysis (ms)
TARGET_EVTS     = {'m1rew2on'};
TARGET_EVT_WINS = {[-250, 250]};   % window relative to TARGET_EVT

% Example for multi-event analysis (kept here for reference)
% TARGET_EVTS = {                                                   % Each cell is a string array of event names in temporal order per trial
%     "m1rew1on","m1rew1off","m1rew2on","m1rew2off",...
%     "m2rew1on","m2rew1off","m2rew2on","m2rew2off"
%     };
% TARGET_EVT_WINS = {
%       [-250 250];   % m1rew1on
%       [-250  250];  % m1rew1off
%       [-250  250];  % m1rew2on
%       [-250  250];  % m1rew2off/m1trialStop
%       [-250  250];  % m2rew1on
%       [-250  250];  % m2rew1off
%       [-250  250];  % m2rew2on
%       [-250  250];  % m2rew2off/m2trialStop
%     };

% Baseline window for z-scoring firing rate (ms)
BASELINE_EVT = 'm1rew1on';
BASELINE_WIN = [-1000, -700];

% ------------------ CONDITION SETS ------------------
% Define condition names exactly how you want them to appear as struct
% fieldnames and logical labels. Any spaces in 'conds' are converted to
% underscores for struct fieldnames (e.g. "OT AI" -> "OT_AI"), so you can
% keep human-readable names on disk if you want.
%
% Example filename mapping:
%   "OT AI"      -> loads "OT AI_condition_pack.mat", stored as Cstruct.OT_AI
%   "Saline AI"  -> loads "Saline AI_condition_pack.mat", stored as Cstruct.Saline_AI

cond_sets = { ...
    % {'AI','Replay','Decoy','Live'}; ...         % Baseline set (optional)
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames = {'OT','Saline'};
% If/when you want baseline included in the loop, uncomment the first line
% in cond_sets above and extend setNames to {'Baseline','OT','Saline'}.

allResults = struct();

% ------------------ MAIN LOOPS: EVENT × CONDITION SET ------------------

for pair = 1:numel(TARGET_EVTS)
    TARGET_EVT     = TARGET_EVTS{pair};
    TARGET_EVT_WIN = TARGET_EVT_WINS{pair};

    for s = 1:numel(cond_sets)
        conds = cond_sets{s};
        fprintf('\n=== Analyzing Set: %s ===\n', setNames{s});

        % ---- Load condition packs ----
        % We build a struct Cstruct with one field per condition, using
        % underscores in the fieldnames so that "OT AI" becomes "OT_AI".
        Cstruct = struct();
        for c = 1:numel(conds)
            cname_on_disk = conds{c};  % e.g. "OT AI", "Saline Replay", etc.
            matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname_on_disk));
            if ~exist(matFile,'file')
                warning('Missing file: %s', matFile);
                continue;
            end
            tmp = load(matFile, 'C');

            % Normalize field name: replace spaces with underscores so we get
            % Cstruct.OT_AI, Cstruct.Saline_AI, etc. instead of "OTAI".
            fieldName = strrep(cname_on_disk, ' ', '_');
            Cstruct.(fieldName) = tmp.C;
        end
        condNames = fieldnames(Cstruct);  % e.g. 'OT_AI','OT_Replay',...

        % ---- RUN GLME ----
        PBuy     = struct();  % per-condition P(BUY) summary
        trialData = table();  % pooled trial-level data for GLME

        for c = 1:numel(condNames)
            cname = condNames{c};     % normalized condition name, e.g. "OT_AI"
            C     = Cstruct.(cname);  % condition-specific data struct

            for f = 1:numel(C.eventTables)
                etab = C.eventTables{f};

                % Skip sessions lacking the "option" column or any option data
                if ~ismember('option', etab.Properties.VariableNames)
                    continue;
                end
                if isempty(etab.option{1})
                    continue;
                end

                % ------------------ NEURAL DATA: Z-SCORED FIRING ------------------
                Z_fr = get_Z_fr(C, f, TARGET_EVT, TARGET_EVT_WIN, BASELINE_EVT, BASELINE_WIN);

                % ------------------ BEHAVIOR: SUBJECT CHOICES ------------------
                % subjChoice is flattened (TR × 1) string array: "BUY","HOLD","SELL"
                subjChoice = mapChoicesToLabels(etab.option{1});
                ntr = numel(subjChoice);

                % We assume each session has 6 blocks × 15 trials = 90 trials.
                % All block-wise lagged variables below use this block structure.
                blockSize = 15;
                if mod(ntr, blockSize) ~= 0
                    error('Session %d: nTrials = %d not divisible by block size %d.', ...
                          f, ntr, blockSize);
                end
                nBlocks = ntr / blockSize;

                % Valid trials (in case some options are missing)
                valid = ~ismissing(subjChoice);

                % Probability of BUY per session (ignores missing)
                pBuy = mean(subjChoice(valid) == "BUY");
                PBuy.(cname)(f,1) = pBuy;  %#ok<AGROW>

                % Subject BUY indicator for valid trials
                subjBuy = subjChoice(valid) == "BUY";

                % ------------------ CONDITION LABEL (AI/Replay/Decoy/Live) ------------------
                % 'cname' includes drug prefix if present (e.g. "OT_AI", "Saline_Replay").
                % We treat the last underscore-delimited token as the task context,
                % so all drug variants share the same Condition levels ("AI","Replay",...).
                tokens   = split(string(cname), '_');
                baseCond = tokens(end);  % "AI","Replay","Decoy","Live"

                condLabel = repmat(baseCond, sum(valid), 1);
                if LIVE_NON_LIVE_ONLY
                    condLabel(condLabel ~= "Live") = "not_Live";
                end

                % ------------------ TRIAL-LEVEL COVARIATES ------------------
                monkey = string(etab.monkey{1}(:));
                session = double(etab.session{1}(:));
                trial   = double(etab.trialNum{1}(:));
                market  = string(etab.marketOrig{1}(:));
                bubble  = string(etab.bubbleMarket{1}(:));

                if string(TARGET_EVT) == "m1rew2on"
                    portfolioVec = double(etab.postPortfolio{1}(:));
                else
                    portfolioVec = double(etab.prePortfolio{1}(:));
                end

                % ------------------ OPPONENT CHOICES & PORTFOLIO ------------------
                % Map opponent's raw choice codes to "BUY","HOLD","SELL"
                oppRaw = mapChoicesToLabels(etab.optionOpp{1});

                % Opponent post-portfolio at current trial
                oppPortVec = etab.postPortfolioOpp{1}(:);

                % ------------------ REWARD PER TRIAL ------------------
                % Market reward experienced by subject on each trial
                choiceRewardVec = etab.postPortfolio{1}(:) .* etab.divPerShare{1}(:);

                % ------------------ BLOCK-WISE LAGGED VARIABLES (KEY FIX) ------------------
                % We now compute all "previous" variables within each 15-trial block
                % *before* linearizing, so there are 6 NaNs per session instead of
                % just 1. This avoids bleeding information across block boundaries.

                % Reshape trial-wise vectors into [blockSize × nBlocks] matrices
                subjChoiceMat    = reshape(subjChoice,      blockSize, nBlocks);
                oppChoiceMat     = reshape(oppRaw,          blockSize, nBlocks);
                oppPortMat       = reshape(oppPortVec,      blockSize, nBlocks);
                choiceRewardMat  = reshape(choiceRewardVec, blockSize, nBlocks);

                % Initialize "previous" matrices with missing/NaN
                prevChoiceMat        = strings(size(subjChoiceMat)); prevChoiceMat(:) = missing;
                opponentChoiceMat    = strings(size(oppChoiceMat));  opponentChoiceMat(:) = missing;
                prevChoiceRewardMat  = nan(size(choiceRewardMat));
                opponentPortfolioMat = nan(size(oppPortMat));

                % Shift within each block: row 2:t gets previous trial's value,
                % row 1 remains missing/NaN for each block.
                prevChoiceMat(2:end,:)        = subjChoiceMat(1:end-1,:);
                opponentChoiceMat(2:end,:)    = oppChoiceMat(1:end-1,:);
                prevChoiceRewardMat(2:end,:)  = choiceRewardMat(1:end-1,:);
                opponentPortfolioMat(2:end,:) = oppPortMat(1:end-1,:);

                % Flatten back to column vectors (TR × 1)
                prevChoice        = prevChoiceMat(:);
                opponentChoice    = opponentChoiceMat(:);
                prevChoiceReward  = prevChoiceRewardMat(:);
                opponentPortfolio = opponentPortfolioMat(:);

                % Derived BUY indicators for previous and opponent choices
                prevBuy     = prevChoice == "BUY";
                opponentBuy = opponentChoice == "BUY";

                % Current-trial portfolio and reward vectors
                portfolio    = portfolioVec;
                choiceReward = choiceRewardVec;

                % ------------------ BUILD TRIAL-LEVEL TABLE ------------------
                % All predictors are filtered by 'valid' (subject choice present).
                tmpT = table( ...
                    subjBuy, prevBuy(valid), opponentBuy(valid), ...
                    condLabel, string(monkey(valid)), session(valid), trial(valid), ...
                    string(market(valid)), string(bubble(valid)), ...
                    string(subjChoice(valid)), string(opponentChoice(valid)), string(prevChoice(valid)), ...
                    choiceReward(valid), prevChoiceReward(valid), ...
                    portfolio(valid), opponentPortfolio(valid), ...
                    Z_fr(valid), ...
                    'VariableNames', {'subjBuy', 'prevBuy', 'opponentBuy', ...
                        'Condition','Monkey','Session','Trial','Market','Bubble', ...
                        'subjChoice','opponentChoice','prevChoice', ...
                        'choiceReward', 'prevChoiceReward', ...
                        'portfolio', 'opponentPortfolio', ...
                        'Z_fr'});

                trialData = [trialData; tmpT]; %#ok<AGROW>
            end
        end

        % ------------------ CATEGORICAL ENCODINGS ------------------
        trialData.Monkey = categorical(trialData.Monkey);
        trialData.Market = categorical(trialData.Market);
        trialData.Bubble = categorical(trialData.Bubble);

        % Condition is the context (AI/Replay/Decoy/Live), regardless of drug.
        if LIVE_NON_LIVE_ONLY
            trialData.Condition = setcats(categorical(trialData.Condition), {'not_Live','Live'});        
        else
            trialData.Condition = setcats(categorical(trialData.Condition), {'AI','Replay','Decoy','Live'});
        end
        trialData.subjChoice      = setcats(categorical(trialData.subjChoice),      {'HOLD','BUY','SELL'});
        trialData.opponentChoice  = setcats(categorical(trialData.opponentChoice),  {'HOLD','BUY','SELL'});
        trialData.prevChoice      = setcats(categorical(trialData.prevChoice),      {'HOLD','BUY','SELL'});

        % ------------------ CONDITION-SPECIFIC SUBSETS ------------------
        trialDataAI      = trialData(trialData.Condition == "AI", :);
        trialDataReplay  = trialData(trialData.Condition == "Replay", :);
        trialDataDecoy   = trialData(trialData.Condition == "Decoy", :);
        trialDataLive    = trialData(trialData.Condition == "Live", :);

        trialDataNonLive   = trialData(trialData.Condition ~= "Live", :);
        trialDataSocial    = trialData(trialData.Condition == "Live" | trialData.Condition == "Decoy", :);
        trialDataNonSocial = trialData(trialData.Condition == "AI"   | trialData.Condition == "Replay", :);

        % By Monkey
        trialDataM1allConds = trialData(trialData.Monkey == "1", :);
        trialDataM2allConds = trialData(trialData.Monkey == "2", :);

        % By Monkey and Condition
        trialDataM1AI      = trialDataM1allConds(trialDataM1allConds.Condition == "AI", :);
        trialDataM2AI      = trialDataM2allConds(trialDataM2allConds.Condition == "AI", :);
        trialDataM1Replay  = trialDataM1allConds(trialDataM1allConds.Condition == "Replay", :);
        trialDataM2Replay  = trialDataM2allConds(trialDataM2allConds.Condition == "Replay", :);
        trialDataM1Decoy   = trialDataM1allConds(trialDataM1allConds.Condition == "Decoy", :);
        trialDataM2Decoy   = trialDataM2allConds(trialDataM2allConds.Condition == "Decoy", :);
        trialDataM1Live    = trialDataM1allConds(trialDataM1allConds.Condition == "Live", :);
        trialDataM2Live    = trialDataM2allConds(trialDataM2allConds.Condition == "Live", :);

        fprintf("Total sessions (all conds) = " + size(trialData, 1)/90 + "\n")
        fprintf("Total sessions (Live)      = " + size(trialDataLive, 1)/90 + "\n")
        fprintf("M1 sessions (all conds)    = " + size(trialDataM1allConds, 1)/90 + "\n")
        fprintf("M2 sessions (all conds)    = " + size(trialDataM2allConds, 1)/90 + "\n")
        fprintf("M1 sessions (Live)         = " + size(trialDataM1Live, 1)/90 + "\n")
        fprintf("M2 sessions (Live)         = " + size(trialDataM2Live, 1)/90 + "\n")

        % ------------------ GLME MODEL SPECS ------------------
        base = "1";
        choicePreds    = {'subjChoice','prevChoice','opponentChoice'};
        rewardPreds    = {'choiceReward','prevChoiceReward'};
        portfolioPreds = {'portfolio','opponentPortfolio'};

        re_all    = {'(1|Monkey)','(1|Monkey:Session)'};  % random intercepts
        re_single = {'(1|Monkey:Session)'};
        rs_all    = {'(buy + prevChoice + prevChoiceReward | Monkey)'}; % random slopes by Monkey (if needed)

        vars = {
            % 1) All conds - choice only
            [{'Condition'}, choicePreds, re_all], ...
            % 2) Live - choice only
            [choicePreds, re_all], ...
            % 3) All conds - reward only
            [{'Condition'}, rewardPreds, re_all], ...
            % 4) Live - reward only
            [rewardPreds, re_all], ...
            % 5) All conds - portfolio only
            [{'Condition'}, portfolioPreds, re_all], ...
            % 6) Live - portfolio only
            [portfolioPreds, re_all], ...
            % 7) All conds - all preds
            [{'Condition'}, choicePreds, rewardPreds, portfolioPreds, re_all], ...
            % 8) Live - all preds
            [choicePreds, rewardPreds, portfolioPreds, re_all], ...
             };

        cats = {'Condition','Monkey','Market','buy','opponentChoice','prevChoice'};    
        tasks = {
            % 'All conds - choice only',   trialData,             1, 'all_choice';
            % 'All conds - reward only',   trialData,             3, 'all_rew';
            % 'All conds - portfolio only',trialData,             5, 'all_pf';
            % 'All conds - all preds',     trialData,             7, 'all_comb';
            % 
            % 'AI - choice only',          trialDataAI,           2, 'ai_choice';
            % 'AI - reward only',          trialDataAI,           4, 'ai_rew';
            % 'AI - portfolio only',       trialDataAI,           6, 'ai_pf';
            % 'AI - all preds',            trialDataAI,           8, 'ai_comb';
            % 
            % 'Replay - choice only',      trialDataReplay,       2, 'rp_choice';
            % 'Replay - reward only',      trialDataReplay,       4, 'rp_rew';
            % 'Replay - portfolio only',   trialDataReplay,       6, 'rp_pf';
            % 'Replay - all preds',        trialDataReplay,       8, 'rp_comb';
            % 
            % 'Decoy - choice only',       trialDataDecoy,        2, 'decoy_choice';
            % 'Decoy - reward only',       trialDataDecoy,        4, 'decoy_rew';
            % 'Decoy - portfolio only',    trialDataDecoy,        6, 'decoy_pf';
            % 'Decoy - all preds',         trialDataDecoy,        8, 'decoy_comb';

            'Live - choice only',         trialDataLive,         2, 'live_choice';
            'Live - reward only',         trialDataLive,         4, 'live_rew';
            'Live - portfolio only',      trialDataLive,         6, 'live_pf';
            % 'Live - all preds',          trialDataLive,         8, 'live_comb';

            'Non-live - choice only',     trialDataNonLive,      2, 'nonlive_choice';
            'Non-live - reward only',     trialDataNonLive,      4, 'nonlive_rew';
            'Non-live - portfolio only',  trialDataNonLive,      6, 'nonlive_pf';
            % 'Non-live - all preds',      trialDataNonLive,      8, 'nonlive_comb';

            % 'Social - choice only',      trialDataSocial,       2, 'social_choice';
            % 'Social - reward only',      trialDataSocial,       4, 'social_rew';
            % 'Social - portfolio only',   trialDataSocial,       6, 'social_pf';
            % 'Social - all preds',        trialDataSocial,       8, 'social_comb';
            % 
            % 'Non-Social - choice only',  trialDataNonSocial,    2, 'nonSocial_choice';
            % 'Non-Social - reward only',  trialDataNonSocial,    4, 'nonSocial_rew';
            % 'Non-Social - portfolio only',trialDataNonSocial,   6, 'nonSocial_pf';
            % 'Non-Social - all preds',    trialDataNonSocial,    8, 'nonSocial_comb';

        };

        supcond = strsplit(conds{1},' '); supcond = char(supcond{1});
        mdl = struct();
        for i = 1:size(tasks,1)
            label   = tasks{i,1};
            data    = tasks{i,2};
            vIdx    = tasks{i,3};
            key     = tasks{i,4};
            fprintf("\n\n%s:\n", label);
            mdl.(key) = runGLMVariants(data, base, vars(vIdx), cats, [], ...
                [supcond ' ' label ' (' num2str(size(data, 1)/90) ' sessions), event = ' TARGET_EVT ' [' num2str(TARGET_EVT_WIN) '] ms']);
        end
    end

end

%% ================== HELPERS ==================

function visualizeLinearGLM(mdl, trialData, keyCatVar, plotName)
% LINEAR GLM(ME) visualization — coefficient forest (95% CI) in z-units
SORT_BY_EFFECT_SIZE = 0;
ADDITIONAL_PLOTS    = 0; % set to 1 to enable partial dependence & grouped means CIs

if nargin < 3 || isempty(keyCatVar)
    preds = mdl.PredictorNames;
    % Prefer 'Condition' if present
    if any(strcmpi(preds,'Condition'))
        keyCatVar = 'Condition';
    else
        keyCatVar = preds{1};
    end
end
if nargin < 4 || isempty(plotName), plotName = 'Linear Mixed Model Summary'; end

coefTbl = mdl.Coefficients;
cnames  = string(mdl.CoefficientNames);
vars = coefTbl.Properties.VarNames; % cell array of char
idx = find(strcmpi(vars,'pValue') | strcmpi(vars,'PValue'), 1, 'first');
if isempty(idx)
    pvals = NaN(height(coefTbl),1);
else
    pvals = coefTbl.(vars{idx});
end

% === Coefficients & 95% CI in β (z-units) ===
b    = coefTbl.Estimate;
ci   = coefCI(mdl);              % [low, high] in β
mask = cnames ~= "(Intercept)";

names = cnames(mask);
bet   = b(mask);
betci = ci(mask,:);
pv    = pvals(mask);

if SORT_BY_EFFECT_SIZE
    [~, idx] = sort(bet);
    names = names(idx); bet = bet(idx); betci = betci(idx,:); pv = pv(idx);
end

fh = figure('Name', char(plotName(1)));
t = tiledlayout(fh, 1, 1, 'TileSpacing','compact', 'Padding','compact');
[line1, line2] = splitAtMiddleSpace(plotName(2));
sgtitle(t, [plotName(1); line1; line2], 'FontSize', 12, 'interpreter', 'none');

% Forest plot on linear scale
ax1 = nexttile; hold(ax1,'on');
x = 1:numel(bet);
errorbar(ax1, x, bet, bet - betci(:,1), betci(:,2) - bet, 'o', 'LineWidth',1.2, 'CapSize',8);
yline(ax1, 0, 'k--');
xlim(ax1, [0, numel(bet)+1])
set(ax1, 'XTick', x, 'XTickLabel', cellstr(names), 'TickLabelInterpreter','none');
ylabel(ax1, 'Effect on z-scored firing rate (β)');
xlabel(ax1, 'Predictor (level vs reference)');
title(ax1, 'Linear mixed model effects (95% CI)'); box(ax1,'on');

% Sig stars
for i = 1:numel(bet)
    si = string(pstars(pv(i)));
    ytxt = max(betci(i,:)) + 0.05*range(betci(:));  % small offset
    if si ~= "ns"
        text(ax1, x(i), ytxt, [' ' si], 'VerticalAlignment','bottom', ...
             'HorizontalAlignment','center','FontWeight','bold','FontSize',16,'Color','r');
    else
        text(ax1, x(i), ytxt, [' ' si], 'VerticalAlignment','bottom', ...
             'HorizontalAlignment','center','FontSize',11,'Color','k');
    end
end
hold(ax1,'off');

if ~ADDITIONAL_PLOTS
    % Save and return
    figdir = fullfile(pwd, 'figs/neur_glm');
    if ~exist(figdir,'dir'), mkdir(figdir); end
    set(fh,'Color','w');
    exportgraphics(fh, fullfile(figdir, [char(plotName(1)) '.svg']), 'Resolution', 300);
    savefig(fh, fullfile(figdir, [char(plotName(1)), '.fig']));
    return
end

% ===== Optional panels (population-level predictions) =====
predVars = mdl.PredictorNames;
baseRow  = trialData(1, predVars);
for v = 1:numel(predVars)
    vn = predVars{v};
    x  = trialData.(vn);
    if iscategorical(x)
        x   = removecats(x);
        cats= categories(x); cnt = countcats(x);
        [~,iMax]   = max(cnt);
        baseRow.(vn) = categorical(cats(iMax), cats);
    else
        baseRow.(vn) = median(double(x), 'omitnan');
    end
end

% Population-average predictions (ignore BLUPs)
phat_all = predict(mdl, trialData, 'Conditional', false);

% Group summaries
vn = char(keyCatVar);
gvar = trialData.(vn);
if ~iscategorical(gvar), gvar = categorical(discretize(double(gvar), 4)); end
[G, grpCats] = findgroups(gvar);

obsMean = splitapply(@(y)mean(double(y),'omitnan'), trialData.Z_fr, G);
obsSE   = splitapply(@(y)std(double(y),'omitnan')/sqrt(sum(~isnan(y))), trialData.Z_fr, G);
predMean= splitapply(@(y)mean(double(y),'omitnan'), phat_all, G);

% Panel: Observed vs Predicted by group
ax2 = nexttile; hold(ax2,'on');
Xcats = categories(grpCats);
M = [obsMean, predMean];
bh = bar(ax2, categorical(Xcats, Xcats), M, 'grouped');
[xtObs, ytObs] = deal(bh(1).XEndPoints(:), bh(1).YData(:));
errorbar(ax2, xtObs, ytObs, obsSE, obsSE, ...
    'LineStyle','none', 'Color','k', 'LineWidth',1.2, 'CapSize',8);
ylabel(ax2, 'z-scored firing rate'); 
title(ax2, ["Observed vs Predicted"; "grouped by " + vn]); 
legend(ax2, {'Observed','Predicted'}, 'Location','northeast'); 
set(ax2, 'TickLabelInterpreter','none'); box(ax2,'on'); hold(ax2,'off');

% Save
figdir = fullfile(pwd, 'figs/neur_glm');
if ~exist(figdir,'dir'), mkdir(figdir); end
set(fh,'Color','w');
exportgraphics(fh, fullfile(figdir, [char(plotName(1)) '.svg']), 'Resolution', 300);
savefig(fh, fullfile(figdir, [char(plotName(1)), '.fig']));

end

function s = pstars(p)
    if isnan(p), s = ''; return; end
    if p < 1e-4, s = '****';
    elseif p < 1e-3, s = '***';
    elseif p < 0.01, s = '**';
    elseif p < 0.05, s = '*';
    else, s = 'ns';
    end
end

function s = pstr(p)
    if isnan(p), s = 'NA';
    elseif p < 1e-4, s = '<1e-4';
    else, s = sprintf('%.3g', p);
    end
end

function pv = coeffPForTerm(mdl, termName)
    coefTbl = mdl.Coefficients;
    rows = string(mdl.CoefficientNames);
    if ismember('pValue', coefTbl.Properties.VarNames)
        pvals = coefTbl.pValue;
    elseif ismember('PValue', coefTbl.Properties.VarNames)
        pvals = coefTbl.PValue;
    else
        pv = NaN; return
    end
    idx = find(rows == string(termName), 1);
    if isempty(idx)
        % try contains (robust fallback)
        idx = find(strcmpi(rows, string(termName)) | (startsWith(rows, string(termName)) & ~contains(rows, "_")), 1);
    end
    if isempty(idx)
        pv = NaN;
    else
        pv = pvals(idx);
    end
end

% Replace your coeffPForLevel(...) with this version (underscore fix)
function pv = coeffPForLevel(mdl, varName, level)
    coefTbl = mdl.Coefficients;
    rows = string(mdl.CoefficientNames);
    if ismember('pValue', coefTbl.Properties.VarNames)
        pvals = coefTbl.pValue;
    elseif ismember('PValue', coefTbl.Properties.VarNames)
        pvals = coefTbl.PValue;
    else
        pv = NaN; return
    end

    vn  = string(varName);
    lvl = string(level);
    lvlSan = string(matlab.lang.makeValidName(char(lvl)));

    % EXACT match: Var_Level  (FIX: include underscore)
    idx = find(rows == (vn + "_" + lvlSan), 1);

    % Fallback: any Var_* candidate that matches this level (exclude interactions)
    if isempty(idx)
        cands = find(startsWith(rows, vn + "_") & ~contains(rows, ":"));
        idx = cands(find(contains(rows(cands), "_" + lvlSan), 1));
    end

    if isempty(idx)
        pv = NaN;
    else
        pv = pvals(idx);
    end
end

function [lo, hi] = wilsonCI(phat, n)
    % 95% Wilson interval for binomial proportion
    z = 1.96;
    A = phat + z.^2./(2*n);
    B = z .* sqrt((phat.*(1-phat)./n) + (z.^2)./(4*n.^2));
    den = 1 + z.^2./n;
    lo = (A - B) ./ den;
    hi = (A + B) ./ den;
    % handle n=0
    lo(n==0) = NaN; hi(n==0) = NaN;
end

function p = termPValue(mdl, varName)
    % Returns a p-value for a named term in a GeneralizedLinearMixedModel.
    % - Fixed effects: uses ANOVA (Wald/Chi^2), falls back to joint coefTest.
    % - Random effects: uses likelihood-ratio test via model comparison after
    %   removing the specified random term.

    vn = string(varName);
    p  = NaN;

    % ===== 1) Try ANOVA for FIXED effects =====
    try
        A = anova(mdl,'summary');  % tests fixed effects only
        pcol = pickcol(A, {'pValue','PValue','pvalue','p'});
        tcol = pickcol(A, {'Term','Source','Name'});
        if pcol ~= "" && tcol ~= ""
            tt = string(A.(tcol));
            row = find(tt==vn | contains(tt, vn), 1);
            if ~isempty(row)
                p = A.(pcol)(row);
                if ~isnan(p), return; end
            end
        end
    catch
        % proceed to next strategy
    end

    % ===== 2) Fallback: joint Wald test of FIXED-effect coefficients =====
    try
        cn = string(mdl.CoefficientNames(:));
        % Build mask for all coefficients belonging to varName
        if contains(vn, ":") % interaction: require all parts
            parts = split(vn, ":");
            mask = true(size(cn));
            for i = 1:numel(parts), mask = mask & contains(cn, parts{i}); end
        else % main effect (including categorical expansions)
            mask = (cn == vn) | startsWith(cn, vn + "_");
            mask = mask & ~contains(cn, ":");
        end
        mask = mask & cn ~= "(Intercept)";

        idx = find(mask);
        if ~isempty(idx)
            R = zeros(numel(idx), numel(cn));
            for r = 1:numel(idx), R(r, idx(r)) = 1; end
            p = coefTest(mdl, R);
            return 
        end
    catch
        % give up gracefully
    end
end

% ---------- helpers ----------
function [isRE, reMatch] = isRandomGrouping(mdl, vn)
% Returns true if VN matches any grouping factor on the RHS of a random term
% in the GLME formula, i.e., '( ... | GROUP )' or '( ... || GROUP )'.

% Normalize inputs
if ~ischar(vn) && ~isstring(vn), error('vn must be char or string'); end
vn = char(strtrim(vn));
vn_nospace = regexprep(vn, '\s+', '');

% Extract RHS of the model formula
fstr = char(mdl.Formula);
tildeIdx = strfind(fstr, '~');
if isempty(tildeIdx), error('Unexpected formula format.'); end
rhs  = strtrim(fstr(tildeIdx+1:end));

% Find random-effect terms and their grouping sides
matches = regexp(rhs, '\([^()]*\|{1,2}[^()]*\)', 'match');
toks    = regexp(rhs, '\(([^()|]*)\|{1,2}([^()]*)\)', 'tokens');

isRE = false;
reMatch = '';

for k = 1:numel(toks)
    grp = strtrim(toks{k}{2});
    grp_nospace = regexprep(grp, '\s+', '');
    if strcmp(grp_nospace, vn_nospace)
        isRE = true;
        reMatch = matches{k}; % exact '( ... | GROUP )' text
        return
    end
end

end

function p = trygetprop(mdl, prop, defaultVal)
% Safe getter for model properties; returns default if missing.
p = defaultVal;
try
    if isprop(mdl, prop)
        p = mdl.(prop);
    end
catch
    % ignore
end
end

function name = pickcol(T, candidates)
% Return the first matching column name in table T from candidates (case-insensitive).
name = "";
tn = string(T.Properties.VariableNames);
for i = 1:numel(candidates)
    c = string(candidates{i});
    hit = tn == c | lower(tn) == lower(c);
    if any(hit)
        name = tn(find(hit,1));
        return;
    end
end
end

function p = localCompareRemoveTerms(mdl, term)
    try
        mdl_red = removeTerms(mdl, term);
        T = compare(mdl_red, mdl, CheckNesting=true);
        pcol = intersect(string(T.Properties.VariableNames), ["pValue","PValue","pValueLR","pValue_F","pValue_CHI2"]);
        if isempty(pcol)
            p = NaN;
        else
            p = T.(pcol(1))(end);
        end
    catch
        p = NaN;
    end
end

function models = runGLMVariants(trialData, basePredictors, variants, categoricalVars, keyCatVar, plotName)
% runGLMVariants - Fit multiple GLMs with different predictor subsets and visualize each.
%
% Inputs:
%   trialData       - table with predictors + outcome 'Z_fr'
%   basePredictors  - predictors always included
%   variants        - cell array of predictor sets per model
%   categoricalVars - list of categorical predictor names
%   keyCatVar       - key categorical var for plotting (default: auto)
%   plotName        - base title string
%
% Output:
%   models          - struct with fields mdl1, mdl2, ... each a GLME

    if nargin < 4 || isempty(categoricalVars), categoricalVars = {}; end
    if nargin < 5, keyCatVar = []; end
    if nargin < 6 || isempty(plotName), plotName = 'Logistic GLM Summary'; end

    models = struct();

    baseExpr = localJoinPredictors(basePredictors);

    for i = 1:numel(variants)
        variantExpr = localJoinPredictors(variants{i});

        if isempty(baseExpr) && isempty(variantExpr)
            rhs = '1';
        elseif isempty(baseExpr)
            rhs = variantExpr;
        elseif isempty(variantExpr)
            rhs = baseExpr;
        else
            rhs = [baseExpr, ' + ', variantExpr];
        end

        formula = sprintf('Z_fr ~ %s', rhs);

        % remove invalide Z_fr
        trialData(isinf(trialData.Z_fr) | isnan(trialData.Z_fr), :) = [];

        mdl = fitglme(trialData, formula, ...
            'Distribution','Normal','Link','identity','DummyVarCoding','reference');
        models.(sprintf('mdl%d', i)) = mdl;

        fprintf('--- Model %d: %s ---\n', i, formula);
        disp(mdl.Coefficients)

        visualizeLinearGLM(mdl, trialData, keyCatVar, [string(plotName); string(formula)]);
    end
end

function s = localJoinPredictors(p)
% Join a predictor spec into 'a + b + c' form.
    if isempty(p)
        s = '';
    elseif ischar(p)
        s = p;
    elseif isstring(p)
        if isscalar(p)
            s = char(p);
        else
            s = strjoin(cellstr(p), ' + ');
        end
    elseif iscellstr(p) || (iscell(p) && all(cellfun(@(x)ischar(x)||isstring(x), p)))
        s = strjoin(cellfun(@char, p, 'UniformOutput', false), ' + ');
    else
        error('Unsupported predictor spec type.');
    end
end

function out = mapChoicesToLabels(x)
% mapChoicesToLabels - map numeric or string codes to "BUY","HOLD","SELL".
% Returns a string array of same length with values "BUY","HOLD","SELL" (or missing)

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

function [left,right,splitIdx] = splitAtMiddleSpace(s)
% Split S into two parts at the space nearest the middle.

if isstring(s)
    s = char(s);
end
if ~ischar(s) || size(s,1)~=1
    error('Input must be a 1xN char or string scalar.');
end

n = length(s);
if n==0
    left=""; right=""; splitIdx=0; return;
end

spaces = find(s==' ');
if isempty(spaces)
    splitIdx = floor(n/2);
else
    mid = (n+1)/2;
    [~,k] = min(abs(spaces - mid));
    splitIdx = spaces(k);
end

l = s(1:max(splitIdx-1,0));
r = s(min(splitIdx+1,n):end);
left  = string(strtrim(l));
right = string(strtrim(r));

end

function fr = get_Z_fr(C, f, target_event, target_win, baseline_event, baseline_win)
% get_Z_fr - Compute z-scored firing rate per trial for a given session.
%
%   fr = get_Z_fr(C, f, target_event, target_win, baseline_event, baseline_win)
%
%   C             - condition struct with fields .dt, .eventTables, .S
%   f             - session index into C.eventTables / C.S
%   target_event  - event name for the analysis window (string/char)
%   target_win    - [t_start t_end] window (ms) relative to target_event
%   baseline_event- event name used to define baseline window
%   baseline_win  - [t_start t_end] window (ms) relative to baseline_event
%
%   Output:
%       fr        - column vector (nTrials × 1) of z-scored firing rates,
%                   averaged across neurons, per trial.

    % Convert ms windows to sample offsets
    target_interval   = (target_win   ./ 1000) ./ C.dt;
    baseline_interval = (baseline_win ./ 1000) ./ C.dt;
    fr = [];
    try
        % Anchor times for target and baseline events
        target_anchor = C.eventTables{f}.(target_event){1,1}(:);
        target_ts = target_anchor + target_interval(1); 
        target_te = target_anchor + target_interval(2);

        baseline_anchor = C.eventTables{f}.(baseline_event){1,1}(:);
        baseline_ts = baseline_anchor + baseline_interval(1); 
        baseline_te = baseline_anchor + baseline_interval(2);
    catch
        % If any event is missing, return NaNs (assumes 90 trials)
        fr = nan(90, 1);
        return
    end

    % Spike trains for all neurons in this session
    Ssel = C.S{f};
    n_neurons = size(Ssel, 1);
    ntr = numel(target_ts);

    % Preallocate matrices: neurons × trials
    S_interval      = nan(n_neurons, ntr);
    baseline_values = nan(n_neurons, ntr);

    % ----- Baseline mean and std per neuron -----
    baseDur_sec = (baseline_win(2) - baseline_win(1)) / 1000;
    baseBins    = baseDur_sec / C.dt;  % number of samples in baseline window

    for t = 1:ntr
        a = baseline_ts(t); 
        b = baseline_te(t);
        if ~isnan(a) && ~isnan(b)
            baseline_values(:, t) = sum(Ssel(:, a:b), 2) ./ baseDur_sec;
        end
    end
    baseline_mean = mean(baseline_values, 2, 'omitnan');
    baseline_std  = std(baseline_values, 0, 2, 'omitnan');

    % ----- Target window firing rate, z-scored with session baseline -----
    targDur_sec = (target_win(2) - target_win(1)) / 1000;

    for t = 1:ntr
        a = target_ts(t); 
        b = target_te(t);
        if ~isnan(a) && ~isnan(b)
            S_interval(:, t) = sum(Ssel(:, a:b), 2) ./ targDur_sec;
            S_interval(:, t) = (S_interval(:, t) - baseline_mean) ./ baseline_std;
        end
    end

    % Average across neurons and return as column vector (nTrials × 1)
    fr = mean(S_interval, 1, 'omitnan');
    fr = fr(:);
end
