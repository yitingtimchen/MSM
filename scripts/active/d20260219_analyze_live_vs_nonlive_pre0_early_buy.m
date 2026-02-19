% File: d20260219_analyze_live_vs_nonlive_pre0_early_buy.m
%% Confirmatory test: Live_First vs NonLive in pre0-early state
% -------------------------------------------------------------------------
% PURPOSE:
%   Run pre-specified confirmatory tests in the state:
%   prePortfolio==0 and trialNum<=5.
%
% PRIMARY ENDPOINT:
%   P(Buy | prePortfolio==0, trialNum<=5), Live_First - NonLive
%
% SECONDARY ENDPOINT:
%   P(Hold | prePortfolio==0, trialNum<=5), Live_First - NonLive
%
% INPUTS:
%   - condition packs in OUTDIR: AI/Replay/Decoy/Live + OT + Saline variants
%   - event table fields: option, trialNum, prePortfolio(optional),
%     year, month, day, SessionOfDay, monkey
%   - exclusion list from d20260219_load_session_exclusions(EXCLUDE_SCOPE)
%
% METHOD OVERVIEW:
%   - subject-only role for valid Live vs NonLive comparison
%   - match Live_First sessions to nearest NonLive by training index
%     within monkey+treatment
%   - paired sign-flip permutation test + bootstrap CI
%
% OUTPUTS:
%   SAVE_DIR = analysis_analyze_live_vs_nonlive_pre0_early_buy
%   - session_endpoints_pre0_early.csv
%   - matched_pairs_pre0_early.csv
%   - confirmatory_results_overall.csv
%   - confirmatory_results_by_monkey.csv
%   - confirmatory_interaction_treatment.csv
%   - confirmatory_treatment_effects.csv
%   - confirmatory_treatment_pairwise_contrasts.csv
% -------------------------------------------------------------------------

clear; clc; close all;
set(groot, 'defaultTextInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter', 'none');
set(groot, 'defaultLegendInterpreter', 'none');

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';
SAVE_DIR = fullfile(OUTDIR, 'analysis_analyze_live_vs_nonlive_pre0_early_buy');
if ~exist(SAVE_DIR, 'dir'), mkdir(SAVE_DIR); end

cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames = {'baseline','OT','Saline'};

EXCLUDE_SCOPE = "behavior";
EXCLUDE_SESSIONS = d20260219_load_session_exclusions(EXCLUDE_SCOPE);
fprintf('Loaded %d excluded sessions from registry (scope=%s)\n', numel(EXCLUDE_SESSIONS), EXCLUDE_SCOPE);

N_PERM = 10000;
N_BOOT = 5000;
RNG_SEED = 19;
rng(RNG_SEED);

%% ------------------------ REQUIRED ETAB FIELD MAP --------------------------
FIELD_option       = 'option';
FIELD_trialNum     = 'trialNum';
FIELD_prePortfolio = 'prePortfolio'; % optional
FIELD_year         = 'year';
FIELD_month        = 'month';
FIELD_day          = 'day';
FIELD_SessionOfDay = 'SessionOfDay';
FIELD_monkey       = 'monkey';

%% --------------------- BUILD SUBJECT-ONLY TRIAL TABLE ----------------------
trialTable = table();

fprintf('\n=== Building trial-level table (subject-only NonLive + Live_First) ===\n');
for s = 1:numel(cond_sets)
    conds = cond_sets{s};
    treatName = setNames{s};
    fprintf('Treatment: %s\n', treatName);

    for c = 1:numel(conds)
        condRaw = conds{c};
        tok = regexp(condRaw, '(AI|Replay|Decoy|Live)$', 'tokens', 'once');
        if isempty(tok), continue; end
        baseCond = string(tok{1});

        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', condRaw));
        if ~exist(matFile, 'file')
            warning('Missing pack: %s', matFile);
            continue;
        end

        S = load(matFile, 'C');
        C = S.C;
        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            needed = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey, FIELD_option, FIELD_trialNum};
            if ~all(ismember(needed, etab.Properties.VariableNames)), continue; end
            if isempty(etab.(FIELD_option){1}) || isempty(etab.(FIELD_trialNum){1}) || isempty(etab.(FIELD_monkey){1}), continue; end

            y   = etab.(FIELD_year){1}(:);
            mth = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(mth) || isempty(d) || isempty(sod), continue; end

            sessionStr = string(sprintf('%04d-%02d-%02d-%02d', y(1), mth(1), d(1), sod(1)));
            if ~isempty(EXCLUDE_SESSIONS) && any(sessionStr == EXCLUDE_SESSIONS), continue; end

            monkeyStr = string(etab.(FIELD_monkey){1}(1));
            if ~isMonkey12(monkeyStr), continue; end

            actionCode = decodeOptionToAction(etab.(FIELD_option){1}(:));
            trialNum = double(etab.(FIELD_trialNum){1}(:));

            if ismember(FIELD_prePortfolio, etab.Properties.VariableNames) && ~isempty(etab.(FIELD_prePortfolio){1})
                prePortfolio = double(etab.(FIELD_prePortfolio){1}(:));
                hasPre = true;
            else
                prePortfolio = nan(size(trialNum));
                hasPre = false;
            end

            nT = min([numel(actionCode), numel(trialNum), numel(prePortfolio)]);
            actionCode = actionCode(1:nT);
            trialNum = trialNum(1:nT);
            prePortfolio = prePortfolio(1:nT);

            keep = isfinite(actionCode) & isfinite(trialNum);
            if ~any(keep), continue; end
            actionCode = actionCode(keep);
            trialNum = trialNum(keep);
            prePortfolio = prePortfolio(keep);

            if ~hasPre || all(~isfinite(prePortfolio))
                prePortfolio = reconstructPrePortfolio(trialNum, actionCode);
            end
            prePortfolio = max(0, round(prePortfolio));

            if baseCond == "Live"
                liveFlag = "Live_First";
            else
                liveFlag = "NonLive";
            end

            n = numel(actionCode);
            rows = table( ...
                repmat(string(treatName), n, 1), ...
                repmat(string(condRaw),   n, 1), ...
                repmat(baseCond,          n, 1), ...
                repmat(liveFlag,          n, 1), ...
                repmat(monkeyStr,         n, 1), ...
                repmat(sessionStr,        n, 1), ...
                trialNum(:), actionCode(:), prePortfolio(:), ...
                'VariableNames', {'Treatment','ConditionRaw','ConditionBase','LiveFlag','Monkey','Session', ...
                                  'TrialNum','ActionCode','PrePortfolio'} ...
            );
            trialTable = [trialTable; rows]; %#ok<AGROW>
        end
    end
end

if isempty(trialTable)
    error('No rows found for confirmatory analysis.');
end

trialTable.Monkey = categorical(trialTable.Monkey);
trialTable.Treatment = categorical(trialTable.Treatment, {'baseline','OT','Saline'});
trialTable.LiveFlag = categorical(trialTable.LiveFlag, {'NonLive','Live_First'});
trialTable.Action = categorical(mapActionCodeToLabel(trialTable.ActionCode), {'Buy','Hold','Sell'});

trialTable.isBuy  = double(trialTable.Action == 'Buy');
trialTable.isHold = double(trialTable.Action == 'Hold');

%% ----------------------- RESTRICT TO TARGET STATE --------------------------
stateMask = (trialTable.PrePortfolio == 0) & (trialTable.TrialNum <= 5);
stateTbl = trialTable(stateMask, :);
if isempty(stateTbl)
    error('No trials in target state prePortfolio==0 & trialNum<=5.');
end

% Session-level endpoint rates in target state
[G, gM, gT, gL, gS] = findgroups(stateTbl.Monkey, stateTbl.Treatment, stateTbl.LiveFlag, stateTbl.Session);
pBuy = splitapply(@meanFinite, stateTbl.isBuy, G);
pHold = splitapply(@meanFinite, stateTbl.isHold, G);
nTrialsState = splitapply(@numel, stateTbl.isBuy, G);

sess = table(gM, gT, gL, string(gS), nTrialsState, pBuy, pHold, ...
    'VariableNames', {'Monkey','Treatment','LiveFlag','Session','NTrialsState','pBuy_pre0_early','pHold_pre0_early'});

% Add training index using all subject trials (not only target state)
[Gall, aM, aT, aL, aS] = findgroups(trialTable.Monkey, trialTable.Treatment, trialTable.LiveFlag, trialTable.Session);
allSess = table(aM, aT, aL, string(aS), 'VariableNames', {'Monkey','Treatment','LiveFlag','Session'});
allSess = unique(allSess, 'rows');
allSess = sortrows(allSess, {'Monkey','Treatment','Session'});
allSess.TrainIdx = zeros(height(allSess),1);
[Gt, ~, ~] = findgroups(allSess.Monkey, allSess.Treatment);
for gi = unique(Gt).'
    idx = find(Gt == gi);
    [~, ord] = sort(allSess.Session(idx));
    idx = idx(ord);
    allSess.TrainIdx(idx) = (1:numel(idx)).';
end

sess = innerjoin(sess, allSess, 'Keys', {'Monkey','Treatment','LiveFlag','Session'});

writetable(sess, fullfile(SAVE_DIR, 'session_endpoints_pre0_early.csv'));

%% --------------------------- MATCHED PAIRS ---------------------------------
pairTbl = table();
for gi = unique(Gt).'
    A = allSess(Gt == gi, :);
    A_nl = A(A.LiveFlag == 'NonLive', :);
    A_lf = A(A.LiveFlag == 'Live_First', :);
    if isempty(A_nl) || isempty(A_lf), continue; end

    [idxN, idxL] = greedyNearestMatch(A_nl.TrainIdx, A_lf.TrainIdx);
    if isempty(idxN), continue; end

    K_nl = A_nl(idxN, {'Monkey','Treatment','LiveFlag','Session','TrainIdx'});
    K_lf = A_lf(idxL, {'Monkey','Treatment','LiveFlag','Session','TrainIdx'});

    S_nl = innerjoin(K_nl(:, {'Monkey','Treatment','LiveFlag','Session','TrainIdx'}), sess, ...
        'Keys', {'Monkey','Treatment','LiveFlag','Session','TrainIdx'});
    S_lf = innerjoin(K_lf(:, {'Monkey','Treatment','LiveFlag','Session','TrainIdx'}), sess, ...
        'Keys', {'Monkey','Treatment','LiveFlag','Session','TrainIdx'});

    nPair = min(height(S_nl), height(S_lf));
    if nPair == 0, continue; end
    S_nl = S_nl(1:nPair,:);
    S_lf = S_lf(1:nPair,:);

    add = table( ...
        string(S_lf.Monkey), string(S_lf.Treatment), ...
        string(S_lf.Session), string(S_nl.Session), ...
        S_lf.TrainIdx, S_nl.TrainIdx, ...
        S_lf.NTrialsState, S_nl.NTrialsState, ...
        S_lf.pBuy_pre0_early - S_nl.pBuy_pre0_early, ...
        S_lf.pHold_pre0_early - S_nl.pHold_pre0_early, ...
        'VariableNames', {'Monkey','Treatment','LiveSession','NonLiveSession','TrainIdxLive','TrainIdxNonLive', ...
                          'NStateTrialsLive','NStateTrialsNonLive','dBuy_pre0_early','dHold_pre0_early'});
    pairTbl = [pairTbl; add]; %#ok<AGROW>
end

if isempty(pairTbl)
    error('No matched pairs with endpoint data found.');
end
writetable(pairTbl, fullfile(SAVE_DIR, 'matched_pairs_pre0_early.csv'));

%% ---------------------- PRIMARY + SECONDARY TESTS --------------------------
[dP, dP_obs, dP_null, nP] = pairedSignPermTest(pairTbl.dBuy_pre0_early, N_PERM);
[dP_ciLo, dP_ciHi] = bootstrapMeanCI(pairTbl.dBuy_pre0_early, N_BOOT);

[dS, dS_obs, dS_null, nS] = pairedSignPermTest(pairTbl.dHold_pre0_early, N_PERM);
[dS_ciLo, dS_ciHi] = bootstrapMeanCI(pairTbl.dHold_pre0_early, N_BOOT);

resultsOverall = table( ...
    ["Primary"; "Secondary"], ...
    ["pBuy_pre0_early"; "pHold_pre0_early"], ...
    [nP; nS], [dP_obs; dS_obs], [dP_ciLo; dS_ciLo], [dP_ciHi; dS_ciHi], ...
    [dP_null; dS_null], [dP; dS], ...
    'VariableNames', {'Type','Endpoint','NPairsUsed','DeltaMean','CI_Lo','CI_Hi','NullMean','pValue'});
writetable(resultsOverall, fullfile(SAVE_DIR, 'confirmatory_results_overall.csv'));

% By monkey
mkys = unique(string(pairTbl.Monkey), 'stable');
resultsByMonkey = table();
for i = 1:numel(mkys)
    Pm = pairTbl(string(pairTbl.Monkey) == mkys(i), :);

    [p1, obs1, null1, n1] = pairedSignPermTest(Pm.dBuy_pre0_early, N_PERM);
    [lo1, hi1] = bootstrapMeanCI(Pm.dBuy_pre0_early, N_BOOT);
    resultsByMonkey = [resultsByMonkey; table(mkys(i), "Primary", "pBuy_pre0_early", n1, obs1, lo1, hi1, null1, p1, ...
        'VariableNames', {'Monkey','Type','Endpoint','NPairsUsed','DeltaMean','CI_Lo','CI_Hi','NullMean','pValue'})]; %#ok<AGROW>

    [p2, obs2, null2, n2] = pairedSignPermTest(Pm.dHold_pre0_early, N_PERM);
    [lo2, hi2] = bootstrapMeanCI(Pm.dHold_pre0_early, N_BOOT);
    resultsByMonkey = [resultsByMonkey; table(mkys(i), "Secondary", "pHold_pre0_early", n2, obs2, lo2, hi2, null2, p2, ...
        'VariableNames', {'Monkey','Type','Endpoint','NPairsUsed','DeltaMean','CI_Lo','CI_Hi','NullMean','pValue'})]; %#ok<AGROW>
end
writetable(resultsByMonkey, fullfile(SAVE_DIR, 'confirmatory_results_by_monkey.csv'));

% Treatment interaction tests (single interaction test; no subgroup inference)
pairTbl.Treatment = categorical(string(pairTbl.Treatment), {'baseline','OT','Saline'});
[pIntBuy, FBuy, df1Buy, df2Buy] = oneWayPermFTest(pairTbl.dBuy_pre0_early, pairTbl.Treatment, N_PERM);
[pIntHold, FHold, df1Hold, df2Hold] = oneWayPermFTest(pairTbl.dHold_pre0_early, pairTbl.Treatment, N_PERM);

interactionResults = table( ...
    ["Primary"; "Secondary"], ...
    ["dBuy_pre0_early"; "dHold_pre0_early"], ...
    [height(pairTbl); height(pairTbl)], ...
    [FBuy; FHold], [df1Buy; df1Hold], [df2Buy; df2Hold], ...
    [pIntBuy; pIntHold], ...
    'VariableNames', {'Type','EndpointDelta','NTotalPairs','FStat','DF1','DF2','pValue'});
writetable(interactionResults, fullfile(SAVE_DIR, 'confirmatory_interaction_treatment.csv'));

% Treatment effect size + direction and pairwise treatment contrasts
treatLevels = {'baseline','OT','Saline'};
endpointDefs = struct( ...
    'Type', {"Primary","Secondary"}, ...
    'Endpoint', {"dBuy_pre0_early","dHold_pre0_early"} ...
);

treatEffects = table();
pairwiseRows = table();
for ei = 1:numel(endpointDefs)
    eType = string(endpointDefs(ei).Type);
    eName = string(endpointDefs(ei).Endpoint);
    dAll = pairTbl.(eName);
    gAll = categorical(string(pairTbl.Treatment), treatLevels);

    for ti = 1:numel(treatLevels)
        tName = string(treatLevels{ti});
        x = dAll(gAll == tName);
        x = x(isfinite(x));
        n = numel(x);
        if n == 0
            obs = NaN; ciLo = NaN; ciHi = NaN;
        else
            obs = mean(x);
            [ciLo, ciHi] = bootstrapMeanCI(x, N_BOOT);
        end
        dirLabel = directionFromCI(obs, ciLo, ciHi);
        treatEffects = [treatEffects; table(eType, eName, tName, n, obs, ciLo, ciHi, dirLabel, ...
            'VariableNames', {'Type','EndpointDelta','Treatment','NPairsUsed','MeanEffect','CI_Lo','CI_Hi','Direction'})]; %#ok<AGROW>
    end

    contrasts = { ...
        'OT - baseline', 'OT', 'baseline'; ...
        'Saline - baseline', 'Saline', 'baseline'; ...
        'OT - Saline', 'OT', 'Saline' ...
    };
    tmpP = nan(size(contrasts,1),1);
    tmpTbl = table();
    for ci = 1:size(contrasts,1)
        cName = string(contrasts{ci,1});
        gA = string(contrasts{ci,2});
        gB = string(contrasts{ci,3});
        xA = dAll(gAll == gA);
        xB = dAll(gAll == gB);
        [pC, obsC, nullC, nA, nB] = twoSamplePermMeanDiff(xA, xB, N_PERM);
        [loC, hiC] = bootstrapDiffCI(xA, xB, N_BOOT);
        tmpP(ci) = pC;
        tmpTbl = [tmpTbl; table(eType, eName, cName, gA, gB, nA, nB, obsC, loC, hiC, nullC, pC, ...
            'VariableNames', {'Type','EndpointDelta','Contrast','GroupA','GroupB','NGroupA','NGroupB', ...
                              'DeltaOfEffects','CI_Lo','CI_Hi','NullMean','pValue'})]; %#ok<AGROW>
    end
    tmpTbl.pValueFDR = fdrBH(tmpP);
    pairwiseRows = [pairwiseRows; tmpTbl]; %#ok<AGROW>
end

writetable(treatEffects, fullfile(SAVE_DIR, 'confirmatory_treatment_effects.csv'));
writetable(pairwiseRows, fullfile(SAVE_DIR, 'confirmatory_treatment_pairwise_contrasts.csv'));

fprintf('\n=== Confirmatory results (overall) ===\n');
disp(resultsOverall);
fprintf('\n=== Confirmatory results (by monkey) ===\n');
disp(resultsByMonkey);
fprintf('\n=== Confirmatory treatment interaction tests ===\n');
disp(interactionResults);
fprintf('\n=== Confirmatory treatment effects (direction/size) ===\n');
disp(treatEffects);
fprintf('\n=== Confirmatory treatment pairwise contrasts ===\n');
disp(pairwiseRows);

fprintf('\nSaved outputs to:\n  %s\n', SAVE_DIR);
fprintf('Done.\n');

%% ================================ HELPERS ==================================
function actionCode = decodeOptionToAction(optionRaw)
    if iscell(optionRaw) || isstring(optionRaw) || ischar(optionRaw)
        s = lower(strtrim(string(optionRaw)));
        actionCode = nan(size(s));
        actionCode(contains(s,'buy')  | s=='b') = 1;
        actionCode(contains(s,'hold') | s=='h') = 2;
        actionCode(contains(s,'sell') | s=='s') = 3;
        return;
    end

    x = double(optionRaw);
    actionCode = nan(size(x));
    u = unique(x(isfinite(x)));
    if isempty(u), return; end

    if all(ismember(u,[1 2 3]))
        actionCode = x;
        return;
    end
    if all(ismember(u,[0 1 2]))
        actionCode(x==0)=1; actionCode(x==1)=2; actionCode(x==2)=3;
        return;
    end
    warning('decodeOptionToAction: unknown numeric option codes: %s', mat2str(u(:)'));
end

function labels = mapActionCodeToLabel(actionCode)
    labels = strings(size(actionCode));
    labels(actionCode==1) = "Buy";
    labels(actionCode==2) = "Hold";
    labels(actionCode==3) = "Sell";
    labels(labels=="") = "Unknown";
end

function prePortfolio = reconstructPrePortfolio(trialNum, actionCode)
    n = numel(actionCode);
    prePortfolio = nan(n,1);
    p = 0;
    for i = 1:n
        if i==1 || trialNum(i)==1
            p = 0;
        end
        prePortfolio(i) = p;
        if actionCode(i)==1
            p = p + 1;
        elseif actionCode(i)==3
            p = max(p - 1, 0);
        end
    end
end

function ok = isMonkey12(monkeyVal)
    s = upper(strtrim(string(monkeyVal)));
    v = str2double(s);
    ok = (isfinite(v) && any(v==[1 2])) || any(s==["M1","M2"]);
end

function m = meanFinite(x)
    x = x(isfinite(x));
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end

function [idxN, idxL] = greedyNearestMatch(trainN, trainL)
    idxN = [];
    idxL = [];
    if isempty(trainN) || isempty(trainL), return; end

    usedN = false(numel(trainN),1);
    [~, orderL] = sort(trainL, 'ascend');
    for ii = 1:numel(orderL)
        j = orderL(ii);
        avail = find(~usedN);
        if isempty(avail), break; end
        [~, kRel] = min(abs(trainN(avail) - trainL(j)));
        k = avail(kRel);
        usedN(k) = true;
        idxN(end+1,1) = k; %#ok<AGROW>
        idxL(end+1,1) = j; %#ok<AGROW>
    end
end

function [pValue, obsMean, nullMean, nUsed] = pairedSignPermTest(d, nperm)
    d = d(isfinite(d));
    nUsed = numel(d);
    if isempty(d)
        pValue = NaN; obsMean = NaN; nullMean = NaN; return;
    end

    obsMean = mean(d);
    nullDist = nan(nperm,1);
    for b = 1:nperm
        sgn = (randi(2,nUsed,1)*2 - 3);
        nullDist(b) = mean(d .* sgn);
    end
    nullMean = meanFinite(nullDist);
    pValue = (sum(abs(nullDist) >= abs(obsMean)) + 1) / (nperm + 1);
end

function [ciLo, ciHi] = bootstrapMeanCI(d, nboot)
    d = d(isfinite(d));
    if isempty(d)
        ciLo = NaN; ciHi = NaN; return;
    end
    n = numel(d);
    bs = nan(nboot,1);
    for b = 1:nboot
        samp = d(randi(n, n, 1));
        bs(b) = mean(samp);
    end
    ci = prctile(bs, [2.5 97.5]);
    ciLo = ci(1);
    ciHi = ci(2);
end

function [pValue, Fobs, df1, df2] = oneWayPermFTest(y, g, nperm)
% One-way permutation ANOVA on y across groups g.
% pValue from permutation null by shuffling group labels.
    keep = isfinite(y) & ~ismissing(g);
    y = y(keep);
    g = categorical(g(keep));
    if isempty(y) || numel(categories(g)) < 2
        pValue = NaN; Fobs = NaN; df1 = NaN; df2 = NaN; return;
    end

    [~, tbl] = anova1(y, g, 'off');
    Fobs = tbl{2,5};
    df1 = tbl{2,3};
    df2 = tbl{3,3};

    nullF = nan(nperm,1);
    n = numel(y);
    for b = 1:nperm
        gp = g(randperm(n));
        [~, tpb] = anova1(y, gp, 'off');
        nullF(b) = tpb{2,5};
    end
    pValue = (sum(nullF >= Fobs) + 1) / (nperm + 1);
end

function [pValue, obsDiff, nullMean, nA, nB] = twoSamplePermMeanDiff(xA, xB, nperm)
% Permutation test for difference in means: mean(xA)-mean(xB).
    xA = xA(isfinite(xA));
    xB = xB(isfinite(xB));
    nA = numel(xA);
    nB = numel(xB);
    if nA == 0 || nB == 0
        pValue = NaN; obsDiff = NaN; nullMean = NaN; return;
    end

    obsDiff = mean(xA) - mean(xB);
    allx = [xA(:); xB(:)];
    n = numel(allx);
    nullDist = nan(nperm,1);
    for b = 1:nperm
        p = randperm(n);
        A = allx(p(1:nA));
        B = allx(p(nA+1:end));
        nullDist(b) = mean(A) - mean(B);
    end
    nullMean = meanFinite(nullDist);
    pValue = (sum(abs(nullDist) >= abs(obsDiff)) + 1) / (nperm + 1);
end

function [ciLo, ciHi] = bootstrapDiffCI(xA, xB, nboot)
% Bootstrap CI for mean(xA)-mean(xB).
    xA = xA(isfinite(xA));
    xB = xB(isfinite(xB));
    if isempty(xA) || isempty(xB)
        ciLo = NaN; ciHi = NaN; return;
    end
    nA = numel(xA);
    nB = numel(xB);
    bs = nan(nboot,1);
    for b = 1:nboot
        A = xA(randi(nA, nA, 1));
        B = xB(randi(nB, nB, 1));
        bs(b) = mean(A) - mean(B);
    end
    ci = prctile(bs, [2.5 97.5]);
    ciLo = ci(1);
    ciHi = ci(2);
end

function label = directionFromCI(meanVal, ciLo, ciHi)
    if ~isfinite(meanVal) || ~isfinite(ciLo) || ~isfinite(ciHi)
        label = "insufficient_data";
    elseif ciLo > 0
        label = "increase";
    elseif ciHi < 0
        label = "decrease";
    else
        label = "no_clear_change";
    end
end

function q = fdrBH(p)
% Benjamini-Hochberg FDR adjustment.
    p = p(:);
    q = nan(size(p));
    keep = isfinite(p);
    pv = p(keep);
    m = numel(pv);
    if m == 0, return; end
    [ps, ord] = sort(pv);
    qs = ps .* m ./ (1:m)';
    for i = m-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end
    qs = min(qs, 1);
    qv = nan(size(pv));
    qv(ord) = qs;
    q(keep) = qv;
end
