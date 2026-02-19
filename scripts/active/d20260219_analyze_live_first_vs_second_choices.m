% File: d20260219_analyze_live_first_vs_second_choices.m
%% Live First vs Live Second Choice Analysis
% -------------------------------------------------------------------------
% PURPOSE:
%   Quantify Live-only role effects: first mover vs second mover choices.
%
% INPUTS:
%   - condition packs in OUTDIR: Live, OT Live, Saline Live
%   - event table fields:
%       option, optionOpp, year, month, day, SessionOfDay, monkey
%   - exclusion list from d20260219_load_session_exclusions(EXCLUDE_SCOPE)
%
% METHOD OVERVIEW:
%   1) Build trial-level table for subject stream (Live_First) and
%      opponent stream (Live_Second).
%   2) Build paired first-vs-second trial table for within-trial follow test.
%   3) Compute raw action probabilities by monkey and role.
%   4) Test following using within-session permutation null:
%      shuffle second actions inside each session while keeping first
%      actions fixed.
%
% OUTPUTS:
%   SAVE_DIR = analysis_live_first_vs_second_choices
%   - trial_table_live_first_vs_second.csv
%   - raw_choice_rates_by_monkey_live_first_vs_second.csv
%   - raw_choice_rates_by_monkey_live_first_vs_second.png
%   - follow_raw_conditionals.csv
%   - follow_rate_summary.csv
%   - follow_permutation_test.csv
%   - follow_treatment_interaction.csv
%   - follow_treatment_effects.csv
%   - follow_treatment_pairwise_contrasts.csv
% -------------------------------------------------------------------------

clear; clc; close all;
set(groot, 'defaultTextInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter', 'none');
set(groot, 'defaultLegendInterpreter', 'none');

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';
SAVE_DIR = fullfile(OUTDIR, 'analysis_live_first_vs_second_choices');
if ~exist(SAVE_DIR, 'dir'), mkdir(SAVE_DIR); end

cond_sets = { ...
    {'Live'}; ...
    {'OT Live'}; ...
    {'Saline Live'} ...
};
setNames = {'baseline','OT','Saline'};

% Central session exclusion registry
EXCLUDE_SCOPE = "behavior";
EXCLUDE_SESSIONS = d20260219_load_session_exclusions(EXCLUDE_SCOPE);
fprintf('Loaded %d excluded sessions from registry (scope=%s)\n', numel(EXCLUDE_SESSIONS), EXCLUDE_SCOPE);

% Permutation settings
NBOOT = 3000;
RNG_SEED = 19;
rng(RNG_SEED);

%% ------------------------ REQUIRED ETAB FIELD MAP --------------------------
FIELD_option       = 'option';
FIELD_optionOpp    = 'optionOpp';
FIELD_year         = 'year';
FIELD_month        = 'month';
FIELD_day          = 'day';
FIELD_SessionOfDay = 'SessionOfDay';
FIELD_monkey       = 'monkey';

%% ----------------------- BUILD LIVE-ONLY TABLES ----------------------------
trialTable = table();
livePairTable = table();

fprintf('\n=== Building trial-level table (Live_First/Live_Second only) ===\n');
for s = 1:numel(cond_sets)
    conds = cond_sets{s};
    treatName = setNames{s};
    fprintf('Treatment: %s\n', treatName);

    for c = 1:numel(conds)
        condRaw = conds{c};
        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', condRaw));
        if ~exist(matFile, 'file')
            warning('Missing pack: %s', matFile);
            continue;
        end

        S = load(matFile, 'C');
        C = S.C;
        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            needed = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey, FIELD_option, FIELD_optionOpp};
            if ~all(ismember(needed, etab.Properties.VariableNames))
                continue;
            end
            if isempty(etab.(FIELD_option){1}) || isempty(etab.(FIELD_optionOpp){1}) || isempty(etab.(FIELD_monkey){1})
                continue;
            end

            y   = etab.(FIELD_year){1}(:);
            mth = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(mth) || isempty(d) || isempty(sod), continue; end

            sessionStr = string(sprintf('%04d-%02d-%02d-%02d', y(1), mth(1), d(1), sod(1)));
            if ~isempty(EXCLUDE_SESSIONS) && any(sessionStr == EXCLUDE_SESSIONS)
                continue;
            end

            subjMonkey = string(etab.(FIELD_monkey){1}(1));
            if ~isMonkey12(subjMonkey), continue; end
            oppMonkey = inferOpponentMonkey(subjMonkey);

            % Subject stream (Live_First)
            subjActionRaw = decodeOptionToAction(etab.(FIELD_option){1}(:));
            keepSubj = isfinite(subjActionRaw);
            if any(keepSubj)
                subjAction = subjActionRaw(keepSubj);
                n = numel(subjAction);
                rows = table( ...
                    repmat(string(treatName), n, 1), ...
                    repmat(string(condRaw),   n, 1), ...
                    repmat("Live",           n, 1), ...
                    repmat("Live_First",     n, 1), ...
                    repmat("First",          n, 1), ...
                    repmat(subjMonkey,        n, 1), ...
                    repmat(sessionStr,        n, 1), ...
                    subjAction(:), ...
                    'VariableNames', {'Treatment','ConditionRaw','ConditionBase','Group2','MoverRole','Monkey','Session','ActionCode'} ...
                );
                trialTable = [trialTable; rows]; %#ok<AGROW>
            end

            % Opponent stream (Live_Second)
            oppActionRaw = decodeOptionToAction(etab.(FIELD_optionOpp){1}(:));
            keepOpp = isfinite(oppActionRaw);
            if any(keepOpp)
                oppAction = oppActionRaw(keepOpp);
                n = numel(oppAction);
                rows = table( ...
                    repmat(string(treatName), n, 1), ...
                    repmat(string(condRaw),   n, 1), ...
                    repmat("Live",           n, 1), ...
                    repmat("Live_Second",    n, 1), ...
                    repmat("Second",         n, 1), ...
                    repmat(oppMonkey,         n, 1), ...
                    repmat(sessionStr,        n, 1), ...
                    oppAction(:), ...
                    'VariableNames', {'Treatment','ConditionRaw','ConditionBase','Group2','MoverRole','Monkey','Session','ActionCode'} ...
                );
                trialTable = [trialTable; rows]; %#ok<AGROW>
            end

            % Paired table (first vs second by trial index)
            nPair = min(numel(subjActionRaw), numel(oppActionRaw));
            if nPair > 0
                a1 = subjActionRaw(1:nPair);
                a2 = oppActionRaw(1:nPair);
                keepPair = isfinite(a1) & isfinite(a2);
                a1 = a1(keepPair);
                a2 = a2(keepPair);
                if ~isempty(a1)
                    np = numel(a1);
                    pairRows = table( ...
                        repmat(string(treatName), np, 1), ...
                        repmat(string(condRaw),   np, 1), ...
                        repmat(subjMonkey,        np, 1), ...
                        repmat(oppMonkey,         np, 1), ...
                        repmat(sessionStr,        np, 1), ...
                        a1(:), a2(:), ...
                        'VariableNames', {'Treatment','ConditionRaw','MonkeyFirst','MonkeySecond','Session','FirstActionCode','SecondActionCode'} ...
                    );
                    livePairTable = [livePairTable; pairRows]; %#ok<AGROW>
                end
            end
        end
    end
end

if isempty(trialTable)
    error('No trial rows found. Check OUTDIR and live condition packs.');
end

trialTable.Group2 = categorical(trialTable.Group2, {'Live_First','Live_Second'});
trialTable.Monkey = categorical(string(trialTable.Monkey));
trialTable.Treatment = categorical(string(trialTable.Treatment), {'baseline','OT','Saline'});
trialTable.SessionKey = categorical(strcat(string(trialTable.Monkey), "|", string(trialTable.Session)));
trialTable.Action = categorical(mapActionCodeToLabel(trialTable.ActionCode), {'Buy','Hold','Sell'});

fprintf('Built %d rows (%d sessions, monkeys=%s)\n', height(trialTable), ...
    numel(unique(trialTable.SessionKey)), strjoin(string(categories(trialTable.Monkey)), ','));
writetable(trialTable, fullfile(SAVE_DIR, 'trial_table_live_first_vs_second.csv'));

%% -------------------- RAW CHOICE-RATE SUMMARY ------------------------------
[G3, kM, kG, kA] = findgroups(trialTable.Monkey, trialTable.Group2, trialTable.Action);
N = splitapply(@numel, trialTable.ActionCode, G3);
rawCounts = table(kM, kG, kA, N, 'VariableNames', {'Monkey','Group2','Action','N'});

[G2, kM2, kG2] = findgroups(trialTable.Monkey, trialTable.Group2);
NGroup = splitapply(@numel, trialTable.ActionCode, G2);
groupN = table(kM2, kG2, NGroup, 'VariableNames', {'Monkey','Group2','NGroup'});

rawSummary = outerjoin(rawCounts, groupN, 'Keys', {'Monkey','Group2'}, 'MergeKeys', true);
rawSummary.P = rawSummary.N ./ rawSummary.NGroup;

writetable(rawSummary, fullfile(SAVE_DIR, 'raw_choice_rates_by_monkey_live_first_vs_second.csv'));

% Treatment-stratified raw summary
[G3t, tT, tM, tG, tA] = findgroups(trialTable.Treatment, trialTable.Monkey, trialTable.Group2, trialTable.Action);
Nt = splitapply(@numel, trialTable.ActionCode, G3t);
rawCountsT = table(tT, tM, tG, tA, Nt, ...
    'VariableNames', {'Treatment','Monkey','Group2','Action','N'});
[G2t, tT2, tM2, tG2] = findgroups(trialTable.Treatment, trialTable.Monkey, trialTable.Group2);
NGroupT = splitapply(@numel, trialTable.ActionCode, G2t);
groupNT = table(tT2, tM2, tG2, NGroupT, ...
    'VariableNames', {'Treatment','Monkey','Group2','NGroup'});
rawSummaryByTreatment = outerjoin(rawCountsT, groupNT, ...
    'Keys', {'Treatment','Monkey','Group2'}, 'MergeKeys', true);
rawSummaryByTreatment.P = rawSummaryByTreatment.N ./ rawSummaryByTreatment.NGroup;
writetable(rawSummaryByTreatment, fullfile(SAVE_DIR, 'raw_choice_rates_by_treatment_monkey_live_first_vs_second.csv'));

mkys = categories(trialTable.Monkey);
groupOrder = {'Live_First','Live_Second'};
actOrder = {'Buy','Hold','Sell'};
actColors = [0.20 0.55 0.85; 0.60 0.60 0.60; 0.85 0.35 0.20];

figRaw = figure('Color','w', 'Name','Raw choice rates: Live First vs Second');
t = tiledlayout(figRaw, numel(mkys), 1, 'TileSpacing','compact','Padding','compact');
for i = 1:numel(mkys)
    ax = nexttile(t); hold(ax,'on');
    M = zeros(numel(groupOrder), numel(actOrder));
    for g = 1:numel(groupOrder)
        for a = 1:numel(actOrder)
            q = rawSummary.Monkey == mkys{i} & rawSummary.Group2 == groupOrder{g} & rawSummary.Action == actOrder{a};
            if any(q), M(g,a) = rawSummary.P(find(q,1,'first')); end %#ok<FNDSB>
        end
    end
    b = bar(ax, M, 'stacked', 'BarWidth', 0.75);
    for a = 1:numel(b), b(a).FaceColor = actColors(a,:); end
    set(ax, 'XTick', 1:numel(groupOrder), 'XTickLabel', groupOrder, 'YLim', [0 1]);
    ylabel(ax, 'Probability');
    title(ax, sprintf('Monkey %s', string(mkys{i})));
    if i == 1, legend(ax, actOrder, 'Location','eastoutside'); end
end
title(t, 'Raw choice rates: Live First vs Live Second');
exportgraphics(figRaw, fullfile(SAVE_DIR, 'raw_choice_rates_by_monkey_live_first_vs_second.png'), 'Resolution', 300);

%% ================== FOLLOWING WHOEVER GOES FIRST ===========================
fprintf('\n=== Live-only following whoever goes first ===\n');
if isempty(livePairTable)
    warning('livePairTable is empty. Following cannot be tested.');
    followSummary = table();
    followPerm = table();
    followRaw = table();
else
    livePairTable.Treatment   = categorical(string(livePairTable.Treatment), {'baseline','OT','Saline'});
    livePairTable.MonkeySecond = categorical(string(livePairTable.MonkeySecond));
    livePairTable.MonkeyFirst  = categorical(string(livePairTable.MonkeyFirst));
    livePairTable.FirstAction  = categorical(mapActionCodeToLabel(livePairTable.FirstActionCode), {'Buy','Hold','Sell'});
    livePairTable.SecondAction = categorical(mapActionCodeToLabel(livePairTable.SecondActionCode), {'Buy','Hold','Sell'});
    livePairTable.SessionKey   = categorical(strcat(string(livePairTable.MonkeySecond), "|", string(livePairTable.Session)));
    livePairTable.follow = double(livePairTable.SecondAction == livePairTable.FirstAction);

    [Ga, aM, aF, aS] = findgroups(livePairTable.MonkeySecond, livePairTable.FirstAction, livePairTable.SecondAction);
    Na = splitapply(@numel, livePairTable.SecondActionCode, Ga);
    followCounts = table(aM, aF, aS, Na, 'VariableNames', {'MonkeySecond','FirstAction','SecondAction','N'});
    [Gb, bM, bF] = findgroups(livePairTable.MonkeySecond, livePairTable.FirstAction);
    Nb = splitapply(@numel, livePairTable.SecondActionCode, Gb);
    followTot = table(bM, bF, Nb, 'VariableNames', {'MonkeySecond','FirstAction','NTotal'});
    followRaw = outerjoin(followCounts, followTot, 'Keys', {'MonkeySecond','FirstAction'}, 'MergeKeys', true);
    followRaw.P = followRaw.N ./ followRaw.NTotal;
    writetable(followRaw, fullfile(SAVE_DIR, 'follow_raw_conditionals.csv'));

    [Gf, fM] = findgroups(livePairTable.MonkeySecond);
    followRate = splitapply(@meanFinite, livePairTable.follow, Gf);
    Nfollow = splitapply(@numel, livePairTable.follow, Gf);
    followSummary = table(fM, Nfollow, followRate, ...
        'VariableNames', {'MonkeySecond','NTrials','FollowRate'});
    followSummary = [followSummary; ...
        table(categorical("All"), height(livePairTable), meanFinite(livePairTable.follow), ...
        'VariableNames', {'MonkeySecond','NTrials','FollowRate'})];
    writetable(followSummary, fullfile(SAVE_DIR, 'follow_rate_summary.csv'));

    followPerm = table();
    [pAll, obsAll, nullMeanAll] = permutationFollowTest(livePairTable, NBOOT);
    followPerm = [followPerm; table(categorical("All"), height(livePairTable), obsAll, nullMeanAll, pAll, ...
        'VariableNames', {'MonkeySecond','NTrials','ObservedFollowRate','NullMeanFollowRate','pValue'})]; %#ok<AGROW>

    mkys2 = categories(livePairTable.MonkeySecond);
    for i = 1:numel(mkys2)
        Tm = livePairTable(livePairTable.MonkeySecond == mkys2{i}, :);
        if isempty(Tm), continue; end
        [pM, obsM, nullMeanM] = permutationFollowTest(Tm, NBOOT);
        followPerm = [followPerm; table(mkys2(i), height(Tm), obsM, nullMeanM, pM, ...
            'VariableNames', {'MonkeySecond','NTrials','ObservedFollowRate','NullMeanFollowRate','pValue'})]; %#ok<AGROW>
    end
    writetable(followPerm, fullfile(SAVE_DIR, 'follow_permutation_test.csv'));

    % Treatment interaction test at session level:
    % endpoint = follow excess = observed follow - within-session null mean.
    [Gsess, kTreat, kMonk, kSess] = findgroups(livePairTable.Treatment, livePairTable.MonkeySecond, livePairTable.Session);
    nSess = max(Gsess);
    sessRows = table();
    npermSessionNull = max(500, round(NBOOT/3));
    for gi = 1:nSess
        Tsi = livePairTable(Gsess == gi, :);
        if isempty(Tsi), continue; end
        [~, obsS, nullS] = permutationFollowTest(Tsi, npermSessionNull);
        sessRows = [sessRows; table(kTreat(gi), kMonk(gi), string(kSess(gi)), height(Tsi), obsS, nullS, obsS-nullS, ...
            'VariableNames', {'Treatment','MonkeySecond','Session','NTrials','ObservedFollowRate','NullMeanFollowRate','FollowExcess'})]; %#ok<AGROW>
    end

    [pIntRate, FRate, df1Rate, df2Rate] = oneWayPermFTest(sessRows.ObservedFollowRate, sessRows.Treatment, NBOOT);
    [pIntEx, FEx, df1Ex, df2Ex] = oneWayPermFTest(sessRows.FollowExcess, sessRows.Treatment, NBOOT);
    interactionResults = table( ...
        ["ObservedFollowRate"; "FollowExcess"], ...
        [height(sessRows); height(sessRows)], ...
        [FRate; FEx], [df1Rate; df1Ex], [df2Rate; df2Ex], ...
        [pIntRate; pIntEx], ...
        'VariableNames', {'Endpoint','NSessions','FStat','DF1','DF2','pValue'});

    % Treatment effect size + direction and pairwise contrasts
    treatLevels = {'baseline','OT','Saline'};
    endpointList = {'ObservedFollowRate','FollowExcess'};
    treatEffects = table();
    pairwiseRows = table();
    for ei = 1:numel(endpointList)
        eName = string(endpointList{ei});
        yAll = sessRows.(eName);
        gAll = categorical(string(sessRows.Treatment), treatLevels);

        for ti = 1:numel(treatLevels)
            tName = string(treatLevels{ti});
            x = yAll(gAll == tName);
            x = x(isfinite(x));
            n = numel(x);
            if n == 0
                obs = NaN; ciLo = NaN; ciHi = NaN;
            else
                obs = mean(x);
                [ciLo, ciHi] = bootstrapMeanCI(x, NBOOT);
            end
            dirLabel = directionFromCI(obs, ciLo, ciHi);
            treatEffects = [treatEffects; table(eName, tName, n, obs, ciLo, ciHi, dirLabel, ...
                'VariableNames', {'Endpoint','Treatment','NSessions','MeanEffect','CI_Lo','CI_Hi','Direction'})]; %#ok<AGROW>
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
            xA = yAll(gAll == gA);
            xB = yAll(gAll == gB);
            [pC, obsC, nullC, nA, nB] = twoSamplePermMeanDiff(xA, xB, NBOOT);
            [loC, hiC] = bootstrapDiffCI(xA, xB, NBOOT);
            tmpP(ci) = pC;
            tmpTbl = [tmpTbl; table(eName, cName, gA, gB, nA, nB, obsC, loC, hiC, nullC, pC, ...
                'VariableNames', {'Endpoint','Contrast','GroupA','GroupB','NGroupA','NGroupB', ...
                                  'DeltaOfEffects','CI_Lo','CI_Hi','NullMean','pValue'})]; %#ok<AGROW>
        end
        tmpTbl.pValueFDR = fdrBH(tmpP);
        pairwiseRows = [pairwiseRows; tmpTbl]; %#ok<AGROW>
    end

    writetable(sessRows, fullfile(SAVE_DIR, 'follow_session_level_summary.csv'));
    writetable(interactionResults, fullfile(SAVE_DIR, 'follow_treatment_interaction.csv'));
    writetable(treatEffects, fullfile(SAVE_DIR, 'follow_treatment_effects.csv'));
    writetable(pairwiseRows, fullfile(SAVE_DIR, 'follow_treatment_pairwise_contrasts.csv'));

    disp(followSummary);
    disp(followPerm);
    disp(interactionResults);
    disp(treatEffects);
    disp(pairwiseRows);
end

fprintf('\nSaved outputs to:\n  %s\n', SAVE_DIR);
fprintf('Done.\n');

%% ================================ HELPERS ==================================
function actionCode = decodeOptionToAction(optionRaw)
% BUY=1, HOLD=2, SELL=3
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

function opp = inferOpponentMonkey(subj)
    s = strtrim(string(subj));
    v = str2double(s);
    if isfinite(v)
        if v==1, opp="2"; return; end
        if v==2, opp="1"; return; end
    end
    su = upper(s);
    if su=="M1", opp="M2"; return; end
    if su=="M2", opp="M1"; return; end
    error("inferOpponentMonkey: cannot infer opponent from '%s'.", s);
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

function [pValue, obsRate, nullMean] = permutationFollowTest(T, nperm)
% One-sided permutation test for follow rate > chance.
% Null: shuffle second actions within session, keep first actions fixed.

    if isempty(T)
        pValue = NaN; obsRate = NaN; nullMean = NaN; return;
    end

    first = double(T.FirstActionCode(:));
    second = double(T.SecondActionCode(:));
    keep = isfinite(first) & isfinite(second);
    first = first(keep);
    second = second(keep);
    if isempty(first)
        pValue = NaN; obsRate = NaN; nullMean = NaN; return;
    end

    obsRate = mean(first == second);

    if ismember('SessionKey', T.Properties.VariableNames)
        sess = categorical(T.SessionKey(keep));
    else
        sess = categorical(T.Session(keep));
    end
    [Gs, ~] = findgroups(sess);
    uniqGs = unique(Gs);

    nullDist = nan(nperm,1);
    for b = 1:nperm
        secP = second;
        for i = 1:numel(uniqGs)
            idx = find(Gs == uniqGs(i));
            if numel(idx) > 1
                secP(idx) = secP(idx(randperm(numel(idx))));
            end
        end
        nullDist(b) = mean(first == secP);
    end

    nullMean = meanFinite(nullDist);
    pValue = (sum(nullDist >= obsRate) + 1) / (nperm + 1);
end

function [pValue, Fobs, df1, df2] = oneWayPermFTest(y, g, nperm)
% One-way permutation ANOVA on y across groups g.
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

    n = numel(y);
    nullF = nan(nperm,1);
    for b = 1:nperm
        gp = g(randperm(n));
        [~, tpb] = anova1(y, gp, 'off');
        nullF(b) = tpb{2,5};
    end
    pValue = (sum(nullF >= Fobs) + 1) / (nperm + 1);
end

function [ciLo, ciHi] = bootstrapMeanCI(x, nboot)
    x = x(isfinite(x));
    if isempty(x)
        ciLo = NaN; ciHi = NaN; return;
    end
    n = numel(x);
    bs = nan(nboot,1);
    for b = 1:nboot
        samp = x(randi(n, n, 1));
        bs(b) = mean(samp);
    end
    ci = prctile(bs, [2.5 97.5]);
    ciLo = ci(1);
    ciHi = ci(2);
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
