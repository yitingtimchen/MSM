%% MSM "Best Choice" Analysis (ROLE-AWARE LIVE): DP-optimal vs Fixed-Fundamental Benchmarks
% -----------------------------------------------------------------------------
% WHAT THIS VERSION CHANGES (vs your current optimal_choice_analysis):
%   1) For baseCond=="Live", it extracts TWO trial streams per file:
%        - Role="Subj": subject's own trials (existing behavior)
%        - Role="Opp" : trials where the subject pack stores opponent choices in *_Opp fields
%   2) For Role="Opp", Monkey is inferred as the counterpart (because etab.monkey is ALWAYS subject)
%   3) Adds a Role column throughout (trialTableAll + sessionSummaryAll)
%   4) Builds the canonical price schedule *from Role=="Subj"* by default (avoids mixing streams)
%   5) Plotting: if cb=="Live" and Role=="Opp", face color is lightened for visibility
%
% NOTE:
%   This script assumes the condition packs store opponent-role vectors in fields like:
%     optionOpp, trialNumOpp, marketOrigOpp, priceBuyOpp, priceSellOpp, prePortfolioOpp
%   If your packs use different names, update the FIELD MAP section only.

clear; clc; close all;

%% ========================== CONFIG (mirror your structure) ===================
OUTDIR = 'C:\Users\plattlab\Tim\Stock_market_experiment\Tim\condition_packs_behavior_only_from_dataTabFULL'; % <-- EDIT

cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames     = {'baseline','OT','Saline'};
baseCondList = {'AI','Replay','Decoy','Live'};

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
EXCLUDE_SESSIONS = { ...
    % '2018-07-07-01', ...
};
EXCLUDE_SESSIONS = string(EXCLUDE_SESSIONS(:));
EXCLUDE_SESSIONS = EXCLUDE_SESSIONS(strlength(EXCLUDE_SESSIONS) > 0);
excludedSessSeen = strings(0,1);

% Treatment -> marker shape
treatOrder  = {'baseline','OT','Saline'};
markerList  = {'o','s','^','d','v','>'}; % only first 3 used
treatMarker = containers.Map(treatOrder, markerList(1:3));

% Base condition -> face color
condFaceColors = containers.Map();
condFaceColors('AI')     = [0.0000 0.4470 0.7410];
condFaceColors('Replay') = [0.8500 0.3250 0.0980];
condFaceColors('Decoy')  = [0.9290 0.6940 0.1250];
condFaceColors('Live')   = [0.4940 0.1840 0.5560];

% Consistent x ticks
xTicksTrials = 0:3:15;

% IMPORTANT: choose how to construct canonical schedule:
%   true  -> schedule built only from Role=="Subj" rows (recommended)
%   false -> schedule built from all rows (will assert invariance across mixed streams)
SCHEDULE_FROM_SUBJECT_ONLY = true;

%% ========================== TASK / MODEL PARAMETERS ==========================
Tmax = 15;
Pmax = 15;
allowance_dollars = 80;
dividend_values_dollars = [0 0.8 2.8 6];
muDiv_DP_dollars = mean(dividend_values_dollars); % 2.4

muDiv_Fund_dollars = muDiv_DP_dollars;
trialIdx = (1:Tmax)';
fundamental_byTrial_dollars = muDiv_Fund_dollars .* (Tmax - trialIdx + 1);

%% ========================== ACTION CODING ===================================
ACTION_BUY  = 1;
ACTION_HOLD = 2;
ACTION_SELL = 3;

%% ========================== FIELD MAP =======================================
% ---- Subject fields (existing) ----
FIELD_trialNum      = 'trialNum';
FIELD_marketOrig    = 'marketOrig';
FIELD_priceBuy      = 'priceBuy';
FIELD_priceSell     = 'priceSell';
FIELD_prePortfolio  = 'prePortfolio';
FIELD_option        = 'option';

% ---- Opponent-role fields (Live only) ----
FIELD_trialNumOpp      = 'trialNumOpp';
FIELD_marketOrigOpp    = 'marketOrigOpp';      % if missing in your data, set = FIELD_marketOrig
FIELD_priceBuyOpp      = 'priceBuyOpp';
FIELD_priceSellOpp     = 'priceSellOpp';
FIELD_prePortfolioOpp  = 'prePortfolioOpp';
FIELD_optionOpp        = 'optionOpp';

% ---- Session identity fields ----
FIELD_year          = 'year';
FIELD_month         = 'month';
FIELD_day           = 'day';
FIELD_SessionOfDay  = 'SessionOfDay';
FIELD_monkey        = 'monkey';

%% ========================== BUILD TRIAL-LEVEL TABLE ==========================
trialTableAll = table();

fprintf('\n=== BUILDING trialTableAll (trial-level rows) ===\n');

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

        % Parse base condition label
        tok = regexp(cname, '(AI|Replay|Decoy|Live)$', 'tokens', 'once');
        if isempty(tok)
            warning('Could not parse base condition from name: %s', cname);
            baseCond = string(cname);
        else
            baseCond = string(tok{1});
        end

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            % ---- Required session fields (role-specific action fields handled below) ----
            neededSession = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey};
            if ~all(ismember(neededSession, etab.Properties.VariableNames))
                warning('Missing required session fields in %s (file %d). Skipping.', cname, f);
                continue;
            end

            % ---- Session string ----
            y   = etab.(FIELD_year){1}(:);
            mth = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(mth) || isempty(d) || isempty(sod)
                continue;
            end
            sessionStr = string(sprintf('%04d-%02d-%02d-%02d', y(1), mth(1), d(1), sod(1)));

            % exclusions
            if ~isempty(EXCLUDE_SESSIONS) && any(sessionStr == EXCLUDE_SESSIONS)
                excludedSessSeen(end+1,1) = sessionStr; %#ok<AGROW>
                continue;
            end

            % subject monkey id (this is ALWAYS the subject in the pack)
            monkeyRaw = etab.(FIELD_monkey){1};
            if isempty(monkeyRaw), continue; end
            subjMonkeyStr = string(monkeyRaw(1));

            % ---- Decide which roles to extract (Live only) ----
            if baseCond == "Live"
                roleList = ["Subj","Opp"];
            else
                roleList = "Subj";
            end

            for role = roleList

                % ---- Role-specific field selection ----
                if role == "Subj"
                    trlField  = FIELD_trialNum;
                    mktField  = FIELD_marketOrig;
                    buyField  = FIELD_priceBuy;
                    sellField = FIELD_priceSell;
                    preField  = FIELD_prePortfolio;
                    optField  = FIELD_option;
                    actorMonkeyStr = subjMonkeyStr;
                else
                    trlField  = FIELD_trialNumOpp;
                    mktField  = FIELD_marketOrigOpp;
                    buyField  = FIELD_priceBuyOpp;
                    sellField = FIELD_priceSellOpp;
                    preField  = FIELD_prePortfolioOpp;
                    optField  = FIELD_optionOpp;

                    % IMPORTANT: for opponent-role rows, infer the ACTOR monkey as counterpart
                    actorMonkeyStr = inferOpponentMonkey(subjMonkeyStr);
                end

                % ---- Required role fields ----
                neededRole = {trlField, buyField, sellField, optField};

                % market field: some packs might not have marketOrigOpp
                if ~ismember(mktField, etab.Properties.VariableNames)
                    if role=="Opp"
                        % fall back: use subject marketOrig if Opp version doesn't exist
                        mktField = FIELD_marketOrig;
                    end
                end
                neededRole = [neededRole, {mktField}];

                if ~all(ismember(neededRole, etab.Properties.VariableNames))
                    % For Live Opp, skip cleanly if Opp fields absent
                    continue;
                end

                if isempty(etab.(trlField){1}) || isempty(etab.(optField){1})
                    continue;
                end

                % prePortfolio optional: reconstruct if missing/empty
                hasPre = ismember(preField, etab.Properties.VariableNames) && ~isempty(etab.(preField){1});

                % ---- Pull vectors ----
                trialNum   = etab.(trlField){1}(:);
                marketOrig = etab.(mktField){1}(:);
                priceBuy   = etab.(buyField){1}(:);
                priceSell  = etab.(sellField){1}(:);
                optionRaw  = etab.(optField){1}(:);

                if hasPre
                    prePortfolio = etab.(preField){1}(:);
                else
                    prePortfolio = nan(size(trialNum));
                end

                % ---- Align lengths ----
                nT = min([numel(trialNum), numel(marketOrig), numel(priceBuy), numel(priceSell), numel(optionRaw), numel(prePortfolio)]);
                trialNum     = trialNum(1:nT);
                marketOrig   = marketOrig(1:nT);
                priceBuy     = priceBuy(1:nT);
                priceSell    = priceSell(1:nT);
                optionRaw    = optionRaw(1:nT);
                prePortfolio = prePortfolio(1:nT);

                % ---- Drop missing essential ----
                keep = ~isnan(trialNum) & ~isnan(marketOrig) & ~isnan(priceBuy) & ~isnan(priceSell);
                if ~any(keep), continue; end

                trialNum     = trialNum(keep);
                marketOrig   = marketOrig(keep);
                priceBuy     = priceBuy(keep);
                priceSell    = priceSell(keep);
                optionRaw    = optionRaw(keep);
                prePortfolio = prePortfolio(keep);

                % ---- Decode actions (robust) ----
                [actionTaken_code, actionMappingReport] = decodeActionsStandalone(optionRaw, trialNum, marketOrig, prePortfolio);

                % ---- Reconstruct prePortfolio if missing/NaN ----
                if ~hasPre || all(isnan(prePortfolio))
                    prePortfolio = reconstructPrePortfolioFromActions(trialNum, marketOrig, actionTaken_code);
                end

                % ---- Emit trial-level rows ----
                nKeep = numel(trialNum);
                newRows = table( ...
                    repmat(string(treatName), nKeep,1), ...
                    repmat(string(cname), nKeep,1), ...
                    repmat(baseCond, nKeep,1), ...
                    repmat(string(role), nKeep,1), ...
                    repmat(actorMonkeyStr, nKeep,1), ...
                    repmat(sessionStr, nKeep,1), ...
                    trialNum, marketOrig, priceBuy, priceSell, prePortfolio, actionTaken_code, ...
                    repmat(actionMappingReport, nKeep,1), ...
                    'VariableNames', { ...
                        'Treatment','ConditionRaw','ConditionBase','Role', ...
                        'Monkey','Session', ...
                        'trialNum','marketOrig','priceBuy','priceSell','prePortfolio','actionTaken', ...
                        'actionMappingReport' ...
                    } ...
                );

                trialTableAll = [trialTableAll; newRows]; %#ok<AGROW>
            end
        end
    end
end

%% -------------------------- Exclusion report ---------------------------------
if ~isempty(EXCLUDE_SESSIONS)
    excludedSessSeenU = unique(excludedSessSeen);
    notFound = setdiff(EXCLUDE_SESSIONS, excludedSessSeenU);

    fprintf('\n=== Exclusion report ===\n');
    fprintf('Requested exclusions: %d\n', numel(EXCLUDE_SESSIONS));
    fprintf('Excluded (found+skipped): %d\n', numel(excludedSessSeenU));
    if ~isempty(excludedSessSeenU)
        disp('Excluded sessions actually seen:'); disp(excludedSessSeenU);
    end
    if ~isempty(notFound)
        disp('Exclusion sessions NOT found in loaded data (check strings):'); disp(notFound);
    end
end

if isempty(trialTableAll)
    error('No trials loaded. Check OUTDIR, condition pack filenames, and field map.');
end

%% ========================== FACTORS + BASIC CLEANUP ==========================
trialTableAll.Monkey = categorical(string(trialTableAll.Monkey));
trialTableAll.Role   = categorical(string(trialTableAll.Role), ["Subj","Opp"]);

trialTableAll.Treatment = categorical(string(trialTableAll.Treatment));
trialTableAll.Treatment = removecats(trialTableAll.Treatment);
trialTableAll.Treatment = reordercats(trialTableAll.Treatment, ...
    treatOrder(ismember(treatOrder, categories(trialTableAll.Treatment))));

trialTableAll.ConditionBase = categorical(string(trialTableAll.ConditionBase), baseCondList);
trialTableAll.LiveFlag = categorical( ...
    string(trialTableAll.ConditionBase) == "Live", ...
    [false true], {'NonLive','Live'} ...
);

% Market templates inferred from *all* loaded rows (safe)
market_templates = unique(trialTableAll.marketOrig(~isnan(trialTableAll.marketOrig)));
market_templates = sort(market_templates(:)');
fprintf('\nMarket templates detected in loaded data: %s\n', mat2str(market_templates));

%% ========================== BUILD CANONICAL PRICE SCHEDULES ==================
% WHY: DP + Fundamental policies require a single canonical schedule per (marketOrig, trialNum).
% ROLE PITFALL: If Opp prices differ (or include NaNs), mixing could break invariance assertion.
if SCHEDULE_FROM_SUBJECT_ONLY
    schedTab = trialTableAll(trialTableAll.Role=="Subj", :);
else
    schedTab = trialTableAll;
end

schedKey = findgroups(schedTab.marketOrig, schedTab.trialNum);

uBuy  = splitapply(@(x) numel(unique(round(x(~isnan(x)), 6))), schedTab.priceBuy,  schedKey);
uSell = splitapply(@(x) numel(unique(round(x(~isnan(x)), 6))), schedTab.priceSell, schedKey);
assert(all(uBuy==1 & uSell==1), 'Price schedule invariance failed within (marketOrig, trialNum).');

marketKey = splitapply(@(x) x(find(~isnan(x),1,'first')), schedTab.marketOrig, schedKey);
trialKey  = splitapply(@(x) x(find(~isnan(x),1,'first')), schedTab.trialNum,   schedKey);
buyKey    = splitapply(@(x) x(find(~isnan(x),1,'first')), schedTab.priceBuy,   schedKey);
sellKey   = splitapply(@(x) x(find(~isnan(x),1,'first')), schedTab.priceSell,  schedKey);

scheduleTable = table(marketKey, trialKey, buyKey, sellKey, ...
    'VariableNames', {'marketOrig','trialNum','priceBuy','priceSell'});

priceBuySched  = nan(numel(market_templates), Tmax);
priceSellSched = nan(numel(market_templates), Tmax);

for mi = 1:numel(market_templates)
    m = market_templates(mi);
    for t = 1:Tmax
        ii = scheduleTable.marketOrig==m & scheduleTable.trialNum==t;
        if nnz(ii)~=1
            error('Missing or duplicate schedule row for marketOrig=%d trial=%d', m, t);
        end
        priceBuySched(mi,t)  = scheduleTable.priceBuy(ii);
        priceSellSched(mi,t) = scheduleTable.priceSell(ii);
    end
end

%% ========================== BENCHMARK #1: DP OPTIMAL POLICY ==================
policyDP = cell(1, numel(market_templates)); % policyDP{mi} is [Tmax x (Pmax+1)]

for mi = 1:numel(market_templates)
    buySched  = priceBuySched(mi,:).';
    sellSched = priceSellSched(mi,:).';

    pol   = nan(Tmax, Pmax+1);
    Vnext = zeros(Pmax+1,1);

    for t = Tmax:-1:1
        V = -inf(Pmax+1,1);
        pMaxFeasible = min(Pmax, t-1);

        for p = 0:pMaxFeasible
            qB = (allowance_dollars - buySched(t)) + muDiv_DP_dollars*(p+1) + Vnext((p+1)+1);
            qH = allowance_dollars + muDiv_DP_dollars*p + Vnext(p+1);

            if p>0
                qS = (allowance_dollars + sellSched(t)) + muDiv_DP_dollars*(p-1) + Vnext((p-1)+1);
            else
                qS = -inf;
            end

            [V(p+1), a] = max([qB qH qS]);
            pol(t,p+1) = a;
        end

        V(pMaxFeasible+2:end) = -inf;
        Vnext = V;
    end

    policyDP{mi} = pol;
end

%% ========================== BENCHMARK #2: FIXED-FUNDAMENTAL ==================
policyFundamental_raw = nan(numel(market_templates), Tmax);

for mi = 1:numel(market_templates)
    for t = 1:Tmax
        if priceSellSched(mi,t) > fundamental_byTrial_dollars(t)
            policyFundamental_raw(mi,t) = ACTION_SELL;
        elseif priceBuySched(mi,t) < fundamental_byTrial_dollars(t)
            policyFundamental_raw(mi,t) = ACTION_BUY;
        else
            policyFundamental_raw(mi,t) = ACTION_HOLD;
        end
    end
end

%% ========================== LABEL BEST ACTIONS PER TRIAL =====================
bestDP_code   = nan(height(trialTableAll),1);
bestFund_code = nan(height(trialTableAll),1);

for i = 1:height(trialTableAll)
    m = trialTableAll.marketOrig(i);
    t = trialTableAll.trialNum(i);
    p = trialTableAll.prePortfolio(i);

    if isnan(m) || isnan(t) || isnan(p) || t<1 || t>Tmax
        continue;
    end

    mi = find(market_templates==m, 1);
    if isempty(mi), continue; end

    p = max(0, min(Pmax, round(p)));

    bestDP_code(i) = policyDP{mi}(t, p+1);

    aFund = policyFundamental_raw(mi,t);
    if (aFund==ACTION_SELL) && (p==0)
        aFund = ACTION_HOLD; % infeasible sell -> hold
    end
    bestFund_code(i) = aFund;
end

trialTableAll.bestDP = bestDP_code;
trialTableAll.bestFundamental = bestFund_code;

trialTableAll.isBestDP = (trialTableAll.actionTaken == trialTableAll.bestDP);
trialTableAll.isBestFundamental = (trialTableAll.actionTaken == trialTableAll.bestFundamental);
trialTableAll.benchmarksDisagree = (trialTableAll.bestDP ~= trialTableAll.bestFundamental);

validRow = ~isnan(trialTableAll.bestDP) & ~isnan(trialTableAll.bestFundamental) & ~isnan(trialTableAll.actionTaken);
trialTableAll = trialTableAll(validRow,:);

%% ========================== SESSION-LEVEL SUMMARY (ROLE-AWARE) ================
% KEY CHANGE: include Role in grouping so Live splits into Subj vs Opp streams.
[G, keyMonkey, keySession, keyTreatment, keyConditionBase, keyLiveFlag, keyRole] = findgroups( ...
    trialTableAll.Monkey, trialTableAll.Session, trialTableAll.Treatment, trialTableAll.ConditionBase, trialTableAll.LiveFlag, trialTableAll.Role);

nTrials   = splitapply(@numel, trialTableAll.isBestDP, G);
pBestDP   = splitapply(@mean,  trialTableAll.isBestDP, G);
pBestFund = splitapply(@mean,  trialTableAll.isBestFundamental, G);
pDisagree = splitapply(@mean,  trialTableAll.benchmarksDisagree, G);

sessionSummaryAll = table(keyMonkey, keySession, keyTreatment, keyConditionBase, keyLiveFlag, keyRole, ...
    nTrials, pBestDP, pBestFund, pDisagree, ...
    'VariableNames', {'Monkey','Session','Treatment','ConditionBase','LiveFlag','Role', ...
                      'nTrials','pBestDP','pBestFV','pBenchmarksDisagree'});

%% -------------------------- sessionIdx per monkey (ROLE-INDEPENDENT ORDER) ----
% You want the x-axis to be session progression within monkey, not duplicated by Role.
sessTab = unique(sessionSummaryAll(:, {'Monkey','Session'}));
[sessTabSorted, sortIdxSess] = sortrows(sessTab, {'Monkey','Session'});

sessionIdxSorted = zeros(height(sessTabSorted),1);
monkeys = categories(sessTabSorted.Monkey);
for im = 1:numel(monkeys)
    mask = sessTabSorted.Monkey == monkeys{im};
    sessionIdxSorted(mask) = (1:nnz(mask)).';
end
sessTabSorted.sessionIdx = sessionIdxSorted;

% join back into sessionSummaryAll
sessionSummaryAll = outerjoin(sessionSummaryAll, sessTabSorted, 'Keys', {'Monkey','Session'}, 'MergeKeys', true);

%% ========================== SAVE TABLES (CSV) =================================
try
    writetable(trialTableAll, fullfile(OUTDIR, 'msm_bestChoice_trialTableAll_ROLEAWARE.csv'));
    writetable(sessionSummaryAll, fullfile(OUTDIR, 'msm_bestChoice_sessionSummaryAll_ROLEAWARE.csv'));
    fprintf('\nSaved CSV outputs to OUTDIR.\n');
catch ME
    warning('Could not write CSV outputs: %s', ME.message);
end

%% ========================== FIG 1: Session trends ============================
figure('Color','w','Name','FIG 1: Session trends (Best-choice rates; role-aware)');
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% pooled copy (3rd column)
pooled = sessionSummaryAll;
pooled.Monkey = categorical(repmat("Pooled", height(pooled),1));

% column tabs
colTabs = cell(1,3);
colNames = cell(1,3);

if numel(monkeys) >= 1
    colTabs{1} = sessionSummaryAll(sessionSummaryAll.Monkey==monkeys{1},:);
    colNames{1} = sprintf('Monkey %s', string(monkeys{1}));
else
    colTabs{1} = sessionSummaryAll([],:); colNames{1} = 'Monkey ?';
end

if numel(monkeys) >= 2
    colTabs{2} = sessionSummaryAll(sessionSummaryAll.Monkey==monkeys{2},:);
    colNames{2} = sprintf('Monkey %s', string(monkeys{2}));
else
    colTabs{2} = sessionSummaryAll([],:); colNames{2} = 'Monkey ?';
end

colTabs{3} = pooled;
colNames{3} = 'Pooled';

treatLabels = treatOrder;
condLabels  = baseCondList;

% Row 1: DP best
for c = 1:3
    ax = nexttile(tl, c);
    plotSessionScatter_ROLE(ax, colTabs{c}, 'pBestDP', sprintf('%s | DP best', colNames{c}), ...
        treatMarker, condFaceColors);
    ylim(ax,[0 1]); yticks(ax, 0:0.25:1);

    hold(ax,'on');
    yline(ax, 1/3, '--', 'HandleVisibility','off');
    hold(ax,'off');
end

% Row 2: Fundamental best
for c = 1:3
    ax = nexttile(tl, 3+c);
    plotSessionScatter_ROLE(ax, colTabs{c}, 'pBestFV', sprintf('%s | Fundamental best', colNames{c}), ...
        treatMarker, condFaceColors);
    ylim(ax,[0 1]); yticks(ax, 0:0.25:1);

    hold(ax,'on');
    yline(ax, 1/3, '--', 'HandleVisibility','off');
    hold(ax,'off');
end

%% ========================== FIG 2: Condition means ============================
figure('Color','w','Name','FIG 2: Condition means (DP vs Fundamental; role-aware)');
tl2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl2,1);
plotConditionMeans_ROLE(ax1, sessionSummaryAll, 'pBestDP', 'Mean best-choice rate (DP)');
ylim(ax1,[0 1]); yticks(ax1, 0:0.25:1);

ax2 = nexttile(tl2,2);
plotConditionMeans_ROLE(ax2, sessionSummaryAll, 'pBestFV', 'Mean best-choice rate (Fundamental)');
ylim(ax2,[0 1]); yticks(ax2, 0:0.25:1);

%% ========================== FIG 3: Trial-wise curves ==========================
figure('Color','w','Name','FIG 3: Trial-wise best-choice curves (pooled; role-aware)');
ax = axes; hold(ax,'on');

Gt = findgroups(trialTableAll.trialNum);
pDP_t   = splitapply(@mean, trialTableAll.isBestDP, Gt);
pFund_t = splitapply(@mean, trialTableAll.isBestFundamental, Gt);
tKeys   = splitapply(@(x) x(1), trialTableAll.trialNum, Gt);

plot(ax, tKeys, pDP_t,   '-',  'LineWidth', 2.2, 'Color', [0 0 0]);
plot(ax, tKeys, pFund_t, '--', 'LineWidth', 2.2, 'Color', [0 0 0]);

xlabel(ax,'Trial'); ylabel(ax,'P(best choice)');
xlim(ax,[0 15]); xticks(ax, xTicksTrials);
ylim(ax,[0 1]); yticks(ax, 0:0.25:1);
legend(ax, {'DP benchmark','Fixed-fundamental benchmark'}, 'Location','southwest'); legend(ax,'boxoff');
title(ax,'Best-choice rate by trial number (pooled; all roles)');

hold(ax,'off');

%% ========================== Helper functions (standalone) =====================
function [actionCode, reportStr] = decodeActionsStandalone(optionRaw, trialNum, marketOrig, prePortfolio)
% Decode action vector into codes: BUY=1, HOLD=2, SELL=3.
% Robust behavior:
%   - If optionRaw is string/cellstr: match on tokens (buy/hold/sell or b/h/s).
%   - If numeric: infer mapping by trying permutations and picking the one
%     that best matches observed prePortfolio (when available).
%   - If prePortfolio missing: fall back to common conventions.

    reportStr = "unknown";

    if iscell(optionRaw) || isstring(optionRaw) || ischar(optionRaw)
        s = lower(strtrim(string(optionRaw)));
        actionCode = nan(size(s));
        actionCode(contains(s,"buy")  | s=="b") = 1;
        actionCode(contains(s,"hold") | s=="h") = 2;
        actionCode(contains(s,"sell") | s=="s") = 3;
        reportStr = "string-mapped";
        return;
    end

    x = double(optionRaw);
    actionCode = nan(size(x));

    u = unique(x(isfinite(x)));
    u = u(:)';

    if isempty(u)
        reportStr = "numeric-empty";
        return;
    end

    if all(ismember(u, [1 2 3]))
        if any(isfinite(prePortfolio))
            [bestMap, bestScore] = inferBestPermutationMapping(x, trialNum, marketOrig, prePortfolio, [1 2 3]);
            actionCode = remapByPairs(x, bestMap);
            reportStr = sprintf("numeric-123 perm=%s score=%.3f", mat2str(bestMap), bestScore);
        else
            actionCode = x;
            reportStr = "numeric-123 assumed";
        end
        return;
    end

    if all(ismember(u, [0 1 2]))
        if any(isfinite(prePortfolio))
            [bestMap, bestScore] = inferBestPermutationMapping(x, trialNum, marketOrig, prePortfolio, [0 1 2]);
            actionCode = remapByPairs(x, bestMap);
            reportStr = sprintf("numeric-012 perm=%s score=%.3f", mat2str(bestMap), bestScore);
        else
            map = [0 1 2; 1 2 3]; % 0=BUY,1=HOLD,2=SELL
            actionCode = remapByPairs(x, map);
            reportStr = "numeric-012 fallback(0->BUY,1->HOLD,2->SELL)";
        end
        return;
    end

    if numel(u) >= 3
        u3 = u(1:3);
        if any(isfinite(prePortfolio))
            [bestMap, bestScore] = inferBestPermutationMapping(x, trialNum, marketOrig, prePortfolio, u3);
            actionCode = remapByPairs(x, bestMap);
            reportStr = sprintf("numeric-generic perm=%s score=%.3f", mat2str(bestMap), bestScore);
        else
            map = [sort(u3); [1 2 3]];
            actionCode = remapByPairs(x, map);
            reportStr = sprintf("numeric-generic fallback(sorted %s)", mat2str(sort(u3)));
        end
        return;
    end

    reportStr = sprintf("numeric-insufficient-codes=%s", mat2str(u));
end

function [bestMap, bestScore] = inferBestPermutationMapping(optionNumeric, trialNum, marketOrig, prePortfolio, codeSet)
% Try all permutations assigning provided codeSet to {BUY,HOLD,SELL}.
% Choose the permutation that best matches observed prePortfolio consistency.
    permsIdx = perms(1:3);
    bestScore = -inf;
    bestMap = [codeSet; [1 2 3]];

    for k = 1:size(permsIdx,1)
        permCodes = codeSet(permsIdx(k,:));
        map = [permCodes; [1 2 3]];

        a = remapByPairs(optionNumeric, map);
        predPre = reconstructPrePortfolioFromActions(trialNum, marketOrig, a);

        valid = isfinite(prePortfolio) & isfinite(predPre);
        if ~any(valid), continue; end
        score = mean(predPre(valid) == round(prePortfolio(valid)));

        if score > bestScore
            bestScore = score;
            bestMap = map;
        end
    end
end

function y = remapByPairs(x, map)
% map: 2xK mapping where map(1,:) are source codes and map(2,:) are target codes
    y = nan(size(x));
    for j = 1:size(map,2)
        y(x == map(1,j)) = map(2,j);
    end
end

function pre = reconstructPrePortfolioFromActions(trialNum, marketOrig, actionCode)
% Reconstruct portfolio BEFORE choice assuming:
%   - portfolio resets to 0 at the start of each market instance
%   - can trade +/-1 per trial; no shorting
% Detect boundaries via trialNum==1.
    n = numel(actionCode);
    pre = nan(n,1);
    p = 0;
    for i = 1:n
        if i==1 || trialNum(i)==1
            p = 0;
        end
        pre(i) = p;

        a = actionCode(i);
        if a==1
            p = p + 1;
        elseif a==3
            p = max(p - 1, 0);
        end
    end
end

function opp = inferOpponentMonkey(subj)
% Infer opponent monkey ID/name from subject monkey.
% Assumes exactly two monkeys with IDs 1 and 2 (or labels 'M1' and 'M2').
% Extend this mapping if your IDs differ.
    s = strtrim(string(subj));

    val = str2double(s);
    if isfinite(val)
        if val == 1, opp = "2"; return; end
        if val == 2, opp = "1"; return; end
    end

    sU = upper(s);
    if sU == "M1", opp = "M2"; return; end
    if sU == "M2", opp = "M1"; return; end

    error("inferOpponentMonkey: cannot infer opponent from subject '%s'. Update mapping.", s);
end

function c2 = lightenColor(c1, frac)
% Move c1 toward white by fraction frac in [0,1].
    c1 = double(c1(:)');
    frac = max(0, min(1, frac));
    c2 = c1 + frac*(1 - c1);
end

function plotSessionScatter_ROLE(ax, tab, yvar, titleStr, treatMarker, condFaceColors)
% Scatter plot of session-level metric vs sessionIdx with:
%   - marker shape = Treatment
%   - marker face color = ConditionBase
%   - ROLE cue (Live only): Opp trials shown as lightened Live color
    hold(ax,'on');
    tab = tab(~isnan(tab.(yvar)),:);

    for i = 1:height(tab)
        tr   = string(tab.Treatment(i));
        cb   = string(tab.ConditionBase(i));
        role = string(tab.Role(i));

        if ~isKey(treatMarker, tr), continue; end
        if ~isKey(condFaceColors, cb), continue; end

        mk = treatMarker(tr);
        fc = condFaceColors(cb);

        if cb=="Live" && role=="Opp"
            fc = lightenColor(fc, 0.55);
        end

        plot(ax, tab.sessionIdx(i), tab.(yvar)(i), mk, ...
            'MarkerFaceColor', fc, 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    end

    xlabel(ax,'Session # (within monkey)');
    ylabel(ax, yvar);
    title(ax, titleStr);
    grid(ax,'on'); box(ax,'off');

    % simple legend note (kept minimal on purpose)
    hold(ax,'off');
end

function plotConditionMeans_ROLE(ax, sessionSummaryAll, yvar, titleStr)
% Bar plot: mean across sessions for each (Treatment x ConditionBase x Role).
% This makes the Live split explicit without you having to filter first.
    hold(ax,'on');

    G = findgroups(sessionSummaryAll.Treatment, sessionSummaryAll.ConditionBase, sessionSummaryAll.Role);
    mu = splitapply(@mean, sessionSummaryAll.(yvar), G);
    se = splitapply(@(x) std(x)/sqrt(numel(x)), sessionSummaryAll.(yvar), G);

    trKey = splitapply(@(x) x(1), sessionSummaryAll.Treatment, G);
    cbKey = splitapply(@(x) x(1), sessionSummaryAll.ConditionBase, G);
    rKey  = splitapply(@(x) x(1), sessionSummaryAll.Role, G);

    T = table(trKey, cbKey, rKey, mu, se, 'VariableNames', {'Treatment','ConditionBase','Role','Mean','SE'});

    xLabels = strings(height(T),1);
    for i = 1:height(T)
        xLabels(i) = string(T.Treatment(i)) + " | " + string(T.ConditionBase(i)) + " | " + string(T.Role(i));
    end

    x = 1:height(T);
    bar(ax, x, T.Mean);
    errorbar(ax, x, T.Mean, T.SE, 'k.', 'LineWidth', 1);

    title(ax, titleStr);
    ylabel(ax, yvar);
    xticks(ax, x);
    xticklabels(ax, xLabels);
    xtickangle(ax, 45);
    grid(ax,'on'); box(ax,'off');
    hold(ax,'off');
end
