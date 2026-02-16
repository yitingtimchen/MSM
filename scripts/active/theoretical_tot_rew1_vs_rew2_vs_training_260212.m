%% Total THEORETICAL Reward over Training Sessions per Monkey (ROLE-AWARE LIVE)
% -------------------------------------------------------------------------
% ROLE-AWARE version of theoretical_tot_rew1_vs_rew2_vs_training_260209
%
% Key additions:
%   - Live condition split into Role="Subj" and Role="Opp" using *_Opp fields
%   - Opp role monkey identity inferred (because etab.monkey is always subject)
%   - Role carried into sessionSummaryAll
%   - sessionIdx shared across roles (same x-position within a session)
%   - Live Opp points lightened in plots
%
% IMPORTANT (per Tim request):
%   - DOES NOT change your treatment/condition ordering.
%       Treatments plotted in setNames order: baseline -> OT -> Saline
%       Conditions plotted in baseCondList order: AI -> Replay -> Decoy -> Live
%
% NOTE ABOUT FIELDNAMES
%   This script assumes your packs use the following *_Opp fieldnames:
%       priceBuyOpp, priceSellOpp, optionOpp, trialNumOpp, prePortfolioOpp, divPerShareOpp
%   If your pack uses different names, edit ONLY the "FIELD MAP" section.

clear; clc;
close all;

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL'; % <-- EDIT

cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames = {'baseline','OT','Saline'};          % <<< YOUR treatment order (do not change)
baseCondList = {'AI','Replay','Decoy','Live'};  % <<< YOUR condition order (do not change)

% Live role handling:
%   "subj_only" -> include only Subj stream in Live
%   "separate"  -> include Subj and Opp as separate role conditions
%   "combined"  -> include both streams but pool them as one Live role;
%                  session indexing uses role-qualified session keys so
%                  effective session count doubles for Live per monkey.
LIVE_ROLE_MODE = "separate";  % "subj_only" | "separate" | "combined"

% Which metric to plot?
%   'choiceGross' : allowance + tradeNet
%   'tradeNet'    : trade-only cashflow (BUY=-buyPrice, SELL=+sellPrice, HOLD=0)
%   'div'         : dividend reward only
%   'totalGross'  : choiceGross + div
%   'totalNet'    : tradeNet + div
%   'fracDiv'     : sumDiv / sumTotalGross
plotMetric = 'choiceGross';  % <-- EDIT

% Allowance constant used to construct choiceGross
allowance_dollars = 80;

% Dividend handling
useRealizedDivIfAvailable = true;
muDiv_fallback_dollars = mean([0 0.8 2.8 6]); % 2.4

% Optional linear scaling for plotting (does NOT change comparisons)
plot_in_pixels = false;
pxPerDollar    = 2.125;

%% ------------------------ REQUIRED ETAB FIELD MAP --------------------------
% If your packs use different names, change these ONLY.

% --- Subject stream fields (standard) ---
FIELD_priceBuy      = 'priceBuy';
FIELD_priceSell     = 'priceSell';
FIELD_option        = 'option';
FIELD_prePortfolio  = 'prePortfolio';   % if missing, we reconstruct
FIELD_trialNum      = 'trialNum';       % used to reset portfolio at trial==1
FIELD_divPerShare   = 'divPerShare';    % optional (realized dividend per share)

% --- Live opponent stream fields (only used when baseCond == "Live") ---
FIELD_priceBuyOpp     = 'priceBuyOpp';
FIELD_priceSellOpp    = 'priceSellOpp';
FIELD_optionOpp       = 'optionOpp';
FIELD_prePortfolioOpp = 'prePortfolioOpp';  % optional
FIELD_trialNumOpp     = 'trialNumOpp';
FIELD_divPerShareOpp  = 'divPerShareOpp';   % optional; many packs do NOT have this

% Session identity fields
FIELD_year          = 'year';
FIELD_month         = 'month';
FIELD_day           = 'day';
FIELD_SessionOfDay  = 'SessionOfDay';
FIELD_monkey        = 'monkey';

%% --------------------- MASTER SESSION SUMMARY TABLE ------------------------
% One row per (session x role):
sessionSummaryAll = table();

%% ======================= MAIN LOOP OVER TREATMENT SETS =====================
for s = 1:numel(cond_sets)
    conds     = cond_sets{s};
    treatName = setNames{s};

    fprintf('\n=== Treatment: %s ===\n', treatName);
    sessLocal = table();

    for c = 1:numel(conds)
        cname   = conds{c};
        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));

        if ~exist(matFile,'file')
            warning('Missing condition pack: %s', matFile);
            continue;
        end

        % Parse base condition from condition name (final token)
        parts = split(string(cname));
        baseCond = string(parts(end));

        tmp = load(matFile, 'C');
        C   = tmp.C;

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            % ---------- Required identity fields (always needed) ----------
            neededID = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey};
            if ~all(ismember(neededID, etab.Properties.VariableNames))
                warning('Missing required identity fields in condition %s, file %d; skipping.', cname, f);
                continue;
            end

            % --- Session ID ---
            y   = etab.(FIELD_year){1}(:);
            m   = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(m) || isempty(d) || isempty(sod)
                continue;
            end

            y0 = y(1); m0 = m(1); d0 = d(1); sod0 = sod(1);
            sessionStr  = string(sprintf('%04d-%02d-%02d-%02d', y0, m0, d0, sod0));
            sessionDate = datetime(y0, m0, d0);

            % --- Subject monkey ID (etab.monkey is ALWAYS the subject) ---
            if isempty(etab.(FIELD_monkey){1})
                warning('Missing monkey field; skipping a file in %s.', cname);
                continue;
            end
            subjMonkeyVal = etab.(FIELD_monkey){1}(1);
            subjMonkeyStr = string(subjMonkeyVal);

            % ---------- Role split logic ----------
            if baseCond == "Live"
                if LIVE_ROLE_MODE == "subj_only"
                    roleList = "Subj";
                else
                    roleList = ["Subj","Opp"];
                end
            else
                roleList = "Subj";
            end

            for role = roleList

                % ---------- Role-specific field selection ----------
                if role == "Subj"
                    f_priceBuy  = FIELD_priceBuy;
                    f_priceSell = FIELD_priceSell;
                    f_option    = FIELD_option;
                    f_trialNum  = FIELD_trialNum;
                    f_prePort   = FIELD_prePortfolio;
                    f_div       = FIELD_divPerShare;

                    actorMonkeyStr = subjMonkeyStr;
                else
                    f_priceBuy  = FIELD_priceBuyOpp;
                    f_priceSell = FIELD_priceSellOpp;
                    f_option    = FIELD_optionOpp;
                    f_trialNum  = FIELD_trialNumOpp;
                    f_prePort   = FIELD_prePortfolioOpp;
                    f_div       = FIELD_divPerShareOpp;

                    actorMonkeyStr = inferOpponentMonkey(subjMonkeyStr);
                end
                if ~isMonkey12(actorMonkeyStr), continue; end

                % ---------- Required fields for THIS role ----------
                neededRole = {f_priceBuy, f_priceSell, f_option, f_trialNum};

                % For Live Opp, if the pack is missing Opp fields, skip cleanly.
                if ~all(ismember(neededRole, etab.Properties.VariableNames))
                    continue;
                end

                if isempty(etab.(f_priceBuy){1}) || isempty(etab.(f_priceSell){1}) || ...
                   isempty(etab.(f_option){1})   || isempty(etab.(f_trialNum){1})
                    continue;
                end

                % ---------- Trial vectors (role stream) ----------
                priceBuy  = double(etab.(f_priceBuy){1}(:));
                priceSell = double(etab.(f_priceSell){1}(:));
                optionRaw = etab.(f_option){1}(:);
                trialNum  = double(etab.(f_trialNum){1}(:));

                % prePortfolio: if missing, we'll reconstruct
                hasPre = ismember(f_prePort, etab.Properties.VariableNames) && ~isempty(etab.(f_prePort){1});
                if hasPre
                    prePortfolio = double(etab.(f_prePort){1}(:));
                else
                    prePortfolio = nan(size(trialNum));
                end

                % dividend per share (realized): optional and often absent for Opp
                hasDiv = ismember(f_div, etab.Properties.VariableNames) && ~isempty(etab.(f_div){1});
                if hasDiv
                    divPerShare = double(etab.(f_div){1}(:));
                else
                    divPerShare = nan(size(trialNum));
                end

                % ---------- Align lengths ----------
                nT = min([numel(priceBuy), numel(priceSell), numel(optionRaw), numel(trialNum), numel(prePortfolio), numel(divPerShare)]);
                priceBuy     = priceBuy(1:nT);
                priceSell    = priceSell(1:nT);
                optionRaw    = optionRaw(1:nT);
                trialNum     = trialNum(1:nT);
                prePortfolio = prePortfolio(1:nT);
                divPerShare  = divPerShare(1:nT);

                % ---------- Basic validity for prices and trialNum ----------
                keep = isfinite(priceBuy) & isfinite(priceSell) & isfinite(trialNum);
                if ~any(keep)
                    continue;
                end

                priceBuy     = priceBuy(keep);
                priceSell    = priceSell(keep);
                optionRaw    = optionRaw(keep);
                trialNum     = trialNum(keep);
                prePortfolio = prePortfolio(keep);
                divPerShare  = divPerShare(keep);

                % ---------- Decode action into BUY/HOLD/SELL codes ----------
                actionCode = decodeOptionToAction(optionRaw);

                % If actionCode has NaNs, drop those trials (unscorable)
                validAction = isfinite(actionCode);
                if ~any(validAction)
                    continue;
                end

                priceBuy     = priceBuy(validAction);
                priceSell    = priceSell(validAction);
                actionCode   = actionCode(validAction);
                trialNum     = trialNum(validAction);
                prePortfolio = prePortfolio(validAction);
                divPerShare  = divPerShare(validAction);

                % ---------- Reconstruct prePortfolio if missing ----------
                if ~hasPre || all(~isfinite(prePortfolio))
                    prePortfolio = reconstructPrePortfolio(trialNum, actionCode);
                end

                % Clamp portfolio (no shorting) + integerize
                prePortfolio = max(0, round(prePortfolio));

                % ---------- Compute postPortfolio (dividend is post-action) ----------
                postPortfolio = prePortfolio;
                postPortfolio(actionCode==1) = prePortfolio(actionCode==1) + 1;               % BUY
                postPortfolio(actionCode==3) = max(prePortfolio(actionCode==3) - 1, 0);       % SELL

                % Diagnostic: SELL at p=0 (should be rare/none)
                nSellAtZero = sum((actionCode==3) & (prePortfolio==0));

                % ---------- Immediate trade-only cashflow (net) ----------
                tradeNet = zeros(size(actionCode));
                tradeNet(actionCode==1) = -priceBuy(actionCode==1);   % BUY
                tradeNet(actionCode==2) = 0;                          % HOLD
                tradeNet(actionCode==3) = +priceSell(actionCode==3);  % SELL

                % ---------- Choice reward (gross = allowance + tradeNet) ----------
                choiceGross = allowance_dollars + tradeNet;

                % ---------- Dividend reward ----------
                if useRealizedDivIfAvailable && any(isfinite(divPerShare))
                    divUsed = divPerShare;
                    divUsed(~isfinite(divUsed)) = muDiv_fallback_dollars;
                else
                    divUsed = repmat(muDiv_fallback_dollars, size(tradeNet));
                end
                divReward = divUsed .* postPortfolio;

                % ---------- Totals ----------
                totalGross = choiceGross + divReward;
                totalNet   = tradeNet   + divReward;

                % Valid trials for session-level sums: require all components finite
                valid = isfinite(tradeNet) & isfinite(choiceGross) & isfinite(divReward) & isfinite(totalGross) & isfinite(totalNet);
                if ~any(valid)
                    continue;
                end

                % ---------- Session-level aggregates ----------
                sumTradeNet    = sum(tradeNet(valid));
                sumChoiceGross = sum(choiceGross(valid));
                sumDiv         = sum(divReward(valid));
                sumTotalGross  = sum(totalGross(valid));
                sumTotalNet    = sum(totalNet(valid));

                nTrialsValid   = sum(valid);

                meanTradeNet    = sumTradeNet    / nTrialsValid;
                meanChoiceGross = sumChoiceGross / nTrialsValid;
                meanDiv         = sumDiv         / nTrialsValid;
                meanTotalGross  = sumTotalGross  / nTrialsValid;
                meanTotalNet    = sumTotalNet    / nTrialsValid;

                if sumTotalGross ~= 0
                    fracDivGross = sumDiv / sumTotalGross;
                else
                    fracDivGross = NaN;
                end

                % Portfolio trend metrics (shares; never scaled)
                meanPrePortfolio  = mean(prePortfolio(valid));
                meanPostPortfolio = mean(postPortfolio(valid));
                endPostPortfolio  = postPortfolio(find(valid, 1, 'last'));

                % Optional scaling for stored reward units
                if plot_in_pixels
                    sumTradeNet    = sumTradeNet    * pxPerDollar;
                    sumChoiceGross = sumChoiceGross * pxPerDollar;
                    sumDiv         = sumDiv         * pxPerDollar;
                    sumTotalGross  = sumTotalGross  * pxPerDollar;
                    sumTotalNet    = sumTotalNet    * pxPerDollar;

                    meanTradeNet    = meanTradeNet    * pxPerDollar;
                    meanChoiceGross = meanChoiceGross * pxPerDollar;
                    meanDiv         = meanDiv         * pxPerDollar;
                    meanTotalGross  = meanTotalGross  * pxPerDollar;
                    meanTotalNet    = meanTotalNet    * pxPerDollar;
                end

                % ---------- Append session-row (role-aware) ----------
                rowRole = string(role);
                analysisSession = sessionStr;
                if baseCond == "Live" && LIVE_ROLE_MODE == "combined"
                    rowRole = "Combined";
                    analysisSession = sessionStr + "-" + string(role);
                end

                newRow = table( ...
                    string(treatName), string(cname), string(baseCond), rowRole, string(actorMonkeyStr), sessionStr, analysisSession, sessionDate, ...
                    sumTradeNet, sumChoiceGross, sumDiv, sumTotalGross, sumTotalNet, ...
                    meanTradeNet, meanChoiceGross, meanDiv, meanTotalGross, meanTotalNet, ...
                    meanPrePortfolio, meanPostPortfolio, endPostPortfolio, ...
                    fracDivGross, nTrialsValid, nSellAtZero, ...
                    'VariableNames', { ...
                        'Treatment','ConditionRaw','ConditionBase','Role','Monkey','Session','AnalysisSession','SessionDate', ...
                        'sumTradeNet','sumChoiceGross','sumDiv','totalGross','totalNet', ...
                        'meanTradeNet','meanChoiceGross','meanDiv','meanTotalGross','meanTotalNet', ...
                        'meanPrePortfolio','meanPostPortfolio','endPostPortfolio', ...
                        'fracDivGross','nTrialsValid','nSellAtZero' ...
                    } ...
                );

                sessLocal = [sessLocal; newRow]; %#ok<AGROW>
            end % role
        end % eventTables
    end % conditions

    if isempty(sessLocal)
        warning('No sessions for treatment %s; skipping.', treatName);
        continue;
    end

    sessionSummaryAll = [sessionSummaryAll; sessLocal]; %#ok<AGROW>
end % treatments

%% ------------------------ CHECK + FACTORS + TRAINING INDEX -------------------
if isempty(sessionSummaryAll)
    warning('No sessions found in any treatment. Nothing to analyze.');
    return;
end

% FORCE orders (do not let MATLAB choose)
treatOrder = setNames;
sessionSummaryAll.Treatment = categorical(string(sessionSummaryAll.Treatment), treatOrder);
sessionSummaryAll.Treatment = removecats(sessionSummaryAll.Treatment);

sessionSummaryAll.ConditionBase = categorical(string(sessionSummaryAll.ConditionBase), baseCondList);
sessionSummaryAll.ConditionBase = removecats(sessionSummaryAll.ConditionBase);

sessionSummaryAll.Monkey = categorical(string(sessionSummaryAll.Monkey));
sessionSummaryAll.Role   = categorical(string(sessionSummaryAll.Role), ["Subj","Opp","Combined"]);

% sessionIdx: chronological order within each monkey.
% In LIVE_ROLE_MODE=="combined", role-qualified session keys double Live density.
sessTab = unique(sessionSummaryAll(:, {'Monkey','AnalysisSession'}));
[sessTabSorted, ~] = sortrows(sessTab, {'Monkey','AnalysisSession'});

sessionIdxSorted = zeros(height(sessTabSorted),1);
monkeys = categories(sessTabSorted.Monkey);

for im = 1:numel(monkeys)
    mID  = monkeys{im};
    mask = sessTabSorted.Monkey == mID;
    sessionIdxSorted(mask) = (1:nnz(mask)).';
end

sessTabSorted.sessionIdx = sessionIdxSorted;

% join back
sessionSummaryAll = outerjoin(sessionSummaryAll, sessTabSorted, 'Keys', {'Monkey','AnalysisSession'}, 'MergeKeys', true);
sessionSummaryAll.Session = sessionSummaryAll.AnalysisSession;

%% ------------------------ PLOT SETTINGS (ORDER-PRESERVING) -------------------
% Use YOUR fixed orders for plotting (and only keep those that exist in data)
treatCats = treatOrder(ismember(treatOrder, categories(sessionSummaryAll.Treatment)));
baseCats  = baseCondList(ismember(baseCondList, categories(sessionSummaryAll.ConditionBase)));

baseColors = lines(numel(baseCats));  % color assignment matches baseCats order
markerList = {'o','s','^','d','v','>'}; % shapes per treatment (in treatCats order)

%% ------------------------ SELECT FIELD/LABELS ------------------------------
if plot_in_pixels
    unitStr = ' (px)';
else
    unitStr = ' ($)';
end

switch lower(plotMetric)
    case 'choicegross'
        yField   = 'sumChoiceGross';
        yLabel   = ['Total choiceGross per session' unitStr];
        fig1Name = 'ChoiceGross vs training (by session order)';
        fig2Name = 'ChoiceGross vs training (by date)';
        sg1Title = 'ChoiceGross over training sessions (sessionIdx, all treatments)';
        sg2Title = 'ChoiceGross over training sessions (session date, all treatments)';
    case 'tradenet'
        yField   = 'sumTradeNet';
        yLabel   = ['Total tradeNet per session' unitStr];
        fig1Name = 'TradeNet vs training (by session order)';
        fig2Name = 'TradeNet vs training (by date)';
        sg1Title = 'TradeNet over training sessions (sessionIdx, all treatments)';
        sg2Title = 'TradeNet over training sessions (session date, all treatments)';
    case 'div'
        yField   = 'sumDiv';
        yLabel   = ['Total dividend reward per session' unitStr];
        fig1Name = 'Dividend reward vs training (by session order)';
        fig2Name = 'Dividend reward vs training (by date)';
        sg1Title = 'Dividend reward over training sessions (sessionIdx, all treatments)';
        sg2Title = 'Dividend reward over training sessions (session date, all treatments)';
    case 'totalnet'
        yField   = 'totalNet';
        yLabel   = ['TotalNet (tradeNet + div) per session' unitStr];
        fig1Name = 'TotalNet vs training (by session order)';
        fig2Name = 'TotalNet vs training (by date)';
        sg1Title = 'TotalNet over training sessions (sessionIdx, all treatments)';
        sg2Title = 'TotalNet over training sessions (session date, all treatments)';
    case 'fracdiv'
        yField   = 'fracDivGross';
        yLabel   = 'Dividend fraction of totalGross per session';
        fig1Name = 'Dividend fraction vs training (by session order)';
        fig2Name = 'Dividend fraction vs training (by date)';
        sg1Title = 'Dividend fraction over training sessions (sessionIdx, all treatments)';
        sg2Title = 'Dividend fraction over training sessions (session date, all treatments)';
    otherwise % totalGross
        yField   = 'totalGross';
        yLabel   = ['TotalGross (choiceGross + div) per session' unitStr];
        fig1Name = 'TotalGross vs training (by session order)';
        fig2Name = 'TotalGross vs training (by date)';
        sg1Title = 'TotalGross over training sessions (sessionIdx, all treatments)';
        sg2Title = 'TotalGross over training sessions (session date, all treatments)';
end

%% ------------------------ PLOTTING ----------------------------------------
fprintf('\nPlotting %s vs training session and date per monkey (ROLE-AWARE; order-preserving)...\n', yField);

figure('Name', fig1Name, 'Color','w');
plot_metric_vs_x_ROLE(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, markerList, ...
    yField, yLabel, 'idx', sg1Title, LIVE_ROLE_MODE);

figure('Name', fig2Name, 'Color','w');
plot_metric_vs_x_ROLE(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, markerList, ...
    yField, yLabel, 'date', sg2Title, LIVE_ROLE_MODE);

%% ------------------------ PORTFOLIO TREND PLOTS ---------------------------
portfolioField = 'meanPostPortfolio';
portfolioLabel = 'Mean post-choice portfolio per session (shares)';

figure('Name', 'Portfolio trend vs training (by session order)', 'Color','w');
plot_metric_vs_x_ROLE(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, markerList, ...
    portfolioField, portfolioLabel, 'idx', ...
    'Mean post-choice portfolio over training sessions (sessionIdx, all treatments)', LIVE_ROLE_MODE);

figure('Name', 'Portfolio trend vs training (by date)', 'Color','w');
plot_metric_vs_x_ROLE(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, markerList, ...
    portfolioField, portfolioLabel, 'date', ...
    'Mean post-choice portfolio over training sessions (session date, all treatments)', LIVE_ROLE_MODE);

%% ========================= LOCAL FUNCTIONS =================================

function actionCode = decodeOptionToAction(optionRaw)
% Map optionRaw -> {BUY,HOLD,SELL} codes:
%   BUY  = 1
%   HOLD = 2
%   SELL = 3
%
% If your numeric encoding differs, edit this function.

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

    if all(ismember(u, [1 2 3]))
        actionCode = x; % assume 1=BUY,2=HOLD,3=SELL
        return;
    end

    if all(ismember(u, [0 1 2]))
        % assume 0=BUY,1=HOLD,2=SELL
        actionCode(x==0) = 1;
        actionCode(x==1) = 2;
        actionCode(x==2) = 3;
        return;
    end

    warning('decodeOptionToAction: unrecognized numeric codes %s. Update mapping.', mat2str(u(:)'));
end

function prePortfolio = reconstructPrePortfolio(trialNum, actionCode)
% Reconstruct prePortfolio from actions:
%   - resets to 0 at trialNum==1
%   - BUY increments, SELL decrements (clipped at 0), HOLD no change

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

function opp = inferOpponentMonkey(subj)
% Infer opponent monkey from subject monkey.
% Assumes exactly two monkeys with IDs 1 and 2 (or labels M1 and M2).
% Update if your IDs differ.

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

function ok = isMonkey12(monkeyVal)
% Keep only rows where monkey is 1/2 (or M1/M2 string equivalent).
    s = upper(strtrim(string(monkeyVal)));
    v = str2double(s);
    ok = (isfinite(v) && any(v == [1 2])) || any(s == ["M1","M2"]);
end

function plot_metric_vs_x_ROLE(sessionSummaryAll, monkeys, treatCats, baseCats, baseColors, ...
                              markerList, yField, yLabel, xMode, sgTitle, LIVE_ROLE_MODE)
% Role-aware scatter + fit stats.
%   - Marker shape encodes Treatment (in treatCats order)
%   - Color encodes ConditionBase (in baseCats order)
%   - For baseCond=="Live": Subj (dark) and Opp (lightened) separately

nMonk = numel(monkeys);

for im = 1:nMonk
    mID   = monkeys{im};
    rowsM = sessionSummaryAll.Monkey == mID;

    subplot(1, nMonk, im); hold on;

    x_all = [];
    y_all = [];
    idx_all = [];

    baseLegendHandles = gobjects(numel(baseCats),1);

    for it = 1:numel(treatCats)
        tName  = treatCats{it};
        rowsT  = sessionSummaryAll.Treatment == tName;
        marker = markerList{min(it, numel(markerList))};

        for ib = 1:numel(baseCats)
            bName = baseCats{ib};
            rowsB = sessionSummaryAll.ConditionBase == bName;

            if strcmp(bName, 'Live') && LIVE_ROLE_MODE == "separate"
                for rr = ["Subj","Opp"]
                    rowsR = sessionSummaryAll.Role == rr;
                    rows = rowsM & rowsT & rowsB & rowsR;
                    if ~any(rows), continue; end

                    switch lower(xMode)
                        case 'idx'
                            x = double(sessionSummaryAll.sessionIdx(rows));
                        case 'date'
                            x = datenum(sessionSummaryAll.SessionDate(rows));
                        otherwise
                            error('Unknown xMode: %s', xMode);
                    end
                    y = sessionSummaryAll.(yField)(rows);

                    fc = baseColors(ib,:);
                    if rr == "Opp"
                        fc = lightenColor(fc, 0.55);
                    end

                    h = scatter(x, y, 50, ...
                        'Marker', marker, ...
                        'MarkerFaceColor', fc, ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerFaceAlpha', 0.7);

                    x_all   = [x_all; x(:)];
                    y_all   = [y_all; y(:)];
                    idx_all = [idx_all; find(rows)];

                    if ~isgraphics(baseLegendHandles(ib))
                        baseLegendHandles(ib) = h;
                    end
                end
            else
                rows = rowsM & rowsT & rowsB;
                if ~any(rows), continue; end

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
                idx_all = [idx_all; find(rows)];

                if ~isgraphics(baseLegendHandles(ib))
                    baseLegendHandles(ib) = h;
                end
            end
        end
    end

    statsTxt = "";

    if numel(x_all) > 1
        xvec = x_all(:);
        yvec = y_all(:);
        mdl  = fitlm(xvec, yvec);

        slope = mdl.Coefficients.Estimate(2);
        R2    = mdl.Rsquared.Ordinary;
        pval  = mdl.Coefficients.pValue(2);

        xfit = linspace(min(xvec), max(xvec), 100)';
        [yfit, yCI] = predict(mdl, xfit, 'Alpha', 0.05);

        fill([xfit; flipud(xfit)], [yCI(:,1); flipud(yCI(:,2))], ...
             [0.7 0.7 0.7], 'EdgeColor','none', 'FaceAlpha', 0.3);

        plot(xfit, yfit, 'k-', 'LineWidth', 1.5);

        residStd    = mdl.Residuals.Standardized;
        outlierMask = abs(residStd) > 2.5;

        if any(outlierMask)
            x_out   = xvec(outlierMask);
            y_out   = yvec(outlierMask);
            idx_out = idx_all(outlierMask);

            plot(x_out, y_out, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, ...
                 'MarkerFaceColor', 'none');

            for k = 1:numel(x_out)
                switch lower(xMode)
                    case 'idx'
                        lbl = sprintf('%d', sessionSummaryAll.sessionIdx(idx_out(k)));
                    case 'date'
                        lbl = char(sessionSummaryAll.Session(idx_out(k)));
                end
                text(x_out(k), y_out(k), ['  ' lbl], ...
                     'HorizontalAlignment','left', ...
                     'VerticalAlignment','bottom', ...
                     'FontSize', 8, 'Color', 'r');
            end
        end

        statsTxt = sprintf('slope = %.3g    R^2 = %.3f    p = %.3g', slope, R2, pval);
    end

    text(0.02, 0.02, statsTxt, ...
         'Units','normalized', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','bottom', ...
         'FontSize', 10, ...
         'Interpreter','none', ...
         'BackgroundColor','none', ...
         'Margin', 2);

    switch lower(xMode)
        case 'idx'
            xlabel('sessionIdx (within monkey, all treatments)');
        case 'date'
            xlabel('Session date');
            if ~isempty(x_all)
                ax = gca;
                ax.XLim = [min(x_all) max(x_all)];
                xt = ax.XTick;
                xt_dt = datetime(xt, 'ConvertFrom', 'datenum');
                ax.XTickLabel = cellstr(datestr(xt_dt, 'yyyy-mm-dd'));
                ax.XTickLabelRotation = 45;
            end
    end

    ylabel(yLabel);
    title(sprintf('Monkey %s', char(mID)), 'Interpreter','none');

    validMask = isgraphics(baseLegendHandles);

    % Treatment legend stubs (in treatCats order)
    treatLegendHandles = gobjects(numel(treatCats),1);
    for it = 1:numel(treatCats)
        mk = markerList{min(it, numel(markerList))};
        treatLegendHandles(it) = scatter(nan, nan, 50, ...
            'Marker', mk, ...
            'MarkerFaceColor', [0.85 0.85 0.85], ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 1.0);
    end

    baseLabels  = "Condition: "  + string(baseCats(validMask));
    treatLabels = "Treatment: " + string(treatCats(:));

    legHandles  = [baseLegendHandles(validMask); treatLegendHandles(:)];
    legLabels   = cellstr([baseLabels(:); treatLabels(:)]);

    legend(legHandles, legLabels, 'Location','northwest', 'Box','off','NumColumns',2);

    hold off;
end

if exist('sgtitle','file')
    sgtitle(sgTitle, 'Interpreter','none');
end
end
