%% MSM Step 1: Action-rate trends and portfolio trends over training
% -------------------------------------------------------------------------
% GOAL (Step 1)
%   Per session, compute:
%     - P(BUY), P(HOLD), P(SELL)
%     - mean (post-pre) portfolio = net share flow per trial (P(BUY)-P(SELL))
%     - mean postPortfolio (post-choice, i.e., after applying action)
%     - end-of-market postPortfolio (mean across markets within session, trialNum==15)
%
%   Then plot these metrics vs sessionIdx (chronological session order within monkey),
%   using the same condition-pack loading structure as your other analyses.
%
% WHY THIS STEP
%   It tells you whether "portfolio increases with training" is driven by:
%     - more BUY, or
%     - less SELL, or
%     - both
%
% INPUTS (from condition packs: C.eventTables)
%   Required (per etab):
%     - option       : subject action code/label (BUY/HOLD/SELL)
%     - trialNum     : trial number within market (1..15)
%     - year, month, day, SessionOfDay
%     - monkey
%   Optional:
%     - prePortfolio : inventory before choice (if missing, reconstructed)
%
% OUTPUTS
%   - sessionSummaryAll (table)
%   - Plots:
%       FIG 1: Action probabilities vs sessionIdx (per monkey)
%       FIG 2: Portfolio summaries vs sessionIdx (per monkey)
%
% NOTE
%   This script does NOT use fb1/fb2 and does NOT use DP labels.
%   It is purely descriptive of what actions are happening over training.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL'; % <-- EDIT

cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames     = {'baseline', 'OT', 'Saline'};
baseCondList = {'AI','Replay','Decoy','Live'};
% cond_sets = { ...
%     {'Live'}; ...
%     {'OT Live'}; ...
%     {'Saline Live'} ...
% };
% setNames     = {'baseline','OT','Saline'};
% baseCondList = {'AI','Replay','Decoy','Live'};

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
EXCLUDE_SESSIONS = { ...
    % '2018-07-07-01', ...
    % '2018-07-20-02', ...
    % '2018-07-21-01', ...
    % '2018-07-06-02', ...
    % '2018-07-06-01', ...
    % '2018-08-01-03', ...
    % '2018-08-22-02', ...
    % '2018-06-07-01', ...
    % '2018-06-07-02', ...
    % '2018-06-08-01', ...
    % '2018-06-08-02', ...
    % '2018-06-08-03', ...
};
EXCLUDE_SESSIONS = string(EXCLUDE_SESSIONS(:));
EXCLUDE_SESSIONS = EXCLUDE_SESSIONS(strlength(EXCLUDE_SESSIONS) > 0);
excludedSessSeen = strings(0,1);

% Live role handling:
%   "subj_only" -> include only Subj stream in Live
%   "separate"  -> include Subj and Opp as separate role conditions
%   "combined"  -> include both streams but pool them as one Live role;
%                  session indexing uses role-qualified session keys so
%                  effective session count doubles for Live per monkey.
LIVE_ROLE_MODE = "separate";  % "subj_only" | "separate" | "combined"

% Plot encodings (match your other scripts)
treatOrder  = {'baseline','OT','Saline'};
markerList  = {'o','s','^','d','v','>'}; % first 3 used

condFaceColors = containers.Map();
condFaceColors('AI')     = [0.0000 0.4470 0.7410];
condFaceColors('Replay') = [0.8500 0.3250 0.0980];
condFaceColors('Decoy')  = [0.9290 0.6940 0.1250];
condFaceColors('Live')   = [0.4940 0.1840 0.5560];

%% ------------------------ REQUIRED ETAB FIELD MAP --------------------------
% If your packs use different names, change these ONLY.
FIELD_option       = 'option';
FIELD_trialNum     = 'trialNum';
FIELD_prePortfolio = 'prePortfolio';   % optional

% Opponent fields (Live only; role = Opp)
FIELD_optionOpp       = 'optionOpp';
FIELD_trialNumOpp     = 'trialNumOpp';
FIELD_prePortfolioOpp = 'prePortfolioOpp'; % optional
FIELD_monkeyOpp       = 'monkeyOpp';

FIELD_year         = 'year';
FIELD_month        = 'month';
FIELD_day          = 'day';
FIELD_SessionOfDay = 'SessionOfDay';
FIELD_monkey       = 'monkey';

%% --------------------- MASTER SESSION SUMMARY TABLE ------------------------
% One row per session
sessionSummaryAll = table();

fprintf('\n=== BUILDING sessionSummaryAll (session-level metrics) ===\n');

%% ======================= MAIN LOOP OVER TREATMENT SETS =====================
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

        % Parse base condition from condition name (AI/Replay/Decoy/Live at end)
        tok = regexp(cname, '(AI|Replay|Decoy|Live)$', 'tokens', 'once');
        if isempty(tok)
            baseCond = string(cname);
        else
            baseCond = string(tok{1});
        end

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            % Required session fields (role-specific action fields handled below)
            neededSession = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey};
            if ~all(ismember(neededSession, etab.Properties.VariableNames))
                warning('Missing required session fields; skipping a file in %s.', cname);
                continue;
            end

            % Session identity
            y   = etab.(FIELD_year){1}(:);
            m   = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(m) || isempty(d) || isempty(sod), continue; end

            sessionStr = string(sprintf('%04d-%02d-%02d-%02d', y(1), m(1), d(1), sod(1)));

            % Exclusions
            if ~isempty(EXCLUDE_SESSIONS) && any(sessionStr == EXCLUDE_SESSIONS)
                excludedSessSeen(end+1,1) = sessionStr; %#ok<AGROW>
                continue;
            end

            % Decide which roles to extract (Live only)
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

                % Role-specific field selection
                if role == "Subj"
                    optField = FIELD_option;
                    trlField = FIELD_trialNum;
                    preField = FIELD_prePortfolio;
                else
                    optField = FIELD_optionOpp;
                    trlField = FIELD_trialNumOpp;
                    preField = FIELD_prePortfolioOpp;
                end

                neededRole = {optField, trlField};
                if ~all(ismember(neededRole, etab.Properties.VariableNames))
                    % For Live Opp, skip cleanly if Opp fields absent
                    continue;
                end

                if isempty(etab.(optField){1}) || isempty(etab.(trlField){1})
                    continue;
                end

                % Monkey identity:
                % etab.monkey is ALWAYS the subject ID (even in *Opp fields*).
                % For opponent-role rows, infer the actor as the counterpart monkey.
                if isempty(etab.(FIELD_monkey){1})
                    warning('Missing monkey field; skipping a file in %s.', cname);
                    continue;
                end
                subjMonkeyStr = string(etab.(FIELD_monkey){1}(1));

                if role == "Subj"
                    monkeyStr = subjMonkeyStr;
                else
                    monkeyStr = inferOpponentMonkey(subjMonkeyStr);
                end
                if ~isMonkey12(monkeyStr), continue; end

                % Trial vectors (role-specific)
                optionRaw = etab.(optField){1}(:);
                trialNum  = double(etab.(trlField){1}(:));

                hasPre = ismember(preField, etab.Properties.VariableNames) && ~isempty(etab.(preField){1});
                if hasPre
                    prePortfolio = double(etab.(preField){1}(:));
                else
                    prePortfolio = nan(size(trialNum));
                end

                % Align lengths
                nT = min([numel(optionRaw), numel(trialNum), numel(prePortfolio)]);
                optionRaw    = optionRaw(1:nT);
                trialNum     = trialNum(1:nT);
                prePortfolio = prePortfolio(1:nT);

                % Decode BUY/HOLD/SELL
                actionCode = decodeOptionToAction(optionRaw); % 1=BUY,2=HOLD,3=SELL

                % Keep only scorable trials
                valid = isfinite(actionCode) & isfinite(trialNum);
                if ~any(valid), continue; end
                actionCode   = actionCode(valid);
                trialNum     = trialNum(valid);
                prePortfolio = prePortfolio(valid);

                % Reconstruct prePortfolio if missing or all-NaN
                if ~hasPre || all(~isfinite(prePortfolio))
                    prePortfolio = reconstructPrePortfolio(trialNum, actionCode);
                end
                prePortfolio = max(0, round(prePortfolio));

                % Post-choice portfolio (postPortfolio)
                postPortfolio = prePortfolio;
                postPortfolio(actionCode==ACTION_BUY)  = prePortfolio(actionCode==ACTION_BUY) + 1;
                postPortfolio(actionCode==ACTION_SELL) = max(prePortfolio(actionCode==ACTION_SELL) - 1, 0);

                % Metrics: action probabilities
                nTrials = numel(actionCode);
                pBuy    = mean(actionCode==ACTION_BUY);
                pHold   = mean(actionCode==ACTION_HOLD);
                pSell   = mean(actionCode==ACTION_SELL);

                % Diagnostics
                nSellAtZero = sum((actionCode==ACTION_SELL) & (prePortfolio==0));

                % Portfolio summary metrics
                deltaPortfolio        = postPortfolio - prePortfolio;  % per-trial net share flow
                meanDeltaPortfolio   = mean(deltaPortfolio);          % = P(BUY) - P(SELL) averaged over trials
                meanPostPortfolio    = mean(postPortfolio);

                % End-of-market postPortfolio: mean across market endings (trialNum==15)
                endMask = (trialNum == 15);
                if any(endMask)
                    endPostPortfolio_mean = mean(postPortfolio(endMask));
                    nMarkets              = sum(endMask);
                else
                    endPostPortfolio_mean = NaN;
                    nMarkets              = NaN;
                end

                rowRole = string(role);
                analysisSession = sessionStr;
                if baseCond == "Live" && LIVE_ROLE_MODE == "combined"
                    rowRole = "Combined";
                    analysisSession = sessionStr + "-" + string(role);
                end

                newRow = table( ...
                    string(treatName), string(cname), baseCond, ...
                    rowRole, ...
                    monkeyStr, sessionStr, ...
                    analysisSession, ...
                    nTrials, nMarkets, ...
                    pBuy, pHold, pSell, ...
                    meanDeltaPortfolio, meanPostPortfolio, endPostPortfolio_mean, ...
                    nSellAtZero, ...
                    'VariableNames', { ...
                        'Treatment','ConditionRaw','ConditionBase', ...
                        'Role', ...
                        'Monkey','Session','AnalysisSession', ...
                        'nTrials','nMarkets', ...
                        'pBuy','pHold','pSell', ...
                        'meanDeltaPortfolio','meanPostPortfolio','endPostPortfolioMean', ...
                        'nSellAtZero' ...
                    } ...
                );

                sessionSummaryAll = [sessionSummaryAll; newRow]; %#ok<AGROW>
            end
        end
    end
end

%% -------------------------- EXCLUSION REPORT -------------------------------
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

%% ------------------------------ CHECK --------------------------------------
if isempty(sessionSummaryAll)
    warning('No sessions found. Nothing to plot.');
    return;
end

% Factorize
sessionSummaryAll.Monkey = categorical(string(sessionSummaryAll.Monkey));

sessionSummaryAll.Treatment = categorical(string(sessionSummaryAll.Treatment));
sessionSummaryAll.Treatment = removecats(sessionSummaryAll.Treatment);
sessionSummaryAll.Treatment = reordercats(sessionSummaryAll.Treatment, ...
    treatOrder(ismember(treatOrder, categories(sessionSummaryAll.Treatment))));

sessionSummaryAll.ConditionBase = categorical(string(sessionSummaryAll.ConditionBase), baseCondList);
sessionSummaryAll.Role = categorical(string(sessionSummaryAll.Role), ["Subj","Opp","Combined"]);

%% ---------------------- SESSION INDEX (training order) ---------------------
% sessionIdx: chronological order within each monkey
sessKeyTab = unique(sessionSummaryAll(:, {'Monkey','AnalysisSession'}));
sessKeyTab = sortrows(sessKeyTab, {'Monkey','AnalysisSession'});

sessKeyTab.sessionIdx = zeros(height(sessKeyTab),1);
monkeys = categories(sessKeyTab.Monkey);
for im = 1:numel(monkeys)
    mask = sessKeyTab.Monkey == monkeys{im};
    sessKeyTab.sessionIdx(mask) = (1:nnz(mask)).';
end

sessionSummaryAll = outerjoin(sessionSummaryAll, sessKeyTab, ...
    'Keys', {'Monkey','AnalysisSession'}, ...
    'MergeKeys', true);
sessionSummaryAll = sortrows(sessionSummaryAll, {'Monkey','sessionIdx'});

%% --------------------------- OPTIONAL SAVE ---------------------------------
% writetable(sessionSummaryAll, fullfile(OUTDIR, 'msm_step1_actionRates_portfolioDelta_sessionSummary.csv'));

%% ============================= FIG 1 ======================================
% Action probabilities vs training index (sessionIdx), per monkey.
figure('Color','w','Name','FIG 1: Action probabilities over time');
tl1 = tiledlayout(3, numel(monkeys), 'TileSpacing','compact', 'Padding','compact');

metricFields = {'pBuy','pHold','pSell'};
metricNames  = {'P(BUY)','P(HOLD)','P(SELL)'};

for r = 1:3
    yField = metricFields{r};
    yName  = metricNames{r};

    for im = 1:numel(monkeys)
        mID  = monkeys{im};
        ax   = nexttile(tl1, (r-1)*numel(monkeys) + im);
        hold(ax,'on');

        rowsM = sessionSummaryAll.Monkey == mID;

        treatCats = categories(sessionSummaryAll.Treatment);

        % scatter by Treatment (shape) and ConditionBase (fill color)
        % scatter by Treatment (shape) and ConditionBase (fill color)
        % Live is role-split only when LIVE_ROLE_MODE=="separate".
        for it = 1:numel(treatCats)
            tName = treatCats{it};
            rowsT = sessionSummaryAll.Treatment == tName;
            marker = markerList{min(it, numel(markerList))};

            for ib = 1:numel(baseCondList)
                bName = baseCondList{ib};
                rowsB = sessionSummaryAll.ConditionBase == bName;

                if strcmp(bName, 'Live') && LIVE_ROLE_MODE == "separate"
                    % Plot Live Subj and Live Opp separately (dark/light)
                    for rr = ["Subj","Opp"]
                        rowsR = sessionSummaryAll.Role == rr;
                        rows  = rowsM & rowsT & rowsB & rowsR;
                        if ~any(rows), continue; end

                        x = sessionSummaryAll.sessionIdx(rows);
                        y = sessionSummaryAll.(yField)(rows);

                        fc = condFaceColors(bName);
                        if rr == "Opp"
                            fc = lightenColor(fc, 0.55);
                        end

                        scatter(ax, x, y, 45, ...
                            'Marker', marker, ...
                            'MarkerFaceColor', fc, ...
                            'MarkerEdgeColor', 'k', ...
                            'MarkerFaceAlpha', 0.75);
                    end
                else
                    rows  = rowsM & rowsT & rowsB;
                    if ~any(rows), continue; end

                    x = sessionSummaryAll.sessionIdx(rows);
                    y = sessionSummaryAll.(yField)(rows);

                    scatter(ax, x, y, 45, ...
                        'Marker', marker, ...
                        'MarkerFaceColor', condFaceColors(bName), ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerFaceAlpha', 0.75);
                end
            end
        end

        % trendline per monkey (pooled)
        xAll = sessionSummaryAll.sessionIdx(rowsM);
        yAll = sessionSummaryAll.(yField)(rowsM);
        ok   = isfinite(xAll) & isfinite(yAll);

        if sum(ok) > 2
            mdl = fitlm(double(xAll(ok)), double(yAll(ok)));
            xFit = linspace(min(xAll(ok)), max(xAll(ok)), 100)';
            [yFit, yCI] = predict(mdl, xFit, 'Alpha', 0.05);
            fill(ax, [xFit; flipud(xFit)], [yCI(:,1); flipud(yCI(:,2))], ...
                [0.7 0.7 0.7], 'EdgeColor','none', 'FaceAlpha', 0.25);
            plot(ax, xFit, yFit, 'k-', 'LineWidth', 1.5);
            txt = sprintf('slope=%.3g, p=%.3g', mdl.Coefficients.Estimate(2), mdl.Coefficients.pValue(2));
            text(ax, 0.02, 0.95, txt, 'Units','normalized', ...
                'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 9);
        end

        ylim(ax, [0 1]);
        xlabel(ax, 'sessionIdx');
        ylabel(ax, yName);
        title(ax, sprintf('Monkey %s', string(mID)));
        hold(ax,'off');
    end
end

if exist('sgtitle','file')
    sgtitle(tl1, 'Action probabilities', 'Interpreter','none');
end

%% --------------------------- LEGEND (shared) -------------------------------
% Put the legend inside the same figure as FIG 1 (instead of a separate figure).
% Create an invisible host axes for dummy legend handles, then attach legend to the tiledlayout.
axLeg = axes('Parent', gcf, 'Position',[0 0 1 1], 'Visible','off');
hold(axLeg,'on');

% Dummy handles for condition colors
% Live is split into Subj/Opp only when LIVE_ROLE_MODE=="separate".
if LIVE_ROLE_MODE == "separate"
    nCondLegend = numel(baseCondList) + 1;
else
    nCondLegend = numel(baseCondList);
end
condHandles = gobjects(nCondLegend,1);
condLabels  = strings(nCondLegend,1);
idx = 1;
for i = 1:numel(baseCondList)
    bName = baseCondList{i};
    if strcmp(bName, 'Live') && LIVE_ROLE_MODE == "separate"
        condHandles(idx) = scatter(axLeg, nan, nan, 60, ...
            'Marker', 'o', ...
            'MarkerFaceColor', condFaceColors(bName), ...
            'MarkerEdgeColor', 'k');
        condLabels(idx) = "Condition: Live (Subj)";
        idx = idx + 1;

        condHandles(idx) = scatter(axLeg, nan, nan, 60, ...
            'Marker', 'o', ...
            'MarkerFaceColor', lightenColor(condFaceColors(bName), 0.55), ...
            'MarkerEdgeColor', 'k');
        condLabels(idx) = "Condition: Live (Opp)";
        idx = idx + 1;
    else
        condHandles(idx) = scatter(axLeg, nan, nan, 60, ...
            'Marker', 'o', ...
            'MarkerFaceColor', condFaceColors(bName), ...
            'MarkerEdgeColor', 'k');
        condLabels(idx) = "Condition: " + string(bName);
        idx = idx + 1;
    end
end

% Dummy handles for treatment shapes
treatCats = categories(sessionSummaryAll.Treatment);
treatHandles = gobjects(numel(treatCats),1);
for i = 1:numel(treatCats)
    treatHandles(i) = scatter(axLeg, nan, nan, 60, ...
        'Marker', markerList{min(i,numel(markerList))}, ...
        'MarkerFaceColor', [0.85 0.85 0.85], ...
        'MarkerEdgeColor', 'k');
end

legHandles = [condHandles; treatHandles];
legLabels  = [condLabels; ("Treatment: " + string(treatCats(:)))];

% NOTE: legend() cannot take a TiledChartLayout handle. Attach the legend
% to an axes within the tiledlayout, then (if supported) dock it into the
% layout.
lgd = legend(axLeg, legHandles, cellstr(legLabels), ...
    'Location','eastoutside', 'Box','on', 'Orientation','vertical');
lgd.AutoUpdate = 'off';

% % For newer MATLAB versions: place the legend into the tiledlayout region.
% if isprop(lgd,'Layout') & isprop(lgd.Layout,'Tile')
%     lgd.Layout.Tile = 'south';
% end

%% ========================= LOCAL FUNCTIONS =================================
function actionCode = decodeOptionToAction(optionRaw)
% Map optionRaw -> {BUY,HOLD,SELL} codes:
%   BUY  = 1
%   HOLD = 2
%   SELL = 3
%
% Robust handling:
%   - If strings: matches tokens 'buy','hold','sell' or 'b','h','s'.
%   - If numeric: assumes common MSM encoding 1/2/3 or 0/1/2.
%     If numeric is different, you must update mapping explicitly here.

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
% Reconstruct prePortfolio from actions under MSM constraints:
%   - starts at 0 at the first trial of each market instance (trialNum==1)
%   - BUY increments, SELL decrements (clipped at 0), HOLD does nothing

    n = numel(actionCode);
    prePortfolio = nan(n,1);

    p = 0;
    for i = 1:n
        if i==1 || trialNum(i)==1
            p = 0;
        end

        prePortfolio(i) = p;

        if actionCode(i)==ACTION_BUY
            p = p + 1;
        elseif actionCode(i)==ACTION_SELL
            p = max(p - 1, 0);
        end
    end
end

function opp = inferOpponentMonkey(subj)
% Infer opponent monkey ID/name from subject monkey.
% Assumes exactly two monkeys with IDs 1 and 2 (or labels 'M1' and 'M2').
% Extend this mapping if your IDs differ.

    s = string(subj);
    s = strtrim(s);

    % Try numeric swap first
    val = str2double(s);
    if isfinite(val)
        if val == 1
            opp = "2";
            return;
        elseif val == 2
            opp = "1";
            return;
        end
    end

    % Try common string labels
    sU = upper(s);
    if sU == "M1"
        opp = "M2";
        return;
    elseif sU == "M2"
        opp = "M1";
        return;
    end

    error("inferOpponentMonkey: cannot infer opponent from subject '%s'. Update mapping.", s);
end

function ok = isMonkey12(monkeyVal)
% Keep only rows where monkey is 1/2 (or M1/M2 string equivalent).
    s = upper(strtrim(string(monkeyVal)));
    v = str2double(s);
    ok = (isfinite(v) && any(v == [1 2])) || any(s == ["M1","M2"]);
end

function c2 = lightenColor(c1, amt)
% amt in [0,1]; higher = lighter
    amt = max(0, min(1, amt));
    c2 = c1 + (1 - c1) * amt;
end

function v = ACTION_BUY
% Named constants (avoid magic numbers)
    v = 1;
end
function v = ACTION_HOLD
    v = 2;
end
function v = ACTION_SELL
    v = 3;
end
