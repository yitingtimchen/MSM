%% MSM Step 4c: Categorical trial/market controls + plot only (prePortfolio, sessionIdx) coefficients
% -------------------------------------------------------------------------
% ROLE-AWARE LIVE (+ toggle)
%
% YOU ASKED FOR (original intent):
%   - Use trialNum and marketOrig as *categorical* predictors (trialCat, marketCat).
%   - Only PLOT coefficients for:
%       zPrePortfolio  and  zSessionIdx
%     for NAIVE vs CONTROLLED models, comparing M1 vs M2 vs pooled.
%   - In the report: do NOT print per-level coefficients for trialCat/marketCat.
%     Instead, report whether they "generally matter" (omnibus LR tests).
%
% NEW IN THIS VERSION:
%   - Live condition can include trials where the actor is the opponent (stored in *_Opp fields)
%   - Toggle: LIVE_ROLE_MODE ("subj_only" / "separate" / "combined")
%   - Adds Role column (Subj/Opp) in trialTableAll
%   - Opp monkey inferred (etab.monkey is always subject)
%   - sessionIdx defined per (Monkey, Session) and shared across roles
%
% MODELS (BUY probability; binomial-logit)
%   NAIVE:
%     isBuy ~ zPrePortfolio + zSessionIdx
%
%   CONTROLLED (categorical controls):
%     isBuy ~ zPrePortfolio + zSessionIdx + trialCat + marketCat
%
% "DO THEY MATTER?" TESTS (likelihood ratio; nested models)
%   - trialCat matters?:  compare NAIVE vs (+trialCat)
%   - marketCat matters?: compare (+trialCat) vs (+trialCat + marketCat)
%
% NOTES
%   - We z-score within each group (M1, M2, pooled) so coefficients are
%     comparable as log-odds per 1 SD change.
%   - trialCat is categorical T1..T15; marketCat is categorical templates (e.g., 1..6).

clear; clc; close all;

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL'; % <-- EDIT

cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames     = {'baseline','OT','Saline'};          % preserve ordering
baseCondList = {'AI','Replay','Decoy','Live'};      % preserve ordering

% Live role handling:
%   "subj_only" -> include only Subj stream in Live
%   "separate"  -> include Subj and Opp as separate role conditions
%   "combined"  -> include both streams but pool them as one Live role;
%                  session indexing uses role-qualified session keys so
%                  effective session count doubles for Live per monkey.
LIVE_ROLE_MODE = "separate";  % "subj_only" | "separate" | "combined"

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
EXCLUDE_SESSIONS = { ...
    % '2018-07-07-01', ...
};
EXCLUDE_SESSIONS = string(EXCLUDE_SESSIONS(:));
EXCLUDE_SESSIONS = EXCLUDE_SESSIONS(strlength(EXCLUDE_SESSIONS) > 0);
excludedSessSeen = strings(0,1);

% Optional filters
FILTER_TREATMENT = "all"; % "all" or "baseline"/"OT"/"Saline"
FILTER_COND_BASE = "all"; % "all" or "AI"/"Replay"/"Decoy"/"Live"

Tmax = 15;

% ---- Adjusted-curve visualization config ----
Pgrid = (0:Tmax)';      % portfolio grid to visualize
minCountPerP = 20;      % for plotting raw P(BUY|p) points with CI

%% ------------------------ REQUIRED ETAB FIELD MAP --------------------------
% --- Subject stream ---
FIELD_option       = 'option';
FIELD_trialNum     = 'trialNum';
FIELD_prePortfolio = 'prePortfolio';     % optional
FIELD_marketOrig   = 'marketOrig';       % required for marketCat test

% --- Live opponent stream (Live only; used unless LIVE_ROLE_MODE="subj_only") ---
FIELD_optionOpp       = 'optionOpp';
FIELD_trialNumOpp     = 'trialNumOpp';
FIELD_prePortfolioOpp = 'prePortfolioOpp';  % optional
FIELD_marketOrigOpp   = 'marketOrigOpp';    % required for marketCat on Opp stream (if absent, Opp rows skipped)

% --- Session identity ---
FIELD_year         = 'year';
FIELD_month        = 'month';
FIELD_day          = 'day';
FIELD_SessionOfDay = 'SessionOfDay';
FIELD_monkey       = 'monkey';

%% --------------------- BUILD MASTER TRIAL-LEVEL TABLE ----------------------
trialTableAll = table();

fprintf('\n=== BUILDING trialTableAll (trial-level rows; role-aware Live) ===\n');

for s = 1:numel(cond_sets)
    conds     = cond_sets{s};
    treatName = setNames{s};
    fprintf('Loading treatment: %s\n', treatName);

    if FILTER_TREATMENT ~= "all" && string(treatName) ~= FILTER_TREATMENT
        continue;
    end

    for c = 1:numel(conds)
        cname = conds{c};

        tok = regexp(cname, '(AI|Replay|Decoy|Live)$', 'tokens', 'once');
        if isempty(tok)
            baseCond = string(cname);
        else
            baseCond = string(tok{1});
        end

        if FILTER_COND_BASE ~= "all" && baseCond ~= FILTER_COND_BASE
            continue;
        end

        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));
        if ~exist(matFile,'file')
            warning('Missing condition pack: %s', matFile);
            continue;
        end

        tmp = load(matFile, 'C');
        C = tmp.C;

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            neededID = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey};
            if ~all(ismember(neededID, etab.Properties.VariableNames))
                warning('Missing required identity fields; skipping a file in %s.', cname);
                continue;
            end

            % Session identity (yyyy-mm-dd-sod)
            y   = etab.(FIELD_year){1}(:);
            mth = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(mth) || isempty(d) || isempty(sod), continue; end

            sessionStr  = string(sprintf('%04d-%02d-%02d-%02d', y(1), mth(1), d(1), sod(1)));
            sessionDate = datetime(y(1), mth(1), d(1));

            if ~isempty(EXCLUDE_SESSIONS) && any(sessionStr == EXCLUDE_SESSIONS)
                excludedSessSeen(end+1,1) = sessionStr; %#ok<AGROW>
                continue;
            end

            % Subject monkey identity (etab.monkey always subject)
            if isempty(etab.(FIELD_monkey){1})
                warning('Missing monkey field; skipping a file in %s.', cname);
                continue;
            end
            subjMonkeyStr = string(etab.(FIELD_monkey){1}(1));

            % Role split: Live only (+toggle)
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
                % Role-specific fieldnames + actor monkey label
                if role == "Subj"
                    f_option   = FIELD_option;
                    f_trialNum = FIELD_trialNum;
                    f_prePort  = FIELD_prePortfolio;
                    f_market   = FIELD_marketOrig;
                    actorMonkeyStr = subjMonkeyStr;
                else
                    f_option   = FIELD_optionOpp;
                    f_trialNum = FIELD_trialNumOpp;
                    f_prePort  = FIELD_prePortfolioOpp;
                    f_market   = FIELD_marketOrigOpp;
                    actorMonkeyStr = inferOpponentMonkey(subjMonkeyStr);
                end
                if ~isMonkey12(actorMonkeyStr), continue; end

                neededRole = {f_option, f_trialNum, f_market};
                if ~all(ismember(neededRole, etab.Properties.VariableNames))
                    % For Live Opp, skip cleanly if Opp fields absent
                    continue;
                end

                if isempty(etab.(f_option){1}) || isempty(etab.(f_trialNum){1}) || isempty(etab.(f_market){1})
                    continue;
                end

                optionRaw = etab.(f_option){1}(:);
                trialNum  = double(etab.(f_trialNum){1}(:));
                marketOrig = double(etab.(f_market){1}(:));

                hasPre = ismember(f_prePort, etab.Properties.VariableNames) && ~isempty(etab.(f_prePort){1});
                if hasPre
                    prePortfolio = double(etab.(f_prePort){1}(:));
                else
                    prePortfolio = nan(size(trialNum));
                end

                nT = min([numel(optionRaw), numel(trialNum), numel(prePortfolio), numel(marketOrig)]);
                optionRaw    = optionRaw(1:nT);
                trialNum     = trialNum(1:nT);
                prePortfolio = prePortfolio(1:nT);
                marketOrig   = marketOrig(1:nT);

                actionCode = decodeOptionToAction(optionRaw); % 1=BUY,2=HOLD,3=SELL

                keep = isfinite(actionCode) & isfinite(trialNum) & isfinite(marketOrig);
                if ~any(keep), continue; end
                actionCode   = actionCode(keep);
                trialNum     = trialNum(keep);
                prePortfolio = prePortfolio(keep);
                marketOrig   = marketOrig(keep);

                if ~hasPre || all(~isfinite(prePortfolio))
                    prePortfolio = reconstructPrePortfolio(trialNum, actionCode);
                end
                prePortfolio = max(0, round(prePortfolio));

                nKeep = numel(actionCode);

                rowRole = string(role);
                analysisSession = sessionStr;
                if baseCond == "Live" && LIVE_ROLE_MODE == "combined"
                    rowRole = "Combined";
                    analysisSession = sessionStr + "-" + string(role);
                end

                newRows = table( ...
                    repmat(string(treatName), nKeep, 1), ...
                    repmat(string(cname),     nKeep, 1), ...
                    repmat(baseCond,          nKeep, 1), ...
                    repmat(rowRole,           nKeep, 1), ...
                    repmat(actorMonkeyStr,    nKeep, 1), ...
                    repmat(sessionStr,        nKeep, 1), ...
                    repmat(analysisSession,   nKeep, 1), ...
                    repmat(sessionDate,       nKeep, 1), ...
                    trialNum(:), ...
                    prePortfolio(:), ...
                    actionCode(:), ...
                    marketOrig(:), ...
                    'VariableNames', { ...
                        'Treatment','ConditionRaw','ConditionBase','Role', ...
                        'Monkey','Session','AnalysisSession','SessionDate', ...
                        'trialNum','prePortfolio','actionCode','marketOrig' ...
                    } ...
                );

                trialTableAll = [trialTableAll; newRows]; %#ok<AGROW>
            end % role
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
if isempty(trialTableAll)
    warning('No trials found. Nothing to analyze.');
    return;
end

% Preserve ordering (don't let MATLAB reorder)
trialTableAll.Monkey        = categorical(string(trialTableAll.Monkey));
trialTableAll.Role          = categorical(string(trialTableAll.Role), ["Subj","Opp","Combined"]);

trialTableAll.Treatment     = categorical(string(trialTableAll.Treatment), setNames);
trialTableAll.Treatment     = removecats(trialTableAll.Treatment);

trialTableAll.ConditionBase = categorical(string(trialTableAll.ConditionBase), baseCondList);
trialTableAll.ConditionBase = removecats(trialTableAll.ConditionBase);

trialTableAll.isBuy = double(trialTableAll.actionCode == ACTION_BUY);

% Categorical controls
trialTableAll.trialCat = categorical(trialTableAll.trialNum, 1:Tmax, compose('T%d', 1:Tmax));
trialTableAll.marketCat = categorical(trialTableAll.marketOrig);

%% ---------------------- SESSION INDEX (training order) ---------------------
% Define sessionIdx per (Monkey, AnalysisSession).
% In LIVE_ROLE_MODE=="combined", Live streams get role-qualified keys and
% therefore double session density for each monkey.
trialTableAll.sessionIdx = nan(height(trialTableAll),1);
monkeys = categories(trialTableAll.Monkey);
if numel(monkeys) < 2
    warning('Expected 2 monkeys for M1/M2 comparison. Found: %d', numel(monkeys));
end

sessKeyTab = unique(trialTableAll(:, {'Monkey','AnalysisSession'}));
sessKeyTab = sortrows(sessKeyTab, {'Monkey','AnalysisSession'});

for im = 1:numel(monkeys)
    mID = monkeys{im};
    rowsM = sessKeyTab.Monkey == mID;
    sessList = sessKeyTab.AnalysisSession(rowsM);
    sessIdxLocal = (1:numel(sessList))';

    for k = 1:numel(sessList)
        mask = (trialTableAll.Monkey == mID) & (trialTableAll.AnalysisSession == sessList(k));
        trialTableAll.sessionIdx(mask) = sessIdxLocal(k);
    end
end

%% ---------------------- FIT GROUPS: M1, M2, pooled -------------------------
% (kept as in your script: per-monkey + pooled; role is included/excluded upstream via toggle)
groupNames = ["M1","M2","Pooled"];

% Use category order as-is; assumes first two categories correspond to M1/M2
if numel(monkeys) >= 2
    groupMasks = { ...
        trialTableAll.Monkey == monkeys{1}, ...
        trialTableAll.Monkey == monkeys{2}, ...
        true(height(trialTableAll),1) ...
    };
else
    groupMasks = {true(height(trialTableAll),1)};
    groupNames = ["Pooled"];
end

% Tables to collect coef of zPrePortfolio and zSessionIdx
coefNaiveAll = table();
coefCtrlAll  = table();

% Omnibus "do trialCat/marketCat matter" tests
omnibusAll = table();

% Adjusted curve storage (raw + marginal-standardized predictions)
curveAll = table();

for gi = 1:numel(groupNames)
    gName = groupNames(gi);
    D = trialTableAll(groupMasks{gi}, :);

    ok = isfinite(D.isBuy) & isfinite(D.prePortfolio) & isfinite(D.sessionIdx) & ...
         isfinite(D.trialNum) & isfinite(D.marketOrig) & ~isundefined(D.trialCat) & ~isundefined(D.marketCat);
    D = D(ok,:);
    if isempty(D)
        warning('No valid rows for group %s; skipping.', gName);
        continue;
    end

    % Z-score within group for comparability
    muPre = mean(double(D.prePortfolio));
    sdPre = std(double(D.prePortfolio));
    muSes = mean(double(D.sessionIdx));
    sdSes = std(double(D.sessionIdx));

    D.zPrePortfolio = (double(D.prePortfolio) - muPre) ./ sdPre;
    D.zSessionIdx   = (double(D.sessionIdx)   - muSes) ./ sdSes;

    % Models for LR tests
    mdlBase  = fitglm(D, 'isBuy ~ zPrePortfolio + zSessionIdx', ...
        'Distribution','binomial','Link','logit');

    mdlTrial = fitglm(D, 'isBuy ~ zPrePortfolio + zSessionIdx + trialCat', ...
        'Distribution','binomial','Link','logit');

    mdlFull  = fitglm(D, 'isBuy ~ zPrePortfolio + zSessionIdx + trialCat + marketCat', ...
        'Distribution','binomial','Link','logit');

    % Extract ONLY the two continuous predictors (what you want plotted)
    coefNaive = extractCoefs(mdlBase, ["zPrePortfolio","zSessionIdx"]);
    coefNaive.Group = repmat(gName, height(coefNaive), 1);
    coefNaiveAll = [coefNaiveAll; coefNaive]; %#ok<AGROW>

    coefCtrl = extractCoefs(mdlFull, ["zPrePortfolio","zSessionIdx"]);
    coefCtrl.Group = repmat(gName, height(coefCtrl), 1);
    coefCtrlAll = [coefCtrlAll; coefCtrl]; %#ok<AGROW>

    % Omnibus LR tests: do the categorical blocks improve fit?
    [pTrial,  dfTrial,  chi2Trial]  = lrBlock(mdlBase,  mdlTrial);
    [pMarket, dfMarket, chi2Market] = lrBlock(mdlTrial, mdlFull);

    omnibusAll = [omnibusAll; table( ...
        gName, height(D), ...
        pTrial,  dfTrial,  chi2Trial, ...
        pMarket, dfMarket, chi2Market, ...
        'VariableNames', {'Group','nTrials', ...
                          'p_trialCat','df_trialCat','chi2_trialCat', ...
                          'p_marketCat','df_marketCat','chi2_marketCat'} ...
    )]; %#ok<AGROW>

    % ------------------ RAW P(BUY|p) + Wilson CI (descriptive) ------------------
    raw_p  = nan(size(Pgrid));
    raw_lo = nan(size(Pgrid));
    raw_hi = nan(size(Pgrid));
    raw_n  = zeros(size(Pgrid));

    for iP = 1:numel(Pgrid)
        p = Pgrid(iP);
        idxP = (D.prePortfolio == p);
        nP = sum(idxP);
        raw_n(iP) = nP;

        if nP >= minCountPerP
            kP = sum(D.isBuy(idxP) == 1);
            [phat, lo, hi] = wilsonCI(kP, nP, 0.05);
            raw_p(iP)  = phat;
            raw_lo(iP) = lo;
            raw_hi(iP) = hi;
        end
    end

    % ------------------ ADJUSTED CURVES (marginal standardization) --------------
    % Fix zSessionIdx = 0 (group mean sessionIdx)
    zSess0 = 0;

    % NAIVE adjusted curve (no trial/market controls)
    pNaive = nan(size(Pgrid));
    for iP = 1:numel(Pgrid)
        zPre = (double(Pgrid(iP)) - muPre) ./ sdPre;
        Tpred = table(zPre, zSess0, 'VariableNames', {'zPrePortfolio','zSessionIdx'});
        pNaive(iP) = predict(mdlBase, Tpred);
    end

    % CONTROLLED adjusted curve: average predictions across observed trialCat/marketCat
    pCtrl = nan(size(Pgrid));
    for iP = 1:numel(Pgrid)
        zPre = (double(Pgrid(iP)) - muPre) ./ sdPre;

        Tpred = table( ...
            repmat(zPre, height(D), 1), ...
            repmat(zSess0, height(D), 1), ...
            D.trialCat, ...
            D.marketCat, ...
            'VariableNames', {'zPrePortfolio','zSessionIdx','trialCat','marketCat'} ...
        );

        pHat = predict(mdlFull, Tpred);
        pCtrl(iP) = mean(pHat, 'omitnan');
    end

    curveAll = [curveAll; table( ...
        repmat(gName, numel(Pgrid), 1), ...
        Pgrid, raw_p, raw_lo, raw_hi, raw_n, pNaive, pCtrl, ...
        repmat(muPre, numel(Pgrid), 1), repmat(sdPre, numel(Pgrid), 1), ...
        'VariableNames', {'Group','prePortfolio','raw_p','raw_lo','raw_hi','raw_n','pNaive','pCtrl','muPre','sdPre'} ...
    )]; %#ok<AGROW>
end

%% ---------------------- REPORT TABLES --------------------------------------
fprintf('\n=== NAIVE: isBuy ~ zPrePortfolio + zSessionIdx (only plotted predictors) ===\n');
disp(coefNaiveAll);

fprintf('\n=== CONTROLLED: isBuy ~ zPrePortfolio + zSessionIdx + trialCat + marketCat (only plotted predictors) ===\n');
disp(coefCtrlAll);

fprintf('\n=== Omnibus tests: do trialCat / marketCat generally matter? (LR tests) ===\n');
disp(omnibusAll);

%% ---------------------- FIG: ADJUSTED CURVES (SHOWS CONTROLS) -----------------
figure('Color','w','Name','Adjusted P(BUY) vs portfolio (controls: trialCat + marketCat)');
tl = tiledlayout(1, numel(groupNames), 'TileSpacing','compact', 'Padding','compact');

for gi = 1:numel(groupNames)
    gName = groupNames(gi);
    ax = nexttile(tl, gi); hold(ax,'on');

    R = curveAll(curveAll.Group == gName, :);

    % Raw points (+Wilson CI) where defined
    okRaw = isfinite(R.raw_p) & isfinite(R.raw_lo) & isfinite(R.raw_hi);
    if any(okRaw)
        errorbar(ax, R.prePortfolio(okRaw), R.raw_p(okRaw), ...
            R.raw_p(okRaw) - R.raw_lo(okRaw), R.raw_hi(okRaw) - R.raw_p(okRaw), ...
            'ko', 'MarkerFaceColor','k', 'MarkerSize', 4, 'LineWidth', 1.0, ...
            'DisplayName', sprintf('raw p(Buy) with CI, n>=%d', minCountPerP));
    end

    % Adjusted curves
    plot(ax, R.prePortfolio, R.pNaive, '--', 'LineWidth', 2.0, ...
        'DisplayName', 'isBuy ~ zPortfolio + zSessionIdx(=0)');
    plot(ax, R.prePortfolio, R.pCtrl,  '-',  'LineWidth', 2.5, ...
        'DisplayName', 'isBuy ~ zPortfolio + zSessionIdx(=0) + trialCat + marketCat');

    xlabel(ax, 'Portfolio size (shares)');
    ylabel(ax, 'P(BUY)');
    title(ax, sprintf('%s', gName));
    ylim(ax, [0 1]);
    xlim(ax, [min(Pgrid) max(Pgrid)]);

    legend(ax, 'Location','southeast', 'Box','off');
    hold(ax,'off');
end

if exist('sgtitle','file')
    sgtitle(tl, 'Adjusted P(BUY) vs inventory: naive vs controlled (trialCat + marketCat controls)', 'Interpreter','none');
end

%% ========================= LOCAL FUNCTIONS =================================

function opp = inferOpponentMonkey(subj)
% Infer opponent monkey ID/name from subject monkey.
% Assumes exactly two monkeys with IDs 1 and 2 (or labels 'M1' and 'M2').
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

function ok = isMonkey12(monkeyVal)
% Keep only rows where monkey is 1/2 (or M1/M2 string equivalent).
    s = upper(strtrim(string(monkeyVal)));
    v = str2double(s);
    ok = (isfinite(v) && any(v == [1 2])) || any(s == ["M1","M2"]);
end

function [phat, lo, hi] = wilsonCI(k, n, alpha)
% Wilson score interval for a binomial proportion.
    if n==0
        phat=NaN; lo=NaN; hi=NaN; return;
    end
    z = norminv(1 - alpha/2);
    phat = k/n;
    denom = 1 + z^2/n;
    center = (phat + z^2/(2*n)) / denom;
    half = (z/denom) * sqrt( (phat*(1-phat) + z^2/(4*n)) / n );
    lo = max(0, center - half);
    hi = min(1, center + half);
end

function Tcoef = extractCoefs(mdl, coefNames)
% Return a tidy table of Estimate/SE/p for selected coefficient names.
    Tcoef = table();
    for i = 1:numel(coefNames)
        nm = coefNames(i);
        idx = strcmp(mdl.CoefficientNames, nm);
        if ~any(idx)
            warning('Coefficient %s not found. Available: %s', nm, strjoin(string(mdl.CoefficientNames),', '));
            est = NaN; se = NaN; p = NaN;
        else
            est = mdl.Coefficients.Estimate(idx);
            se  = mdl.Coefficients.SE(idx);
            p   = mdl.Coefficients.pValue(idx);
        end
        Tcoef = [Tcoef; table(nm, est, se, p, 'VariableNames', {'Predictor','Estimate','SE','pValue'})]; %#ok<AGROW>
    end
end

function [p, df, chi2] = lrBlock(mdl0, mdl1)
% Likelihood ratio test for nested models mdl0 (reduced) vs mdl1 (full).
% Returns p-value, df difference, and chi-square statistic.
    ll0 = mdl0.LogLikelihood;
    ll1 = mdl1.LogLikelihood;
    df  = mdl1.NumEstimatedCoefficients - mdl0.NumEstimatedCoefficients;
    chi2 = 2*(ll1 - ll0);
    if df <= 0 || ~isfinite(chi2)
        p = NaN; return;
    end
    p = 1 - chi2cdf(chi2, df);
end

function actionCode = decodeOptionToAction(optionRaw)
% Map optionRaw -> {BUY,HOLD,SELL} codes:
%   BUY  = 1
%   HOLD = 2
%   SELL = 3
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
% Reconstruct prePortfolio given only trialNum and actionCode.
% Resets to 0 whenever trialNum==1 (new market).
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

function v = ACTION_BUY
    v = 1;
end
function v = ACTION_HOLD
    v = 2;
end
function v = ACTION_SELL
    v = 3;
end
