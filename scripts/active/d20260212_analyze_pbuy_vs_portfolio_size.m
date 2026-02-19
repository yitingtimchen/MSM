% File: d20260212_analyze_pbuy_vs_portfolio_size.m
%% MSM Step 3: Test for an inventory setpoint ("buy-until-K") strategy over training (ROLE-AWARE LIVE)
% -------------------------------------------------------------------------
% ROLE-AWARE version of pBuy_vs_portfolio_size_260209
%
% Key additions (mirrors your "Live includes opponent-role trials" pattern):
%   - Live condition split into Role="Subj" and Role="Opp" using *_Opp fields
%   - Opp role monkey identity inferred (because etab.monkey is always subject)
%   - Role carried into trialTableAll
%   - sessionIdx is defined per (Monkey, Session) and shared across roles
%   - Everything else (filters, early/late split, CI computation, plots) preserved
%
% IMPORTANT:
%   - Does NOT change your treatment/condition ordering:
%       setNames order: baseline -> OT -> Saline
%       baseCondList order: AI -> Replay -> Decoy -> Live
%
% NOTE ABOUT FIELDNAMES
%   This script assumes your packs use the following *_Opp fieldnames:
%       optionOpp, trialNumOpp, prePortfolioOpp
%   If your pack uses different names, edit ONLY the "FIELD MAP" section.

clear; clc; close all;

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL'; % <-- EDIT

cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames     = {'baseline','OT','Saline'};          % <<< preserve order
baseCondList = {'AI','Replay','Decoy','Live'};      % <<< preserve order

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
EXCLUDE_SESSIONS = { ...
    % '2018-07-07-01', ...
};
EXCLUDE_SESSIONS = string(EXCLUDE_SESSIONS(:));
EXCLUDE_SESSIONS = EXCLUDE_SESSIONS(strlength(EXCLUDE_SESSIONS) > 0);
excludedSessSeen = strings(0,1);

% Optional filters (leave as "all" to include everything)
FILTER_TREATMENT   = "all";  % "all" or "baseline"/"OT"/"Saline"
FILTER_COND_BASE   = "all";  % "all" or "AI"/"Replay"/"Decoy"/"Live"

% Live role handling:
%   "subj_only" -> include only Subj stream in Live
%   "separate"  -> include Subj and Opp as separate role conditions
%   "combined"  -> include both streams but pool them as one Live role;
%                  session indexing uses role-qualified session keys so
%                  effective session count doubles for Live per monkey.
LIVE_ROLE_MODE = "subj_only";  % "subj_only" | "separate" | "combined"


% Training split within monkey
earlyQuant = 0.33;
lateQuant  = 0.67;

% Inventory range to analyze/plot
PmaxPlot = 15;     % 15 trials => p never exceeds 15
minCountPerP = 20; % minimum trials at that p within early/late group to show CI/point

%% ------------------------ REQUIRED ETAB FIELD MAP --------------------------
% --- Subject stream fields ---
FIELD_option       = 'option';
FIELD_trialNum     = 'trialNum';
FIELD_prePortfolio = 'prePortfolio';  % optional

% --- Live opponent stream fields (Live only) ---
FIELD_optionOpp       = 'optionOpp';
FIELD_trialNumOpp     = 'trialNumOpp';
FIELD_prePortfolioOpp = 'prePortfolioOpp'; % optional

% --- Session identity fields ---
FIELD_year         = 'year';
FIELD_month        = 'month';
FIELD_day          = 'day';
FIELD_SessionOfDay = 'SessionOfDay';
FIELD_monkey       = 'monkey';

%% --------------------- BUILD MASTER TRIAL-LEVEL TABLE ----------------------
trialTableAll = table();

fprintf('\n=== BUILDING trialTableAll (trial-level rows; role-aware) ===\n');

for s = 1:numel(cond_sets)
    conds     = cond_sets{s};
    treatName = setNames{s};
    fprintf('Loading treatment: %s\n', treatName);

    if FILTER_TREATMENT ~= "all" && string(treatName) ~= FILTER_TREATMENT
        continue;
    end

    for c = 1:numel(conds)
        cname   = conds{c};

        % Parse base condition token
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
        C   = tmp.C;

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};

            % Identity fields always required
            neededID = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey};
            if ~all(ismember(neededID, etab.Properties.VariableNames))
                warning('Missing required identity fields; skipping a file in %s.', cname);
                continue;
            end

            % Session identity
            y   = etab.(FIELD_year){1}(:);
            m   = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(m) || isempty(d) || isempty(sod), continue; end

            sessionStr  = string(sprintf('%04d-%02d-%02d-%02d', y(1), m(1), d(1), sod(1)));
            sessionDate = datetime(y(1), m(1), d(1));

            % Exclusions
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

            % Role split: Live only (toggle controls whether/how to include Opp)
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
                    actorMonkeyStr = subjMonkeyStr;
                else
                    f_option   = FIELD_optionOpp;
                    f_trialNum = FIELD_trialNumOpp;
                    f_prePort  = FIELD_prePortfolioOpp;
                    actorMonkeyStr = inferOpponentMonkey(subjMonkeyStr);
                end
                if ~isMonkey12(actorMonkeyStr), continue; end

                % Required role fields
                neededRole = {f_option, f_trialNum};
                if ~all(ismember(neededRole, etab.Properties.VariableNames))
                    % For Live Opp, skip cleanly if Opp fields absent
                    continue;
                end

                if isempty(etab.(f_option){1}) || isempty(etab.(f_trialNum){1})
                    continue;
                end

                optionRaw = etab.(f_option){1}(:);
                trialNum  = double(etab.(f_trialNum){1}(:));

                hasPre = ismember(f_prePort, etab.Properties.VariableNames) && ~isempty(etab.(f_prePort){1});
                if hasPre
                    prePortfolio = double(etab.(f_prePort){1}(:));
                else
                    prePortfolio = nan(size(trialNum));
                end

                nT = min([numel(optionRaw), numel(trialNum), numel(prePortfolio)]);
                optionRaw    = optionRaw(1:nT);
                trialNum     = trialNum(1:nT);
                prePortfolio = prePortfolio(1:nT);

                actionCode = decodeOptionToAction(optionRaw); % 1=BUY,2=HOLD,3=SELL

                keep = isfinite(actionCode) & isfinite(trialNum);
                if ~any(keep), continue; end
                actionCode   = actionCode(keep);
                trialNum     = trialNum(keep);
                prePortfolio = prePortfolio(keep);

                % Reconstruct prePortfolio if missing or all-NaN
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
                    'VariableNames', { ...
                        'Treatment','ConditionRaw','ConditionBase','Role', ...
                        'Monkey','Session','AnalysisSession','SessionDate', ...
                        'trialNum','prePortfolio','actionCode' ...
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

% Preserve category ordering
trialTableAll.Monkey        = categorical(string(trialTableAll.Monkey));
trialTableAll.Role          = categorical(string(trialTableAll.Role), ["Subj","Opp","Combined"]);

trialTableAll.Treatment     = categorical(string(trialTableAll.Treatment), setNames);
trialTableAll.Treatment     = removecats(trialTableAll.Treatment);

trialTableAll.ConditionBase = categorical(string(trialTableAll.ConditionBase), baseCondList);
trialTableAll.ConditionBase = removecats(trialTableAll.ConditionBase);

monkeys = categories(trialTableAll.Monkey);

%% ---------------------- SESSION INDEX (training order) ---------------------
% Define sessionIdx per (Monkey, AnalysisSession).
% In LIVE_ROLE_MODE=="combined", Live streams get role-qualified keys and
% therefore double session density for each monkey.
trialTableAll.sessionIdx = nan(height(trialTableAll),1);

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

%% ======================= COMPUTE P(BUY|p), P(SELL|p) =======================
results = struct([]);
resIdx = 0;

for im = 1:numel(monkeys)
    mID = monkeys{im};
    Dm = trialTableAll(trialTableAll.Monkey == mID, :);
    if isempty(Dm), continue; end

    % Keep early/late split defined at the monkey level (shared across roles).
    sessLo = quantile(Dm.sessionIdx, earlyQuant);
    sessHi = quantile(Dm.sessionIdx, lateQuant);

    if LIVE_ROLE_MODE == "separate"
        roleListAnalysis = ["Subj","Opp"];
    else
        roleListAnalysis = "All";
    end

    for r = 1:numel(roleListAnalysis)
        roleName = roleListAnalysis(r);
        if roleName == "All"
            D = Dm;
        else
            D = Dm(Dm.Role == roleName, :);
        end
        if isempty(D), continue; end

        isEarly = D.sessionIdx <= sessLo;
        isLate  = D.sessionIdx >= sessHi;

        pGrid = (0:PmaxPlot)';

        % Initialize
        pBuyEarly = nan(size(pGrid)); loBuyEarly = nan(size(pGrid)); hiBuyEarly = nan(size(pGrid)); nEarly = zeros(size(pGrid));
        pBuyLate  = nan(size(pGrid)); loBuyLate  = nan(size(pGrid)); hiBuyLate  = nan(size(pGrid)); nLate  = zeros(size(pGrid));
        pSellEarly= nan(size(pGrid)); loSellEarly= nan(size(pGrid)); hiSellEarly= nan(size(pGrid));
        pSellLate = nan(size(pGrid)); loSellLate = nan(size(pGrid)); hiSellLate = nan(size(pGrid));

        for iP = 1:numel(pGrid)
            p = pGrid(iP);
            inP = (D.prePortfolio == p);

            idxE = inP & isEarly;
            idxL = inP & isLate;

            nE = sum(idxE);
            nL = sum(idxL);

            nEarly(iP) = nE;
            nLate(iP)  = nL;

            if nE >= minCountPerP
                kBuyE = sum(D.actionCode(idxE) == ACTION_BUY);
                [phat, lo, hi] = wilsonCI(kBuyE, nE, 0.05);
                pBuyEarly(iP) = phat; loBuyEarly(iP) = lo; hiBuyEarly(iP) = hi;

                kSellE = sum(D.actionCode(idxE) == ACTION_SELL);
                [phat, lo, hi] = wilsonCI(kSellE, nE, 0.05);
                pSellEarly(iP) = phat; loSellEarly(iP) = lo; hiSellEarly(iP) = hi;
            end

            if nL >= minCountPerP
                kBuyL = sum(D.actionCode(idxL) == ACTION_BUY);
                [phat, lo, hi] = wilsonCI(kBuyL, nL, 0.05);
                pBuyLate(iP) = phat; loBuyLate(iP) = lo; hiBuyLate(iP) = hi;

                kSellL = sum(D.actionCode(idxL) == ACTION_SELL);
                [phat, lo, hi] = wilsonCI(kSellL, nL, 0.05);
                pSellLate(iP) = phat; loSellLate(iP) = lo; hiSellLate(iP) = hi;
            end

            if p==0 % cannot Sell at 0
                pSellEarly(iP) = NaN; loSellEarly(iP)=NaN; hiSellEarly(iP)=NaN;
                pSellLate(iP)  = NaN; loSellLate(iP) =NaN; hiSellLate(iP) =NaN;
                continue;
            end
        end

        resIdx = resIdx + 1;
        results(resIdx).Monkey = mID;
        results(resIdx).Role   = roleName;
        results(resIdx).pGrid = pGrid;
        results(resIdx).sessLo = sessLo;
        results(resIdx).sessHi = sessHi;

        results(resIdx).pBuyEarly = pBuyEarly;
        results(resIdx).loBuyEarly = loBuyEarly;
        results(resIdx).hiBuyEarly = hiBuyEarly;
        results(resIdx).pBuyLate = pBuyLate;
        results(resIdx).loBuyLate = loBuyLate;
        results(resIdx).hiBuyLate = hiBuyLate;

        results(resIdx).pSellEarly = pSellEarly;
        results(resIdx).loSellEarly = loSellEarly;
        results(resIdx).hiSellEarly = hiSellEarly;
        results(resIdx).pSellLate = pSellLate;
        results(resIdx).loSellLate = loSellLate;
        results(resIdx).hiSellLate = hiSellLate;

        results(resIdx).nEarlyByP = nEarly;
        results(resIdx).nLateByP  = nLate;
    end
end

if isempty(results)
    warning('No analyzable rows after role/session filtering. Nothing to plot.');
    return;
end

% --- FIG 1 ---
figure('Color','w','Name','FIG 1: P(BUY | portfolio size) + occupancy (top/bottom)');
tl1 = tiledlayout(2, numel(results), 'TileSpacing','compact', 'Padding','compact');

for ir = 1:numel(results)
    R = results(ir);

    % TOP
    axTop = nexttile(tl1, ir); hold(axTop,'on');
    he = plotWithCI(axTop, R.pGrid, R.pBuyEarly, R.loBuyEarly, R.hiBuyEarly, '-');
    hl = plotWithCI(axTop, R.pGrid, R.pBuyLate,  R.loBuyLate,  R.hiBuyLate,  '--');

    ylabel(axTop, 'P(BUY | portfolio size)');
    if LIVE_ROLE_MODE == "separate"
        roleTag = string(R.Role);
    else
        roleTag = "All";
    end
    title(axTop, sprintf('Monkey %s | Role %s', string(R.Monkey), roleTag));
    ylim(axTop, [0 1]); xlim(axTop, [0 PmaxPlot]);
    axTop.XTickLabel = [];

    legend(axTop, [he hl], ...
        {sprintf('Early sessions (<=%d%%)', earlyQuant*100), sprintf('Late sessions (>=%d%%)', lateQuant*100)}, ...
        'Location','southeast', 'Box','off');
    hold(axTop,'off');

    % BOTTOM
    axBot = nexttile(tl1, numel(results) + ir); hold(axBot,'on');

    x      = R.pGrid(:) + 0.5;
    yEarly = R.nEarlyByP(:);
    yLate  = R.nLateByP(:);

    b = bar(axBot, x, [yEarly yLate], 1.0, 'stacked');
    b(1).EdgeColor = 'none'; b(2).EdgeColor = 'none';
    yline(axBot, minCountPerP, '--k', 'LineWidth',2);

    xlim(axBot, [0 PmaxPlot]);
    xlabel(axBot, 'Portfolio size (shares)');
    ylabel(axBot, 'Occupancy');
    hold(axBot,'off');
end

% Make bottom row shorter by manually resizing axes (older MATLAB compatible)
drawnow;
for ir = 1:numel(results)
    axTop = nexttile(tl1, ir);
    axBot = nexttile(tl1, numel(results) + ir);

    pTop = axTop.Position;
    pBot = axBot.Position;

    gap = 0.01 * pTop(4);     % small vertical gap (relative)
    botFrac = 0.18;           % bottom height as fraction of combined height

    y0 = pBot(2);
    y1 = pTop(2) + pTop(4);

    H  = y1 - y0;
    Hb = botFrac * H;
    Ht = H - Hb - gap;

    axBot.Position = [pTop(1), y0,     pTop(3), Hb];
    axTop.Position = [pTop(1), y0+Hb+gap, pTop(3), Ht];
end

if exist('sgtitle','file')
    sgtitle(tl1, 'Setpoint test: BUY probability vs inventory (top) + occupancy (bottom)', 'Interpreter','none');
end

% --- FIG 2 ---
figure('Color','w','Name','FIG 2: P(SELL | portfolio size) + occupancy (top/bottom)');
tl2 = tiledlayout(2, numel(results), 'TileSpacing','compact', 'Padding','compact');

for ir = 1:numel(results)
    R = results(ir);

    % TOP
    axTop = nexttile(tl2, ir); hold(axTop,'on');
    he = plotWithCI(axTop, R.pGrid, R.pSellEarly, R.loSellEarly, R.hiSellEarly, '-');
    hl = plotWithCI(axTop, R.pGrid, R.pSellLate,  R.loSellLate,  R.hiSellLate,  '--');

    ylabel(axTop, 'P(SELL | portfolio size)');
    if LIVE_ROLE_MODE == "separate"
        roleTag = string(R.Role);
    else
        roleTag = "All";
    end
    title(axTop, sprintf('Monkey %s | Role %s', string(R.Monkey), roleTag));
    ylim(axTop, [0 1]); xlim(axTop, [0 PmaxPlot]);
    axTop.XTickLabel = [];

    legend(axTop, [he hl], ...
        {sprintf('Early sessions (<=%d%%)', earlyQuant*100), sprintf('Late sessions (>=%d%%)', lateQuant*100)}, ...
        'Location','northeast', 'Box','off');
    hold(axTop,'off');

    % BOTTOM
    axBot = nexttile(tl2, numel(results) + ir); hold(axBot,'on');

    x      = R.pGrid(:) + 0.5;
    yEarly = R.nEarlyByP(:);
    yLate  = R.nLateByP(:);

    b = bar(axBot, x, [yEarly yLate], 1.0, 'stacked');
    b(1).EdgeColor = 'none'; b(2).EdgeColor = 'none';

    xlim(axBot, [0 PmaxPlot]);
    xlabel(axBot, 'Portfolio size (shares)');
    ylabel(axBot, 'Occupancy');
    hold(axBot,'off');
end

% Make bottom row shorter (older MATLAB compatible)
drawnow;
for ir = 1:numel(results)
    axTop = nexttile(tl2, ir);
    axBot = nexttile(tl2, numel(results) + ir);

    pTop = axTop.Position;
    pBot = axBot.Position;

    gap = 0.01 * pTop(4);
    botFrac = 0.18;

    y0 = pBot(2);
    y1 = pTop(2) + pTop(4);

    H  = y1 - y0;
    Hb = botFrac * H;
    Ht = H - Hb - gap;

    axBot.Position = [pTop(1), y0,     pTop(3), Hb];
    axTop.Position = [pTop(1), y0+Hb+gap, pTop(3), Ht];
end

if exist('sgtitle','file')
    sgtitle(tl2, 'Setpoint test: SELL probability vs inventory (top) + occupancy (bottom)', 'Interpreter','none');
end

%% ========================= LOCAL FUNCTIONS =================================
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

function hLine = plotWithCI(ax, x, y, lo, hi, lineStyle)
% Plot y vs x with CI bands, but DO NOT bridge gaps where data is missing.
% CI band is drawn in contiguous x segments only.

    ok = isfinite(x) & isfinite(y) & isfinite(lo) & isfinite(hi);
    x = x(ok); y = y(ok); lo = lo(ok); hi = hi(ok);
    if isempty(x), return; end

    % Plot the mean line first to grab its color from MATLAB's color order
    hLine = plot(ax, x, y, lineStyle, 'LineWidth', 1.8);
    c = hLine.Color;

    % Find contiguous segments (for integer pGrid, contiguity is diff==1)
    breaks = [1; find(diff(x) > 1) + 1; numel(x) + 1];

    for s = 1:numel(breaks)-1
        idx = breaks(s):(breaks(s+1)-1);
        if numel(idx) < 2, continue; end  % need at least 2 points to fill

        hFill = fill(ax, ...
            [x(idx); flipud(x(idx))], ...
            [lo(idx); flipud(hi(idx))], ...
            c, 'EdgeColor','none', 'FaceAlpha', 0.20);

        uistack(hFill, 'bottom'); % keep fill behind the line
    end

    uistack(hLine, 'top'); % ensure the mean line stays on top
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
