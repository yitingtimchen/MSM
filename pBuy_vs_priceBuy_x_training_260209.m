%% MSM: Is BUY becoming more price-sensitive over training? (with equation on plot; beta_price highlighted)
% -------------------------------------------------------------------------
% QUESTION
%   Are monkeys simply buying more over sessions (action bias),
%   OR are they buying more selectively when BUY price is low (price sensitivity),
%   and does that price sensitivity strengthen with training?
%
% DATA SOURCE
%   - Uses your condition-pack structure: load(..._condition_pack.mat) -> C.eventTables{f}
%   - Reads subject choices from etab.option
%   - Reads trialNum and priceBuy from etab.trialNum, etab.priceBuy
%
% MODEL (per monkey; logistic regression / binomial GLM)
%   Let Y_i = 1 if BUY on trial i else 0.
%   We fit:
%     p(BUY)_i = sigmoid( β0 + βP * z(priceBuy)_i + βS * z(sessionIdx)_i + βPS * z(priceBuy)_i*z(sessionIdx)_i )
%   where sigmoid(x)=1/(1+exp(-x)), and z(.) is z-scored WITHIN MONKEY.
%
% KEY COEFFICIENTS
%   βP  (price-only term):  main effect of z(priceBuy)
%       - negative => higher priceBuy reduces BUY probability (price sensitivity)
%   βS  (training-only):    main effect of z(sessionIdx)
%       - positive => BUY increases over training (action bias / increased BUY tendency)
%   βPS (interaction):      does price sensitivity change over training?
%       - negative => sensitivity strengthens (price effect becomes more negative) with training
%
% WHAT THE PLOT SHOWS (per monkey)
%   - Points: binned observed P(BUY) vs priceBuy in EARLY vs LATE sessions
%   - Lines: model-predicted P(BUY) vs priceBuy at representative early vs late sessionIdx
%   - Text: equation in probability form (not logit), with βP highlighted + its estimate/p-value
%
% NOTES / ASSUMPTIONS
%   - option encoding is decoded by decodeOptionToAction(); update if your encoding differs.
%   - sessionIdx is chronological order within each monkey based on Session string (yyyy-mm-dd-sod).
%   - This analysis ignores opponent rows; it uses the subject option field only.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\Tim\\condition_packs_behavior_only_from_dataTabFULL'; % <-- EDIT

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

% Optional filters (set to "all" to include everything)
FILTER_TREATMENT   = "all";  % "all" or "baseline"/"OT"/"Saline"
FILTER_COND_BASE   = "all";  % "all" or "AI"/"Replay"/"Decoy"/"Live"

% Binning config for descriptive plot
nPriceBins   = 10;    % quantile bins
earlyQuant   = 0.33;  % early sessions = bottom 33% of sessionIdx
lateQuant    = 0.67;  % late  sessions = top    33% of sessionIdx
minBinCount  = 15;    % minimum trials per bin to plot (per early/late group)

% Plot styling
betaP_color = [0 0.6 0];      % green highlight for beta_price line
lw_curve    = 2.0;
ms_points   = 6;

%% ------------------------ REQUIRED ETAB FIELD MAP --------------------------
FIELD_option       = 'option';
FIELD_trialNum     = 'trialNum';
FIELD_priceBuy     = 'priceBuy';

FIELD_year         = 'year';
FIELD_month        = 'month';
FIELD_day          = 'day';
FIELD_SessionOfDay = 'SessionOfDay';
FIELD_monkey       = 'monkey';

%% --------------------- BUILD MASTER TRIAL-LEVEL TABLE ----------------------
trialTableAll = table();

fprintf('\n=== BUILDING trialTableAll (trial-level rows) ===\n');

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

            needed = {FIELD_option, FIELD_trialNum, FIELD_priceBuy, ...
                      FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay, FIELD_monkey};
            if ~all(ismember(needed, etab.Properties.VariableNames))
                warning('Missing required fields; skipping a file in %s.', cname);
                continue;
            end

            if isempty(etab.(FIELD_option){1}) || isempty(etab.(FIELD_trialNum){1}) || isempty(etab.(FIELD_priceBuy){1})
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

            % Monkey identity
            if isempty(etab.(FIELD_monkey){1})
                warning('Missing monkey field; skipping a file in %s.', cname);
                continue;
            end
            monkeyStr = string(etab.(FIELD_monkey){1}(1));

            % Trial vectors
            optionRaw = etab.(FIELD_option){1}(:);
            trialNum  = double(etab.(FIELD_trialNum){1}(:));
            priceBuy  = double(etab.(FIELD_priceBuy){1}(:));

            nT = min([numel(optionRaw), numel(trialNum), numel(priceBuy)]);
            optionRaw = optionRaw(1:nT);
            trialNum  = trialNum(1:nT);
            priceBuy  = priceBuy(1:nT);

            actionCode = decodeOptionToAction(optionRaw); % 1=BUY,2=HOLD,3=SELL
            isBuy      = (actionCode == ACTION_BUY);

            % Keep scorable trials with finite price and trialNum
            keep = isfinite(actionCode) & isfinite(trialNum) & isfinite(priceBuy);
            if ~any(keep), continue; end

            isBuy      = isBuy(keep);
            trialNum   = trialNum(keep);
            priceBuy   = priceBuy(keep);

            nKeep = numel(priceBuy);

            newRows = table( ...
                repmat(string(treatName), nKeep, 1), ...
                repmat(string(cname),     nKeep, 1), ...
                repmat(baseCond,          nKeep, 1), ...
                repmat(monkeyStr,         nKeep, 1), ...
                repmat(sessionStr,        nKeep, 1), ...
                repmat(sessionDate,       nKeep, 1), ...
                trialNum(:), ...
                priceBuy(:), ...
                double(isBuy(:)), ...
                'VariableNames', { ...
                    'Treatment','ConditionRaw','ConditionBase', ...
                    'Monkey','Session','SessionDate', ...
                    'trialNum','priceBuy','isBuy' ...
                } ...
            );

            trialTableAll = [trialTableAll; newRows]; %#ok<AGROW>
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

trialTableAll.Monkey        = categorical(string(trialTableAll.Monkey));
trialTableAll.Treatment     = categorical(string(trialTableAll.Treatment));
trialTableAll.ConditionBase = categorical(string(trialTableAll.ConditionBase), baseCondList);

monkeys = categories(trialTableAll.Monkey);

%% ---------------------- SESSION INDEX (training order) ---------------------
% Assign sessionIdx within each monkey using unique sessions sorted by Session string.
trialTableAll.sessionIdx = nan(height(trialTableAll),1);

sessKeyTab = unique(trialTableAll(:, {'Monkey','Session'}));
sessKeyTab = sortrows(sessKeyTab, {'Monkey','Session'});

for im = 1:numel(monkeys)
    mID = monkeys{im};
    rowsM = sessKeyTab.Monkey == mID;
    sessList = sessKeyTab.Session(rowsM);
    sessIdxLocal = (1:numel(sessList))';

    for k = 1:numel(sessList)
        mask = (trialTableAll.Monkey == mID) & (trialTableAll.Session == sessList(k));
        trialTableAll.sessionIdx(mask) = sessIdxLocal(k);
    end
end

%% ======================= MODEL + PLOTS PER MONKEY ==========================
modelSummary = table();

figure('Color','w','Name','BUY vs priceBuy with training (per monkey)');
tl = tiledlayout(1, numel(monkeys), 'TileSpacing','compact', 'Padding','compact');

for im = 1:numel(monkeys)
    mID = monkeys{im};
    D = trialTableAll(trialTableAll.Monkey == mID, :);

    % Z-score predictors within monkey
    muPrice = mean(D.priceBuy);
    sdPrice = std(D.priceBuy);
    muSess  = mean(D.sessionIdx);
    sdSess  = std(D.sessionIdx);

    D.zPriceBuy   = (D.priceBuy   - muPrice) / sdPrice;
    D.zSessionIdx = (D.sessionIdx - muSess)  / sdSess;

    % Logistic regression: isBuy ~ zPriceBuy * zSessionIdx
    mdl = fitglm(D, 'isBuy ~ zPriceBuy * zSessionIdx', 'Distribution','binomial', 'Link','logit');

    % Extract coefficients
    coef = mdl.Coefficients;
    bP   = coef.Estimate(strcmp(mdl.CoefficientNames,'zPriceBuy'));
    bS   = coef.Estimate(strcmp(mdl.CoefficientNames,'zSessionIdx'));
    bPS  = coef.Estimate(strcmp(mdl.CoefficientNames,'zPriceBuy:zSessionIdx'));

    pP   = coef.pValue(strcmp(mdl.CoefficientNames,'zPriceBuy'));
    pS   = coef.pValue(strcmp(mdl.CoefficientNames,'zSessionIdx'));
    pPS  = coef.pValue(strcmp(mdl.CoefficientNames,'zPriceBuy:zSessionIdx'));

    modelSummary = [modelSummary; table( ...
        string(mID), height(D), bP, pP, bS, pS, bPS, pPS, ...
        'VariableNames', {'Monkey','nTrials','b_zPriceBuy','p_zPriceBuy','b_zSessionIdx','p_zSessionIdx','b_interaction','p_interaction'})]; %#ok<AGROW>

    % EARLY vs LATE session split
    sessLo = quantile(D.sessionIdx, earlyQuant);
    sessHi = quantile(D.sessionIdx, lateQuant);
    isEarly = D.sessionIdx <= sessLo;
    isLate  = D.sessionIdx >= sessHi;

    % Price bins (quantiles)
    edges = quantile(D.priceBuy, linspace(0,1,nPriceBins+1));
    edges(1)   = edges(1) - 1e-6;
    edges(end) = edges(end) + 1e-6;

    [~, binID] = histc(D.priceBuy, edges); %#ok<HISTC>

    binCenters = nan(nPriceBins,1);
    pEarly = nan(nPriceBins,1);
    pLate  = nan(nPriceBins,1);

    for b = 1:nPriceBins
        inBin = (binID == b);
        if ~any(inBin), continue; end
        binCenters(b) = mean(D.priceBuy(inBin));

        ee = inBin & isEarly;
        ll = inBin & isLate;

        if sum(ee) >= minBinCount
            pEarly(b) = mean(D.isBuy(ee));
        end
        if sum(ll) >= minBinCount
            pLate(b) = mean(D.isBuy(ll));
        end
    end

    % Model-predicted curves at representative early vs late sessionIdx
    priceGrid  = linspace(min(D.priceBuy), max(D.priceBuy), 200)';
    zPriceGrid = (priceGrid - muPrice) / sdPrice;

    zSessEarly = (sessLo - muSess) / sdSess;
    zSessLate  = (sessHi - muSess) / sdSess;

    TpredEarly = table(zPriceGrid, repmat(zSessEarly, numel(zPriceGrid),1), 'VariableNames', {'zPriceBuy','zSessionIdx'});
    TpredLate  = table(zPriceGrid, repmat(zSessLate,  numel(zPriceGrid),1), 'VariableNames', {'zPriceBuy','zSessionIdx'});

    pHatEarly = predict(mdl, TpredEarly);
    pHatLate  = predict(mdl, TpredLate);

    % Plot
    ax = nexttile(tl, im); hold(ax,'on');

    hCurveEarly = plot(ax, priceGrid, pHatEarly, '-',  'LineWidth', lw_curve);
    hCurveLate  = plot(ax, priceGrid, pHatLate,  '--', 'LineWidth', lw_curve);

    plot(ax, binCenters, pEarly, 'o', 'MarkerSize', ms_points, 'LineWidth', 1.2, ...
        'Color', hCurveEarly.Color, 'MarkerEdgeColor', hCurveEarly.Color, 'MarkerFaceColor', hCurveEarly.Color);

    plot(ax, binCenters, pLate,  's', 'MarkerSize', ms_points, 'LineWidth', 1.2, ...
        'Color', hCurveLate.Color,  'MarkerEdgeColor', hCurveLate.Color,  'MarkerFaceColor', hCurveLate.Color);

    xlabel(ax, 'priceBuy ($)');
    ylabel(ax, 'P(BUY)');
    title(ax, sprintf('Monkey %s', string(mID)));
    ylim(ax, [0 1]);

    % Force TeX for this annotation (most reliable in MATLAB figures)
    % Split into two lines using eta:
    %   eta = ...
    %   p(BUY)=1/(1+exp(-eta))
    
    eqn1 = "\eta = \beta_0 + \beta_P z(priceBuy) + \beta_S z(sessionIdx) + \beta_{PS} z(priceBuy)z(sessionIdx)";
    eqn2 = "p(BUY)=1/(1+exp(-\eta))";
    
    eqnP  = "\beta_P="  + sprintf("%.3g", bP)  + " (p=" + sprintf("%.3g", pP)  + ")";  % price-only term (highlight)
    eqnS  = "\beta_S="  + sprintf("%.3g", bS)  + " (p=" + sprintf("%.3g", pS)  + ")";
    eqnPS = "\beta_{PS}="+ sprintf("%.3g", bPS) + " (p=" + sprintf("%.3g", pPS) + ")";
    
    text(ax, 0.02, 0.98, eqn1, 'Units','normalized', 'Interpreter','tex', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 10);
    text(ax, 0.02, 0.92, eqn2, 'Units','normalized', 'Interpreter','tex', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 10);
    
    text(ax, 0.02, 0.84, eqnP, 'Units','normalized', 'Interpreter','tex', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 10, 'Color', betaP_color);
    text(ax, 0.02, 0.78, eqnS, 'Units','normalized', 'Interpreter','tex', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 10);
    text(ax, 0.02, 0.72, eqnPS, 'Units','normalized', 'Interpreter','tex', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 10);
    
    legend(ax, ...
        {sprintf('Early sessions (<=%.0f%%)', earlyQuant*100), sprintf('Late sessions (>=%.0f%%)', lateQuant*100)}, ...
        'Location','southeast', 'Box','off');

    hold(ax,'off');
end

if exist('sgtitle','file')
    sgtitle(tl, 'Does BUY become more price-sensitive with training?', 'Interpreter','none');
end

%% Print model summary
disp('=== Logistic model summary per monkey: isBuy ~ zPriceBuy * zSessionIdx ===');
disp(modelSummary);

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
%     If numeric is different, update mapping explicitly here.

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

function v = ACTION_BUY
    v = 1;
end
function v = ACTION_HOLD
    v = 2;
end
function v = ACTION_SELL
    v = 3;
end
