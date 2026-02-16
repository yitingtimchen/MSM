%% MSM: Is BUY becoming more price-sensitive over training?
% (UPDATED: include Live sessions as subject + as opponent using *Opp fields)

clear; clc; close all;

%% ------------------------------- CONFIG -----------------------------------
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL'; % <-- EDIT

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

% Live role handling:
%   "subj_only" -> include only Subj stream in Live
%   "separate"  -> include Subj and Opp as separate role conditions
%   "combined"  -> include both streams but pool them as one Live role;
%                  session indexing uses role-qualified session keys so
%                  effective session count doubles for Live per monkey.
LIVE_ROLE_MODE = "separate";  % "subj_only" | "separate" | "combined"

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

% Opponent fields (Live only; role = Opp)
FIELD_optionOpp    = 'optionOpp';
FIELD_trialNumOpp  = 'trialNumOpp';
FIELD_priceBuyOpp  = 'priceBuyOpp';

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

            % Always-needed session fields
            neededSession = {FIELD_year, FIELD_month, FIELD_day, FIELD_SessionOfDay};

            if ~all(ismember(neededSession, etab.Properties.VariableNames))
                warning('Missing required session fields; skipping a file in %s.', cname);
                continue;
            end

            % Session identity
            y   = etab.(FIELD_year){1}(:);
            mth = etab.(FIELD_month){1}(:);
            d   = etab.(FIELD_day){1}(:);
            sod = etab.(FIELD_SessionOfDay){1}(:);
            if isempty(y) || isempty(mth) || isempty(d) || isempty(sod), continue; end

            sessionStr  = string(sprintf('%04d-%02d-%02d-%02d', y(1), mth(1), d(1), sod(1)));
            sessionDate = datetime(y(1), mth(1), d(1));

            % Exclusions
            if ~isempty(EXCLUDE_SESSIONS) && any(sessionStr == EXCLUDE_SESSIONS)
                excludedSessSeen(end+1,1) = sessionStr; %#ok<AGROW>
                continue;
            end

            % Decide which roles to extract
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

                if role == "Subj"
                    optField   = FIELD_option;
                    trlField   = FIELD_trialNum;
                    prcField   = FIELD_priceBuy;
                    monkeyStr  = string(etab.(FIELD_monkey){1}(1));
                else
                    optField   = FIELD_optionOpp;
                    trlField   = FIELD_trialNumOpp;
                    prcField   = FIELD_priceBuyOpp;
                    monkeyStr  = inferOpponentMonkey(string(etab.(FIELD_monkey){1}(1)));
                end

                neededRole = {optField, trlField, prcField, FIELD_monkey};
                if ~all(ismember(neededRole, etab.Properties.VariableNames))
                    % For non-Live this won't happen; for Live Opp just skip if absent.
                    continue;
                end

                if isempty(etab.(optField){1}) || isempty(etab.(trlField){1}) || isempty(etab.(prcField){1})
                    continue;
                end

                % Monkey identity (actor = role-specific)
                if isempty(etab.(FIELD_monkey){1})
                    continue;
                end
                if ~isMonkey12(monkeyStr), continue; end

                % Trial vectors (role-specific)
                optionRaw = etab.(optField){1}(:);
                trialNum  = double(etab.(trlField){1}(:));
                priceBuy  = double(etab.(prcField){1}(:));

                nT = min([numel(optionRaw), numel(trialNum), numel(priceBuy)]);
                optionRaw = optionRaw(1:nT);
                trialNum  = trialNum(1:nT);
                priceBuy  = priceBuy(1:nT);

                actionCode = decodeOptionToAction(optionRaw); % 1=BUY,2=HOLD,3=SELL
                isBuy      = (actionCode == ACTION_BUY);

                % Keep scorable trials with finite price and trialNum
                keep = isfinite(actionCode) & isfinite(trialNum) & isfinite(priceBuy);
                if ~any(keep), continue; end

                isBuy    = isBuy(keep);
                trialNum = trialNum(keep);
                priceBuy = priceBuy(keep);

                nKeep = numel(priceBuy);

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
                    repmat(monkeyStr,         nKeep, 1), ...
                    repmat(sessionStr,        nKeep, 1), ...
                    repmat(analysisSession,   nKeep, 1), ...
                    repmat(sessionDate,       nKeep, 1), ...
                    trialNum(:), ...
                    priceBuy(:), ...
                    double(isBuy(:)), ...
                    'VariableNames', { ...
                        'Treatment','ConditionRaw','ConditionBase','Role', ...
                        'Monkey','Session','AnalysisSession','SessionDate', ...
                        'trialNum','priceBuy','isBuy' ...
                    } ...
                );

                trialTableAll = [trialTableAll; newRows]; %#ok<AGROW>
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
if isempty(trialTableAll)
    warning('No trials found. Nothing to analyze.');
    return;
end

trialTableAll.Monkey        = categorical(string(trialTableAll.Monkey));
trialTableAll.Treatment     = categorical(string(trialTableAll.Treatment));
trialTableAll.ConditionBase = categorical(string(trialTableAll.ConditionBase), baseCondList);
trialTableAll.Role          = categorical(string(trialTableAll.Role), ["Subj","Opp","Combined"]);

monkeys = categories(trialTableAll.Monkey);

%% ---------------------- SESSION INDEX (training order) ---------------------
% Assign sessionIdx within each monkey using AnalysisSession.
% In LIVE_ROLE_MODE=="combined", role-qualified keys double Live density.
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

%% ======================= MODEL + PLOTS PER MONKEY ==========================
modelSummary = table();

figure('Color','w','Name','BUY vs priceBuy with training (per monkey; Live Subj vs Opp)');
tl = tiledlayout(1, numel(monkeys), 'TileSpacing','compact', 'Padding','compact');

for im = 1:numel(monkeys)
    mID = monkeys{im};
    Dall = trialTableAll(trialTableAll.Monkey == mID, :);

    % Early vs late thresholds computed over ALL sessions for this monkey
    sessLo = quantile(Dall.sessionIdx, earlyQuant);
    sessHi = quantile(Dall.sessionIdx, lateQuant);

    ax = nexttile(tl, im); hold(ax,'on');

    % Plot per role (if separated) or pooled (if combined/subj_only)
    if LIVE_ROLE_MODE == "separate"
        roleOrder = ["Subj","Opp"];
    elseif LIVE_ROLE_MODE == "combined"
        roleOrder = "Combined";
    else
        roleOrder = "Subj";
    end
    hForLegend = gobjects(0,1);
    legStr     = strings(0,1);

    % Grab base color from axes ColorOrder in a stable way
    co = ax.ColorOrder;
    baseColor = co( mod(im-1, size(co,1)) + 1, : );

    for r = 1:numel(roleOrder)
        role = roleOrder(r);
        D = Dall(Dall.Role == role, :);
        if isempty(D), continue; end

        if role == "Subj"
            thisColor = baseColor;
            roleTag   = "Live as subject";
        elseif role == "Opp"
            thisColor = lightenColor(baseColor, 0.55); % lighter
            roleTag   = "Live as opponent";
        else
            thisColor = baseColor;
            roleTag   = "Live pooled (Subj+Opp)";
        end

        % Z-score predictors within monkey+role
        muPrice = mean(D.priceBuy);
        sdPrice = std(D.priceBuy);
        muSess  = mean(D.sessionIdx);
        sdSess  = std(D.sessionIdx);

        if sdPrice == 0 || sdSess == 0
            continue;
        end

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
            string(mID), string(role), height(D), bP, pP, bS, pS, bPS, pPS, ...
            'VariableNames', {'Monkey','Role','nTrials','b_zPriceBuy','p_zPriceBuy','b_zSessionIdx','p_zSessionIdx','b_interaction','p_interaction'})]; %#ok<AGROW>

        % Early vs late masks (based on monkey-wide thresholds)
        isEarly = D.sessionIdx <= sessLo;
        isLate  = D.sessionIdx >= sessHi;

        % Price bins (quantiles within role)
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

        % Plot (force role color; use linestyle for early/late)
        hCurveEarly = plot(ax, priceGrid, pHatEarly, '-',  'LineWidth', lw_curve, 'Color', thisColor);
        hCurveLate  = plot(ax, priceGrid, pHatLate,  '--', 'LineWidth', lw_curve, 'Color', thisColor);

        plot(ax, binCenters, pEarly, 'o', 'MarkerSize', ms_points, 'LineWidth', 1.2, ...
            'Color', thisColor, 'MarkerEdgeColor', thisColor, 'MarkerFaceColor', thisColor);

        plot(ax, binCenters, pLate,  's', 'MarkerSize', ms_points, 'LineWidth', 1.2, ...
            'Color', thisColor, 'MarkerEdgeColor', thisColor, 'MarkerFaceColor', thisColor);

        % Legend handles (one per role for early/late)
        hForLegend(end+1,1) = hCurveEarly; %#ok<AGROW>
        legStr(end+1,1) = sprintf('%s (early)', roleTag); %#ok<AGROW>
        hForLegend(end+1,1) = hCurveLate; %#ok<AGROW>
        legStr(end+1,1) = sprintf('%s (late)', roleTag); %#ok<AGROW>

        % Add equation text ONLY for Subj (keeps plot readable)
        if role == "Subj"
            eqn1 = "\eta = \beta_0 + \beta_P z(priceBuy) + \beta_S z(sessionIdx) + \beta_{PS} z(priceBuy)z(sessionIdx)";
            eqn2 = "p(BUY)=1/(1+exp(-\eta))";

            eqnP  = "\beta_P="  + sprintf("%.3g", bP)  + " (p=" + sprintf("%.3g", pP)  + ")";
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
        end
    end

    xlabel(ax, 'priceBuy ($)');
    ylabel(ax, 'P(BUY)');
    title(ax, sprintf('Monkey %s', string(mID)));
    ylim(ax, [0 1]);

    if ~isempty(hForLegend)
        legend(ax, hForLegend, cellstr(legStr), 'Location','southeast', 'Box','off');
    end

    hold(ax,'off');
end

if exist('sgtitle','file')
    sgtitle(tl, 'Does BUY become more price-sensitive with training? (Live: subject vs opponent)', 'Interpreter','none');
end

%% Print model summary
disp('=== Logistic model summary per monkey+role: isBuy ~ zPriceBuy * zSessionIdx ===');
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

function c2 = lightenColor(c1, amt)
% amt in [0,1]; higher = lighter
    amt = max(0, min(1, amt));
    c2 = c1 + (1 - c1) * amt;
end

function opp = inferOpponentMonkey(subj)
% Infer opponent monkey ID/name from subject monkey.
    s = upper(strtrim(string(subj)));
    v = str2double(s);
    if isfinite(v)
        if v == 1, opp = "2"; return; end
        if v == 2, opp = "1"; return; end
    end
    if s == "M1", opp = "M2"; return; end
    if s == "M2", opp = "M1"; return; end
    error("inferOpponentMonkey: cannot infer opponent from subject '%s'. Update mapping.", s);
end

function ok = isMonkey12(monkeyVal)
% Keep only rows where monkey is 1/2 (or M1/M2 string equivalent).
    s = upper(strtrim(string(monkeyVal)));
    v = str2double(s);
    ok = (isfinite(v) && any(v == [1 2])) || any(s == ["M1","M2"]);
end
