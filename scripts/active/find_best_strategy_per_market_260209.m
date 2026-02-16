%% MSM stock market task: TWO benchmark policies + publication-ready plots
% -------------------------------------------------------------------------
% WHAT THIS SCRIPT DOES
%   This script computes and visualizes TWO different benchmarks for the MSM task:
%
%   (A) DP OPTIMAL POLICY (STATE-DEPENDENT, NORMATIVE)
%       - Full dynamic programming (finite-horizon, inventory-constrained).
%       - Policy depends on (Trial, Portfolio-before-choice).
%       - Output: policy heatmaps + implied optimal trajectories from p0=0.
%
%   (B) FIXED-FUNDAMENTAL HEURISTIC (STATE-INDEPENDENT, FIRST-ORDER)
%       - "Fundamental value" per trial is defined as remaining expected dividends:
%           fundamentalValue(t) = muDiv_fundamental * (Tmax - t + 1)
%         (inclusive of the current trial because dividends are post-action).
%       - Heuristic best action per trial compares price schedules to that scalar:
%           If priceSell(t) > fundamentalValue(t)  -> SELL signal (overpriced)
%           Else if priceBuy(t) < fundamentalValue(t) -> BUY signal (underpriced)
%           Else -> HOLD
%       - This heuristic is NOT DP-optimal under inventory constraints; it is a
%         first-order reference used in some prior literature.
%
% IMPORTANT TERMINOLOGY
%   - DP benchmark: "optimal" (given the task mechanics and constraints).
%   - Fundamental benchmark: "heuristic" or "trial-wise reference" (not optimal).
%
% CONFIRMED TASK ASSUMPTIONS (used for DP)
%   - 15 trials per market (Tmax=15).
%   - Each trial has 2 rows: Subject first, Opponent second (always).
%   - priceBuy/priceSell invariant within each marketOrig template.
%   - Can trade +/-1 share per trial; no shorting (SELL not allowed at p=0).
%   - Starting portfolio always 0 and resets each market.
%   - No liquidation / cash-out after trial 15.
%   - Dividend paid post-action each trial, including trial 15.
%   - Dividend per share distribution in {0, 0.8, 2.8, 6} dollars; DP uses mean.
%
% DISPLAY / PIXEL CALIBRATION
%   - 1 dollar corresponds to 2.125 pixels (linear scaling).
%   - Scaling does NOT change argmax decisions; DP is computed in dollars.
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% --------------------------- User settings --------------------------------
data_xlsx_path  = 'C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx'; % <-- EDIT
data_sheet_name = 'Sheet1';                                                               % <-- EDIT

market_templates = 1:6;  % marketOrig templates
Tmax = 15;               % trials per market
Pmax = 15;               % portfolio states tracked (0..15 sufficient for 15 trials)

% DP reward parameters (dollars)
allowance_dollars = 80;  % constant per trial (adds to all actions equally; does not affect argmax)
dividend_values_dollars = [0 0.8 2.8 6];
muDiv_dollars = mean(dividend_values_dollars);  % = 2.4

% Fixed-fundamental heuristic parameter (set equal to muDiv by default)
% If you want to match a paper's expected dividend, override this value.
muDiv_fundamental_dollars = muDiv_dollars;

% Optional plotting in pixels (DP is still computed in dollars)
pxPerDollar = 2.125;
plot_in_pixels = false;

% Invariance check tolerance (rounding decimals)
priceTolDec = 6;

% Publication-consistent x ticks
xTicksTrials = 0:3:15;

%% --------------------------- Action codes ---------------------------------
ACTION_BUY  = 1;
ACTION_HOLD = 2;
ACTION_SELL = 3;
ACTION_AMBIG = 4; % only used for the heuristic if BUY and SELL signals conflict

actionNames = {'BUY','HOLD','SELL','AMBIG'};

% Plot colors (consistent across figures)
colorPriceBuy  = [0.15 0.55 0.95];   % blue-ish
colorPriceSell = [0.90 0.25 0.25];   % red-ish
colorFundamental = [0 0 0];          % black
colorTrajectory = [0 0.6 0];         % green (for DP-derived trajectories)

% Line styles: solid=Subject, dashed=Opponent (consistent everywhere)
lineStyleSubject  = '-';
lineStyleOpponent = '--';

PERSPECTIVE_SUBJECT  = 0;
PERSPECTIVE_OPPONENT = 1;
perspectiveNames = {'Subject','Opponent'};

%% --------------------------- Load data ------------------------------------
T = readtable(data_xlsx_path, 'Sheet', data_sheet_name);

% Keep only monkeys 1 and 2 (safety guard).
if ismember('monkey', T.Properties.VariableNames)
    T = T(ismember(double(T.monkey), [1 2]), :);
end

% Force numeric types (robust to Excel import quirks)
numVars = {'trialNum','trialBlock','marketOrig','priceBuy','priceSell','prePortfolio','option', ...
           'session','year','month','day','SessionOfDay'};
for i = 1:numel(numVars)
    if ismember(numVars{i}, T.Properties.VariableNames)
        T.(numVars{i}) = double(T.(numVars{i}));
    end
end

% Required columns
req = {'trialNum','marketOrig','priceBuy','priceSell'};
for i = 1:numel(req)
    assert(ismember(req{i}, T.Properties.VariableNames), 'Missing required column: %s', req{i});
end

%% --------------------------- Sort rows ------------------------------------
% Preserve original within-trial order using origRow (important for pairing)
T.origRow = (1:height(T))';
sortKeys = intersect({'year','month','day','SessionOfDay','session','marketOrig','trialBlock','trialNum','origRow'}, ...
                     T.Properties.VariableNames, 'stable');
if ~isempty(sortKeys)
    T = sortrows(T, sortKeys);
end

%% --------------------------- Define market blocks --------------------------
% Each market block = 15 trials x 2 players = 30 rows.
if all(ismember({'year','month','day','SessionOfDay','session'}, T.Properties.VariableNames))
    marketBlockID = findgroups(T.year, T.month, T.day, T.SessionOfDay, T.session);
elseif ismember('session', T.Properties.VariableNames)
    marketBlockID = findgroups(T.session);
elseif ismember('trialBlock', T.Properties.VariableNames)
    marketBlockID = ceil(T.trialBlock/30);
else
    marketBlockID = floor(((1:height(T))' - 1)/30) + 1;
end

%% --------------------------- Pairing: Subject vs Opponent -----------------
% Confirmed: within each (marketBlockID, trialNum), first row = Subject, second = Opponent.
pairGroup = findgroups(marketBlockID, T.trialNum);
nPerPair = splitapply(@numel, T.trialNum, pairGroup);
assert(all(nPerPair==2), ...
    'Pairing check failed: some (marketBlockID, trialNum) groups are not exactly 2 rows.');

T.perspectiveCode = nan(height(T),1); % 0=Subject, 1=Opponent
for g = 1:max(pairGroup)
    idx = find(pairGroup==g);
    T.perspectiveCode(idx) = (0:numel(idx)-1)';
end
assert(all(ismember(T.perspectiveCode, [0 1])), 'Perspective coding contains values other than 0/1.');

%% --------------------------- Sanity: one marketOrig per block --------------
blkGroup = findgroups(marketBlockID);
nUniqueMarket = splitapply(@(x) numel(unique(x(~isnan(x)))), T.marketOrig, blkGroup);
assert(all(nUniqueMarket==1), ...
    'Block sanity check failed: at least one market block contains >1 marketOrig.');

%% --------------------------- Assert invariance of price schedules ----------
% Invariant within (marketOrig, trialNum, perspectiveCode).
schedGroup = findgroups(T.marketOrig, T.trialNum, T.perspectiveCode);

uBuy  = splitapply(@(x) numel(unique(round(x(~isnan(x)), priceTolDec))), T.priceBuy,  schedGroup);
uSell = splitapply(@(x) numel(unique(round(x(~isnan(x)), priceTolDec))), T.priceSell, schedGroup);

assert(all(uBuy==1 & uSell==1), ...
    'Price invariance failed: priceBuy/priceSell vary within a fixed (marketOrig, trialNum, perspective).');

%% --------------------------- Build canonical schedules (ROBUST) ------------
% IMPORTANT: derive keys via splitapply to avoid any group/key misalignment.
marketOrigKey = splitapply(@(x) x(find(~isnan(x),1,'first')), T.marketOrig,      schedGroup);
trialNumKey   = splitapply(@(x) x(find(~isnan(x),1,'first')), T.trialNum,        schedGroup);
perspKey      = splitapply(@(x) x(find(~isnan(x),1,'first')), T.perspectiveCode, schedGroup);

priceBuyKey   = splitapply(@(x) x(find(~isnan(x),1,'first')), T.priceBuy,  schedGroup);
priceSellKey  = splitapply(@(x) x(find(~isnan(x),1,'first')), T.priceSell, schedGroup);

priceBuyStd   = splitapply(@(x) std(x(~isnan(x))),  T.priceBuy,  schedGroup);
priceSellStd  = splitapply(@(x) std(x(~isnan(x))),  T.priceSell, schedGroup);

priceScheduleTable = table(marketOrigKey, trialNumKey, perspKey, priceBuyKey, priceSellKey, priceBuyStd, priceSellStd, ...
    'VariableNames', {'marketOrig','trialNum','perspectiveCode','priceBuy','priceSell','stdBuy','stdSell'});

% Schedule arrays (dollars):
%   priceBuy_dollars(mi,pp,t)  and priceSell_dollars(mi,pp,t)
% where pp = 1 (Subject) or 2 (Opponent)
priceBuy_dollars  = nan(numel(market_templates), 2, Tmax);
priceSell_dollars = nan(numel(market_templates), 2, Tmax);

for mi = 1:numel(market_templates)
    m = market_templates(mi);
    for pp = 1:2
        perspCode = pp-1; % 0=Subject, 1=Opponent
        for t = 1:Tmax
            ii = (priceScheduleTable.marketOrig==m) & ...
                 (priceScheduleTable.trialNum==t) & ...
                 (priceScheduleTable.perspectiveCode==perspCode);
            assert(nnz(ii)==1, 'Missing/duplicate schedule entry for marketOrig=%d trial=%d perspective=%d', m, t, perspCode);
            priceBuy_dollars(mi,pp,t)  = priceScheduleTable.priceBuy(ii);
            priceSell_dollars(mi,pp,t) = priceScheduleTable.priceSell(ii);
        end
    end
end

%% --------------------------- Fixed-fundamental values ----------------------
% Fundamental value per share at trial t (dollars) for the heuristic benchmark.
% Because dividend is post-action and includes trial 15, remaining payments from trial t are:
%   (Tmax - t + 1) = number of remaining trials INCLUDING the current trial.
trialIdx = (1:Tmax)';
fundamental_dollars_byTrial = muDiv_fundamental_dollars .* (Tmax - trialIdx + 1);

%% ==========================================================================
%                           BENCHMARK (A): DP OPTIMAL
% ==========================================================================
% DP optimal policy depends on state (Trial, Portfolio-before-choice).
% policyDP_action{pp,mi} is [Tmax x (Pmax+1)] with actions coded as 1/2/3.
policyDP_action = cell(2, numel(market_templates));

for mi = 1:numel(market_templates)
    for pp = 1:2
        buySched  = squeeze(priceBuy_dollars(mi,pp,:));   % [Tmax x 1]
        sellSched = squeeze(priceSell_dollars(mi,pp,:));  % [Tmax x 1]

        pol   = nan(Tmax, Pmax+1);
        Vnext = zeros(Pmax+1,1);  % V_{t+1}(p) ; terminal value is 0 (no liquidation)

        for t = Tmax:-1:1
            V = -inf(Pmax+1,1);

            % Starting portfolio always 0; max feasible p at trial t is (t-1).
            pMaxFeasible = min(Pmax, t-1);

            for p = 0:pMaxFeasible
                % BUY -> p+1
                pB = p+1;
                qB = (allowance_dollars - buySched(t)) + muDiv_dollars*pB + Vnext(pB+1);

                % HOLD -> p
                pH = p;
                qH = allowance_dollars + muDiv_dollars*pH + Vnext(pH+1);

                % SELL -> p-1 (no shorting)
                if p>0
                    pS = p-1;
                    qS = (allowance_dollars + sellSched(t)) + muDiv_dollars*pS + Vnext(pS+1);
                else
                    qS = -inf;
                end

                [V(p+1), a] = max([qB, qH, qS]);
                pol(t,p+1) = a;
            end

            % Mark unreachable states as -inf to prevent accidental use
            V(pMaxFeasible+2:end) = -inf;
            Vnext = V;
        end

        policyDP_action{pp,mi} = pol;
    end
end

%% ==========================================================================
%                    BENCHMARK (B): FIXED-FUNDAMENTAL HEURISTIC
% ==========================================================================
% Heuristic "best action" depends ONLY on Trial (state-independent signal).
% We compute a per-trial action signal for each (market, perspective).
%
% Rule (priority = SELL first):
%   If priceSell(t) > fundamental(t) -> SELL
%   Else if priceBuy(t) < fundamental(t) -> BUY
%   Else -> HOLD
%
% If both SELL and BUY conditions are true simultaneously, label as AMBIG.
% (This is rare; it means sell price is above fundamental while buy price is below.)
policyFundamental_action_trialOnly = nan(2, numel(market_templates), Tmax); % [pp x mi x t]

for mi = 1:numel(market_templates)
    for pp = 1:2
        buySched  = squeeze(priceBuy_dollars(mi,pp,:));
        sellSched = squeeze(priceSell_dollars(mi,pp,:));

        for t = 1:Tmax
            isSell = sellSched(t) > fundamental_dollars_byTrial(t);
            isBuy  = buySched(t)  < fundamental_dollars_byTrial(t);

            if isSell && ~isBuy
                policyFundamental_action_trialOnly(pp,mi,t) = ACTION_SELL;
            elseif isBuy && ~isSell
                policyFundamental_action_trialOnly(pp,mi,t) = ACTION_BUY;
            elseif ~isBuy && ~isSell
                policyFundamental_action_trialOnly(pp,mi,t) = ACTION_HOLD;
            else
                policyFundamental_action_trialOnly(pp,mi,t) = ACTION_AMBIG;
            end
        end
    end
end

%% =========================== FIG A ========================================
% PRICE SCHEDULE INVARIANCE CHECK
% WHAT IT SHOWS
%   - The canonical BUY and SELL price schedules by trial for each market template.
%   - Subject vs Opponent shown via line style (solid vs dashed).
%   - Same color = same price type (BUY blue, SELL red).
% WHY IT EXISTS
%   - Sanity check: schedules should be invariant across sessions for a given marketOrig.
figure('Color','w','Name','FIG A: Price schedule invariance');
tlA = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for mi = 1:numel(market_templates)
    m = market_templates(mi);
    ax = nexttile(tlA, mi);

    subjBuy  = squeeze(priceBuy_dollars(mi,1,:));
    oppBuy   = squeeze(priceBuy_dollars(mi,2,:));
    subjSell = squeeze(priceSell_dollars(mi,1,:));
    oppSell  = squeeze(priceSell_dollars(mi,2,:));

    if plot_in_pixels
        subjBuy  = subjBuy*pxPerDollar;   oppBuy  = oppBuy*pxPerDollar;
        subjSell = subjSell*pxPerDollar;  oppSell = oppSell*pxPerDollar;
        ylab = 'Price (px)';
    else
        ylab = 'Price ($)';
    end

    plot(ax, trialIdx, subjBuy,  lineStyleSubject,  'Color', colorPriceBuy,  'LineWidth', 1.8); hold(ax,'on');
    plot(ax, trialIdx, oppBuy,   lineStyleOpponent, 'Color', colorPriceBuy,  'LineWidth', 1.8);
    plot(ax, trialIdx, subjSell, lineStyleSubject,  'Color', colorPriceSell, 'LineWidth', 1.8);
    plot(ax, trialIdx, oppSell,  lineStyleOpponent, 'Color', colorPriceSell, 'LineWidth', 1.8);
    hold(ax,'off');

    title(ax, sprintf('Market %d', m));
    xlabel(ax, 'Trial'); ylabel(ax, ylab);
    xlim(ax, [0 15]); xticks(ax, xTicksTrials);

    legend(ax, {'BUY (Subject)','BUY (Opponent)','SELL (Subject)','SELL (Opponent)'}, 'Location','northeast');
    legend(ax,'boxoff');
end

%% =========================== FIG B ========================================
% DP OPTIMAL POLICY HEATMAPS (STATE-DEPENDENT)
% WHAT IT SHOWS
%   - Heatmap = DP-optimal action at each (Trial, Portfolio-before-choice).
%   - Each column is a trial; each row is portfolio size p before choosing.
%   - Green line overlay = deterministic DP-optimal trajectory starting at p0=0.
% WHY IT LOOKS "SPLIT"
%   - Because DP optimal action depends on current portfolio state.
%
% IMPLEMENTATION NOTES
%   - Unreachable states (p > t-1) masked as NaN (white background).
%   - Two separate figures: Subject and Opponent (3x2 layout each).
cmapActions = [0.15 0.55 0.95;   % BUY
              0.60 0.60 0.60;   % HOLD
              0.90 0.25 0.25];  % SELL
alphaHeatmap = 0.6;

for pp = 1:2
    figName = sprintf('FIG B: DP optimal policy (%s)', perspectiveNames{pp});
    figure('Color','w','Name',figName);
    tlB = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    trajStyle = lineStyleSubject;
    if pp==2, trajStyle = lineStyleOpponent; end

    for mi = 1:numel(market_templates)
        m = market_templates(mi);
        ax = nexttile(tlB, mi);

        A = policyDP_action{pp,mi}'; % (Pmax+1) x Tmax

        % Mask unreachable states p > t-1
        [Pgrid, Tgrid] = ndgrid(0:Pmax, 1:Tmax);
        A(Pgrid > (Tgrid-1)) = NaN;

        hImg = imagesc(ax, 1:Tmax, 0:Pmax, A);
        set(ax, 'YDir','normal', 'Color', [1 1 1]);
        set(hImg, 'AlphaData', alphaHeatmap * ~isnan(A));
        colormap(ax, cmapActions);
        caxis(ax, [1 3]);

        xlabel(ax, 'Trial');
        ylabel(ax, {'Portfolio size','(before choice)'});
        xlim(ax, [0 15]); xticks(ax, xTicksTrials);
        title(ax, sprintf('Market %d', m));

        % Overlay DP-optimal trajectory from p0=0 (prePortfolio each trial)
        pol = policyDP_action{pp,mi};
        p = 0;
        prePortfolioPath = zeros(Tmax,1);
        for t = 1:Tmax
            prePortfolioPath(t) = p;
            a = pol(t, p+1);
            if a==ACTION_BUY
                p = p+1;
            elseif a==ACTION_SELL
                p = max(p-1,0);
            end
        end

        hold(ax,'on');
        plot(ax, 1:Tmax, prePortfolioPath, trajStyle, 'Color', colorTrajectory, 'LineWidth', 2.2);
        hold(ax,'off');

        cb = colorbar(ax);
        cb.Ticks = [1 2 3];
        cb.TickLabels = {'BUY','HOLD','SELL'}; % colorbar is the action legend
    end

    sgtitle(sprintf('DP optimal policy – %s', perspectiveNames{pp}));
end

%% =========================== FIG C ========================================
% DP OPTIMAL INVENTORY TRAJECTORY (STATE-DEPENDENT, IMPLIED BEHAVIOR)
% WHAT IT SHOWS
%   - Deterministic portfolio trajectory if the agent follows the DP policy,
%     starting from p0=0.
%   - x-axis is "Trial" = 0..15 (0 = before any action).
%   - Both Subject and Opponent trajectories shown in GREEN; solid vs dashed.
figure('Color','w','Name','FIG C: DP optimal inventory trajectory from p0=0');
tlC = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for mi = 1:numel(market_templates)
    m = market_templates(mi);
    ax = nexttile(tlC, mi);

    % Subject trajectory
    polSubj = policyDP_action{1,mi};
    p = 0;
    trajSubj = zeros(Tmax+1,1);
    for t = 1:Tmax
        a = polSubj(t, p+1);
        if a==ACTION_BUY
            p = p+1;
        elseif a==ACTION_SELL
            p = max(p-1,0);
        end
        trajSubj(t+1) = p;
    end

    % Opponent trajectory
    polOpp = policyDP_action{2,mi};
    p = 0;
    trajOpp = zeros(Tmax+1,1);
    for t = 1:Tmax
        a = polOpp(t, p+1);
        if a==ACTION_BUY
            p = p+1;
        elseif a==ACTION_SELL
            p = max(p-1,0);
        end
        trajOpp(t+1) = p;
    end

    stairs(ax, 0:Tmax, trajSubj, lineStyleSubject,  'Color', colorTrajectory, 'LineWidth', 2.2); hold(ax,'on');
    stairs(ax, 0:Tmax, trajOpp,  lineStyleOpponent, 'Color', colorTrajectory, 'LineWidth', 2.2); hold(ax,'off');

    title(ax, sprintf('Market %d', m));
    xlabel(ax, 'Trial');
    ylabel(ax, {'Portfolio size','(after choice)'});
    xlim(ax, [0 15]); xticks(ax, xTicksTrials);
    ylim(ax, [0 15]);

    legend(ax, {'Subject','Opponent'}, 'Location','northeast');
    legend(ax,'boxoff');
end

%% =========================== FIG D ========================================
% FIXED-FUNDAMENTAL HEURISTIC: DETERMINISTIC POLICY DESCRIPTOR PLOTS
% WHAT IT SHOWS
%   - BUY price schedule (blue) and SELL price schedule (red) for a given perspective.
%   - Fundamental value curve (black) computed from remaining expected dividends.
%   - Markers at the top of each subplot encode the heuristic "best action per trial".
% IMPORTANT STYLE RULE (consistent across all figures)
%   - Subject = solid lines, Opponent = dashed lines.

markerSize = 6;
colorHoldMarker = [0.5 0.5 0.5];

for pp = 1:2
    figName = sprintf('FIG D: Fixed-fundamental heuristic descriptor (%s)', perspectiveNames{pp});
    figure('Color','w','Name',figName);
    tlD = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % Line style by perspective (Subject solid; Opponent dashed)
    if pp==1
        ls = lineStyleSubject;
    else
        ls = lineStyleOpponent;
    end

    for mi = 1:numel(market_templates)
        m = market_templates(mi);
        ax = nexttile(tlD, mi);

        buySched  = squeeze(priceBuy_dollars(mi,pp,:));
        sellSched = squeeze(priceSell_dollars(mi,pp,:));

        % Optionally plot in pixels (policy comparisons are invariant to scaling)
        if plot_in_pixels
            buyPlot  = buySched*pxPerDollar;
            sellPlot = sellSched*pxPerDollar;
            fundPlot = fundamental_dollars_byTrial*pxPerDollar;
            ylab = 'Value (px)';
        else
            buyPlot  = buySched;
            sellPlot = sellSched;
            fundPlot = fundamental_dollars_byTrial;
            ylab = 'Value ($)';
        end

        % Curves: style encodes perspective; color encodes quantity
        plot(ax, trialIdx, buyPlot,  ls, 'Color', colorPriceBuy,      'LineWidth', 1.8); hold(ax,'on');
        plot(ax, trialIdx, sellPlot, ls, 'Color', colorPriceSell,     'LineWidth', 1.8);
        plot(ax, trialIdx, fundPlot, ls, 'Color', colorFundamental,  'LineWidth', 1.8); % fundamental is shared, always solid
        hold(ax,'off');

        title(ax, sprintf('Market %d', m));
        xlabel(ax, 'Trial'); ylabel(ax, ylab);
        xlim(ax, [0 15]); xticks(ax, xTicksTrials);

        % Marker row encoding the heuristic action signal (trial-only action)
        yL = ylim(ax);
        yMarker = yL(2) - 0.06*(yL(2)-yL(1)); % near the top

        actions_t = squeeze(policyFundamental_action_trialOnly(pp,mi,:)); % [Tmax x 1]

        hold(ax,'on');
        for t = 1:Tmax
            a = actions_t(t);
            if a==ACTION_BUY
                plot(ax, t, yMarker, '^', 'MarkerSize', markerSize, ...
                    'MarkerFaceColor', colorPriceBuy, 'MarkerEdgeColor', colorPriceBuy);
            elseif a==ACTION_SELL
                plot(ax, t, yMarker, 'v', 'MarkerSize', markerSize, ...
                    'MarkerFaceColor', colorPriceSell, 'MarkerEdgeColor', colorPriceSell);
            elseif a==ACTION_HOLD
                plot(ax, t, yMarker, 'o', 'MarkerSize', markerSize, ...
                    'MarkerFaceColor', colorHoldMarker, 'MarkerEdgeColor', colorHoldMarker);
            else
                plot(ax, t, yMarker, 's', 'MarkerSize', markerSize, ...
                    'MarkerFaceColor', colorFundamental, 'MarkerEdgeColor', colorFundamental);
            end
        end
        hold(ax,'off');

        % Legend describes curves only (markers described in caption)
        legend(ax, {'priceBuy','priceSell','fundamental'}, 'Location','northeast');
        legend(ax,'boxoff');
    end

    sgtitle(sprintf('Fixed-fundamental heuristic descriptor – %s', perspectiveNames{pp}));
end

%% =========================== FIG E ========================================
% FIXED-FUNDAMENTAL HEURISTIC: IMPLIED INVENTORY TRAJECTORY FROM p0=0
% WHAT IT SHOWS
%   - Deterministic portfolio trajectory if the agent follows the FIXED-FUNDAMENTAL
%     heuristic action signal (trial-only), starting from p0=0.
%   - x-axis is "Trial" = 0..15 (0 = before any action).
%   - Subject and Opponent shown in GREEN; solid vs dashed (consistent with Fig C).
%
% IMPORTANT
%   - This is NOT DP-optimal. It is the implied trajectory under the heuristic.
%   - Because the heuristic is trial-only, feasibility constraints are enforced here:
%       * If the heuristic says SELL at p=0 -> we force HOLD (cannot short).
%
figure('Color','w','Name','FIG E: Fixed-fundamental heuristic inventory trajectory from p0=0');
tlE = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for mi = 1:numel(market_templates)
    m = market_templates(mi);
    ax = nexttile(tlE, mi);

    % Subject trajectory under heuristic
    trajSubj = zeros(Tmax+1,1);
    p = 0;
    for t = 1:Tmax
        a = policyFundamental_action_trialOnly(1,mi,t); % trial-only signal
        if a==ACTION_BUY
            p = p + 1;
        elseif a==ACTION_SELL
            if p>0
                p = p - 1;
            else
                % infeasible SELL -> HOLD
                p = p;
            end
        else
            % HOLD or AMBIG -> HOLD (conservative)
            p = p;
        end
        trajSubj(t+1) = p;
    end

    % Opponent trajectory under heuristic
    trajOpp = zeros(Tmax+1,1);
    p = 0;
    for t = 1:Tmax
        a = policyFundamental_action_trialOnly(2,mi,t);
        if a==ACTION_BUY
            p = p + 1;
        elseif a==ACTION_SELL
            if p>0
                p = p - 1;
            else
                p = p;
            end
        else
            p = p;
        end
        trajOpp(t+1) = p;
    end

    stairs(ax, 0:Tmax, trajSubj, lineStyleSubject,  'Color', colorTrajectory, 'LineWidth', 2.2); hold(ax,'on');
    stairs(ax, 0:Tmax, trajOpp,  lineStyleOpponent, 'Color', colorTrajectory, 'LineWidth', 2.2); hold(ax,'off');

    title(ax, sprintf('Market %d', m));
    xlabel(ax, 'Trial');
    ylabel(ax, {'Portfolio size','(after choice)'});
    xlim(ax, [0 15]); xticks(ax, xTicksTrials);
    ylim(ax, [0 15]);

    legend(ax, {'Subject','Opponent'}, 'Location','northeast');
    legend(ax,'boxoff');
end

%% =========================== Optional: Save figures ------------------------
% saveas(findobj('Name','FIG A: Price schedule invariance'), 'FIG_A_prices.png');
% saveas(findobj('Name','FIG B: DP optimal policy (Subject)'), 'FIG_B_DP_Subject.png');
% saveas(findobj('Name','FIG B: DP optimal policy (Opponent)'), 'FIG_B_DP_Opponent.png');
% saveas(findobj('Name','FIG C: DP optimal inventory trajectory from p0=0'), 'FIG_C_DP_traj.png');
% saveas(findobj('Name','FIG D: Fixed-fundamental heuristic descriptor (Subject)'), 'FIG_D_fund_Subject.png');
% saveas(findobj('Name','FIG D: Fixed-fundamental heuristic descriptor (Opponent)'), 'FIG_D_fund_Opponent.png');
