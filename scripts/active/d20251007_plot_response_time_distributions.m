% File: d20251007_plot_response_time_distributions.m
% croak
clear; close all; clc;

% ------------------- Toggles -------------------
LIVE_NONLIVE         = 0;   % 1: collapse to Non-Live vs Live
FILTER_TRESP         = 0;   % 1: exclude invalid tResp by thresholds
TRESP_MIN            = 0;
TRESP_MAX            = inf;
ADJUST_TRESP_BASELINE= 0;
AGGREGATE_BY_SESSION = 0;   % 1: aggregate per session-of-day using unique T.time

% ------------------- Load & prep -------------------
T = readtable('C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx','Sheet','Sheet1');
% % if we want sessions included in the behavior-only analysis 
% T = readtable('C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL\dataTabFULL_filtered_from_recoveryFile.xlsx','Sheet','Sheet1');
T.supCond = floor(T.condition);
T = T(mod(T.trialBlock,2) == 1, :);   % subject monkey only

if FILTER_TRESP
    T = T(T.tResp < TRESP_MAX & T.tResp > TRESP_MIN, :);
end

if ADJUST_TRESP_BASELINE
    T.tResp = T.tResp + 0.1;  % Adjust tResp by adding 0.1 s
end

% ------------------- Per-monkey loop -------------------
for mky = 1:2
    Tmonkey = T(T.monkey == mky, :);

    % Raw for histogram (before any remap/aggregation)
    y_raw = Tmonkey.tResp;
    if isduration(y_raw), y_raw = seconds(y_raw); end
    y_raw = double(y_raw(:));

    % Working vectors
    y = Tmonkey.tResp;
    if isduration(y), y = seconds(y); end
    y = double(y(:));
    g = double(Tmonkey.supCond(:));
    tSess = Tmonkey.time .* Tmonkey.session;  % session-of-day unique ID

    if LIVE_NONLIVE
        g(g ~= 4) = 1; g(g == 4) = 2;
        labels = ["Non-Live","Live"];
    else
        labels = ["AI","Replay","Decoy","Live"];
    end
    K_all = numel(labels);

    % valid & recode
    valid = isfinite(y) & y > 0 & ismember(g,1:K_all);
    y = y(valid); g = g(valid); tSess = tSess(valid);

    % -------- Optional aggregation by session-of-day --------
    if AGGREGATE_BY_SESSION
        [G, keySess, keyCond] = findgroups(tSess, g); %#ok<ASGLU>
        y = splitapply(@median, y, G);       % one value per (session, condition)
        g = keyCond;                          % keep condition codes with gaps possible
    end

    % -------- Densify group ids (handles missing conditions) --------
    present = intersect(1:K_all, unique(g(:))','stable');
    map = zeros(1, K_all); map(present) = 1:numel(present);
    g_dense = arrayfun(@(x) map(x), g);      % remap to 1..K_present
    labels_present = labels(present);
    K = numel(labels_present);

    % ------------------- Plots -------------------
    % Histogram (raw distribution; independent of aggregation)
    figure
    histogram(y_raw, 'BinEdges', logspace(-4.5, 2, 50), 'Normalization','probability')
    xscale('log'); xticks(logspace(-4, 2, 7)); xlim([0.5e-4, 2e2])
    xlabel("tResp (s)"); ylabel("Frequency")
    title("M" + mky + " combined conditions tResp distribution")

    % Swarm of the analysis sample (aggregated if AGGREGATE_BY_SESSION=1)
    figure; hold on
    if ~isempty(y)
        swarmchart(g_dense, y, 6, 'filled', 'MarkerFaceAlpha',0.5, 'XJitter','rand');
        set(gca,'XTick',1:K,'XTickLabel',labels_present);
        xlabel('Condition'); ylabel('Response time (s)');
        % ylim([max(min(y)*0.5, realmin), max([mean(y)*2, prctile(y,99)])]); 
        ylim([max(min(y)*0.5, realmin), 1]); 
        yscale('log');
    else
        text(0.5,0.5,'No valid data', 'Units','normalized','HorizontalAlignment','center');
        set(gca,'XTick',1:K,'XTickLabel',labels_present); yscale('log'); 
    end

    % ------------------- Stats: KW + Dunn + Effect sizes -------------------
    if K >= 2 && numel(y) >= 2
        [pKW, tblKW, statsKW] = kruskalwallis(y, g_dense, 'off');
        eps2 = kw_epsilon2(tblKW, numel(y), K);  % overall effect size

        % Dunn–Šidák pairwise (if >=2 groups)
        cmp = [];
        try
            cmp = multcompare(statsKW, 'CType','dunn-sidak', 'Display','off');  % [i j lo est hi p]
        catch
            % leave cmp empty if not computable
        end

        % Pairwise nonparametric effect sizes (Cliff's delta, AUC, median diff)
        eff = pairwise_effects(y, g_dense, labels_present);

        % Title with KW + epsilon^2
        title(sprintf('M%d Response time by condition (KW p = %.3g, \\epsilon^2 = %.3f)', mky, pKW, eps2));

        % ---- Annotate significant pairs with brackets and stars ----
        if ~isempty(cmp)
            alpha = 0.05;
            sigIdx = find(cmp(:,6) < alpha);
            pairs  = cmp(sigIdx,1:2);
            pvals  = cmp(sigIdx,6);

            % stack overlapping brackets into levels
            levels = zeros(size(sigIdx));
            level_end = [];
            for m = 1:numel(sigIdx)
                i = pairs(m,1); j = pairs(m,2);
                placed = false;
                for L = 1:numel(level_end)
                    if i > level_end(L)
                        levels(m) = L; level_end(L) = j; placed = true; break
                    end
                end
                if ~placed
                    levels(m) = numel(level_end)+1; level_end(end+1) = j; %#ok<SAGROW>
                end
            end

            fac = 2; 
            baseTop = 10; 
            yMaxNeeded = baseTop * fac^(max([levels;0])+1);
            yl = ylim; if yMaxNeeded > yl(2), ylim([yl(1) yMaxNeeded]); yl = ylim; end

            for m = 1:numel(sigIdx)
                i = pairs(m,1); j = pairs(m,2);
                L = levels(m);
                ypairTop = max([y(g_dense==i); y(g_dense==j)]);
                yPos = max(baseTop, ypairTop) * fac^(L);
                yTick = yPos / 1.05;  % Adjusted for better spacing
                plot([i i j j],[yTick yPos yPos yTick], 'k', 'LineWidth', 1.2);
                text(mean([i,j]), yPos*1.02, p2stars(pvals(m)), ...  % Adjusted for better positioning
                    'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
            end

            % ------- Print results -------
            res = table(labels_present(cmp(:,1))', labels_present(cmp(:,2))', cmp(:,4), cmp(:,6), ...
                'VariableNames', {'Group1','Group2','MeanRankDiff','PValue'});
            fprintf('\nM%d — Kruskal–Wallis: p=%.3g, epsilon^2=%.3f\n', mky, pKW, eps2);
            disp(res);
        else
            fprintf('\nM%d — Kruskal–Wallis: p=%.3g, epsilon^2=%.3f (no pairwise cmp)\n', mky, pKW, eps2);
        end

        % Effect sizes table
        disp(eff);
    else
        fprintf('\nM%d — Not enough data/groups for KW.\n', mky);
    end
end

% ------------------- Helpers -------------------
function s = p2stars(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end

function eps2 = kw_epsilon2(tblKW, N, K)
    % Epsilon-squared effect size for Kruskal–Wallis:
    % eps^2 = (H - K + 1) / (N - K)
    H = NaN;
    try
        H = tblKW{2,5};  % typical location for Chi-square in MATLAB's KW table
        if ~isnumeric(H) || ~isfinite(H), H = NaN; end
    catch
        % fallback: scan numeric
        for r = 1:size(tblKW,1)
            for c = 1:size(tblKW,2)
                if isnumeric(tblKW{r,c}) && isfinite(tblKW{r,c}), H = tblKW{r,c}; break; end
            end
            if ~isnan(H), break; end
        end
    end
    if isnan(H), eps2 = NaN; else, eps2 = (H - K + 1) / max(1, (N - K)); end
end

function effTab = pairwise_effects(y, g, labels_present)
    % g must be dense 1..K codes; labels_present is Kx1 string
    K = numel(labels_present);
    G1 = strings(0,1); G2 = strings(0,1);
    n1v = zeros(0,1); n2v = zeros(0,1);
    aucv = zeros(0,1); deltav = zeros(0,1);
    medDiffv = zeros(0,1); pRSv = zeros(0,1);

    for i = 1:K
        for j = i+1:K
            xi = y(g==i); xj = y(g==j);
            if isempty(xi) || isempty(xj), continue; end
            [pRS,~,st] = ranksum(xi, xj, 'method','approx');
            ni = numel(xi); nj = numel(xj);
            U  = st.ranksum - ni*(ni+1)/2;
            AUC = U / (ni*nj);
            delta = 2*AUC - 1;           % Cliff's delta
            medDiff = median(xi) - median(xj);

            G1(end+1,1) = labels_present(i);
            G2(end+1,1) = labels_present(j);
            n1v(end+1,1) = ni; n2v(end+1,1) = nj;
            aucv(end+1,1) = AUC; deltav(end+1,1) = delta;
            medDiffv(end+1,1) = medDiff; pRSv(end+1,1) = pRS;
        end
    end

    % Build table (handles 0-row case cleanly)
    effTab = table(G1, G2, n1v, n2v, aucv, deltav, medDiffv, pRSv, ...
        'VariableNames', {'Group1','Group2','N1','N2','AUC','CliffsDelta','MedianDiff','PValue_RankSum'});
end
