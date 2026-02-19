% File: d20251009_plot_buy_vs_bubble.m
% This code is written to compare the number of BUY's monkey chooses in
% bubble vs non-bubble markets in each condition. Tim 251103

% ---------- Buy-count distributions by Bubble vs Non-Bubble ----------
clear; clc;

% ------------------- Toggles -------------------
LIVE_NONLIVE = 0; % 1: collapse to {Non-Live, Live}
SUBJECT_ONLY = 1; % 1: subject monkey only (odd trialBlock)
PER_MONKEY = 1; % 1: also plot per-monkey distributions

% ------------------- Load -------------------
try
T = readtable('C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx','Sheet','Sheet1');
% % if we want sessions included in the behavior-only analysis
% T = readtable('C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx','Sheet','Sheet1');
catch
if ~exist('T','var'), error('Could not read table. Provide T or update the path.'); end
end

% ------------------- Required columns -------------------
need = {'condition','option','marketOrig','year','month','day','SessionOfDay','bubbleMarket'};
miss = need(~ismember(need, T.Properties.VariableNames));
if ~isempty(miss), error('Missing columns: %s', strjoin(miss, ', ')); end

% ------------------- Prep -------------------
T.supCond = floor(T.condition);
if SUBJECT_ONLY && ismember('trialBlock', T.Properties.VariableNames)
T = T(mod(T.trialBlock,2) == 1, :);
end

% Bubble flag
raw = T.bubbleMarket;
if isnumeric(raw) || islogical(raw)
isBubble = logical(raw);
else
s = string(raw);
isBubble = strcmpi(s,'bubble') | strcmpi(s,'bubbles') | strcmpi(s,'true') | strcmpi(s,'yes') | strcmpi(s,'1');
end

% Conditions (respect LIVE_NONLIVE)
if LIVE_NONLIVE
sup = T.supCond; sup(sup ~= 4) = 1; sup(sup == 4) = 2;
labels = ["Non-Live","Live"]; K_all = 2;
else
sup = T.supCond;
labels = ["AI","Replay","Decoy","Live"]; K_all = 4;
end

% Buy flag
isBuy = (T.option == 1);

% Densify condition ids to handle missing conditions
present = intersect(1:K_all, unique(sup(:))','stable');
if isempty(present), error('No valid conditions present.'); end
maxSup = max([K_all; sup(:)]);
map = zeros(1, maxSup); map(present) = 1:numel(present);
sup_d = arrayfun(@(x) map(x), sup);
labels_present = labels(present);
K = numel(labels_present);

% ------------------- Build observations (marketOrig × Year × Month × Day × SessionOfDay) -------------------
ObsTbl = table( ...
sup_d(:), isBubble(:), string(T.marketOrig(:)), T.year(:), T.month(:), T.day(:), T.SessionOfDay(:), isBuy(:), ...
'VariableNames', {'sup_d','isBubble','marketOrig','Year','Month','Day','SessionOfDay','isBuy'});

Gobs = groupsummary(ObsTbl, {'sup_d','isBubble','marketOrig','Year','Month','Day','SessionOfDay'}, 'sum','isBuy');
Gobs.sum_isBuy = double(Gobs.sum_isBuy);
Gobs.isBubble = logical(Gobs.isBubble);

% Global bin edges (integer counts)
maxN = max([0; Gobs.sum_isBuy]);
edges = (-0.5):(1):(maxN+0.5);

% ------------------- Overall distribution plots -------------------
ncols = min(4, K); nrows = ceil(K / ncols);
figure('Name','Overall: Buy-count distributions per condition (obs = marketOrig × Y-M-D × SessionOfDay)');
tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');
for k = 1:K
nexttile; hold on
yNB = Gobs.sum_isBuy(Gobs.sup_d==k & ~Gobs.isBubble);
yB = Gobs.sum_isBuy(Gobs.sup_d==k & Gobs.isBubble);
if ~isempty(yNB)
histogram(yNB, 'BinEdges', edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',1.2);
end
if ~isempty(yB)
histogram(yB, 'BinEdges', edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',1.2);
end
title(sprintf('%s (n_{NB}=%d, n_{B}=%d)', labels_present(k), numel(yNB), numel(yB)));
xlabel('Buy count per [marketOrig × Y-M-D × SessionOfDay]'); ylabel('Probability');
if ~isempty(yNB) || ~isempty(yB)
lg = legend({'Non-Bubble','Bubble'}, 'Location','northeast'); set(lg,'Box','off');
else
text(0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
end
box on
end

% ------------------- Overall: statistical tests (Bubble vs Non-Bubble, per condition) -------------------
RS_p = nan(K,1); KS_p = nan(K,1); AUC = nan(K,1); Delta = nan(K,1);
nNB = zeros(K,1); nB = zeros(K,1); medNB = nan(K,1); medB = nan(K,1); meanNB = nan(K,1); meanB = nan(K,1);
for k = 1:K
yNB = Gobs.sum_isBuy(Gobs.sup_d==k & ~Gobs.isBubble);
yB = Gobs.sum_isBuy(Gobs.sup_d==k & Gobs.isBubble);
nNB(k) = numel(yNB); nB(k) = numel(yB);
if nNB(k) >= 1 && nB(k) >= 1
medNB(k) = median(yNB); medB(k) = median(yB);
meanNB(k) = mean(yNB); meanB(k) = mean(yB);
try
[RS_p(k), ~, st] = ranksum(yNB, yB, 'method','approx');
U = st.ranksum - nNB(k)*(nNB(k)+1)/2;
AUC(k) = U / (nNB(k)*nB(k));
Delta(k) = 2*AUC(k) - 1; % Cliff's delta (NB vs B)
catch
RS_p(k) = NaN; AUC(k) = NaN; Delta(k) = NaN;
end
try
[~,KS_p(k)] = kstest2(yNB, yB);
catch
KS_p(k) = NaN;
end
end
end
Q_fdr = fdr_bh(RS_p);

ResOverall = table(labels_present(:), nNB, nB, medNB, medB, meanNB, meanB, AUC, Delta, RS_p, KS_p, Q_fdr, ...
'VariableNames', {'Condition','n_NonBubble','n_Bubble','Median_NonBubble','Median_Bubble','Mean_NonBubble','Mean_Bubble','AUC_NBgtB','CliffsDelta_NBvsB','RankSum_p','KS_p','FDR_q'});
disp('=== Overall tests (obs = marketOrig × Y-M-D × SessionOfDay) ===');
disp(ResOverall);

% ------------------- Optional per-monkey plots + tests -------------------
if PER_MONKEY && ismember('monkey', T.Properties.VariableNames)
um = [1 2];
for mky = um
msel = T.monkey == mky;
if ~any(msel), continue; end

    ObsM = table( ...
        sup_d(msel), isBubble(msel), string(T.marketOrig(msel)), T.year(msel), T.month(msel), T.day(msel), T.SessionOfDay(msel), isBuy(msel), ...
        'VariableNames', {'sup_d','isBubble','marketOrig','Year','Month','Day','SessionOfDay','isBuy'});

    Gm = groupsummary(ObsM, {'sup_d','isBubble','marketOrig','Year','Month','Day','SessionOfDay'}, 'sum','isBuy');
    Gm.sum_isBuy = double(Gm.sum_isBuy); Gm.isBubble = logical(Gm.isBubble);

    maxNm = max([0; Gm.sum_isBuy]); edgesM = (-0.5):(1):(maxNm+0.5);

    figure('Name', sprintf('M%d: Buy-count distributions per condition (obs = marketOrig × Y-M-D × SessionOfDay)', mky));
    tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');
    for k = 1:K
        nexttile; hold on
        yNB = Gm.sum_isBuy(Gm.sup_d==k & ~Gm.isBubble);
        yB  = Gm.sum_isBuy(Gm.sup_d==k &  Gm.isBubble);
        if ~isempty(yNB)
            histogram(yNB, 'BinEdges', edgesM, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',1.2);
        end
        if ~isempty(yB)
            histogram(yB,  'BinEdges', edgesM, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',1.2);
        end
        title(sprintf('M%d — %s (n_{NB}=%d, n_{B}=%d)', mky, labels_present(k), numel(yNB), numel(yB)));
        xlabel('Buy count per [marketOrig × Y-M-D × SessionOfDay]'); ylabel('Probability');
        if ~isempty(yNB) || ~isempty(yB)
            lg = legend({'Non-Bubble','Bubble'}, 'Location','northeast'); set(lg,'Box','off');
        else
            text(0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
        end
        box on
    end

    RS_pM = nan(K,1); KS_pM = nan(K,1); AUCM = nan(K,1); DeltaM = nan(K,1);
    nNBM = zeros(K,1); nBM = zeros(K,1); medNBM = nan(K,1); medBM = nan(K,1); meanNBM = nan(K,1); meanBM = nan(K,1);
    for k = 1:K
        yNB = Gm.sum_isBuy(Gm.sup_d==k & ~Gm.isBubble);
        yB  = Gm.sum_isBuy(Gm.sup_d==k &  Gm.isBubble);
        nNBM(k) = numel(yNB); nBM(k) = numel(yB);
        if nNBM(k) >= 1 && nBM(k) >= 1
            medNBM(k) = median(yNB); medBM(k) = median(yB);
            meanNBM(k) = mean(yNB);  meanBM(k) = mean(yB);
            try
                [RS_pM(k), ~, st] = ranksum(yNB, yB, 'method','approx');
                U = st.ranksum - nNBM(k)*(nNBM(k)+1)/2;
                AUCM(k) = U / (nNBM(k)*nBM(k));
                DeltaM(k) = 2*AUCM(k) - 1;
            catch
                RS_pM(k) = NaN; AUCM(k) = NaN; DeltaM(k) = NaN;
            end
            try
                [~,KS_pM(k)] = kstest2(yNB, yB);
            catch
                KS_pM(k) = NaN;
            end
        end
    end
    Q_fdrM = fdr_bh(RS_pM);

    ResM = table(labels_present(:), nNBM, nBM, medNBM, medBM, meanNBM, meanBM, AUCM, DeltaM, RS_pM, KS_pM, Q_fdrM, ...
        'VariableNames', {'Condition','n_NonBubble','n_Bubble','Median_NonBubble','Median_Bubble','Mean_NonBubble','Mean_Bubble','AUC_NBgtB','CliffsDelta_NBvsB','RankSum_p','KS_p','FDR_q'});
    fprintf('=== M%d tests (obs = marketOrig × Y-M-D × SessionOfDay) ===\n', mky);
    disp(ResM);
end


end

% ------------------- Helpers -------------------
function q = fdr_bh(p)
q = nan(size(p));
x = p; x(~isfinite(x)) = 1;
[sp, idx] = sort(x(:), 'ascend');
m = sum(isfinite(p));
if m == 0, return; end
q_sorted = sp .* (m ./ (1:numel(sp))');
q_sorted = cummin(flipud(q_sorted));
q_sorted = flipud(q_sorted);
q(idx) = q_sorted;
end
