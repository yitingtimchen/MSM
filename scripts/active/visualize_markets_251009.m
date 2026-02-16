% This code was made to visualize the markets played by the monkeys and how
% bubbly they are. Tim, 251104

T = readtable('C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx','Sheet','Sheet1');  % only Sheet1

T.fundamental = 2.35*(15-T.trialNum+1);
T.bubbleness = (T.priceBuy - T.fundamental);



% plot_event_relationships(T(T.month == 8, :)); 
plot_markets(T);
plot_bubbleness(T);



function plot_event_relationships(T) 
% Extract relevant data from T for plotting
buyReward = T.sizeB; 
buyPrice = T.priceBuy;   
fb1 = T.fb1;
fb2 = T.fb2;
reward1 = zeros(size(T, 1), 1);
possibleReward1 = [T.sizeB, T.sizeH, T.sizeS];
for ii = 1:size(possibleReward1, 1)
    reward1(ii) = possibleReward1(ii, T.option(ii));
end
reward2 = sqrt(T.divPerShare) .* T.postPortfolio;

% jitterx = rand(size(T,1), 1);
% jittery = rand(size(T,1), 1);
% reward2 = reward2 + jitterx;
% fb2 = fb2 + jittery;

monk = categorical(T.postPortfolio);
cats = categories(monk);

figure; hold on
for k = 1:numel(cats)
    idx = monk == cats{k};
    scatter(reward2(idx), fb2(idx), 150, 'Marker','.','DisplayName', char(cats{k}));
end
xlabel('reward 2'); ylabel('fb2');
legend('Location','best'); hold off

% figure; hold on
% for k = 1:numel(cats)
%     idx = monk == cats{k};
%     scatter(reward1(idx), fb1(idx), 50, 'Marker','.','DisplayName', char(cats{k}));
% end
% xlabel('reward 1'); ylabel('fb1');
% legend('Location','best'); hold off

end

function plot_bubbleness(T)
figure
tiledlayout(3, 2)
% bubbleness plots
for market = 1:6
    % Extract bubbleness data for the current market
    marketData = T(T.marketOrig == market, :);
    bubbleValues = marketData.bubbleness;

    % Plot the bubbleness distribution
    nexttile
    hold on
    histogram(bubbleValues, 'BinEdges', -30:5:30, 'Normalization', 'pdf');
    xline(mean(bubbleValues), '--r', 'LineWidth',3)
    xline(0, '--k', 'LineWidth',3)
    hold off
    xlabel('Bubbleness');
    ylabel('pdf');
    title(sprintf('Market %d Bubbleness Distribution\nmean bubbleness = %.2f, std = %.2f', market, mean(bubbleValues), std(bubbleValues)));
end
end

% croak
function plot_markets(T)
% Overlay fundamental vs trialBlock and color priceBuy by relation to FV.

if ~ismember('fundamental', T.Properties.VariableNames)
    if ismember('trialNum', T.Properties.VariableNames)
        T.fundamental = 2.35*(15 - T.trialNum + 1);
    else
        error('T must contain fundamental or trialNum.');
    end
end

markets = sort(unique(T.marketOrig,'stable'));
if numel(markets) < 6, error('Expected 6 markets, found %d.', numel(markets)); end
markets = markets(1:6);

figure('Name','priceBuy & priceSell vs trialBlock by Market');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for i = 1:6
    nexttile; hold on
    idx   = ismember(T.marketOrig, markets(i));
    x     = T.trialBlock(idx);
    yBuy  = T.priceBuy(idx);
    ySell = T.priceSell(idx);
    F     = T.fundamental(idx);

    valid = isfinite(x) & isfinite(yBuy) & isfinite(F);
    xv = x(valid); yv = yBuy(valid); Fv = F(valid);

    % Color-code priceBuy vs fundamental
    lo = yv <= Fv; hi = ~lo;
    hBuyLo = scatter(xv(lo), yv(lo), 16, 'filled', 'Marker','o');   % green when <= FV
    hBuyHi = scatter(xv(hi), yv(hi), 16, 'filled', 'Marker','o');   % red when >  FV
    set(hBuyLo, 'MarkerFaceColor',[0 0.6 0], 'MarkerFaceAlpha',0.6, 'MarkerEdgeAlpha',0.6);
    set(hBuyHi, 'MarkerFaceColor',[0.85 0 0], 'MarkerFaceAlpha',0.6, 'MarkerEdgeAlpha',0.6);

    % priceSell (keep as separate series)
    hSell = scatter(x, ySell, 16, 'filled', 'Marker','^');
    set(hSell,'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6,'MarkerFaceColor',[0 0 0]);

    % Fundamental line (dashed) as mean FV per trialBlock
    [xu,~,ic] = unique(xv);
    Fmean = accumarray(ic, Fv, [], @(v) mean(v,'omitnan'));
    [xu, order] = sort(xu);
    Fmean = Fmean(order);
    hF = plot(xu, Fmean, 'b--', 'LineWidth', 0.5);

    xlabel('trialBlock'); ylabel('Price');
    title(sprintf('Market %s', string(markets(i))));
    legend([hBuyLo hBuyHi hSell hF], {'priceBuy\leqFV','priceBuy > FV','priceSell','fundamental'}, 'Location','eastoutside');
    hold off
end
end
