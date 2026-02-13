clear all
% plotStuff

% Plot psth from Nex data
% buy data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1buy.mat')
aibuyUnits = spikes;
aibuyNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1buy.mat')
rebuyUnits = spikes(:,1:92);
rebuyNames = units(:,1:92);
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1buy.mat')
debuyUnits = spikes;
debuyNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1buy.mat')
livebuyUnits = spikes;
livebuyNames = units;

% hold data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1hold.mat')
holdUnits = spikes;
holdNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1hold.mat')
holdUnits = [holdUnits, spikes(:,1:92)];
holdNames = [holdNames, units(:,1:92)];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1hold.mat')
holdUnits = [holdUnits, spikes];
holdNames = [holdNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1hold.mat')
holdUnits = [holdUnits, spikes];
holdNames = [holdNames, units];

% sell data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1sell.mat')
sellUnits = spikes;
sellNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1sell.mat')
sellUnits = [sellUnits, spikes(:,1:92)];
sellNames = [sellNames, units(:,1:92)];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1sell.mat')
sellUnits = [sellUnits, spikes];
sellNames = [sellNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1sell.mat')
sellUnits = [sellUnits, spikes];
sellNames = [sellNames, units];

clear units2 units3 spikes2

aibuyMean = mean(aibuyUnits,2);

SEM_aibuy = std(aibuyUnits, [], 2)./ sqrt(size(aibuyUnits,2));    % Calculate Standard Error Of The Mean
CI95_aibuy = bsxfun(@plus, mean(aibuyUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_aibuy));   % 95% Confidence Intervals

rebuyMean = mean(rebuyUnits,2);

SEM_rebuy = std(rebuyUnits, [], 2)./ sqrt(size(rebuyUnits,2));    % Calculate Standard Error Of The Mean
CI95_rebuy = bsxfun(@plus, mean(rebuyUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_rebuy));   % 95% Confidence Intervals

debuyMean = mean(debuyUnits,2);

SEM_debuy = std(debuyUnits, [], 2)./ sqrt(size(debuyUnits,2));    % Calculate Standard Error Of The Mean
CI95_debuy = bsxfun(@plus, mean(debuyUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_debuy));   % 95% Confidence Intervals

livebuyMean = mean(livebuyUnits,2);

SEM_livebuy = std(livebuyUnits, [], 2)./ sqrt(size(livebuyUnits,2));    % Calculate Standard Error Of The Mean
CI95_livebuy = bsxfun(@plus, mean(livebuyUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_livebuy));   % 95% Confidence Intervals

% holdMean = mean(holdUnits,2);
% 
% SEM_hold = std(holdUnits, [], 2)./ sqrt(size(holdUnits,2));    % Calculate Standard Error Of The Mean
% CI95_hold = bsxfun(@plus, mean(holdUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_hold)); 
% 
% sellMean = mean(sellUnits,2);
% 
% SEM_sell = std(sellUnits, [], 2)./ sqrt(size(sellUnits,2));    % Calculate Standard Error Of The Mean
% CI95_sell = bsxfun(@plus, mean(sellUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_sell)); 


e1 = errorbar(-0.5:0.005:1.495,aibuyMean,SEM_aibuy,'CapSize',0);
hold on
e2 = errorbar(-0.5:0.005:1.495,rebuyMean,SEM_rebuy,'CapSize',0);
e3 = errorbar(-0.5:0.005:1.495,debuyMean,SEM_debuy,'CapSize',0);
e4 = errorbar(-0.5:0.005:1.495,livebuyMean,SEM_livebuy,'CapSize',0);
% e2 = errorbar(-0.5:0.005:1.495,holdMean, SEM_hold,'CapSize',0);
% e3 = errorbar(-0.5:0.005:1.495,sellMean, SEM_sell,'CapSize',0);
hold off
e1.Bar.LineStyle = 'dotted';
e2.Bar.LineStyle = 'dotted';
e3.Bar.LineStyle = 'dotted';
e4.Bar.LineStyle = 'dotted';
e1.CapSize = 0.1;

% add background box from internet example:
% figure
% %Plot something
% plot(1:10)
% Add lines
% h1 = line([-0.5 -0.5],[-1.5 2.5]);
% h2 = line([-0.1 -0.1],[-1.5 2.5]);
% % Set properties of lines
% set([h1 h2],'Color','k','LineWidth',2)
% % Add a patch
% gray = [0.4 0.4 0.4];
% patch([-0.5 -0.01 -0.01 -0.5],[-1.5 -1.5 2.5 2.5],gray)
% % The order of the "children" of the plot determines which one appears on top.
% % I need to flip it here.
% set(gca,'children',flipud(get(gca,'children')))

% numUnits = numel(nexColumnNames);
% 
% figure
% hold on
% 
% for i = 2:numUnits
%     sig = nex(:,i)';
%     plot(-0.5:0.005:2.495,sig)
% end
% 
% plot(-2.5:0.005:3.995,mNex)
% 
% hold off
%     

