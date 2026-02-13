clear all
% plotStuff

% Plot psth from Nex data
% larger Portfolio data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1largerPstart.mat')
largerPstartUnits = spikes;
largerPstartNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1largerPstart.mat')
largerPstartUnits = [largerPstartUnits, spikes(:,1:92)];
largerPstartNames = [largerPstartNames, units(:,1:92)];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1largerPstart.mat')
largerPstartUnits = [largerPstartUnits, spikes];
largerPstartNames = [largerPstartNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1largerPstart.mat')
largerPstartUnits = [largerPstartUnits, spikes];
largerPstartNames = [largerPstartNames, units];

% evenPstart data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1evenPstart.mat')
evenPstartUnits = spikes;
evenPstartNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1evenPstart.mat')
evenPstartUnits = [evenPstartUnits, spikes(:,1:92)];
evenPstartNames = [evenPstartNames, units(:,1:92)];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1evenPstart.mat')
evenPstartUnits = [evenPstartUnits, spikes];
evenPstartNames = [evenPstartNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1evenPstart.mat')
evenPstartUnits = [evenPstartUnits, spikes];
evenPstartNames = [evenPstartNames, units];

% smallerPstart data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1smallerPstart.mat')
smallerPstartUnits = spikes;
smallerPstartNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1smallerPstart.mat')
smallerPstartUnits = [smallerPstartUnits, spikes(:,1:92)];
smallerPstartNames = [smallerPstartNames, units(:,1:92)];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1smallerPstart.mat')
smallerPstartUnits = [smallerPstartUnits, spikes];
smallerPstartNames = [smallerPstartNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1smallerPstart.mat')
smallerPstartUnits = [smallerPstartUnits, spikes];
smallerPstartNames = [smallerPstartNames, units];

clear units2 units3 spikes2

largerPstartMean = mean(largerPstartUnits,2);

SEM_largerPstart = std(largerPstartUnits, [], 2)./ sqrt(size(largerPstartUnits,2));    % Calculate Standard Error Of The Mean
CI95_largerPstart = bsxfun(@plus, mean(largerPstartUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_largerPstart));   % 95% Confidence Intervals

evenPstartMean = mean(evenPstartUnits,2);

SEM_evenPstart = std(evenPstartUnits, [], 2)./ sqrt(size(evenPstartUnits,2));    % Calculate Standard Error Of The Mean
CI95_evenPstart = bsxfun(@plus, mean(evenPstartUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_evenPstart)); 

smallerPstartMean = mean(smallerPstartUnits,2);

SEM_smallerPstart = std(smallerPstartUnits, [], 2)./ sqrt(size(smallerPstartUnits,2));    % Calculate Standard Error Of The Mean
CI95_smallerPstart = bsxfun(@plus, mean(smallerPstartUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_smallerPstart)); 


e1 = errorbar(-0.5:0.005:1.495,largerPstartMean,SEM_largerPstart,'CapSize',0,'DisplayName','Larger Portfolio');
hold on
e2 = errorbar(-0.5:0.005:1.495,evenPstartMean, SEM_evenPstart,'CapSize',0,'DisplayName','Even Portfolio');
e3 = errorbar(-0.5:0.005:1.495,smallerPstartMean, SEM_smallerPstart,'CapSize',0,'DisplayName','Smaller Portfolio');
% Add lines
h1 = line([0 0],[-1.5 2.5]);
% h2 = line([-0.1 -0.1],[-1.5 2.5]);
% % Set properties of lines
set(h1,'Color',[0.5 0.5 0.5],'LineWidth',2)
% % Add a patch
% gray = [0.4 0.4 0.4];

% e1.Bar.LineStyle = 'dotted';
% e2.Bar.LineStyle = 'dotted';
% e3.Bar.LineStyle = 'dotted';
% e1.CapSize = 0.1;

% Create ylabel
ylabel({'Mean Normalized Firing Rate'});

% Create xlabel
xlabel({'Time (secs)'});


hold off


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

