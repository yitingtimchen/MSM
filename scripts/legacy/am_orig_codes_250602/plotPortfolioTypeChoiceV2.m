clear all
% plotStuff

% Plot psth from Nex data
% larger Portfolio data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1largerP.mat')
largerPUnits = spikes;
largerPNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1largerP.mat')
largerPUnits = [largerPUnits, spikes(:,1:92)];
largerPNames = [largerPNames, units(:,1:92)];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1largerP.mat')
largerPUnits = [largerPUnits, spikes];
largerPNames = [largerPNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1largerP.mat')
largerPUnits = [largerPUnits, spikes];
largerPNames = [largerPNames, units];

% evenP data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1evenP.mat')
evenPUnits = spikes;
evenPNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1evenP.mat')
evenPUnits = [evenPUnits, spikes(:,1:92)];
evenPNames = [evenPNames, units(:,1:92)];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1evenP.mat')
evenPUnits = [evenPUnits, spikes];
evenPNames = [evenPNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1evenP.mat')
evenPUnits = [evenPUnits, spikes];
evenPNames = [evenPNames, units];

% smallerP data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1smallerP.mat')
smallerPUnits = spikes;
smallerPNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1smallerP.mat')
smallerPUnits = [smallerPUnits, spikes(:,1:92)];
smallerPNames = [smallerPNames, units(:,1:92)];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1smallerP.mat')
smallerPUnits = [smallerPUnits, spikes];
smallerPNames = [smallerPNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1smallerP.mat')
smallerPUnits = [smallerPUnits, spikes];
smallerPNames = [smallerPNames, units];

clear units2 units3 spikes2

largerPMean = mean(largerPUnits,2);

SEM_largerP = std(largerPUnits, [], 2)./ sqrt(size(largerPUnits,2));    % Calculate Standard Error Of The Mean
CI95_largerP = bsxfun(@plus, mean(largerPUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_largerP));   % 95% Confidence Intervals

evenPMean = mean(evenPUnits,2);

SEM_evenP = std(evenPUnits, [], 2)./ sqrt(size(evenPUnits,2));    % Calculate Standard Error Of The Mean
CI95_evenP = bsxfun(@plus, mean(evenPUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_evenP)); 

smallerPMean = mean(smallerPUnits,2);

SEM_smallerP = std(smallerPUnits, [], 2)./ sqrt(size(smallerPUnits,2));    % Calculate Standard Error Of The Mean
CI95_smallerP = bsxfun(@plus, mean(smallerPUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_smallerP)); 


e1 = errorbar(-0.5:0.005:1.495,largerPMean,SEM_largerP,'CapSize',0,'DisplayName','Larger Portfolio');
hold on
e2 = errorbar(-0.5:0.005:1.495,evenPMean, SEM_evenP,'CapSize',0,'DisplayName','Even Portfolio');
e3 = errorbar(-0.5:0.005:1.495,smallerPMean, SEM_smallerP,'CapSize',0,'DisplayName','Smaller Portfolio');
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

