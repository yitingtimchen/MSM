clear all
% plotStuff

% Plot psth from Nex data
% bubble data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1bubble.mat')
bubbleUnits_ai = spikes;
bubbleNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1bubble.mat')
bubbleUnits_replay = spikes;
bubbleNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1bubble.mat')
bubbleUnits_decoy = spikes;
bubbleNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1bubble.mat')
bubbleUnits_live = spikes;
bubbleNames_live = units;

% Non-Bubble data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1nonbubble.mat')
nonbubbleUnits_ai = spikes;
nonbubbleNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1nonbubble.mat')
nonbubbleUnits_replay = spikes;
nonbubbleNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1nonbubble.mat')
nonbubbleUnits_decoy = spikes;
nonbubbleNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1nonbubble.mat')
nonbubbleUnits_live = spikes;
nonbubbleNames_live = units;

clear units2 units3 spikes2 data

% reduce units to make equal across conditions
bubbleUnits_ai = bubbleUnits_ai(:,1:160);
bubbleUnits_decoy = [bubbleUnits_decoy(:,1:114) bubbleUnits_decoy(:,129:139) bubbleUnits_decoy(:,end-34:end)];
bubbleUnits_live = [bubbleUnits_live(:,1:62) bubbleUnits_live(:,end-97:end)];
nonbubbleUnits_ai = nonbubbleUnits_ai(:,1:160);
nonbubbleUnits_decoy = [nonbubbleUnits_decoy(:,1:114) nonbubbleUnits_decoy(:,129:139) nonbubbleUnits_decoy(:,end-34:end)];
nonbubbleUnits_live = [nonbubbleUnits_live(:,1:62) nonbubbleUnits_live(:,end-97:end)];


bubbleMean_ai = mean(bubbleUnits_ai,2);

SEM_bubble_ai = std(bubbleUnits_ai, [], 2)./ sqrt(size(bubbleUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_bubble_ai = bsxfun(@plus, mean(bubbleUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_bubble_ai));   % 95% Confidence Intervals

nonbubbleMean_ai = mean(nonbubbleUnits_ai,2);

SEM_nonbubble_ai = std(nonbubbleUnits_ai, [], 2)./ sqrt(size(nonbubbleUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_nonbubble_ai = bsxfun(@plus, mean(nonbubbleUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_nonbubble_ai)); 


% replay
bubbleMean_replay = mean(bubbleUnits_replay,2);

SEM_bubble_replay = std(bubbleUnits_replay, [], 2)./ sqrt(size(bubbleUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_bubble_replay = bsxfun(@plus, mean(bubbleUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_bubble_replay));   % 95% Confidence Intervals

nonbubbleMean_replay = mean(nonbubbleUnits_replay,2);

SEM_nonbubble_replay = std(nonbubbleUnits_replay, [], 2)./ sqrt(size(nonbubbleUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_nonbubble_replay = bsxfun(@plus, mean(nonbubbleUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_nonbubble_replay)); 

% decoy
bubbleMean_decoy = mean(bubbleUnits_decoy,2);

SEM_bubble_decoy = std(bubbleUnits_decoy, [], 2)./ sqrt(size(bubbleUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_bubble_decoy = bsxfun(@plus, mean(bubbleUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_bubble_decoy));   % 95% Confidence Intervals

nonbubbleMean_decoy = mean(nonbubbleUnits_decoy,2);

SEM_nonbubble_decoy = std(nonbubbleUnits_decoy, [], 2)./ sqrt(size(nonbubbleUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_nonbubble_decoy = bsxfun(@plus, mean(nonbubbleUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_nonbubble_decoy)); 


% live
bubbleMean_live = mean(bubbleUnits_live,2);

SEM_bubble_live = std(bubbleUnits_live, [], 2)./ sqrt(size(bubbleUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_bubble_live = bsxfun(@plus, mean(bubbleUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_bubble_live));   % 95% Confidence Intervals

nonbubbleMean_live = mean(nonbubbleUnits_live,2);

SEM_nonbubble_live = std(nonbubbleUnits_live, [], 2)./ sqrt(size(nonbubbleUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_nonbubble_live = bsxfun(@plus, mean(nonbubbleUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_nonbubble_live)); 


data(:,:,1) = bubbleUnits_ai;
data(:,:,2) = nonbubbleUnits_ai;
data(:,:,3) = bubbleUnits_replay;
data(:,:,4) = nonbubbleUnits_replay;
data(:,:,5) = bubbleUnits_decoy;
data(:,:,6) = nonbubbleUnits_decoy;
data(:,:,7) = bubbleUnits_live;
data(:,:,8) = nonbubbleUnits_live;



time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 8);
temp = colors(1,:);
colors(1,:) = colors(3,:);
colors(3,:) = temp;

nl = @(s) strrep(s,'\n',char(10));

subplot(4,1,1);
bl1 = boundedline(time, bubbleMean_ai, SEM_bubble_ai, ...
    time, nonbubbleMean_ai, SEM_nonbubble_ai, ...
    'cmap', colors);
hold on
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color','k','LineWidth',0.75)
title('AI');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 1.5]); ylabel(nl('Normalized Firing Rate'));
hold off

subplot(4,1,2);
bl2 = boundedline(time, bubbleMean_replay, SEM_bubble_replay, ...
    time, nonbubbleMean_replay, SEM_nonbubble_replay, ...
    'cmap', colors);
hold on
% Add lines
h4 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h4,'Color','k','LineWidth',0.75)
title('Replay');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 1.5]); ylabel(nl('Normalized Firing Rate'));
hold off

subplot(4,1,3);
bl3 = boundedline(time, bubbleMean_decoy, SEM_bubble_decoy, ...
    time, nonbubbleMean_decoy, SEM_nonbubble_decoy, ...
    'cmap', colors);
hold on
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color','k','LineWidth',0.75)
title('Decoy');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 1.5]); ylabel(nl('Normalized Firing Rate'));
hold off

subplot(4,1,4);
bl4 = boundedline(time, bubbleMean_live, SEM_bubble_live, ...
    time, nonbubbleMean_live, SEM_nonbubble_live, ...
    'cmap', colors);
hold on
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color','k','LineWidth',0.75,'LineStyle','--')
title('Live');

    
 
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 1.5]); ylabel(nl('Normalized Firing Rate'));

 
% instead of a legend, show colored text
lh = legend(bl1);
legnames = {'Bubble','Non-Bubble'};
for i = 1:length(legnames),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(1) = lpos(1) + 0.15;
lh.Position = lpos;

hold off

% hold on
% % Add lines
% 
% h1 = line([0 0],[-1.5 2.5]);
% h3 = line([1.495 1.495],[-1.5 2.5]);
% set(h1,'Color',[0.4 0.4 0.4],'LineWidth',1)
% 
% set(h3,'Color','k','LineWidth',1)
% 
% hold off  

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData
print('-fillpage','bubblePanelbyCondition','-dpdf');

