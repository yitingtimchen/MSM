clear all
% plotStuff

% Plot psth from Nex data
% largerP data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1largerPrew1.mat')
largerPUnits_ai = spikes;
largerPNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1largerPrew1.mat')
largerPUnits_replay = spikes;
largerPNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1largerPrew1.mat')
largerPUnits_decoy = spikes;
largerPNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1largerPrew1.mat')
largerPUnits_live = spikes;
largerPNames_live = units;

% EvenP data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1evenPrew1.mat')
evenPUnits_ai = spikes;
evenPNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1evenPrew1.mat')
evenPUnits_replay = spikes;
evenPNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1evenPrew1.mat')
evenPUnits_decoy = spikes;
evenPNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1evenPrew1.mat')
evenPUnits_live = spikes;
evenPNames_live = units;

% smallerP data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1smallerPrew1.mat')
smallerPUnits_ai = spikes;
smallerPNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1smallerPrew1.mat')
smallerPUnits_replay = spikes;
smallerPNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1smallerPrew1.mat')
smallerPUnits_decoy = spikes;
smallerPNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1smallerPrew1.mat')
smallerPUnits_live = spikes;
smallerPNames_live = units;


clear units2 units3 spikes2

% reduce units to make equal across conditions
% smallerPUnits_ai = smallerPUnits_ai(:,1:160);
% smallerPUnits_decoy = [smallerPUnits_decoy(:,1:114) smallerPUnits_decoy(:,129:139) smallerPUnits_decoy(:,end-34:end)];
% smallerPUnits_live = [smallerPUnits_live(:,1:62) smallerPUnits_live(:,end-97:end)];
% evenPUnits_ai = evenPUnits_ai(:,1:160);
% evenPUnits_decoy = [evenPUnits_decoy(:,1:114) evenPUnits_decoy(:,129:139) evenPUnits_decoy(:,end-34:end)];
% evenPUnits_live = [evenPUnits_live(:,1:62) evenPUnits_live(:,end-97:end)];
% largerPUnits_ai = largerPUnits_ai(:,1:160);
% largerPUnits_decoy = [largerPUnits_decoy(:,1:114) largerPUnits_decoy(:,129:139) largerPUnits_decoy(:,end-34:end)];
% largerPUnits_live = [largerPUnits_live(:,1:62) largerPUnits_live(:,end-97:end)];



largerPMean_ai = mean(largerPUnits_ai,2);

SEM_largerP_ai = std(largerPUnits_ai, [], 2)./ sqrt(size(largerPUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_largerP_ai = bsxfun(@plus, mean(largerPUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_largerP_ai));   % 95% Confidence Intervals

evenPMean_ai = mean(evenPUnits_ai,2);

SEM_evenP_ai = std(evenPUnits_ai, [], 2)./ sqrt(size(evenPUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_evenP_ai = bsxfun(@plus, mean(evenPUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_evenP_ai)); 

smallerPMean_ai = mean(smallerPUnits_ai,2);

SEM_smallerP_ai = std(smallerPUnits_ai, [], 2)./ sqrt(size(smallerPUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_smallerP_ai = bsxfun(@plus, mean(smallerPUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_smallerP_ai)); 

% replay
largerPMean_replay = mean(largerPUnits_replay,2);

SEM_largerP_replay = std(largerPUnits_replay, [], 2)./ sqrt(size(largerPUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_largerP_replay = bsxfun(@plus, mean(largerPUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_largerP_replay));   % 95% Confidence Intervals

evenPMean_replay = mean(evenPUnits_replay,2);

SEM_evenP_replay = std(evenPUnits_replay, [], 2)./ sqrt(size(evenPUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_evenP_replay = bsxfun(@plus, mean(evenPUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_evenP_replay)); 

smallerPMean_replay = mean(smallerPUnits_replay,2);

SEM_smallerP_replay = std(smallerPUnits_replay, [], 2)./ sqrt(size(smallerPUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_smallerP_replay = bsxfun(@plus, mean(smallerPUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_smallerP_replay)); 

% decoy
largerPMean_decoy = mean(largerPUnits_decoy,2);

SEM_largerP_decoy = std(largerPUnits_decoy, [], 2)./ sqrt(size(largerPUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_largerP_decoy = bsxfun(@plus, mean(largerPUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_largerP_decoy));   % 95% Confidence Intervals

evenPMean_decoy = mean(evenPUnits_decoy,2);

SEM_evenP_decoy = std(evenPUnits_decoy, [], 2)./ sqrt(size(evenPUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_evenP_decoy = bsxfun(@plus, mean(evenPUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_evenP_decoy)); 

smallerPMean_decoy = mean(smallerPUnits_decoy,2);

SEM_smallerP_decoy = std(smallerPUnits_decoy, [], 2)./ sqrt(size(smallerPUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_smallerP_decoy = bsxfun(@plus, mean(smallerPUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_smallerP_decoy)); 

% live
largerPMean_live = mean(largerPUnits_live,2);

SEM_largerP_live = std(largerPUnits_live, [], 2)./ sqrt(size(largerPUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_largerP_live = bsxfun(@plus, mean(largerPUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_largerP_live));   % 95% Confidence Intervals

evenPMean_live = mean(evenPUnits_live,2);

SEM_evenP_live = std(evenPUnits_live, [], 2)./ sqrt(size(evenPUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_evenP_live = bsxfun(@plus, mean(evenPUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_evenP_live)); 

smallerPMean_live = mean(smallerPUnits_live,2);

SEM_smallerP_live = std(smallerPUnits_live, [], 2)./ sqrt(size(smallerPUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_smallerP_live = bsxfun(@plus, mean(smallerPUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_smallerP_live)); 

% USE this as guide:
% aiUnits = aiUnits;
% decoyUnits = [decoyUnits(:,1:114) decoyUnits(:,129:139) decoyUnits(:,end-34:end)];
% liveUnits = [liveUnits(:,1:62) liveUnits(:,end-97:end)];

% data(:,:,1) = largerPUnits_ai;
% data(:,:,2) = evenPUnits_ai;
% data(:,:,3) = smallerPUnits_ai;
% data(:,:,4) = largerPUnits_replay;
% data(:,:,5) = evenPUnits_replay;
% data(:,:,6) = smallerPUnits_replay;
% data(:,:,7) = largerPUnits_decoy;
% data(:,:,8) = evenPUnits_decoy;
% data(:,:,9) = smallerPUnits_decoy;
% data(:,:,10) = largerPUnits_live;
% data(:,:,11) = evenPUnits_live;
% data(:,:,12) = smallerPUnits_live;



time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 8);
temp = colors(1,:);
colors(1,:) = colors(3,:);
colors(3,:) = temp;

nl = @(s) strrep(s,'\n',char(10));

subplot(4,1,1);
bl1 = boundedline(time, largerPMean_ai, SEM_largerP_ai, ...
    time, evenPMean_ai, SEM_evenP_ai, ...
    time, smallerPMean_ai, SEM_smallerP_ai, ...
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
bl2 = boundedline(time, largerPMean_replay, SEM_largerP_replay, ...
    time, evenPMean_replay, SEM_evenP_replay, ...
    time, smallerPMean_replay, SEM_smallerP_replay, ...
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
bl3 = boundedline(time, largerPMean_decoy, SEM_largerP_decoy, ...
    time, evenPMean_decoy, SEM_evenP_decoy, ...
    time, smallerPMean_decoy, SEM_smallerP_decoy, ...
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
bl4 = boundedline(time, largerPMean_live, SEM_largerP_live, ...
    time, evenPMean_live, SEM_evenP_live, ...
    time, smallerPMean_live, SEM_smallerP_live, ...
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
legnames = {'Larger Portfolio','Even Portfolio','Smaller Portfolio'};
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

hold on
% Add lines

h1 = line([0 0],[-1.5 2.5]);
h3 = line([1.495 1.495],[-1.5 2.5]);
set(h1,'Color',[0.4 0.4 0.4],'LineWidth',1)

set(h3,'Color','k','LineWidth',1)

hold off  

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData
print('-fillpage','portfolioDiffPanelbyConditionFM','-dpdf');

