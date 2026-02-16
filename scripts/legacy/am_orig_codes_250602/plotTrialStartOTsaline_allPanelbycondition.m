clear all
% plotStuff

% Plot psth from Nex data
% Saline
% start data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m1start.mat')
startUnits_aiS = spikes;
startNames_aiS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m1start.mat')
startUnits_replayS = spikes;
startNames_replayS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m1start.mat')
startUnits_decoyS = spikes;
startNames_decoyS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m1start.mat')
startUnits_liveS = spikes;
startNames_liveS = units;



% OT
% start data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1start.mat')
startUnits_aiOT = spikes;
startNames_aiOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1start.mat')
startUnits_replayOT = spikes;
startNames_replayOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1start.mat')
startUnits_decoyOT = spikes;
startNames_decoyOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1start.mat')
startUnits_liveOT = spikes;
startNames_liveOT = units;



clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(startUnits_aiS);
n = n-1;


    startUnits_aiS = startUnits_aiS(:,end-n:end);
    startUnits_replayS = startUnits_replayS(:,end-n:end);
    startUnits_decoyS = startUnits_decoyS(:,end-n:end);
    startUnits_liveS = startUnits_liveS(:,end-n:end);
    startUnits_aiOT = startUnits_aiOT(:,end-n:end);
    startUnits_replayOT = startUnits_replayOT(:,end-n:end);
    startUnits_decoyOT = startUnits_decoyOT(:,end-n:end);
    startUnits_liveOT = [startUnits_liveOT(:,32:93) startUnits_liveOT(:,103:end)];


% Saline means
startMean_aiS = mean(startUnits_aiS,2);

SEM_start_aiS = std(startUnits_aiS, [], 2)./ sqrt(size(startUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_start_aiS = bsxfun(@plus, mean(startUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_start_aiS));   % 95% Confidence Intervals

% replay
startMean_replayS = mean(startUnits_replayS,2);

SEM_start_replayS = std(startUnits_replayS, [], 2)./ sqrt(size(startUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_start_replayS = bsxfun(@plus, mean(startUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_start_replayS));   % 95% Confidence Intervals


% decoy
startMean_decoyS = mean(startUnits_decoyS,2);

SEM_start_decoyS = std(startUnits_decoyS, [], 2)./ sqrt(size(startUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_start_decoyS = bsxfun(@plus, mean(startUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_start_decoyS));   % 95% Confidence Intervals

% live
startMean_liveS = mean(startUnits_liveS,2);

SEM_start_liveS = std(startUnits_liveS, [], 2)./ sqrt(size(startUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_start_liveS = bsxfun(@plus, mean(startUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_start_liveS));   % 95% Confidence Intervals


% OT means
startMean_aiOT = mean(startUnits_aiOT,2);

SEM_start_aiOT = std(startUnits_aiOT, [], 2)./ sqrt(size(startUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_start_aiOT = bsxfun(@plus, mean(startUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_start_aiOT));   % 95% Confidence Intervals


% replay
startMean_replayOT = mean(startUnits_replayOT,2);

SEM_start_replayOT = std(startUnits_replayOT, [], 2)./ sqrt(size(startUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_start_replayOT = bsxfun(@plus, mean(startUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_start_replayOT));   % 95% Confidence Intervals


% decoy
startMean_decoyOT = mean(startUnits_decoyOT,2);

SEM_start_decoyOT = std(startUnits_decoyOT, [], 2)./ sqrt(size(startUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_start_decoyOT = bsxfun(@plus, mean(startUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_start_decoyOT));   % 95% Confidence Intervals


% live
startMean_liveOT = mean(startUnits_liveOT,2);

SEM_start_liveOT = std(startUnits_liveOT, [], 2)./ sqrt(size(startUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_start_liveOT = bsxfun(@plus, mean(startUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_start_liveOT));   % 95% Confidence Intervals


data(:,:,1) = startUnits_aiS;
data(:,:,2) = startUnits_replayS(:,:);
data(:,:,3) = startUnits_decoyS(:,:);
data(:,:,4) = startUnits_liveS;
data(:,:,5) = startUnits_aiOT;
data(:,:,6) = startUnits_replayOT;
data(:,:,7) = startUnits_decoyOT;
data(:,:,8) = startUnits_liveOT;




time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 12);
% temp = colors(1,:);
colors(4,:) = colors(5,:);
colors(2,:) = colors(7,:);

nl = @(s) strrep(s,'\n',char(10));



subplot(1,2,1);
bl1 = boundedline(time, startMean_aiS, SEM_start_aiS, ...
    time, startMean_replayS, SEM_start_replayS, ...
    time, startMean_decoyS, SEM_start_decoyS, ...
    time, startMean_liveS, SEM_start_liveS, ...
    'cmap', colors);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-2.0 3.0]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Saline');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-2.0 3.0]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(1,2,2);
bl1 = boundedline(time, startMean_aiOT, SEM_start_aiOT, ...
    time, startMean_replayOT, SEM_start_replayOT, ...
    time, startMean_decoyOT, SEM_start_decoyOT, ...
    time, startMean_liveOT, SEM_start_liveOT, ...
    'cmap', colors);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-2.0 3.0]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('OT');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-2.0 3.0]); 
% ylabel(nl('Normalized Firing Rate'));

 
% instead of a legend, show colored text
lh = legend(bl1);
legnames = {'AI','Replay','Decoy','Live'};
for i = 1:length(legnames),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(1) = lpos(1) + 0.14;
lpos(2) = lpos(2) + 0.125;
lh.Position = lpos;
lh.FontSize = 16;



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

orient('landscape')
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData
print('OTsalinestartPanelbyCondition','-dpdf');

