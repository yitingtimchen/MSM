clear all
% plotStuff

% Plot psth from Nex data
% Saline
% rew2off data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m1rew2off.mat')
rew2offUnits_aiS = spikes;
rew2offNames_aiS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m1rew2off.mat')
rew2offUnits_replayS = spikes;
rew2offNames_replayS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m1rew2off.mat')
rew2offUnits_decoyS = spikes;
rew2offNames_decoyS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m1rew2off.mat')
rew2offUnits_liveS = spikes;
rew2offNames_liveS = units;



% OT
% rew2off data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1rew2off.mat')
rew2offUnits_aiOT = spikes;
rew2offNames_aiOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1rew2off.mat')
rew2offUnits_replayOT = spikes;
rew2offNames_replayOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1rew2off.mat')
rew2offUnits_decoyOT = spikes;
rew2offNames_decoyOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1rew2off.mat')
rew2offUnits_liveOT = spikes;
rew2offNames_liveOT = units;



clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(rew2offUnits_aiS);
n = n-1;


    rew2offUnits_aiS = rew2offUnits_aiS(:,end-n:end);
    rew2offUnits_replayS = rew2offUnits_replayS(:,end-n:end);
    rew2offUnits_decoyS = rew2offUnits_decoyS(:,end-n:end);
    rew2offUnits_liveS = rew2offUnits_liveS(:,end-n:end);
    rew2offUnits_aiOT = rew2offUnits_aiOT(:,end-n:end);
    rew2offUnits_replayOT = rew2offUnits_replayOT(:,end-n:end);
    rew2offUnits_decoyOT = rew2offUnits_decoyOT(:,end-n:end);
    rew2offUnits_liveOT = [rew2offUnits_liveOT(:,32:93) rew2offUnits_liveOT(:,103:end)];


% Saline means
rew2offMean_aiS = mean(rew2offUnits_aiS,2);

SEM_rew2off_aiS = std(rew2offUnits_aiS, [], 2)./ sqrt(size(rew2offUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_aiS = bsxfun(@plus, mean(rew2offUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_aiS));   % 95% Confidence Intervals

% replay
rew2offMean_replayS = mean(rew2offUnits_replayS,2);

SEM_rew2off_replayS = std(rew2offUnits_replayS, [], 2)./ sqrt(size(rew2offUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_replayS = bsxfun(@plus, mean(rew2offUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_replayS));   % 95% Confidence Intervals


% decoy
rew2offMean_decoyS = mean(rew2offUnits_decoyS,2);

SEM_rew2off_decoyS = std(rew2offUnits_decoyS, [], 2)./ sqrt(size(rew2offUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_decoyS = bsxfun(@plus, mean(rew2offUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_decoyS));   % 95% Confidence Intervals

% live
rew2offMean_liveS = mean(rew2offUnits_liveS,2);

SEM_rew2off_liveS = std(rew2offUnits_liveS, [], 2)./ sqrt(size(rew2offUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_liveS = bsxfun(@plus, mean(rew2offUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_liveS));   % 95% Confidence Intervals


% OT means
rew2offMean_aiOT = mean(rew2offUnits_aiOT,2);

SEM_rew2off_aiOT = std(rew2offUnits_aiOT, [], 2)./ sqrt(size(rew2offUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_aiOT = bsxfun(@plus, mean(rew2offUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_aiOT));   % 95% Confidence Intervals


% replay
rew2offMean_replayOT = mean(rew2offUnits_replayOT,2);

SEM_rew2off_replayOT = std(rew2offUnits_replayOT, [], 2)./ sqrt(size(rew2offUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_replayOT = bsxfun(@plus, mean(rew2offUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_replayOT));   % 95% Confidence Intervals


% decoy
rew2offMean_decoyOT = mean(rew2offUnits_decoyOT,2);

SEM_rew2off_decoyOT = std(rew2offUnits_decoyOT, [], 2)./ sqrt(size(rew2offUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_decoyOT = bsxfun(@plus, mean(rew2offUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_decoyOT));   % 95% Confidence Intervals


% live
rew2offMean_liveOT = mean(rew2offUnits_liveOT,2);

SEM_rew2off_liveOT = std(rew2offUnits_liveOT, [], 2)./ sqrt(size(rew2offUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2off_liveOT = bsxfun(@plus, mean(rew2offUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2off_liveOT));   % 95% Confidence Intervals


data(:,:,1) = rew2offUnits_aiS;
data(:,:,2) = rew2offUnits_replayS(:,:);
data(:,:,3) = rew2offUnits_decoyS(:,:);
data(:,:,4) = rew2offUnits_liveS;
data(:,:,5) = rew2offUnits_aiOT;
data(:,:,6) = rew2offUnits_replayOT;
data(:,:,7) = rew2offUnits_decoyOT;
data(:,:,8) = rew2offUnits_liveOT;




time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 12);
% temp = colors(1,:);
colors(4,:) = colors(5,:);
colors(2,:) = colors(7,:);

nl = @(s) strrep(s,'\n',char(10));



subplot(1,2,1);
bl1 = boundedline(time, rew2offMean_aiS, SEM_rew2off_aiS, ...
    time, rew2offMean_replayS, SEM_rew2off_replayS, ...
    time, rew2offMean_decoyS, SEM_rew2off_decoyS, ...
    time, rew2offMean_liveS, SEM_rew2off_liveS, ...
    'cmap', colors);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Saline');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 2.5]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(1,2,2);
bl1 = boundedline(time, rew2offMean_aiOT, SEM_rew2off_aiOT, ...
    time, rew2offMean_replayOT, SEM_rew2off_replayOT, ...
    time, rew2offMean_decoyOT, SEM_rew2off_decoyOT, ...
    time, rew2offMean_liveOT, SEM_rew2off_liveOT, ...
    'cmap', colors);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('OT');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 2.5]); 
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
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData
print('OTsalinerew2offPanelbyCondition','-dpdf');

