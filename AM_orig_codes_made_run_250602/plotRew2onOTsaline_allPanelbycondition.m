clear all
% plotStuff

% Plot psth from Nex data
% Saline
% rew2on data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m1rew2on.mat')
rew2onUnits_aiS = spikes;
rew2onNames_aiS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m1rew2on.mat')
rew2onUnits_replayS = spikes;
rew2onNames_replayS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m1rew2on.mat')
rew2onUnits_decoyS = spikes;
rew2onNames_decoyS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m1rew2on.mat')
rew2onUnits_liveS = spikes;
rew2onNames_liveS = units;



% OT
% rew2on data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1rew2on.mat')
rew2onUnits_aiOT = spikes;
rew2onNames_aiOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1rew2on.mat')
rew2onUnits_replayOT = spikes;
rew2onNames_replayOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1rew2on.mat')
rew2onUnits_decoyOT = spikes;
rew2onNames_decoyOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1rew2on.mat')
rew2onUnits_liveOT = spikes;
rew2onNames_liveOT = units;



clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(rew2onUnits_aiS);
n = n-1;


    rew2onUnits_aiS = rew2onUnits_aiS(:,end-n:end);
    rew2onUnits_replayS = rew2onUnits_replayS(:,end-n:end);
    rew2onUnits_decoyS = rew2onUnits_decoyS(:,end-n:end);
    rew2onUnits_liveS = rew2onUnits_liveS(:,end-n:end);
    rew2onUnits_aiOT = rew2onUnits_aiOT(:,end-n:end);
    rew2onUnits_replayOT = rew2onUnits_replayOT(:,end-n:end);
    rew2onUnits_decoyOT = rew2onUnits_decoyOT(:,end-n:end);
    rew2onUnits_liveOT = [rew2onUnits_liveOT(:,32:93) rew2onUnits_liveOT(:,103:end)];


% Saline means
rew2onMean_aiS = mean(rew2onUnits_aiS,2);

SEM_rew2on_aiS = std(rew2onUnits_aiS, [], 2)./ sqrt(size(rew2onUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_aiS = bsxfun(@plus, mean(rew2onUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_aiS));   % 95% Confidence Intervals

% replay
rew2onMean_replayS = mean(rew2onUnits_replayS,2);

SEM_rew2on_replayS = std(rew2onUnits_replayS, [], 2)./ sqrt(size(rew2onUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_replayS = bsxfun(@plus, mean(rew2onUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_replayS));   % 95% Confidence Intervals


% decoy
rew2onMean_decoyS = mean(rew2onUnits_decoyS,2);

SEM_rew2on_decoyS = std(rew2onUnits_decoyS, [], 2)./ sqrt(size(rew2onUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_decoyS = bsxfun(@plus, mean(rew2onUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_decoyS));   % 95% Confidence Intervals

% live
rew2onMean_liveS = mean(rew2onUnits_liveS,2);

SEM_rew2on_liveS = std(rew2onUnits_liveS, [], 2)./ sqrt(size(rew2onUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_liveS = bsxfun(@plus, mean(rew2onUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_liveS));   % 95% Confidence Intervals


% OT means
rew2onMean_aiOT = mean(rew2onUnits_aiOT,2);

SEM_rew2on_aiOT = std(rew2onUnits_aiOT, [], 2)./ sqrt(size(rew2onUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_aiOT = bsxfun(@plus, mean(rew2onUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_aiOT));   % 95% Confidence Intervals


% replay
rew2onMean_replayOT = mean(rew2onUnits_replayOT,2);

SEM_rew2on_replayOT = std(rew2onUnits_replayOT, [], 2)./ sqrt(size(rew2onUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_replayOT = bsxfun(@plus, mean(rew2onUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_replayOT));   % 95% Confidence Intervals


% decoy
rew2onMean_decoyOT = mean(rew2onUnits_decoyOT,2);

SEM_rew2on_decoyOT = std(rew2onUnits_decoyOT, [], 2)./ sqrt(size(rew2onUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_decoyOT = bsxfun(@plus, mean(rew2onUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_decoyOT));   % 95% Confidence Intervals


% live
rew2onMean_liveOT = mean(rew2onUnits_liveOT,2);

SEM_rew2on_liveOT = std(rew2onUnits_liveOT, [], 2)./ sqrt(size(rew2onUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_rew2on_liveOT = bsxfun(@plus, mean(rew2onUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew2on_liveOT));   % 95% Confidence Intervals


data(:,:,1) = rew2onUnits_aiS;
data(:,:,2) = rew2onUnits_replayS(:,:);
data(:,:,3) = rew2onUnits_decoyS(:,:);
data(:,:,4) = rew2onUnits_liveS;
data(:,:,5) = rew2onUnits_aiOT;
data(:,:,6) = rew2onUnits_replayOT;
data(:,:,7) = rew2onUnits_decoyOT;
data(:,:,8) = rew2onUnits_liveOT;




time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 12);
% temp = colors(1,:);
colors(4,:) = colors(5,:);
colors(2,:) = colors(7,:);

nl = @(s) strrep(s,'\n',char(10));



subplot(1,2,1);
bl1 = boundedline(time, rew2onMean_aiS, SEM_rew2on_aiS, ...
    time, rew2onMean_replayS, SEM_rew2on_replayS, ...
    time, rew2onMean_decoyS, SEM_rew2on_decoyS, ...
    time, rew2onMean_liveS, SEM_rew2on_liveS, ...
    'cmap', colors);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Saline');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-0.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(1,2,2);
bl1 = boundedline(time, rew2onMean_aiOT, SEM_rew2on_aiOT, ...
    time, rew2onMean_replayOT, SEM_rew2on_replayOT, ...
    time, rew2onMean_decoyOT, SEM_rew2on_decoyOT, ...
    time, rew2onMean_liveOT, SEM_rew2on_liveOT, ...
    'cmap', colors);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('OT');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-0.5 2]); 
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
lpos(1) = lpos(1) + 0.11;
lpos(2) = lpos(2) + 0.11;
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

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData
print('-fillpage','OTsalineRew2OnPanelbyCondition','-dpdf');

