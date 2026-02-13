clear all
% plotStuff

% Plot psth from Nex data
% Saline
% rew1on data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m1rew1on.mat')
rew1onUnits_aiS = spikes;
rew1onNames_aiS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m1rew1on.mat')
rew1onUnits_replayS = spikes;
rew1onNames_replayS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m1rew1on.mat')
rew1onUnits_decoyS = spikes;
rew1onNames_decoyS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m1rew1on.mat')
rew1onUnits_liveS = spikes;
rew1onNames_liveS = units;



% OT
% rew1on data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1rew1on.mat')
rew1onUnits_aiOT = spikes;
rew1onNames_aiOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1rew1on.mat')
rew1onUnits_replayOT = spikes;
rew1onNames_replayOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1rew1on.mat')
rew1onUnits_decoyOT = spikes;
rew1onNames_decoyOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1rew1on.mat')
rew1onUnits_liveOT = spikes;
rew1onNames_liveOT = units;



clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(rew1onUnits_aiS);
n = n-1;


    rew1onUnits_aiS = rew1onUnits_aiS(:,end-n:end);
    rew1onUnits_replayS = rew1onUnits_replayS(:,end-n:end);
    rew1onUnits_decoyS = rew1onUnits_decoyS(:,end-n:end);
    rew1onUnits_liveS = rew1onUnits_liveS(:,end-n:end);
    rew1onUnits_aiOT = rew1onUnits_aiOT(:,end-n:end);
    rew1onUnits_replayOT = rew1onUnits_replayOT(:,end-n:end);
    rew1onUnits_decoyOT = rew1onUnits_decoyOT(:,end-n:end);
    rew1onUnits_liveOT = [rew1onUnits_liveOT(:,32:93) rew1onUnits_liveOT(:,103:end)];


% Saline means
rew1onMean_aiS = mean(rew1onUnits_aiS,2);

SEM_rew1on_aiS = std(rew1onUnits_aiS, [], 2)./ sqrt(size(rew1onUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_aiS = bsxfun(@plus, mean(rew1onUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_aiS));   % 95% Confidence Intervals

% replay
rew1onMean_replayS = mean(rew1onUnits_replayS,2);

SEM_rew1on_replayS = std(rew1onUnits_replayS, [], 2)./ sqrt(size(rew1onUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_replayS = bsxfun(@plus, mean(rew1onUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_replayS));   % 95% Confidence Intervals


% decoy
rew1onMean_decoyS = mean(rew1onUnits_decoyS,2);

SEM_rew1on_decoyS = std(rew1onUnits_decoyS, [], 2)./ sqrt(size(rew1onUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_decoyS = bsxfun(@plus, mean(rew1onUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_decoyS));   % 95% Confidence Intervals

% live
rew1onMean_liveS = mean(rew1onUnits_liveS,2);

SEM_rew1on_liveS = std(rew1onUnits_liveS, [], 2)./ sqrt(size(rew1onUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_liveS = bsxfun(@plus, mean(rew1onUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_liveS));   % 95% Confidence Intervals


% OT means
rew1onMean_aiOT = mean(rew1onUnits_aiOT,2);

SEM_rew1on_aiOT = std(rew1onUnits_aiOT, [], 2)./ sqrt(size(rew1onUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_aiOT = bsxfun(@plus, mean(rew1onUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_aiOT));   % 95% Confidence Intervals


% replay
rew1onMean_replayOT = mean(rew1onUnits_replayOT,2);

SEM_rew1on_replayOT = std(rew1onUnits_replayOT, [], 2)./ sqrt(size(rew1onUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_replayOT = bsxfun(@plus, mean(rew1onUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_replayOT));   % 95% Confidence Intervals


% decoy
rew1onMean_decoyOT = mean(rew1onUnits_decoyOT,2);

SEM_rew1on_decoyOT = std(rew1onUnits_decoyOT, [], 2)./ sqrt(size(rew1onUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_decoyOT = bsxfun(@plus, mean(rew1onUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_decoyOT));   % 95% Confidence Intervals


% live
rew1onMean_liveOT = mean(rew1onUnits_liveOT,2);

SEM_rew1on_liveOT = std(rew1onUnits_liveOT, [], 2)./ sqrt(size(rew1onUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_rew1on_liveOT = bsxfun(@plus, mean(rew1onUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_rew1on_liveOT));   % 95% Confidence Intervals


data(:,:,1) = rew1onUnits_aiS;
data(:,:,2) = rew1onUnits_replayS(:,:);
data(:,:,3) = rew1onUnits_decoyS(:,:);
data(:,:,4) = rew1onUnits_liveS;
data(:,:,5) = rew1onUnits_aiOT;
data(:,:,6) = rew1onUnits_replayOT;
data(:,:,7) = rew1onUnits_decoyOT;
data(:,:,8) = rew1onUnits_liveOT;




time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 12);
% temp = colors(1,:);
colors(4,:) = colors(5,:);
colors(2,:) = colors(7,:);

nl = @(s) strrep(s,'\n',char(10));



subplot(1,2,1);
bl1 = boundedline(time, rew1onMean_aiS, SEM_rew1on_aiS, ...
    time, rew1onMean_replayS, SEM_rew1on_replayS, ...
    time, rew1onMean_decoyS, SEM_rew1on_decoyS, ...
    time, rew1onMean_liveS, SEM_rew1on_liveS, ...
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
bl1 = boundedline(time, rew1onMean_aiOT, SEM_rew1on_aiOT, ...
    time, rew1onMean_replayOT, SEM_rew1on_replayOT, ...
    time, rew1onMean_decoyOT, SEM_rew1on_decoyOT, ...
    time, rew1onMean_liveOT, SEM_rew1on_liveOT, ...
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
print('OTsalineRew1OnPanelbyCondition','-dpdf');

