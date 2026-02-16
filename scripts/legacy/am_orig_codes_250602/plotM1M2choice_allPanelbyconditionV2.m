clear all
% plotStuff

% Plot psth from Nex data
% Saline
% m1rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m1rew1on.mat')
m1rew1onUnits_aiS = spikes;
m1rew1onNames_aiS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m1rew1on.mat')
m1rew1onUnits_replayS = spikes;
m1rew1onNames_replayS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m1rew1on.mat')
m1rew1onUnits_decoyS = spikes;
m1rew1onNames_decoyS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m1rew1on.mat')
m1rew1onUnits_liveS = spikes;
m1rew1onNames_liveS = units;

% m2rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m2rew1on.mat')
m2rew1onUnits_aiS = spikes;
m2rew1onNames_aiS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m2rew1on.mat')
m2rew1onUnits_replayS = spikes;
m2rew1onNames_replayS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m2rew1on.mat')
m2rew1onUnits_decoyS = spikes;
m2rew1onNames_decoyS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m2rew1on.mat')
m2rew1onUnits_liveS = spikes;
m2rew1onNames_liveS = units;

% OT
% m1rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1rew1on.mat')
m1rew1onUnits_aiOT = spikes;
m1rew1onNames_aiOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1rew1on.mat')
m1rew1onUnits_replayOT = spikes;
m1rew1onNames_replayOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1rew1on.mat')
m1rew1onUnits_decoyOT = spikes;
m1rew1onNames_decoyOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1rew1on.mat')
m1rew1onUnits_liveOT = spikes;
m1rew1onNames_liveOT = units;

% m2rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m2rew1on.mat')
m2rew1onUnits_aiOT = spikes;
m2rew1onNames_aiOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m2rew1on.mat')
m2rew1onUnits_replayOT = spikes;
m2rew1onNames_replayOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m2rew1on.mat')
m2rew1onUnits_decoyOT = spikes;
m2rew1onNames_decoyOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m2rew1on.mat')
m2rew1onUnits_liveOT = spikes;
m2rew1onNames_liveOT = units;

clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(m1rew1onUnits_aiS);
n = n-1;


    m1rew1onUnits_aiS = m1rew1onUnits_aiS(:,end-n:end);
    m2rew1onUnits_aiS = m2rew1onUnits_aiS(:,end-n:end);
    m1rew1onUnits_replayS = m1rew1onUnits_replayS(:,end-n:end);
    m2rew1onUnits_replayS = m2rew1onUnits_replayS(:,end-n:end);
    m1rew1onUnits_decoyS = m1rew1onUnits_decoyS(:,end-n:end);
    m2rew1onUnits_decoyS = m2rew1onUnits_decoyS(:,end-n:end);
    m1rew1onUnits_liveS = m1rew1onUnits_liveS(:,end-n:end);
    m2rew1onUnits_liveS = m2rew1onUnits_liveS(:,end-n:end);
    m1rew1onUnits_aiOT = m1rew1onUnits_aiOT(:,end-n:end);
    m2rew1onUnits_aiOT = m2rew1onUnits_aiOT(:,end-n:end);
    m1rew1onUnits_replayOT = m1rew1onUnits_replayOT(:,end-n:end);
    m2rew1onUnits_replayOT = m2rew1onUnits_replayOT(:,end-n:end);
    m1rew1onUnits_decoyOT = m1rew1onUnits_decoyOT(:,end-n:end);
    m2rew1onUnits_decoyOT = m2rew1onUnits_decoyOT(:,end-n:end);
    m1rew1onUnits_liveOT = [m1rew1onUnits_liveOT(:,32:93) m1rew1onUnits_liveOT(:,103:end)];
    m2rew1onUnits_liveOT = [m2rew1onUnits_liveOT(:,32:93) m2rew1onUnits_liveOT(:,103:end)];

% Saline means
m1rew1onMean_aiS = mean(m1rew1onUnits_aiS,2);

SEM_m1rew1on_aiS = std(m1rew1onUnits_aiS, [], 2)./ sqrt(size(m1rew1onUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_aiS = bsxfun(@plus, mean(m1rew1onUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_aiS));   % 95% Confidence Intervals

m2rew1onMean_aiS = mean(m2rew1onUnits_aiS,2);

SEM_m2rew1on_aiS = std(m2rew1onUnits_aiS, [], 2)./ sqrt(size(m2rew1onUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_aiS = bsxfun(@plus, mean(m2rew1onUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_aiS)); 

% replay
m1rew1onMean_replayS = mean(m1rew1onUnits_replayS,2);

SEM_m1rew1on_replayS = std(m1rew1onUnits_replayS, [], 2)./ sqrt(size(m1rew1onUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_replayS = bsxfun(@plus, mean(m1rew1onUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_replayS));   % 95% Confidence Intervals

m2rew1onMean_replayS = mean(m2rew1onUnits_replayS,2);

SEM_m2rew1on_replayS = std(m2rew1onUnits_replayS, [], 2)./ sqrt(size(m2rew1onUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_replayS = bsxfun(@plus, mean(m2rew1onUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_replayS)); 

% decoy
m1rew1onMean_decoyS = mean(m1rew1onUnits_decoyS,2);

SEM_m1rew1on_decoyS = std(m1rew1onUnits_decoyS, [], 2)./ sqrt(size(m1rew1onUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_decoyS = bsxfun(@plus, mean(m1rew1onUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_decoyS));   % 95% Confidence Intervals

m2rew1onMean_decoyS = mean(m2rew1onUnits_decoyS,2);

SEM_m2rew1on_decoyS = std(m2rew1onUnits_decoyS, [], 2)./ sqrt(size(m2rew1onUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_decoyS = bsxfun(@plus, mean(m2rew1onUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_decoyS)); 

% live
m1rew1onMean_liveS = mean(m1rew1onUnits_liveS,2);

SEM_m1rew1on_liveS = std(m1rew1onUnits_liveS, [], 2)./ sqrt(size(m1rew1onUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_liveS = bsxfun(@plus, mean(m1rew1onUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_liveS));   % 95% Confidence Intervals

m2rew1onMean_liveS = mean(m2rew1onUnits_liveS,2);

SEM_m2rew1on_liveS = std(m2rew1onUnits_liveS, [], 2)./ sqrt(size(m2rew1onUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_liveS = bsxfun(@plus, mean(m2rew1onUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_liveS)); 

% OT means
m1rew1onMean_aiOT = mean(m1rew1onUnits_aiOT,2);

SEM_m1rew1on_aiOT = std(m1rew1onUnits_aiOT, [], 2)./ sqrt(size(m1rew1onUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_aiOT = bsxfun(@plus, mean(m1rew1onUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_aiOT));   % 95% Confidence Intervals

m2rew1onMean_aiOT = mean(m2rew1onUnits_aiOT,2);

SEM_m2rew1on_aiOT = std(m2rew1onUnits_aiOT, [], 2)./ sqrt(size(m2rew1onUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_aiOT = bsxfun(@plus, mean(m2rew1onUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_aiOT)); 

% replay
m1rew1onMean_replayOT = mean(m1rew1onUnits_replayOT,2);

SEM_m1rew1on_replayOT = std(m1rew1onUnits_replayOT, [], 2)./ sqrt(size(m1rew1onUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_replayOT = bsxfun(@plus, mean(m1rew1onUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_replayOT));   % 95% Confidence Intervals

m2rew1onMean_replayOT = mean(m2rew1onUnits_replayOT,2);

SEM_m2rew1on_replayOT = std(m2rew1onUnits_replayOT, [], 2)./ sqrt(size(m2rew1onUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_replayOT = bsxfun(@plus, mean(m2rew1onUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_replayOT)); 

% decoy
m1rew1onMean_decoyOT = mean(m1rew1onUnits_decoyOT,2);

SEM_m1rew1on_decoyOT = std(m1rew1onUnits_decoyOT, [], 2)./ sqrt(size(m1rew1onUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_decoyOT = bsxfun(@plus, mean(m1rew1onUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_decoyOT));   % 95% Confidence Intervals

m2rew1onMean_decoyOT = mean(m2rew1onUnits_decoyOT,2);

SEM_m2rew1on_decoyOT = std(m2rew1onUnits_decoyOT, [], 2)./ sqrt(size(m2rew1onUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_decoyOT = bsxfun(@plus, mean(m2rew1onUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_decoyOT)); 

% live
m1rew1onMean_liveOT = mean(m1rew1onUnits_liveOT,2);

SEM_m1rew1on_liveOT = std(m1rew1onUnits_liveOT, [], 2)./ sqrt(size(m1rew1onUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_liveOT = bsxfun(@plus, mean(m1rew1onUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_liveOT));   % 95% Confidence Intervals

m2rew1onMean_liveOT = mean(m2rew1onUnits_liveOT,2);

SEM_m2rew1on_liveOT = std(m2rew1onUnits_liveOT, [], 2)./ sqrt(size(m2rew1onUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_liveOT = bsxfun(@plus, mean(m2rew1onUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_liveOT)); 


data(:,:,1) = m1rew1onUnits_aiS;
data(:,:,2) = m2rew1onUnits_aiS;
data(:,:,3) = m1rew1onUnits_replayS(:,:);
data(:,:,4) = m2rew1onUnits_replayS(:,:);
data(:,:,5) = m1rew1onUnits_decoyS(:,:);
data(:,:,6) = m2rew1onUnits_decoyS(:,:);
data(:,:,7) = m1rew1onUnits_liveS;
data(:,:,8) = m2rew1onUnits_liveS;
data(:,:,9) = m1rew1onUnits_aiOT;
data(:,:,10) = m2rew1onUnits_aiOT;
data(:,:,11) = m1rew1onUnits_replayOT;
data(:,:,12) = m2rew1onUnits_replayOT;
data(:,:,13) = m1rew1onUnits_decoyOT;
data(:,:,14) = m2rew1onUnits_decoyOT;
data(:,:,15) = m1rew1onUnits_liveOT;
data(:,:,16) = m2rew1onUnits_liveOT;



time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('seq', 'YlGnBu', 9);
colors = [colors(5:7,:);colors(9,:)];

nl = @(s) strrep(s,'\n',char(10));



subplot(1,2,1);
bl1 = boundedline(time, m2rew1onMean_aiS, SEM_m2rew1on_aiS, ...
    time, m2rew1onMean_replayS, SEM_m2rew1on_replayS, ...
    time, m2rew1onMean_decoyS, SEM_m2rew1on_decoyS, ...
    time, m2rew1onMean_liveS, SEM_m2rew1on_liveS, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1 -1 1 1],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1 1]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('AI');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1 1]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(1,2,2);
bl1 = boundedline(time, m2rew1onMean_aiOT, SEM_m2rew1on_aiOT, ...
    time, m2rew1onMean_replayOT, SEM_m2rew1on_replayOT, ...
    time, m2rew1onMean_decoyOT, SEM_m2rew1on_decoyOT, ...
    time, m2rew1onMean_liveOT, SEM_m2rew1on_liveOT, ...
    'cmap', colors,'transparency',0.5);
hold on
rect2 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1 -1 1 1],'yellow','LineStyle','none'); 
rect2.FaceAlpha=0.2;
uistack(rect2,'bottom');
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1 1]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')

xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1 1]); 
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
lpos(1) = lpos(1) + 0.1;
lpos(2) = lpos(2) + 0.05;
lh.Position = lpos;
lh.FontSize = 16;

hold off

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData
print('OTsalineM1M2rew1onPanelbyConditionV2','-dpdf');

