clear all

% plotStuff

% Plot psth from Nex data
% Saline
% m1rew2on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1rew2off.mat')
m1rew2offUnits_aiS = spikes;
m1rew2offNames_aiS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1rew2off.mat')
m1rew2offUnits_replayS = spikes;
m1rew2offNames_replayS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1rew2off.mat')
m1rew2offUnits_decoyS = spikes;
m1rew2offNames_decoyS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1rew2off.mat')
m1rew2offUnits_liveS = spikes;
m1rew2offNames_liveS = units;

% m2rew2off data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m2rew2off.mat')
m2rew2offUnits_aiS = spikes;
m2rew2offNames_aiS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m2rew2off.mat')
m2rew2offUnits_replayS = spikes;
m2rew2offNames_replayS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m2rew2off.mat')
m2rew2offUnits_decoyS = spikes;
m2rew2offNames_decoyS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m2rew2off.mat')
m2rew2offUnits_liveS = spikes;
m2rew2offNames_liveS = units;

% OT
% m1rew2off data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1rew2off.mat')
m1rew2offUnits_aiOT = spikes;
m1rew2offNames_aiOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1rew2off.mat')
m1rew2offUnits_replayOT = spikes;
m1rew2offNames_replayOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1rew2off.mat')
m1rew2offUnits_decoyOT = spikes;
m1rew2offNames_decoyOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1rew2off.mat')
m1rew2offUnits_liveOT = spikes;
m1rew2offNames_liveOT = units;

% m2rew2off data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m2rew2off.mat')
m2rew2offUnits_aiOT = spikes;
m2rew2offNames_aiOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m2rew2off.mat')
m2rew2offUnits_replayOT = spikes;
m2rew2offNames_replayOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m2rew2off.mat')
m2rew2offUnits_decoyOT = spikes;
m2rew2offNames_decoyOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m2rew2off.mat')
m2rew2offUnits_liveOT = spikes;
m2rew2offNames_liveOT = units;

clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(m2rew2offUnits_replayS);
n = n-50;


    m2rew2offUnits_aiS = m2rew2offUnits_aiS(:,end-n:end-9);
    m2rew2offUnits_replayS = m2rew2offUnits_replayS(:,end-n:end-9);
    m2rew2offUnits_decoyS = m2rew2offUnits_decoyS(:,end-n:end-9);
    m2rew2offUnits_liveS = m2rew2offUnits_liveS(:,end-n:end-9);
   
% Saline means
m2rew2offMean_aiS = mean(m2rew2offUnits_aiS,2);

SEM_m2rew2off_aiS = std(m2rew2offUnits_aiS, [], 2)./ sqrt(size(m2rew2offUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew2off_aiS = bsxfun(@plus, mean(m2rew2offUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew2off_aiS)); 

% replay
m2rew2offMean_replayS = mean(m2rew2offUnits_replayS,2);

SEM_m2rew2off_replayS = std(m2rew2offUnits_replayS, [], 2)./ sqrt(size(m2rew2offUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew2off_replayS = bsxfun(@plus, mean(m2rew2offUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew2off_replayS)); 

% decoy
m2rew2offMean_decoyS = mean(m2rew2offUnits_decoyS,2);

SEM_m2rew2off_decoyS = std(m2rew2offUnits_decoyS, [], 2)./ sqrt(size(m2rew2offUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew2off_decoyS = bsxfun(@plus, mean(m2rew2offUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew2off_decoyS)); 

% live
m2rew2offMean_liveS = mean(m2rew2offUnits_liveS,2);

SEM_m2rew2off_liveS = std(m2rew2offUnits_liveS, [], 2)./ sqrt(size(m2rew2offUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_m2rew2off_liveS = bsxfun(@plus, mean(m2rew2offUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew2off_liveS)); 



data(:,:,1) = m2rew2offUnits_aiS;
data(:,:,2) = m2rew2offUnits_replayS(:,:);
data(:,:,3) = m2rew2offUnits_decoyS(:,:);
data(:,:,4) = m2rew2offUnits_liveS;




time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('seq', 'YlGnBu', 9);
colors = [colors(5:7,:);colors(9,:)];

nl = @(s) strrep(s,'\n',char(10));



bl1 = boundedline(time, m2rew2offMean_aiS, SEM_m2rew2off_aiS, ...
    time, m2rew2offMean_replayS, SEM_m2rew2off_replayS, ...
    time, m2rew2offMean_decoyS, SEM_m2rew2off_decoyS, ...
    time, m2rew2offMean_liveS, SEM_m2rew2off_liveS, ...
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
ylabel(nl('Normalized Firing Rate'));
 
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
print('OTsalineM2rew2offV2','-dpdf');

