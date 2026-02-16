clear all
% plotStuff

% Plot psth from Nex data
% Saline
% repeat data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m1repeat.mat')
repeatUnits_aiS = spikes;
repeatNames_aiS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m1repeat.mat')
repeatUnits_replayS = spikes;
repeatNames_replayS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m1repeat.mat')
repeatUnits_decoyS = spikes;
repeatNames_decoyS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m1repeat.mat')
repeatUnits_liveS = spikes;
repeatNames_liveS = units;

% nonrepeat data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m1nonrepeat.mat')
nonrepeatUnits_aiS = spikes;
nonrepeatNames_aiS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m1nonrepeat.mat')
nonrepeatUnits_replayS = spikes;
nonrepeatNames_replayS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m1nonrepeat.mat')
nonrepeatUnits_decoyS = spikes;
nonrepeatNames_decoyS = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m1nonrepeat.mat')
nonrepeatUnits_liveS = spikes;
nonrepeatNames_liveS = units;

% OT
% repeat data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1repeat.mat')
repeatUnits_aiOT = spikes;
repeatNames_aiOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1repeat.mat')
repeatUnits_replayOT = spikes;
repeatNames_replayOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1repeat.mat')
repeatUnits_decoyOT = spikes;
repeatNames_decoyOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1repeat.mat')
repeatUnits_liveOT = spikes;
repeatNames_liveOT = units;

% nonrepeat data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m1nonrepeat.mat')
nonrepeatUnits_aiOT = spikes;
nonrepeatNames_aiOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m1nonrepeat.mat')
nonrepeatUnits_replayOT = spikes;
nonrepeatNames_replayOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m1nonrepeat.mat')
nonrepeatUnits_decoyOT = spikes;
nonrepeatNames_decoyOT = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m1nonrepeat.mat')
nonrepeatUnits_liveOT = spikes;
nonrepeatNames_liveOT = units;

clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(repeatUnits_aiS);
n = n-1;


    repeatUnits_aiS = repeatUnits_aiS(:,end-n:end);
    nonrepeatUnits_aiS = nonrepeatUnits_aiS(:,end-n:end);
    repeatUnits_replayS = repeatUnits_replayS(:,end-n:end);
    nonrepeatUnits_replayS = nonrepeatUnits_replayS(:,end-n:end);
    repeatUnits_decoyS = repeatUnits_decoyS(:,end-n:end);
    nonrepeatUnits_decoyS = nonrepeatUnits_decoyS(:,end-n:end);
    repeatUnits_liveS = repeatUnits_liveS(:,end-n:end);
    nonrepeatUnits_liveS = nonrepeatUnits_liveS(:,end-n:end);
    repeatUnits_aiOT = repeatUnits_aiOT(:,end-n:end);
    nonrepeatUnits_aiOT = nonrepeatUnits_aiOT(:,end-n:end);
    repeatUnits_replayOT = repeatUnits_replayOT(:,end-n:end);
    nonrepeatUnits_replayOT = nonrepeatUnits_replayOT(:,end-n:end);
    repeatUnits_decoyOT = repeatUnits_decoyOT(:,end-n:end);
    nonrepeatUnits_decoyOT = nonrepeatUnits_decoyOT(:,end-n:end);
    repeatUnits_liveOT = [repeatUnits_liveOT(:,32:93) repeatUnits_liveOT(:,103:end)];
    nonrepeatUnits_liveOT = [nonrepeatUnits_liveOT(:,32:93) nonrepeatUnits_liveOT(:,103:end)];

% Saline means
repeatMean_aiS = mean(repeatUnits_aiS,2);

SEM_repeat_aiS = std(repeatUnits_aiS, [], 2)./ sqrt(size(repeatUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_repeat_aiS = bsxfun(@plus, mean(repeatUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_aiS));   % 95% Confidence Intervals

nonrepeatMean_aiS = mean(nonrepeatUnits_aiS,2);

SEM_nonrepeat_aiS = std(nonrepeatUnits_aiS, [], 2)./ sqrt(size(nonrepeatUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_aiS = bsxfun(@plus, mean(nonrepeatUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_aiS)); 

% replay
repeatMean_replayS = mean(repeatUnits_replayS,2);

SEM_repeat_replayS = std(repeatUnits_replayS, [], 2)./ sqrt(size(repeatUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_repeat_replayS = bsxfun(@plus, mean(repeatUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_replayS));   % 95% Confidence Intervals

nonrepeatMean_replayS = mean(nonrepeatUnits_replayS,2);

SEM_nonrepeat_replayS = std(nonrepeatUnits_replayS, [], 2)./ sqrt(size(nonrepeatUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_replayS = bsxfun(@plus, mean(nonrepeatUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_replayS)); 

% decoy
repeatMean_decoyS = mean(repeatUnits_decoyS,2);

SEM_repeat_decoyS = std(repeatUnits_decoyS, [], 2)./ sqrt(size(repeatUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_repeat_decoyS = bsxfun(@plus, mean(repeatUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_decoyS));   % 95% Confidence Intervals

nonrepeatMean_decoyS = mean(nonrepeatUnits_decoyS,2);

SEM_nonrepeat_decoyS = std(nonrepeatUnits_decoyS, [], 2)./ sqrt(size(nonrepeatUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_decoyS = bsxfun(@plus, mean(nonrepeatUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_decoyS)); 

% live
repeatMean_liveS = mean(repeatUnits_liveS,2);

SEM_repeat_liveS = std(repeatUnits_liveS, [], 2)./ sqrt(size(repeatUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_repeat_liveS = bsxfun(@plus, mean(repeatUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_liveS));   % 95% Confidence Intervals

nonrepeatMean_liveS = mean(nonrepeatUnits_liveS,2);

SEM_nonrepeat_liveS = std(nonrepeatUnits_liveS, [], 2)./ sqrt(size(nonrepeatUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_liveS = bsxfun(@plus, mean(nonrepeatUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_liveS)); 

% OT means
repeatMean_aiOT = mean(repeatUnits_aiOT,2);

SEM_repeat_aiOT = std(repeatUnits_aiOT, [], 2)./ sqrt(size(repeatUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_repeat_aiOT = bsxfun(@plus, mean(repeatUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_aiOT));   % 95% Confidence Intervals

nonrepeatMean_aiOT = mean(nonrepeatUnits_aiOT,2);

SEM_nonrepeat_aiOT = std(nonrepeatUnits_aiOT, [], 2)./ sqrt(size(nonrepeatUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_aiOT = bsxfun(@plus, mean(nonrepeatUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_aiOT)); 

% replay
repeatMean_replayOT = mean(repeatUnits_replayOT,2);

SEM_repeat_replayOT = std(repeatUnits_replayOT, [], 2)./ sqrt(size(repeatUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_repeat_replayOT = bsxfun(@plus, mean(repeatUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_replayOT));   % 95% Confidence Intervals

nonrepeatMean_replayOT = mean(nonrepeatUnits_replayOT,2);

SEM_nonrepeat_replayOT = std(nonrepeatUnits_replayOT, [], 2)./ sqrt(size(nonrepeatUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_replayOT = bsxfun(@plus, mean(nonrepeatUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_replayOT)); 

% decoy
repeatMean_decoyOT = mean(repeatUnits_decoyOT,2);

SEM_repeat_decoyOT = std(repeatUnits_decoyOT, [], 2)./ sqrt(size(repeatUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_repeat_decoyOT = bsxfun(@plus, mean(repeatUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_decoyOT));   % 95% Confidence Intervals

nonrepeatMean_decoyOT = mean(nonrepeatUnits_decoyOT,2);

SEM_nonrepeat_decoyOT = std(nonrepeatUnits_decoyOT, [], 2)./ sqrt(size(nonrepeatUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_decoyOT = bsxfun(@plus, mean(nonrepeatUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_decoyOT)); 

% live
repeatMean_liveOT = mean(repeatUnits_liveOT,2);

SEM_repeat_liveOT = std(repeatUnits_liveOT, [], 2)./ sqrt(size(repeatUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_repeat_liveOT = bsxfun(@plus, mean(repeatUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_liveOT));   % 95% Confidence Intervals

nonrepeatMean_liveOT = mean(nonrepeatUnits_liveOT,2);

SEM_nonrepeat_liveOT = std(nonrepeatUnits_liveOT, [], 2)./ sqrt(size(nonrepeatUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_liveOT = bsxfun(@plus, mean(nonrepeatUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_liveOT)); 


data(:,:,1) = repeatUnits_aiS;
data(:,:,2) = nonrepeatUnits_aiS;
data(:,:,3) = repeatUnits_replayS(:,:);
data(:,:,4) = nonrepeatUnits_replayS(:,:);
data(:,:,5) = repeatUnits_decoyS(:,:);
data(:,:,6) = nonrepeatUnits_decoyS(:,:);
data(:,:,7) = repeatUnits_liveS;
data(:,:,8) = nonrepeatUnits_liveS;
data(:,:,9) = repeatUnits_aiOT;
data(:,:,10) = nonrepeatUnits_aiOT;
data(:,:,11) = repeatUnits_replayOT;
data(:,:,12) = nonrepeatUnits_replayOT;
data(:,:,13) = repeatUnits_decoyOT;
data(:,:,14) = nonrepeatUnits_decoyOT;
data(:,:,15) = repeatUnits_liveOT;
data(:,:,16) = nonrepeatUnits_liveOT;



time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set2', 8);
colors = [colors(1,:);colors(8,:)];

nl = @(s) strrep(s,'\n',char(10));



subplot(2,4,1);
bl1 = boundedline(time, repeatMean_aiS, SEM_repeat_aiS, ...
    time, nonrepeatMean_aiS, SEM_nonrepeat_aiS, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('AI');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,5);
bl1 = boundedline(time, repeatMean_aiOT, SEM_repeat_aiOT, ...
    time, nonrepeatMean_aiOT, SEM_nonrepeat_aiOT, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('AI');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,2);
bl2 = boundedline(time, repeatMean_replayS, SEM_repeat_replayS, ...
    time, nonrepeatMean_replayS, SEM_nonrepeat_replayS, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h4 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h4,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Replay');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,6);
bl2 = boundedline(time, repeatMean_replayOT, SEM_repeat_replayOT, ...
    time, nonrepeatMean_replayOT, SEM_nonrepeat_replayOT, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h4 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h4,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Replay');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,3);
bl3 = boundedline(time, repeatMean_decoyS, SEM_repeat_decoyS, ...
    time, nonrepeatMean_decoyS, SEM_nonrepeat_decoyS, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Decoy');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,7);
bl3 = boundedline(time, repeatMean_decoyOT, SEM_repeat_decoyOT, ...
    time, nonrepeatMean_decoyOT, SEM_nonrepeat_decoyOT, ...
    'cmap', colors, 'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Decoy');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,4);
bl4 = boundedline(time, repeatMean_liveS, SEM_repeat_liveS, ...
    time, nonrepeatMean_liveS, SEM_nonrepeat_liveS, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Live');
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,8);
bl4 = boundedline(time, repeatMean_liveOT, SEM_repeat_liveOT, ...
    time, nonrepeatMean_liveOT, SEM_nonrepeat_liveOT, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Live');

    
 
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));

 
% instead of a legend, show colored text
lh = legend(bl1);
legnames = {'Repeat','Non-Repeat'};
for i = 1:length(legnames),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(1) = lpos(1) -0.1;
lpos(2) = lpos(2) + 0.10;
lh.Position = lpos;
lh.FontSize = 16;



hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 6];

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
print('OTsalineRepeatPanelbyCondition','-dpdf');

