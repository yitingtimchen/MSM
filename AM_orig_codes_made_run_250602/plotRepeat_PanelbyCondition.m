clear all
% plotStuff

% Plot psth from Nex data
% repeat data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1repeat.mat')
repeatUnits_ai = spikes;
repeatNames_ai = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1repeat.mat')
repeatUnits_replay = spikes;
repeatNames_replay = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1repeat.mat')
repeatUnits_decoy = spikes;
repeatNames_decoy = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1repeat.mat')
repeatUnits_live = spikes;
repeatNames_live = units;

% Non-repeat data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1nonrepeat.mat')
nonrepeatUnits_ai = spikes;
nonrepeatNames_ai = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1nonrepeat.mat')
nonrepeatUnits_replay = spikes;
nonrepeatNames_replay = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1nonrepeat.mat')
nonrepeatUnits_decoy = spikes;
nonrepeatNames_decoy = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1nonrepeat.mat')
nonrepeatUnits_live = spikes;
nonrepeatNames_live = units;

clear units2 units3 spikes2

repeatUnits_ai = repeatUnits_ai(:,1:160);
repeatUnits_decoy = [repeatUnits_decoy(:,1:114) repeatUnits_decoy(:,129:139) repeatUnits_decoy(:,end-34:end)];
repeatUnits_live = [repeatUnits_live(:,1:62) repeatUnits_live(:,end-97:end)];
nonrepeatUnits_ai = nonrepeatUnits_ai(:,1:160);
nonrepeatUnits_decoy = [nonrepeatUnits_decoy(:,1:114) nonrepeatUnits_decoy(:,129:139) nonrepeatUnits_decoy(:,end-34:end)];
nonrepeatUnits_live = [nonrepeatUnits_live(:,1:62) nonrepeatUnits_live(:,end-97:end)];

repeatMean_ai = mean(repeatUnits_ai,2);

SEM_repeat_ai = std(repeatUnits_ai, [], 2)./ sqrt(size(repeatUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_repeat_ai = bsxfun(@plus, mean(repeatUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_ai));   % 95% Confidence Intervals

nonrepeatMean_ai = mean(nonrepeatUnits_ai,2);

SEM_nonrepeat_ai = std(nonrepeatUnits_ai, [], 2)./ sqrt(size(nonrepeatUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_ai = bsxfun(@plus, mean(nonrepeatUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_ai)); 


% replay
repeatMean_replay = mean(repeatUnits_replay,2);

SEM_repeat_replay = std(repeatUnits_replay, [], 2)./ sqrt(size(repeatUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_repeat_replay = bsxfun(@plus, mean(repeatUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_replay));   % 95% Confidence Intervals

nonrepeatMean_replay = mean(nonrepeatUnits_replay,2);

SEM_nonrepeat_replay = std(nonrepeatUnits_replay, [], 2)./ sqrt(size(nonrepeatUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_replay = bsxfun(@plus, mean(nonrepeatUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_replay)); 

% decoy
repeatMean_decoy = mean(repeatUnits_decoy,2);

SEM_repeat_decoy = std(repeatUnits_decoy, [], 2)./ sqrt(size(repeatUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_repeat_decoy = bsxfun(@plus, mean(repeatUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_decoy));   % 95% Confidence Intervals

nonrepeatMean_decoy = mean(nonrepeatUnits_decoy,2);

SEM_nonrepeat_decoy = std(nonrepeatUnits_decoy, [], 2)./ sqrt(size(nonrepeatUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_decoy = bsxfun(@plus, mean(nonrepeatUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_decoy)); 


% live
repeatMean_live = mean(repeatUnits_live,2);

SEM_repeat_live = std(repeatUnits_live, [], 2)./ sqrt(size(repeatUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_repeat_live = bsxfun(@plus, mean(repeatUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_repeat_live));   % 95% Confidence Intervals

nonrepeatMean_live = mean(nonrepeatUnits_live,2);

SEM_nonrepeat_live = std(nonrepeatUnits_live, [], 2)./ sqrt(size(nonrepeatUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_nonrepeat_live = bsxfun(@plus, mean(nonrepeatUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_nonrepeat_live)); 


data(:,:,1) = repeatUnits_ai(:,1:94);
data(:,:,2) = nonrepeatUnits_ai(:,1:94);
data(:,:,4) = repeatUnits_replay(:,1:94);
data(:,:,5) = nonrepeatUnits_replay(:,1:94);
data(:,:,7) = repeatUnits_decoy(:,1:94);
data(:,:,8) = nonrepeatUnits_decoy(:,1:94);
data(:,:,10) = repeatUnits_live(:,1:94);
data(:,:,11) = nonrepeatUnits_live(:,1:94);



time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set2', 8);
colors = colors(3:4,:);

nl = @(s) strrep(s,'\n',char(10));

subplot(1,4,1);
bl1 = boundedline(time, repeatMean_ai, SEM_repeat_ai, ...
    time, nonrepeatMean_ai, SEM_nonrepeat_ai, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2 2],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('AI');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1 2]); ylabel(nl('Normalized Firing Rate'));
hold off

subplot(1,4,2);
bl2 = boundedline(time, repeatMean_replay, SEM_repeat_replay, ...
    time, nonrepeatMean_replay, SEM_nonrepeat_replay, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2 2],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h1 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h1,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Replay');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1 2]);
hold off

subplot(1,4,3);
bl3 = boundedline(time, repeatMean_decoy, SEM_repeat_decoy, ...
    time, nonrepeatMean_decoy, SEM_nonrepeat_decoy, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2 2],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Decoy');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1 2]);
hold off

subplot(1,4,4);
bl4 = boundedline(time, repeatMean_live, SEM_repeat_live, ...
    time, nonrepeatMean_live, SEM_nonrepeat_live, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2 2],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h6 = line([0 0],[-1.5 2.5]);
set(h6,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% Set properties of lines
title('Live');

    
 
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1 2]);

 
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
lpos(1) = lpos(1) + 0.87;
lpos(2) = lpos(2) - 0.5;
lh.Position = lpos;
lh.FontSize = 16;
lh.Interruptible = 'on';
hold off

fig = gcf;
orient(fig,'landscape');
fig.PaperType = 'uslegal';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 14 3.5];

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData
print('repeatPanelbyCondition','-dpdf');

