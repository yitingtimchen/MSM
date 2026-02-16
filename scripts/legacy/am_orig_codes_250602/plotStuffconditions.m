clear all
% plotStuff

% Plot psth from Nex data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1rew1on.mat')
aiUnits = spikes;
aiNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1rew1on.mat')
replayUnits = spikes;
replayNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1rew1on.mat')
decoyUnits = spikes;
decoyNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1rew1on.mat')
liveUnits = spikes;
liveNames = units;

aiUnits = aiUnits(:,1:160);
decoyUnits = [decoyUnits(:,1:114) decoyUnits(:,129:139) decoyUnits(:,end-34:end)];
liveUnits = [liveUnits(:,1:62) liveUnits(:,end-97:end)];

clear units2 units3 spikes2

aiUnits = aiUnits(:,1:160);
decoyUnits = [decoyUnits(:,1:114) decoyUnits(:,129:139) decoyUnits(:,end-34:end)];
liveUnits = [liveUnits(:,1:62) liveUnits(:,end-97:end)];

aiMean = mean(aiUnits,2);

SEM_ai = std(aiUnits, [], 2)./ sqrt(size(aiUnits,2));    % Calculate Standard Error Of The Mean
CI95_ai = bsxfun(@plus, mean(aiUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_ai));   % 95% Confidence Intervals

replayMean = mean(replayUnits,2);

SEM_replay = std(replayUnits, [], 2)./ sqrt(size(replayUnits,2));    % Calculate Standard Error Of The Mean
CI95_replay = bsxfun(@plus, mean(replayUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_replay)); 

decoyMean = mean(decoyUnits,2);

SEM_decoy = std(decoyUnits, [], 2)./ sqrt(size(decoyUnits,2));    % Calculate Standard Error Of The Mean
CI95_decoy = bsxfun(@plus, mean(decoyUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_decoy)); 

liveMean = mean(liveUnits,2);

SEM_live = std(liveUnits, [], 2)./ sqrt(size(liveUnits,2));    % Calculate Standard Error Of The Mean
CI95_live = bsxfun(@plus, mean(liveUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_live)); 


data(:,:,1) = aiUnits;
data(:,:,2) = replayUnits;
data(:,:,3) = decoyUnits;
data(:,:,4) = liveUnits;

time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('seq', 'YlGnBu', 9);
colors = [colors(5:8,:);colors(9,:)];
 
% subplot(4,4,[13 14]);  % plot across two subplots

bl = boundedline(time, aiMean, SEM_ai, ...
    time, replayMean, SEM_replay, ...
    time, decoyMean, SEM_decoy,...
    time, liveMean, SEM_live,...
    'cmap', colors,'transparency',0.5);
hold on
rect2 = fill([-0.25 -0.05 -0.05 -0.25 ],[-2 -2 2.5 2.5],'yellow','LineStyle','none'); 
rect2.FaceAlpha=0.2;
uistack(rect2,'bottom');
set(gca,'FontSize',16);


h2 = line([0 0],[-2 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')

xlim([-0.5 1]); xlabel('Time (s)'); ylim([-2 2.5]); ylabel('Normalized Firing Rate');

 
% instead of a legend, show colored text
lh = legend(bl);
legnames = {'AI', 'Replay', 'Decoy', 'Live','Choice//Reward 1 Onset','Reward 2 Onset'};
for i = 1:length(legnames)-2,
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(2) = lpos(2) - 0.4;
lh.Position = lpos;
lh.FontSize = 16;
 

hold off

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData
print('rew1on_allConditions','-dpdf');

