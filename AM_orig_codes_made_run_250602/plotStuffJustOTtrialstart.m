clear all
% plotStuff

% Plot psth from Nex data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\AI_OT
load('m1start.mat')
aiUnits_OT = spikes;
aiNames_OT = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Replay_OT
load('m1start.mat')
replayUnits_OT = spikes;
replayNames_OT = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy_OT
load('m1start.mat')
decoyUnits_OT = spikes;
decoyNames_OT = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Live_OT
load('m1start.mat')
liveUnits_OT = spikes;
liveNames_OT = units;

clear units2 units3 spikes2

aiMean_OT = mean(aiUnits_OT,2);

SEM_ai_OT = std(aiUnits_OT, [], 2)./ sqrt(size(aiUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_ai_OT = bsxfun(@plus, mean(aiUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_ai_OT));   % 95% Confidence Intervals

replayMean_OT = mean(replayUnits_OT,2);

SEM_replay_OT = std(replayUnits_OT, [], 2)./ sqrt(size(replayUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_replay_OT = bsxfun(@plus, mean(replayUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_replay_OT)); 

decoyMean_OT = mean(decoyUnits_OT,2);

SEM_decoy_OT = std(decoyUnits_OT, [], 2)./ sqrt(size(decoyUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_decoy_OT = bsxfun(@plus, mean(decoyUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_decoy_OT)); 

liveMean_OT = mean(liveUnits_OT,2);

SEM_live_OT = std(liveUnits_OT, [], 2)./ sqrt(size(liveUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_live_OT = bsxfun(@plus, mean(liveUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_live_OT)); 


data(:,:,1) = aiUnits_OT(:,1:83);
data(:,:,2) = replayUnits_OT(:,1:83);
data(:,:,3) = decoyUnits_OT(:,1:83);
data(:,:,4) = liveUnits_OT(:,1:83);



time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('div', 'RdBu', 4);

 
hold on;
bl = boundedline(time, aiMean_OT, SEM_ai_OT, ...
    time, replayMean_OT, SEM_replay_OT, ...
    time, decoyMean_OT, SEM_decoy_OT, ...
    time, liveMean_OT, SEM_live_OT, ...
    'cmap', colors);
% boundedline has an 'alpha' option, which makes the errorbars transparent
% (so it's nice when they overlap). However, when saving to pdf this makes
% the files HUGE, so better to keep your hands off alpha and make the final
% figure transparant in illustrator
 
xlim([-0.5 1.5]); xlabel('Time (s)'); ylabel('Normalized Firing Rate');

% Add lines
h2 = line([-0.5 -0.5],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color','k','LineWidth',0.75)
 
% instead of a legend, show colored text
lh = legend(bl);
legnames = {'AI (OT)','Replay (OT)','Decoy (OT)','Live (OT)'};
for i = 1:length(legnames),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(1) = lpos(1) + 0.15;
lh.Position = lpos;



 
% you'll still have the lines indicating the data. So far I haven't been
% able to find a good way to remove those, so you can either remove those
% in Illustrator, or use the text command to plot the legend (but then
% you'll have to specify the right x and y position for the text to go,
% which can take a bit of fiddling).
 
% we might want to add significance indicators, to show when the time
% courses are different from each other. In this case, use an uncorrected
% t-test
% for t = 1:length(time)
%     [~, pval1(t)] = ttest(data(:, t, 1), data(:, t, 2));
%     [~, pval2(t)] = ttest(data(:, t, 1), data(:, t, 3));
% end
% % convert to logical
% signific1 = nan(1, length(time)); signific1(pval1<0.001) = 1;
% signific2 = nan(1, length(time)); signific2(pval2<0.001) = 1;
% plot(time, signific1 * -3, '.k');
% plot(time, signific2 * -3, '.k');
% % indicate what we're showing
% text(10.2, -3, 'p < 0.001');
% 


% Add a patch
% gray = [0.4 0.4 0.4];
% patch([-0.5 -0.01 -0.01 -0.5],[-1.5 -1.5 2.5 2.5],gray)

hold off

hold on
% Add lines

h1 = line([0 0],[-1.5 2.5]);
h3 = line([1.495 1.495],[-1.5 2.5]);
set(h1,'Color',[0.4 0.4 0.4],'LineWidth',1)

set(h3,'Color','k','LineWidth',1)

hold off

