clear all
% plotStuff

% Plot psth from Nex data

%----------------------- Saline --------------------------------------

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\AI_Saline
load('m1rew1on.mat')
aiUnits_Saline = spikes;
aiNames_Saline = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay_Saline
load('m1rew1on.mat')
replayUnits_Saline = spikes;
replayNames_Saline = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy_Saline
load('m1rew1on.mat')
decoyUnits_Saline = spikes;
decoyNames_Saline = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Live_Saline
load('m1rew1on.mat')
liveUnits_Saline = spikes;
liveNames_Saline = units;

clear units2 units3 spikes2

aiMean_Saline = mean(aiUnits_Saline,2);

SEM_ai_Saline = std(aiUnits_Saline, [], 2)./ sqrt(size(aiUnits_Saline,2));    % Calculate Standard Error Of The Mean
CI95_ai_Saline = bsxfun(@plus, mean(aiUnits_Saline,2), bsxfun(@times, [-1  1]*1.96, SEM_ai_Saline));   % 95% Confidence Intervals

replayMean_Saline = mean(replayUnits_Saline,2);

SEM_replay_Saline = std(replayUnits_Saline, [], 2)./ sqrt(size(replayUnits_Saline,2));    % Calculate Standard Error Of The Mean
CI95_replay_Saline = bsxfun(@plus, mean(replayUnits_Saline,2), bsxfun(@times, [-1  1]*1.96, SEM_replay_Saline)); 

decoyMean_Saline = mean(decoyUnits_Saline,2);

SEM_decoy_Saline = std(decoyUnits_Saline, [], 2)./ sqrt(size(decoyUnits_Saline,2));    % Calculate Standard Error Of The Mean
CI95_decoy_Saline = bsxfun(@plus, mean(decoyUnits_Saline,2), bsxfun(@times, [-1  1]*1.96, SEM_decoy_Saline)); 

liveMean_Saline = mean(liveUnits_Saline,2);

SEM_live_Saline = std(liveUnits_Saline, [], 2)./ sqrt(size(liveUnits_Saline,2));    % Calculate Standard Error Of The Mean
CI95_live_Saline = bsxfun(@plus, mean(liveUnits_Saline,2), bsxfun(@times, [-1  1]*1.96, SEM_live_Saline)); 


data(:,:,1) = aiUnits_Saline(:,1:57);
data(:,:,2) = replayUnits_Saline(:,1:57);
data(:,:,3) = decoyUnits_Saline(:,1:57);
data(:,:,4) = liveUnits_Saline(:,1:57);


time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('div', 'RdYlBu', 4);


hold on;
bl = boundedline(time, aiMean_Saline, SEM_ai_Saline, ...
    time, replayMean_Saline, SEM_replay_Saline, ...
    time, decoyMean_Saline, SEM_decoy_Saline, ...
    time, liveMean_Saline, SEM_live_Saline, ...
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
legnames = {'AI (Saline)','Replay (Saline)','Decoy (Saline)','Live (Saline)'};
for i = 1:length(legnames)
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

