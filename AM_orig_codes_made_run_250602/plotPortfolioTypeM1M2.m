clear all
% plotStuff

% Plot psth from Nex data
% m1rew1on data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1rew1on.mat')
m1rew1onUnits = spikes;
m1rew1onNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1rew1on.mat')
m1rew1onUnits = [m1rew1onUnits, spikes(:,1:92)];
m1rew1onNames = [m1rew1onNames, units(:,1:92)];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1rew1on.mat')
m1rew1onUnits = [m1rew1onUnits, spikes];
m1rew1onNames = [m1rew1onNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1rew1on.mat')
m1rew1onUnits = [m1rew1onUnits, spikes];
m1rew1onNames = [m1rew1onNames, units];

% m2rew1on data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m2rew1on.mat')
m2rew1onUnits = spikes;
m2rew1onNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m2rew1on.mat')
m2rew1onUnits = [m2rew1onUnits, spikes(:,1:76)];
m2rew1onNames = [m2rew1onNames, units(:,1:76)];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m2rew1on.mat')
m2rew1onUnits = [m2rew1onUnits, spikes];
m2rew1onNames = [m2rew1onNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m2rew1on.mat')
m2rew1onUnits = [m2rew1onUnits, spikes];
m2rew1onNames = [m2rew1onNames, units];


clear units2 units3 spikes2


m1rew1onMean = mean(m1rew1onUnits,2);

SEM_m1rew1on = std(m1rew1onUnits, [], 2)./ sqrt(size(m1rew1onUnits,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on = bsxfun(@plus, mean(m1rew1onUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on));   % 95% Confidence Intervals

m2rew1onMean = mean(m2rew1onUnits,2);

SEM_m2rew1on = std(m2rew1onUnits, [], 2)./ sqrt(size(m2rew1onUnits,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on = bsxfun(@plus, mean(m2rew1onUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on)); 


data(:,:,1) = m1rew1onUnits(:,1:398);
data(:,:,2) = m2rew1onUnits;


time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set2', 8);
 
% subplot(4,4,[13 14]);  % plot across two subplots
hold on;
bl = boundedline(time, m1rew1onMean, SEM_m1rew1on, ...
    time, m2rew1onMean, SEM_m2rew1on, ...
    'cmap', colors);
% boundedline has an 'alpha' option, which makes the errorbars transparent
% (so it's nice when they overlap). However, when saving to pdf this makes
% the files HUGE, so better to keep your hands off alpha and make the final
% figure transparant in illustrator
 
xlim([-0.5 max(time)]); xlabel('Time (s)'); ylabel('Normalized Firing Rate');

% Add lines
h2 = line([-0.5 -0.5],[-1.5 2.5]);
h1 = line([0 0],[-1.5 2.5]);
h3 = line([0.5 0.5],[-1.5 2.5]);
% Set properties of lines
set(h1,'Color',[0.4 0.4 0.4],'LineWidth',1)
set(h3,'Color','k','LineWidth',1)
set(h2,'Color','k','LineWidth',0.75)

 
% instead of a legend, show colored text
lh = legend(bl);
legnames = {'M1 Trial', 'M2 Trial','Reward 1 Onset','Reward 2 Onset'};
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

