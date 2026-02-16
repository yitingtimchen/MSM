clear all
% plotStuff

% Plot psth from Nex data
% m1rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1rew1on.mat')
m1rew1onUnits_ai = spikes;
m1rew1onNames_ai = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1rew1on.mat')
m1rew1onUnits_replay = spikes;
m1rew1onNames_replay = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1rew1on.mat')
m1rew1onUnits_decoy = spikes;
m1rew1onNames_decoy = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1rew1on.mat')
m1rew1onUnits_live = spikes;
m1rew1onNames_live = units;

% m2rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m2rew1on.mat')
m2rew1onUnits_ai = spikes;
m2rew1onNames_ai = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m2rew1on.mat')
m2rew1onUnits_replay = spikes;
m2rew1onNames_replay = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m2rew1on.mat')
m2rew1onUnits_decoy = spikes;
m2rew1onNames_decoy = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m2rew1on.mat')
m2rew1onUnits_live = spikes;
m2rew1onNames_live = units;


clear units2 units3 spikes2

% M1

% ai
m1rew1onMean_ai = mean(m1rew1onUnits_ai,2);

SEM_m1rew1on_ai = std(m1rew1onUnits_ai, [], 2)./ sqrt(size(m1rew1onUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_ai = bsxfun(@plus, mean(m1rew1onUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_ai));   % 95% Confidence Intervals
% replay
m1rew1onMean_replay = mean(m1rew1onUnits_replay,2);

SEM_m1rew1on_replay = std(m1rew1onUnits_replay, [], 2)./ sqrt(size(m1rew1onUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_replay = bsxfun(@plus, mean(m1rew1onUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_replay));   % 95% Confidence Intervals
% decoy
m1rew1onMean_decoy = mean(m1rew1onUnits_decoy,2);

SEM_m1rew1on_decoy = std(m1rew1onUnits_decoy, [], 2)./ sqrt(size(m1rew1onUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_decoy = bsxfun(@plus, mean(m1rew1onUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_decoy));   % 95% Confidence Intervals
% live
m1rew1onMean_live = mean(m1rew1onUnits_live,2);

SEM_m1rew1on_live = std(m1rew1onUnits_live, [], 2)./ sqrt(size(m1rew1onUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_m1rew1on_live = bsxfun(@plus, mean(m1rew1onUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_m1rew1on_live));   % 95% Confidence Intervals



% M2

% ai
m2rew1onMean_ai = mean(m2rew1onUnits_ai,2);

SEM_m2rew1on_ai = std(m2rew1onUnits_ai, [], 2)./ sqrt(size(m2rew1onUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_ai = bsxfun(@plus, mean(m2rew1onUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_ai));   % 95% Confidence Intervals
% replay
m2rew1onMean_replay = mean(m2rew1onUnits_replay,2);

SEM_m2rew1on_replay = std(m2rew1onUnits_replay, [], 2)./ sqrt(size(m2rew1onUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_replay = bsxfun(@plus, mean(m2rew1onUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_replay));   % 95% Confidence Intervals
% decoy
m2rew1onMean_decoy = mean(m2rew1onUnits_decoy,2);

SEM_m2rew1on_decoy = std(m2rew1onUnits_decoy, [], 2)./ sqrt(size(m2rew1onUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_decoy = bsxfun(@plus, mean(m2rew1onUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_decoy));   % 95% Confidence Intervals
% live
m2rew1onMean_live = mean(m2rew1onUnits_live,2);

SEM_m2rew1on_live = std(m2rew1onUnits_live, [], 2)./ sqrt(size(m2rew1onUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_m2rew1on_live = bsxfun(@plus, mean(m2rew1onUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_m2rew1on_live));   % 95% Confidence Intervals


data(:,:,1) = m1rew1onUnits_ai(:,1:76);
data(:,:,2) = m1rew1onUnits_replay(:,1:76);
data(:,:,3) = m1rew1onUnits_decoy(:,1:76);
data(:,:,4) = m1rew1onUnits_live(:,1:76);
data(:,:,5) = m2rew1onUnits_ai(:,1:76);
data(:,:,6) = m2rew1onUnits_replay(:,1:76);
data(:,:,7) = m2rew1onUnits_decoy(:,1:76);
data(:,:,8) = m2rew1onUnits_live(:,1:76);


time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set2', 8);
 
% subplot(4,4,[13 14]);  % plot across two subplots
hold on;
bl = boundedline(time, m1rew1onMean_ai, SEM_m1rew1on_ai, ...
    time, m1rew1onMean_replay, SEM_m1rew1on_replay, ...
    time, m1rew1onMean_decoy, SEM_m1rew1on_decoy, ...
    time, m1rew1onMean_live, SEM_m1rew1on_live, ...
    time, m2rew1onMean_ai, SEM_m2rew1on_ai, ...
    time, m2rew1onMean_replay, SEM_m2rew1on_replay, ...
    time, m2rew1onMean_decoy, SEM_m2rew1on_decoy, ...
    time, m2rew1onMean_live, SEM_m2rew1on_live, ...
    'cmap', colors);
% boundedline has an 'alpha' option, which makes the errorbars transparent
% (so it's nice when they overlap). However, when saving to pdf this makes
% the files HUGE, so better to keep your hands off alpha and make the final
% figure transparant in illustrator
 
 
xlim([-0.5 1.5]); xlabel('Time (s)'); ylabel('Normalized Firing Rate');

% Add lines
h2 = line([-0.5 -0.5],[-1.5 1.0]);
% Set properties of lines
set(h2,'Color','k','LineWidth',0.75)

 
% instead of a legend, show colored text
lh = legend(bl);
legnames = {'M1 AI', 'M1 Replay','M1 Decoy','M1 Live', 'M2 AI', 'M2 Replay', 'M2 Decoy', 'M2 Live'};
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

h1 = line([0 0],[-1.5 1.0]);
h3 = line([1.495 1.495],[-1.5 1.0]);
set(h1,'Color',[0.4 0.4 0.4],'LineWidth',1)
set(h3,'Color','k','LineWidth',1)

hold off

