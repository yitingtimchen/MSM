clear all
% plotStuff

% Plot psth from Nex data
% larger Portfolio data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1largerPrew1.mat')
largerPUnits = spikes;
largerPNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1largerPrew1.mat')
largerPUnits = [largerPUnits, spikes];
largerPNames = [largerPNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1largerPrew1.mat')
largerPUnits = [largerPUnits, spikes];
largerPNames = [largerPNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1largerPrew1.mat')
largerPUnits = [largerPUnits, spikes];
largerPNames = [largerPNames, units];

% evenP data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1evenPrew1.mat')
evenPUnits = spikes;
evenPNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1evenPrew1.mat')
evenPUnits = [evenPUnits, spikes];
evenPNames = [evenPNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1evenPrew1.mat')
evenPUnits = [evenPUnits, spikes];
evenPNames = [evenPNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1evenPrew1.mat')
evenPUnits = [evenPUnits, spikes];
evenPNames = [evenPNames, units];

% smallerP data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1smallerPrew1.mat')
smallerPUnits = spikes;
smallerPNames = units;
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1smallerPrew1.mat')
smallerPUnits = [smallerPUnits, spikes];
smallerPNames = [smallerPNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1smallerPrew1.mat')
smallerPUnits = [smallerPUnits, spikes];
smallerPNames = [smallerPNames, units];
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1smallerPrew1.mat')
smallerPUnits = [smallerPUnits, spikes];
smallerPNames = [smallerPNames, units];

clear units2 units3 spikes2

largerPMean = mean(largerPUnits,2);

SEM_largerP = std(largerPUnits, [], 2)./ sqrt(size(largerPUnits,2));    % Calculate Standard Error Of The Mean
CI95_largerP = bsxfun(@plus, mean(largerPUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_largerP));   % 95% Confidence Intervals

evenPMean = mean(evenPUnits,2);

SEM_evenP = std(evenPUnits, [], 2)./ sqrt(size(evenPUnits,2));    % Calculate Standard Error Of The Mean
CI95_evenP = bsxfun(@plus, mean(evenPUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_evenP)); 

smallerPMean = mean(smallerPUnits,2);

SEM_smallerP = std(smallerPUnits, [], 2)./ sqrt(size(smallerPUnits,2));    % Calculate Standard Error Of The Mean
CI95_smallerP = bsxfun(@plus, mean(smallerPUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_smallerP)); 




data(:,:,1) = largerPUnits;
data(:,:,2) = evenPUnits;
data(:,:,3) = smallerPUnits;

time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Dark2', 8);
 
% subplot(4,4,[13 14]);  % plot across two subplots
hold on;
bl = boundedline(time, largerPMean, SEM_largerP, ...
    time, evenPMean, SEM_evenP, ...
    time, smallerPMean, SEM_smallerP,...
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
legnames = {'M1 > M2', 'M1 = M2', 'M1 < M2'};
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

