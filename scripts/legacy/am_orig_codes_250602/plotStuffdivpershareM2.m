clear all
% plotStuff

% Plot psth from Nex data
% m2rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m2div00.mat')
m2div00Units = spikes;
m2div00Names = units;
load('m2div08.mat')
m2div08Units = spikes;
m2div08Names = units;
load('m2div28.mat')
m2div28Units = spikes;
m2div28Names = units;
load('m2div60.mat')
m2div60Units = spikes;
m2div60Names = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m2div00.mat')
m2div00Units = [m2div00Units, spikes(:,1:92)];
m2div00Names = [m2div00Names, units(:,1:92)];
load('m2div08.mat')
m2div08Units = [m2div08Units, spikes(:,1:92)];
m2div08Names = [m2div08Names, units(:,1:92)];
load('m2div28.mat')
m2div28Units = [m2div28Units, spikes(:,1:92)];
m2div28Names = [m2div28Names, units(:,1:92)];
load('m2div60.mat')
m2div60Units = [m2div60Units, spikes(:,1:92)];
m2div60Names = [m2div60Names, units(:,1:92)];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m2div00.mat')
m2div00Units = [m2div00Units, spikes];
m2div00Names = [m2div00Names, units];
load('m2div08.mat')
m2div08Units = [m2div08Units, spikes];
m2div08Names = [m2div08Names, units];
load('m2div28.mat')
m2div28Units = [m2div28Units, spikes];
m2div28Names = [m2div28Names, units];
load('m2div60.mat')
m2div60Units = [m2div60Units, spikes];
m2div60Names = [m2div60Names, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m2div00.mat')
m2div00Units = [m2div00Units, spikes];
m2div00Names = [m2div00Names, units];
load('m2div08.mat')
m2div08Units = [m2div08Units, spikes];
m2div08Names = [m2div08Names, units];
load('m2div28.mat')
m2div28Units = [m2div28Units, spikes];
m2div28Names = [m2div28Names, units];
load('m2div60.mat')
m2div60Units = [m2div60Units, spikes];
m2div60Names = [m2div60Names, units];



clear units2 units3 spikes2


m2div00Mean = mean(m2div00Units,2);

SEM_m2div00 = std(m2div00Units, [], 2)./ sqrt(size(m2div00Units,2));    % Calculate Standard Error Of The Mean
CI95_m2div00 = bsxfun(@plus, mean(m2div00Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m2div00));   % 95% Confidence Intervals

m2div08Mean = mean(m2div08Units,2);

SEM_m2div08 = std(m2div08Units, [], 2)./ sqrt(size(m2div08Units,2));    % Calculate Standard Error Of The Mean
CI95_m2div08 = bsxfun(@plus, mean(m2div08Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m2div08)); 

m2div28Mean = mean(m2div28Units,2);

SEM_m2div28 = std(m2div28Units, [], 2)./ sqrt(size(m2div28Units,2));    % Calculate Standard Error Of The Mean
CI95_m2div28 = bsxfun(@plus, mean(m2div28Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m2div28)); 

m2div60Mean = mean(m2div60Units,2);

SEM_m2div60 = std(m2div60Units, [], 2)./ sqrt(size(m2div60Units,2));    % Calculate Standard Error Of The Mean
CI95_m2div60 = bsxfun(@plus, mean(m2div60Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m2div60)); 

data(:,:,1) = m2div00Units;
data(:,:,2) = m2div08Units;
data(:,:,3) = m2div28Units;
data(:,:,4) = m2div60Units;


time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('seq', 'YlGn', 5);
colors = colors(2:5,:);
 
% subplot(4,4,[13 14]);  % plot across two subplots
hold on;
bl = boundedline(time, m2div00Mean, SEM_m2div00, ...
    time, m2div08Mean, SEM_m2div08, ...
    time, m2div28Mean, SEM_m2div28, ...
    time, m2div60Mean, SEM_m2div60, ...
    'cmap', colors);
% boundedline has an 'alpha' option, which makes the errorbars transparent
% (so it's nice when they overlap). However, when saving to pdf this makes
% the files HUGE, so better to keep your hands off alpha and make the final
% figure transparant in illustrator
 
xlim([-0.5 max(time)]); xlabel('Time (s)'); ylabel('Normalized Firing Rate');
ylim([-1.5 1.0]);
% Add lines
h2 = line([-0.5 -0.5],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color','k','LineWidth',0.75)
 
% instead of a legend, show colored text
lh = legend(bl);
legnames = {'0.0', '0.8','2.8','6.0'};
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

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData
print('divPerSharePanelbyConditionFM','-dpdf');

