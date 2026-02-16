clear all
% plotStuff

% Plot psth from Nex data
% m1rew1on data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1preP0.mat')
m1preP0Units = spikes;
m1preP0Names = units;
load('m1preP5.mat')
m1preP5Units = spikes;
m1preP5Names = units;
load('m1preP10.mat')
m1preP10Units = spikes;
m1preP10Names = units;
load('m1preP15.mat')
m1preP15Units = spikes;
m1preP15Names = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1preP0.mat')
m1preP0Units = [m1preP0Units, spikes(:,1:92)];
m1preP0Names = [m1preP0Names, units(:,1:92)];
load('m1preP5.mat')
m1preP5Units = [m1preP5Units, spikes(:,1:92)];
m1preP5Names = [m1preP5Names, units(:,1:92)];
load('m1preP10.mat')
m1preP10Units = [m1preP10Units, spikes(:,1:92)];
m1preP10Names = [m1preP10Names, units(:,1:92)];
load('m1preP15.mat')
m1preP15Units = [m1preP15Units, spikes(:,1:92)];
m1preP15Names = [m1preP15Names, units(:,1:92)];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1preP0.mat')
m1preP0Units = [m1preP0Units, spikes];
m1preP0Names = [m1preP0Names, units];
load('m1preP5.mat')
m1preP5Units = [m1preP5Units, spikes];
m1preP5Names = [m1preP5Names, units];
load('m1preP10.mat')
m1preP10Units = [m1preP10Units, spikes];
m1preP10Names = [m1preP10Names, units];
load('m1preP15.mat')
m1preP15Units = [m1preP15Units, spikes];
m1preP15Names = [m1preP15Names, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1preP0.mat')
m1preP0Units = [m1preP0Units, spikes];
m1preP0Names = [m1preP0Names, units];
load('m1preP5.mat')
m1preP5Units = [m1preP5Units, spikes];
m1preP5Names = [m1preP5Names, units];
load('m1preP10.mat')
m1preP10Units = [m1preP10Units, spikes];
m1preP10Names = [m1preP10Names, units];
load('m1preP15.mat')
m1preP15Units = [m1preP15Units, spikes];
m1preP15Names = [m1preP15Names, units];



clear units2 units3 spikes2 data


m1preP0Mean = mean(m1preP0Units,2);

SEM_m1preP0 = std(m1preP0Units, [], 2)./ sqrt(size(m1preP0Units,2));    % Calculate Standard Error Of The Mean
CI95_m1preP0 = bsxfun(@plus, mean(m1preP0Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m1preP0));   % 95% Confidence Intervals

m1preP5Mean = mean(m1preP5Units,2);

SEM_m1preP5 = std(m1preP5Units, [], 2)./ sqrt(size(m1preP5Units,2));    % Calculate Standard Error Of The Mean
CI95_m1preP5 = bsxfun(@plus, mean(m1preP5Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m1preP5)); 

m1preP10Mean = mean(m1preP10Units,2);

SEM_m1preP10 = std(m1preP10Units, [], 2)./ sqrt(size(m1preP10Units,2));    % Calculate Standard Error Of The Mean
CI95_m1preP10 = bsxfun(@plus, mean(m1preP10Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m1preP10)); 

m1preP15Mean = mean(m1preP15Units,2);

SEM_m1preP15 = std(m1preP15Units, [], 2)./ sqrt(size(m1preP15Units,2));    % Calculate Standard Error Of The Mean
CI95_m1preP15 = bsxfun(@plus, mean(m1preP15Units,2), bsxfun(@times, [-1  1]*1.96, SEM_m1preP15)); 

% data(:,:,1) = m1preP0Units;
% data(:,:,2) = m1preP5Units;
% data(:,:,3) = m1preP10Units;
% data(:,:,4) = m1preP15Units;


time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('seq', 'YlGn', 5);
colors = colors(2:5,:);
 
% subplot(4,4,[13 14]);  % plot across two subplots
hold on;
bl = boundedline(time, m1preP0Mean, SEM_m1preP0, ...
    time, m1preP5Mean, SEM_m1preP5, ...
    time, m1preP10Mean, SEM_m1preP10, ...
    time, m1preP15Mean, SEM_m1preP15, ...
    'cmap', colors);
% boundedline has an 'alpha' option, which makes the errorbars transparent
% (so it's nice when they overlap). However, when saving to pdf this makes
% the files HUGE, so better to keep your hands off alpha and make the final
% figure transparant in illustrator
 
xlim([-0.5 max(time)]); xlabel('Time (s)'); ylabel('Normalized Firing Rate');

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

