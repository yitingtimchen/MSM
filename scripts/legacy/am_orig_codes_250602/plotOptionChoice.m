clear all
% plotStuff

% Plot psth from Nex data
% buy data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1buy.mat')
buyUnits = spikes;
buyNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1buy.mat')
buyUnits = [buyUnits, spikes];
buyNames = [buyNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1buy.mat')
buyUnits = [buyUnits, spikes];
buyNames = [buyNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1buy.mat')
buyUnits = [buyUnits, spikes];
buyNames = [buyNames, units];

% hold data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1hold.mat')
holdUnits = spikes;
holdNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1hold.mat')
holdUnits = [holdUnits, spikes];
holdNames = [holdNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1hold.mat')
holdUnits = [holdUnits, spikes];
holdNames = [holdNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1hold.mat')
holdUnits = [holdUnits, spikes];
holdNames = [holdNames, units];

% sell data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1sell.mat')
sellUnits = spikes;
sellNames = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay
load('m1sell.mat')
sellUnits = [sellUnits, spikes];
sellNames = [sellNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy
load('m1sell.mat')
sellUnits = [sellUnits, spikes];
sellNames = [sellNames, units];
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1sell.mat')
sellUnits = [sellUnits, spikes];
sellNames = [sellNames, units];

clear units2 units3 spikes2

buyMean = mean(buyUnits,2);

SEM_buy = std(buyUnits, [], 2)./ sqrt(size(buyUnits,2));    % Calculate Standard Error Of The Mean
CI95_buy = bsxfun(@plus, mean(buyUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_buy));   % 95% Confidence Intervals

holdMean = mean(holdUnits,2);

SEM_hold = std(holdUnits, [], 2)./ sqrt(size(holdUnits,2));    % Calculate Standard Error Of The Mean
CI95_hold = bsxfun(@plus, mean(holdUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_hold)); 

sellMean = mean(sellUnits,2);

SEM_sell = std(sellUnits, [], 2)./ sqrt(size(sellUnits,2));    % Calculate Standard Error Of The Mean
CI95_sell = bsxfun(@plus, mean(sellUnits,2), bsxfun(@times, [-1  1]*1.96, SEM_sell)); 


data(:,:,1) = buyUnits;
data(:,:,2) = holdUnits;
data(:,:,3) = sellUnits;




time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 8);
temp = colors(1,:);
colors(1,:) = colors(3,:);
colors(3,:) = temp;

 
bl = boundedline(time, buyMean, SEM_buy, ...
    time, holdMean, SEM_hold, ...
    time, sellMean, SEM_sell, ...
    'cmap', colors,'transparency',0.5);

hold on
rect2 = fill([-0.25 -0.05 -0.05 -0.25 ],[-2 -2 2.5 2.5],'yellow','LineStyle','none'); 
rect2.FaceAlpha=0.2;
uistack(rect2,'bottom');
set(gca,'FontSize',16)
 
xlim([-0.5 1]); xlabel('Time (s)'); ylim([-2 2.5]); ylabel('Normalized Firing Rate');


% Add lines
h2 = line([-0.5 -0.5],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color','k','LineWidth',0.75)
 
% instead of a legend, show colored text
lh = legend(bl);
legnames = {'Buy)','Hold','Sell'};
for i = 1:length(legnames),
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
print('optionChoice_allConditions','-dpdf');


