clear all
% plotStuff

% Plot psth from Nex data
% buy data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1buy.mat')
buyUnits_ai = spikes;
buyNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1buy.mat')
buyUnits_live = spikes;
buyNames_live = units;

% hold data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1hold.mat')
holdUnits_ai = spikes;
holdNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1hold.mat')
holdUnits_live = spikes;
holdNames_live = units;

% sell data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\AI
load('m1sell.mat')
sellUnits_ai = spikes;
sellNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Live
load('m1sell.mat')
sellUnits_live = spikes;
sellNames_live = units;

% Tim 250603, Fill in hold and sell for replay and decoy
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay\
load('m1buy.mat')
buyUnits_replay = spikes;
buyNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy\
load('m1buy.mat')
buyUnits_decoy = spikes;
buyNames_decoy = units;

% Replay hold data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay\
load('m1hold.mat')
holdUnits_replay = spikes;
holdNames_replay = units;

% Decoy hold data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy\
load('m1hold.mat')
holdUnits_decoy = spikes;
holdNames_decoy = units;

% Replay sell data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Replay\
load('m1sell.mat')
sellUnits_replay = spikes;
sellNames_replay = units;

% Decoy sell data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\Decoy\
load('m1sell.mat')
sellUnits_decoy = spikes;
sellNames_decoy = units;

clear units2 units3 spikes2

buyMean_ai = mean(buyUnits_ai,2);

SEM_buy_ai = std(buyUnits_ai, [], 2)./ sqrt(size(buyUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_buy_ai = bsxfun(@plus, mean(buyUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_ai));   % 95% Confidence Intervals

holdMean_ai = mean(holdUnits_ai,2);

SEM_hold_ai = std(holdUnits_ai, [], 2)./ sqrt(size(holdUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_hold_ai = bsxfun(@plus, mean(holdUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_ai)); 

sellMean_ai = mean(sellUnits_ai,2);

SEM_sell_ai = std(sellUnits_ai, [], 2)./ sqrt(size(sellUnits_ai,2));    % Calculate Standard Error Of The Mean
CI95_sell_ai = bsxfun(@plus, mean(sellUnits_ai,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_ai)); 

% replay
buyMean_replay = mean(buyUnits_replay,2);

SEM_buy_replay = std(buyUnits_replay, [], 2)./ sqrt(size(buyUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_buy_replay = bsxfun(@plus, mean(buyUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_replay));   % 95% Confidence Intervals

holdMean_replay = mean(holdUnits_replay,2);

SEM_hold_replay = std(holdUnits_replay, [], 2)./ sqrt(size(holdUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_hold_replay = bsxfun(@plus, mean(holdUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_replay)); 

sellMean_replay = mean(sellUnits_replay,2);

SEM_sell_replay = std(sellUnits_replay, [], 2)./ sqrt(size(sellUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_sell_replay = bsxfun(@plus, mean(sellUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_replay)); 

% decoy
buyMean_decoy = mean(buyUnits_decoy,2);

SEM_buy_decoy = std(buyUnits_decoy, [], 2)./ sqrt(size(buyUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_buy_decoy = bsxfun(@plus, mean(buyUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_decoy));   % 95% Confidence Intervals

holdMean_decoy = mean(holdUnits_decoy,2);

SEM_hold_decoy = std(holdUnits_decoy, [], 2)./ sqrt(size(holdUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_hold_decoy = bsxfun(@plus, mean(holdUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_decoy)); 

sellMean_decoy = mean(sellUnits_decoy,2);

SEM_sell_decoy = std(sellUnits_decoy, [], 2)./ sqrt(size(sellUnits_decoy,2));    % Calculate Standard Error Of The Mean
CI95_sell_decoy = bsxfun(@plus, mean(sellUnits_decoy,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_decoy)); 

% live
buyMean_replay = mean(buyUnits_replay,2);

SEM_buy_replay = std(buyUnits_replay, [], 2)./ sqrt(size(buyUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_buy_replay = bsxfun(@plus, mean(buyUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_replay));   % 95% Confidence Intervals

holdMean_replay = mean(holdUnits_replay,2);

SEM_hold_replay = std(holdUnits_replay, [], 2)./ sqrt(size(holdUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_hold_replay = bsxfun(@plus, mean(holdUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_replay)); 

sellMean_replay = mean(sellUnits_replay,2);

SEM_sell_replay = std(sellUnits_replay, [], 2)./ sqrt(size(sellUnits_replay,2));    % Calculate Standard Error Of The Mean
CI95_sell_replay = bsxfun(@plus, mean(sellUnits_replay,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_replay)); 

% Tim 250603, commented out, seems unused
% data(:,:,1) = buyUnits_ai;
% data(:,:,2) = holdUnits_ai;
% data(:,:,3) = sellUnits_ai;
% data(:,:,4) = buyUnits_replay;
% data(:,:,5) = holdUnits_replay;
% data(:,:,6) = sellUnits_replay;




time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 8);
temp = colors(1,:);
colors(1,:) = colors(3,:);
colors(3,:) = temp;

 
hold on;
bl = boundedline(time, buyMean, SEM_buy, ...
    time, holdMean, SEM_hold, ...
    time, sellMean, SEM_sell, ...
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
legnames = {'Buy)','Hold','Sell'};
for i = 1:length(legnames),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(1) = lpos(1) + 0.15;
lh.Position = lpos;

hold off

hold on
% Add lines

h1 = line([0 0],[-1.5 2.5]);
h3 = line([1.495 1.495],[-1.5 2.5]);
set(h1,'Color',[0.4 0.4 0.4],'LineWidth',1)

set(h3,'Color','k','LineWidth',1)

hold off   

