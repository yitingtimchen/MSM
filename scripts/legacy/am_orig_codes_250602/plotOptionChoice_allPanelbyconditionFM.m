clear all
% plotStuff

% Plot psth from Nex data
% buy data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1buy.mat')
buyUnits_ai = spikes;
buyNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1buy.mat')
buyUnits_replay = spikes;
buyNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1buy.mat')
buyUnits_decoy = spikes;
buyNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1buy.mat')
buyUnits_live = spikes;
buyNames_live = units;

% hold data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1hold.mat')
holdUnits_ai = spikes;
holdNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1hold.mat')
holdUnits_replay = spikes;
holdNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1hold.mat')
holdUnits_decoy = spikes;
holdNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1hold.mat')
holdUnits_live = spikes;
holdNames_live = units;

% sell data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\AI
load('m1sell.mat')
sellUnits_ai = spikes;
sellNames_ai = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Replay
load('m1sell.mat')
sellUnits_replay = spikes;
sellNames_replay = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Decoy
load('m1sell.mat')
sellUnits_decoy = spikes;
sellNames_decoy = units;

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\FreeMarket\Live
load('m1sell.mat')
sellUnits_live = spikes;
sellNames_live = units;

clear units2 units3 spikes2

% Cut units so that each condition has the same number 
% [m,n] = size(sellUnits_replay);
% n = n-50;
% 
% 
%     buyUnits_ai = buyUnits_ai(:,end-n:end-9);
%     buyUnits_replay = buyUnits_replay(:,end-n:end-9);
%     buyUnits_decoy = buyUnits_decoy(:,end-n:end-9);
%     buyUnits_live = buyUnits_live(:,end-n:end-9);
%     holdUnits_ai = holdUnits_ai(:,end-n:end-9);
%     holdUnits_replay = holdUnits_replay(:,end-n:end-9);
%     holdUnits_decoy = holdUnits_decoy(:,end-n:end-9);
%     holdUnits_live = holdUnits_live(:,end-n:end-9);
%     sellUnits_ai = sellUnits_ai(:,end-n:end-9);
%     sellUnits_replay = sellUnits_replay(:,end-n:end-9);
%     sellUnits_decoy = sellUnits_decoy(:,end-n:end-9);
%     sellUnits_live = sellUnits_live(:,end-n:end-9);


% Means
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
buyMean_live = mean(buyUnits_live,2);

SEM_buy_live = std(buyUnits_live, [], 2)./ sqrt(size(buyUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_buy_live = bsxfun(@plus, mean(buyUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_live));   % 95% Confidence Intervals

holdMean_live = mean(holdUnits_live,2);

SEM_hold_live = std(holdUnits_live, [], 2)./ sqrt(size(holdUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_hold_live = bsxfun(@plus, mean(holdUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_live)); 

sellMean_live = mean(sellUnits_live,2);

SEM_sell_live = std(sellUnits_live, [], 2)./ sqrt(size(sellUnits_live,2));    % Calculate Standard Error Of The Mean
CI95_sell_live = bsxfun(@plus, mean(sellUnits_live,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_live)); 


% data(:,:,1) = buyUnits_ai;
% data(:,:,2) = holdUnits_ai;
% data(:,:,3) = sellUnits_ai;
% data(:,:,4) = buyUnits_replay;
% data(:,:,5) = holdUnits_replay;
% data(:,:,6) = sellUnits_replay;
% data(:,:,7) = buyUnits_decoy;
% data(:,:,8) = holdUnits_decoy;
% data(:,:,9) = sellUnits_decoy;
% data(:,:,10) = buyUnits_live;
% data(:,:,11) = holdUnits_live;
% data(:,:,12) = sellUnits_live;
% 


time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 8);
temp = colors(1,:);
colors(1,:) = colors(3,:);
colors(3,:) = temp;

nl = @(s) strrep(s,'\n',char(10));

subplot(1,4,1);
bl1 = boundedline(time, buyMean_ai, SEM_buy_ai, ...
    time, holdMean_ai, SEM_hold_ai, ...
    time, sellMean_ai, SEM_sell_ai, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('AI');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 2]); ylabel(nl('Normalized Firing Rate'));
hold off

subplot(1,4,2);
bl2 = boundedline(time, buyMean_replay, SEM_buy_replay, ...
    time, holdMean_replay, SEM_hold_replay, ...
    time, sellMean_replay, SEM_sell_replay, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h4 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h4,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Replay');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 2]);
hold off

subplot(1,4,3);
bl3 = boundedline(time, buyMean_decoy, SEM_buy_decoy, ...
    time, holdMean_decoy, SEM_hold_decoy, ...
    time, sellMean_decoy, SEM_sell_decoy, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Decoy');
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 2]);
hold off

subplot(1,4,4);
bl4 = boundedline(time, buyMean_live, SEM_buy_live, ...
    time, holdMean_live, SEM_hold_live, ...
    time, sellMean_live, SEM_sell_live, ...
    'cmap', colors,'transparency',0.5);
hold on
rect1 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect1.FaceAlpha=0.2;
uistack(rect1,'bottom');
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
title('Live');

    
 
xlim([-0.5 1.5]); xlabel('Time (s)'); ylim([-1.5 2]);

 
% instead of a legend, show colored text
% lh = legend(bl1);
% legnames = {'Buy','Hold','Sell'};
% for i = 1:length(legnames),
%     str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
% end
% lh.String = str;
% lh.Box = 'off';
 
% move a bit closer
% move a bit closer
% lpos = lh.Position;
% lpos(1) = lpos(1) + 0.65;
% lpos(2) = lpos(2) - 0.5;
% lh.Position = lpos;
% lh.FontSize = 16;
% uistack(lh,'top');

hold off
fig = gcf;
orient(fig,'landscape');
fig.PaperType = 'uslegal';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 14 4.5];

cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData
print('optionPanelbyConditionFM','-dpdf');

