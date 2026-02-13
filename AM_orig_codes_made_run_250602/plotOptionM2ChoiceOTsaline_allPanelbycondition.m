clear all
% plotStuff

% Plot psth from Nex data
% Saline
% buy data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m2buy.mat')
buyUnits_aiS = spikes;
buyNames_aiS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m2buy.mat')
buyUnits_replayS = spikes;
buyNames_replayS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m2buy.mat')
buyUnits_decoyS = spikes;
buyNames_decoyS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m2buy.mat')
buyUnits_liveS = spikes;
buyNames_liveS = units;

% hold data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m2hold.mat')
holdUnits_aiS = spikes;
holdNames_aiS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m2hold.mat')
holdUnits_replayS = spikes;
holdNames_replayS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m2hold.mat')
holdUnits_decoyS = spikes;
holdNames_decoyS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m2hold.mat')
holdUnits_liveS = spikes;
holdNames_liveS = units;

% sell data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\AI
load('m2sell.mat')
sellUnits_aiS = spikes;
sellNames_aiS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Replay
load('m2sell.mat')
sellUnits_replayS = spikes;
sellNames_replayS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Decoy
load('m2sell.mat')
sellUnits_decoyS = spikes;
sellNames_decoyS = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\Saline\Live
load('m2sell.mat')
sellUnits_liveS = spikes;
sellNames_liveS = units;

% OT
% buy data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m2buy.mat')
buyUnits_aiOT = spikes;
buyNames_aiOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m2buy.mat')
buyUnits_replayOT = spikes;
buyNames_replayOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m2buy.mat')
buyUnits_decoyOT = spikes;
buyNames_decoyOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m2buy.mat')
buyUnits_liveOT = spikes;
buyNames_liveOT = units;

% hold data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m2hold.mat')
holdUnits_aiOT = spikes;
holdNames_aiOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m2hold.mat')
holdUnits_replayOT = spikes;
holdNames_replayOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m2hold.mat')
holdUnits_decoyOT = spikes;
holdNames_decoyOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m2hold.mat')
holdUnits_liveOT = spikes;
holdNames_liveOT = units;

% sell data
cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\AI
load('m2sell.mat')
sellUnits_aiOT = spikes;
sellNames_aiOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Replay
load('m2sell.mat')
sellUnits_replayOT = spikes;
sellNames_replayOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy
load('m2sell.mat')
sellUnits_decoyOT = spikes;
sellNames_decoyOT = units;

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData\OT\Live
load('m2sell.mat')
sellUnits_liveOT = spikes;
sellNames_liveOT = units;

clear units2 units3 spikes2

% Cut units so that each condition has the same number 
[m,n] = size(buyUnits_aiS);
n = n-1;


    buyUnits_aiS = buyUnits_aiS(:,end-n:end);
    holdUnits_aiS = holdUnits_aiS(:,end-n:end);
    sellUnits_aiS = sellUnits_aiS(:,end-n:end);
    buyUnits_replayS = buyUnits_replayS(:,end-n:end);
    holdUnits_replayS = holdUnits_replayS(:,end-n:end);
    sellUnits_replayS = sellUnits_replayS(:,end-n:end);
    buyUnits_decoyS = buyUnits_decoyS(:,end-n:end);
    holdUnits_decoyS = holdUnits_decoyS(:,end-n:end);
    sellUnits_decoyS = sellUnits_decoyS(:,end-n:end);
    buyUnits_liveS = buyUnits_liveS(:,end-n:end);
    holdUnits_liveS = holdUnits_liveS(:,end-n:end);
    sellUnits_liveS = sellUnits_liveS(:,end-n:end);
    buyUnits_aiOT = buyUnits_aiOT(:,end-n:end);
    holdUnits_aiOT = holdUnits_aiOT(:,end-n:end);
    sellUnits_aiOT = sellUnits_aiOT(:,end-n:end);
    buyUnits_replayOT = buyUnits_replayOT(:,end-n:end);
    holdUnits_replayOT = holdUnits_replayOT(:,end-n:end);
    sellUnits_replayOT = sellUnits_replayOT(:,end-n:end);
    buyUnits_decoyOT = buyUnits_decoyOT(:,end-n:end);
    holdUnits_decoyOT = holdUnits_decoyOT(:,end-n:end);
    sellUnits_decoyOT = sellUnits_decoyOT(:,end-n:end);
    buyUnits_liveOT = [buyUnits_liveOT(:,32:93) buyUnits_liveOT(:,end-45:end)];
    holdUnits_liveOT = [holdUnits_liveOT(:,32:93) holdUnits_liveOT(:,end-45:end)];
    sellUnits_liveOT = [sellUnits_liveOT(:,32:93) sellUnits_liveOT(:,end-45:end)];

% Saline means
buyMean_aiS = mean(buyUnits_aiS,2);

SEM_buy_aiS = std(buyUnits_aiS, [], 2)./ sqrt(size(buyUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_buy_aiS = bsxfun(@plus, mean(buyUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_aiS));   % 95% Confidence Intervals

holdMean_aiS = mean(holdUnits_aiS,2);

SEM_hold_aiS = std(holdUnits_aiS, [], 2)./ sqrt(size(holdUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_hold_aiS = bsxfun(@plus, mean(holdUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_aiS)); 

sellMean_aiS = mean(sellUnits_aiS,2);

SEM_sell_aiS = std(sellUnits_aiS, [], 2)./ sqrt(size(sellUnits_aiS,2));    % Calculate Standard Error Of The Mean
CI95_sell_aiS = bsxfun(@plus, mean(sellUnits_aiS,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_aiS)); 

% replay
buyMean_replayS = mean(buyUnits_replayS,2);

SEM_buy_replayS = std(buyUnits_replayS, [], 2)./ sqrt(size(buyUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_buy_replayS = bsxfun(@plus, mean(buyUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_replayS));   % 95% Confidence Intervals

holdMean_replayS = mean(holdUnits_replayS,2);

SEM_hold_replayS = std(holdUnits_replayS, [], 2)./ sqrt(size(holdUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_hold_replayS = bsxfun(@plus, mean(holdUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_replayS)); 

sellMean_replayS = mean(sellUnits_replayS,2);

SEM_sell_replayS = std(sellUnits_replayS, [], 2)./ sqrt(size(sellUnits_replayS,2));    % Calculate Standard Error Of The Mean
CI95_sell_replayS = bsxfun(@plus, mean(sellUnits_replayS,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_replayS)); 

% decoy
buyMean_decoyS = mean(buyUnits_decoyS,2);

SEM_buy_decoyS = std(buyUnits_decoyS, [], 2)./ sqrt(size(buyUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_buy_decoyS = bsxfun(@plus, mean(buyUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_decoyS));   % 95% Confidence Intervals

holdMean_decoyS = mean(holdUnits_decoyS,2);

SEM_hold_decoyS = std(holdUnits_decoyS, [], 2)./ sqrt(size(holdUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_hold_decoyS = bsxfun(@plus, mean(holdUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_decoyS)); 

sellMean_decoyS = mean(sellUnits_decoyS,2);

SEM_sell_decoyS = std(sellUnits_decoyS, [], 2)./ sqrt(size(sellUnits_decoyS,2));    % Calculate Standard Error Of The Mean
CI95_sell_decoyS = bsxfun(@plus, mean(sellUnits_decoyS,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_decoyS)); 

% live
buyMean_liveS = mean(buyUnits_liveS,2);

SEM_buy_liveS = std(buyUnits_liveS, [], 2)./ sqrt(size(buyUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_buy_liveS = bsxfun(@plus, mean(buyUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_liveS));   % 95% Confidence Intervals

holdMean_liveS = mean(holdUnits_liveS,2);

SEM_hold_liveS = std(holdUnits_liveS, [], 2)./ sqrt(size(holdUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_hold_liveS = bsxfun(@plus, mean(holdUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_liveS)); 

sellMean_liveS = mean(sellUnits_liveS,2);

SEM_sell_liveS = std(sellUnits_liveS, [], 2)./ sqrt(size(sellUnits_liveS,2));    % Calculate Standard Error Of The Mean
CI95_sell_liveS = bsxfun(@plus, mean(sellUnits_liveS,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_liveS)); 

% OT means
buyMean_aiOT = mean(buyUnits_aiOT,2);

SEM_buy_aiOT = std(buyUnits_aiOT, [], 2)./ sqrt(size(buyUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_buy_aiOT = bsxfun(@plus, mean(buyUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_aiOT));   % 95% Confidence Intervals

holdMean_aiOT = mean(holdUnits_aiOT,2);

SEM_hold_aiOT = std(holdUnits_aiOT, [], 2)./ sqrt(size(holdUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_hold_aiOT = bsxfun(@plus, mean(holdUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_aiOT)); 

sellMean_aiOT = mean(sellUnits_aiOT,2);

SEM_sell_aiOT = std(sellUnits_aiOT, [], 2)./ sqrt(size(sellUnits_aiOT,2));    % Calculate Standard Error Of The Mean
CI95_sell_aiOT = bsxfun(@plus, mean(sellUnits_aiOT,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_aiOT)); 

% replay
buyMean_replayOT = mean(buyUnits_replayOT,2);

SEM_buy_replayOT = std(buyUnits_replayOT, [], 2)./ sqrt(size(buyUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_buy_replayOT = bsxfun(@plus, mean(buyUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_replayOT));   % 95% Confidence Intervals

holdMean_replayOT = mean(holdUnits_replayOT,2);

SEM_hold_replayOT = std(holdUnits_replayOT, [], 2)./ sqrt(size(holdUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_hold_replayOT = bsxfun(@plus, mean(holdUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_replayOT)); 

sellMean_replayOT = mean(sellUnits_replayOT,2);

SEM_sell_replayOT = std(sellUnits_replayOT, [], 2)./ sqrt(size(sellUnits_replayOT,2));    % Calculate Standard Error Of The Mean
CI95_sell_replayOT = bsxfun(@plus, mean(sellUnits_replayOT,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_replayOT)); 

% decoy
buyMean_decoyOT = mean(buyUnits_decoyOT,2);

SEM_buy_decoyOT = std(buyUnits_decoyOT, [], 2)./ sqrt(size(buyUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_buy_decoyOT = bsxfun(@plus, mean(buyUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_decoyOT));   % 95% Confidence Intervals

holdMean_decoyOT = mean(holdUnits_decoyOT,2);

SEM_hold_decoyOT = std(holdUnits_decoyOT, [], 2)./ sqrt(size(holdUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_hold_decoyOT = bsxfun(@plus, mean(holdUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_decoyOT)); 

sellMean_decoyOT = mean(sellUnits_decoyOT,2);

SEM_sell_decoyOT = std(sellUnits_decoyOT, [], 2)./ sqrt(size(sellUnits_decoyOT,2));    % Calculate Standard Error Of The Mean
CI95_sell_decoyOT = bsxfun(@plus, mean(sellUnits_decoyOT,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_decoyOT)); 

% live
buyMean_liveOT = mean(buyUnits_liveOT,2);

SEM_buy_liveOT = std(buyUnits_liveOT, [], 2)./ sqrt(size(buyUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_buy_liveOT = bsxfun(@plus, mean(buyUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_buy_liveOT));   % 95% Confidence Intervals

holdMean_liveOT = mean(holdUnits_liveOT,2);

SEM_hold_liveOT = std(holdUnits_liveOT, [], 2)./ sqrt(size(holdUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_hold_liveOT = bsxfun(@plus, mean(holdUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_hold_liveOT)); 

sellMean_liveOT = mean(sellUnits_liveOT,2);

SEM_sell_liveOT = std(sellUnits_liveOT, [], 2)./ sqrt(size(sellUnits_liveOT,2));    % Calculate Standard Error Of The Mean
CI95_sell_liveOT = bsxfun(@plus, mean(sellUnits_liveOT,2), bsxfun(@times, [-1  1]*1.96, SEM_sell_liveOT)); 


data(:,:,1) = buyUnits_aiS;
data(:,:,2) = holdUnits_aiS;
data(:,:,3) = sellUnits_aiS;
data(:,:,4) = buyUnits_replayS(:,:);
data(:,:,5) = holdUnits_replayS(:,:);
data(:,:,6) = sellUnits_replayS(:,:);
data(:,:,7) = buyUnits_decoyS(:,:);
data(:,:,8) = holdUnits_decoyS(:,:);
data(:,:,9) = sellUnits_decoyS(:,:);
data(:,:,10) = buyUnits_liveS;
data(:,:,11) = holdUnits_liveS;
data(:,:,12) = sellUnits_liveS;
data(:,:,13) = buyUnits_aiOT;
data(:,:,14) = holdUnits_aiOT;
data(:,:,15) = sellUnits_aiOT;
data(:,:,16) = buyUnits_replayOT;
data(:,:,17) = holdUnits_replayOT;
data(:,:,18) = sellUnits_replayOT;
data(:,:,19) = buyUnits_decoyOT;
data(:,:,20) = holdUnits_decoyOT;
data(:,:,21) = sellUnits_decoyOT;
data(:,:,22) = buyUnits_liveOT;
data(:,:,23) = holdUnits_liveOT;
data(:,:,24) = sellUnits_liveOT;



time = -0.5:0.005:1.495; % seconds, sampled at 5 Hz

colors = cbrewer2('qual', 'Set1', 9);
colors = [colors(3,:);colors(2,:);colors(1,:)];

nl = @(s) strrep(s,'\n',char(10));



subplot(2,4,1);
bl1 = boundedline(time, buyMean_aiS, SEM_buy_aiS, ...
    time, holdMean_aiS, SEM_hold_aiS, ...
    time, sellMean_aiS, SEM_sell_aiS, ...
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
% title('AI');
xlim([-0.5 1]); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,5);
bl2 = boundedline(time, buyMean_aiOT, SEM_buy_aiOT, ...
    time, holdMean_aiOT, SEM_hold_aiOT, ...
    time, sellMean_aiOT, SEM_sell_aiOT, ...
    'cmap', colors,'transparency',0.5);
hold on
set(gca,'FontSize',16);
rect2 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect2.FaceAlpha=0.2;
uistack(rect2,'bottom');
% Add lines
h2 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h2,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('AI');
xlim([-0.5 1]);  ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,2);
bl3 = boundedline(time, buyMean_replayS, SEM_buy_replayS, ...
    time, holdMean_replayS, SEM_hold_replayS, ...
    time, sellMean_replayS, SEM_sell_replayS, ...
    'cmap', colors,'transparency',0.5);
hold on
rect3 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect3.FaceAlpha=0.2;
uistack(rect3,'bottom');
set(gca,'FontSize',16);
% Add lines
h4 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h4,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Replay');
xlim([-0.5 1]); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,6);
bl4 = boundedline(time, buyMean_replayOT, SEM_buy_replayOT, ...
    time, holdMean_replayOT, SEM_hold_replayOT, ...
    time, sellMean_replayOT, SEM_sell_replayOT, ...
    'cmap', colors,'transparency',0.5);
hold on
rect4 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect4.FaceAlpha=0.2;
uistack(rect4,'bottom');
set(gca,'FontSize',16);
% Add lines
h4 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h4,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Replay');
xlim([-0.5 1]); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,3);
bl5 = boundedline(time, buyMean_decoyS, SEM_buy_decoyS, ...
    time, holdMean_decoyS, SEM_hold_decoyS, ...
    time, sellMean_decoyS, SEM_sell_decoyS, ...
    'cmap', colors,'transparency',0.5);
hold on
rect5 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect5.FaceAlpha=0.2;
uistack(rect5,'bottom');
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Decoy');
xlim([-0.5 1]); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,7);
bl6 = boundedline(time, buyMean_decoyOT, SEM_buy_decoyOT, ...
    time, holdMean_decoyOT, SEM_hold_decoyOT, ...
    time, sellMean_decoyOT, SEM_sell_decoyOT, ...
    'cmap', colors,'transparency',0.5);
hold on
rect6 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect6.FaceAlpha=0.2;
uistack(rect6,'bottom');
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Decoy');
xlim([-0.5 1]); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,4);
bl7 = boundedline(time, buyMean_liveS, SEM_buy_liveS, ...
    time, holdMean_liveS, SEM_hold_liveS, ...
    time, sellMean_liveS, SEM_sell_liveS, ...
    'cmap', colors,'transparency',0.5);
hold on
rect7 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect7.FaceAlpha=0.2;
uistack(rect7,'bottom');
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Live');
xlim([-0.5 1]); ylim([-1.5 2]); 
% ylabel(nl('Normalized Firing Rate'));
hold off

subplot(2,4,8);
bl8 = boundedline(time, buyMean_liveOT, SEM_buy_liveOT, ...
    time, holdMean_liveOT, SEM_hold_liveOT, ...
    time, sellMean_liveOT, SEM_sell_liveOT, ...
    'cmap', colors,'transparency',0.5);
hold on
rect8 = fill([-0.25 -0.05 -0.05 -0.25 ],[-1.5 -1.5 2.5 2.5],'yellow','LineStyle','none'); 
rect8.FaceAlpha=0.2;
uistack(rect8,'bottom');
set(gca,'FontSize',16);
% Add lines
h5 = line([0 0],[-1.5 2.5]);
% Set properties of lines
set(h5,'Color',[0.4 0.4 0.4],'LineWidth',0.75,'LineStyle','--')
% title('Live');
xlim([-0.5 1]); ylim([-1.5 2]); 
% xlabel('Time (s)');
% ylabel(nl('Normalized Firing Rate'));

 
% instead of a legend, show colored text
lh = legend(bl8);
legnames = {'Buy','Hold','Sell'};
for i = 1:length(legnames),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(1) = lpos(1) + 0.09;
lpos(2) = lpos(2) + 0.81;
lh.Position = lpos;
lh.FontSize = 16;

hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 6];

cd C:\\Users\\plattlab\\Tim\\Stock_market_experiment\\MATLAB_Copy\\MSM_SpikeData
print('OTsalineM2choicePanelbyConditionV2','-dpdf');

