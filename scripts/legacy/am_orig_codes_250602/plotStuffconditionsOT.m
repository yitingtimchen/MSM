clear all
% plotStuff

% Plot psth from Nex data
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\AI_OT
load('m1rew1on.mat')
aiUnits_OT = spikes;
aiNames_OT = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Replay_OT
load('m1rew1on.mat')
replayUnits_OT = spikes;
replayNames_OT = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Decoy_OT
load('m1rew1on.mat')
decoyUnits_OT = spikes;
decoyNames_OT = units;
cd C:\\Users\\plattlab\\MSM\\data\\tooling\\MATLAB_Copy\\MSM_SpikeData\OT\Live_OT
load('m1rew1on.mat')
liveUnits_OT = spikes;
liveNames_OT = units;

clear units2 units3 spikes2

aiMean_OT = mean(aiUnits_OT,2);

SEM_ai_OT = std(aiUnits_OT, [], 2)./ sqrt(size(aiUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_ai_OT = bsxfun(@plus, mean(aiUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_ai_OT));   % 95% Confidence Intervals

replayMean_OT = mean(replayUnits_OT,2);

SEM_replay_OT = std(replayUnits_OT, [], 2)./ sqrt(size(replayUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_replay_OT = bsxfun(@plus, mean(replayUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_replay_OT)); 

decoyMean_OT = mean(decoyUnits_OT,2);

SEM_decoy_OT = std(decoyUnits_OT, [], 2)./ sqrt(size(decoyUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_decoy_OT = bsxfun(@plus, mean(decoyUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_decoy_OT)); 

liveMean_OT = mean(liveUnits_OT,2);

SEM_live_OT = std(liveUnits_OT, [], 2)./ sqrt(size(liveUnits_OT,2));    % Calculate Standard Error Of The Mean
CI95_live_OT = bsxfun(@plus, mean(liveUnits_OT,2), bsxfun(@times, [-1  1]*1.96, SEM_live_OT)); 

e1 = errorbar(-0.5:0.005:1.495,aiMean_OT,SEM_ai_OT,'CapSize',0);
hold on
e2 = errorbar(-0.5:0.005:1.495,replayMean_OT, SEM_replay_OT,'CapSize',0);
e3 = errorbar(-0.5:0.005:1.495,decoyMean_OT, SEM_decoy_OT,'CapSize',0);
e4 = errorbar(-0.5:0.005:1.495,liveMean_OT, SEM_live_OT,'CapSize',0);
hold off
e1.Bar.LineStyle = 'dotted';
e2.Bar.LineStyle = 'dotted';
e3.Bar.LineStyle = 'dotted';
e4.Bar.LineStyle = 'dotted';
e1.CapSize = 0.1;

% add background box from internet example:
% figure
% %Plot something
% plot(1:10)
% Add lines
% h1 = line([-0.5 -0.5],[-1.5 2.5]);
% h2 = line([-0.1 -0.1],[-1.5 2.5]);
% % Set properties of lines
% set([h1 h2],'Color','k','LineWidth',2)
% % Add a patch
% gray = [0.4 0.4 0.4];
% patch([-0.5 -0.01 -0.01 -0.5],[-1.5 -1.5 2.5 2.5],gray)
% % The order of the "children" of the plot determines which one appears on top.
% % I need to flip it here.
% set(gca,'children',flipud(get(gca,'children')))

% numUnits = numel(nexColumnNames);
% 
% figure
% hold on
% 
% for i = 2:numUnits
%     sig = nex(:,i)';
%     plot(-0.5:0.005:2.495,sig)
% end
% 
% plot(-2.5:0.005:3.995,mNex)
% 
% hold off
%     

