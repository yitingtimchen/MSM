clear all
% plotStuff

% Plot psth from Nex data


numUnits = numel(nexColumnNames);

figure
hold on

for i = 2:numUnits
    sig = nex(:,i)';
    plot(-0.5:0.005:2.495,sig)
end

plot(-2.5:0.005:3.995,mNex)

hold off
    

