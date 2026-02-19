% File: d20251216_plot_event_relationships.m
% This function was written to validate the association between the box
% size on the display, the cost/gain it represents, and the recorded juice
% delivery time. Tim 251103

function d20251216_plot_event_relationships(etab) 
close all;
% Extract relevant data from etab for plotting
buyReward = etab.sizeB{1}(:); 
buyPrice = etab.priceBuy{1}(:);   
fb1 = etab.fb1{1}(:);
fb2 = etab.fb2{1}(:);

reward1 = zeros(90, 1);
possibleReward1 = [etab.sizeB{1}(:), etab.sizeH{1}(:), etab.sizeS{1}(:)];
for ii = 1:size(possibleReward1, 1)
    reward1(ii) = possibleReward1(ii, etab.option{1}(ii));
end

% reward2 = etab.divPerShare{1}(:) .* etab.postPortfolio{1}(:);
reward2 = zeros(90, 1);
for ii = 1:90
    currDivPerShare = etab.divPerShare{1}(ii);
    currPortfolio = etab.postPortfolio{1}(ii);
    divSpeed = 1;
    % if currDivPerShare == 0
    %     divSpeed = 788;
    % elseif currDivPerShare == 0.8
    %     divSpeed = 394;
    % elseif currDivPerShare == 2.8
    %     divSpeed = 230;
    % elseif currDivPerShare == 6
    %     divSpeed = 150;
    % end
    reward2(ii) = currPortfolio / divSpeed; % 
end

figure
scatter(buyPrice, buyReward);
xlabel('stock price');
ylabel('reward for choosing BUY');

figure
scatter(reward1, fb1);
xlabel('reward 1');
ylabel('fb1');

figure
scatter(reward2, fb2);
xlabel('reward 2');
ylabel('fb2');
end
