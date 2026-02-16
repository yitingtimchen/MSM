%% Buy Bias Analysis: Live vs Other Contexts
clear; clc; close all;

% ---- Step toggles ----
RUN_STEP3 = 0;
RUN_STEP5 = 1;

BUY_THRES = 9;


% Path to condition packs
% OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs';
% OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only';
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';

% Condition sets
cond_sets = { ...
{'AI','Replay','Decoy','Live'}; ...
{'OT AI','OT Replay','OT Decoy','OT Live'}; ...
{'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
};
setNames = {'baseline','OT','Saline'};
% cond_sets = { ...
% {'AI','Replay','Decoy','Live'} ...
% };
% setNames = {'baseline'};

allResults = struct();

for s = 1:numel(cond_sets)
conds = cond_sets{s};
fprintf('\n=== Analyzing Set: %s ===\n', setNames{s});

% ---- Load condition packs ----
Cstruct = struct();
for c = 1:numel(conds)
    cname = conds{c};
    matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));
    if ~exist(matFile,'file')
        warning('Missing file: %s', matFile);
        continue;
    end
    tmp = load(matFile, 'C');
    Cstruct.(matlab.lang.makeValidName(cname)) = tmp.C;
end
condNames = fieldnames(Cstruct);


%% ---- Step 3. Plot distributions ----
if RUN_STEP3
    % Ensure independence: compute PBuy if missing
    if ~exist('PBuy','var') || isempty(PBuy)
        PBuy = struct();
        trialData = table();
        for c = 1:numel(condNames)
            cname = condNames{c};
            C = Cstruct.(cname);
            for f = 1:numel(C.eventTables)
                etab = C.eventTables{f};
                if ~ismember('option', etab.Properties.VariableNames), continue; end
                if isempty(etab.option{1}), continue; end
                opts = etab.option{1};
                if isnumeric(opts)
                    map = {'Buy','Hold','Sell'};
                    opts = arrayfun(@(x) map{x}, opts, 'UniformOutput', false);
                end
                choices = opts(:);
                valid = ~cellfun(@isempty, choices);
                pBuy = mean(strcmp(choices(valid), 'Buy'));
                PBuy.(cname)(f,1) = pBuy;
            end
        end
    end

    liveName = condNames(contains(condNames,'Live','IgnoreCase',true));
    otherNames = setdiff(condNames, liveName);
    nonsocialNames = condNames(contains(condNames,{'AI','Replay'},'IgnoreCase',true));

    liveVec = PBuy.(liveName{1});

    % Per-session averages
    otherMat = cell2mat(cellfun(@(f) PBuy.(f), otherNames, 'UniformOutput', false));
    allOthersMean = mean(otherMat, 2); % average AI, Replay, Decoy

    nonMat = cell2mat(cellfun(@(f) PBuy.(f), nonsocialNames, 'UniformOutput', false));
    nonMean = mean(nonMat, 2); % average AI + Replay

    % Plot
    G = {liveVec(:), nonMean(:), allOthersMean(:)};
    labels = {'Live','Non-social','All Others'};
    
    figure; hold on
    w = 0.2; n = 200;
    allvals = vertcat(G{:}); allvals = allvals(~isnan(allvals));
    
    for i = 1:numel(G)
        g = G{i}; g = g(~isnan(g));
        if numel(g) < 2
            scatter(i*ones(size(g)), g, 25, 'filled', 'MarkerFaceAlpha',0.6);
            continue
        end
        [f, xi] = ksdensity(g, 'NumPoints', n);
        f = f / max(f) * w;
        x = [i - f, fliplr(i + f)];
        y = [xi, fliplr(xi)];
        patch(x, y, [0.4 0.4 1], 'FaceAlpha', 0.1, 'EdgeColor', [0.1 0.1 0.1]);
        plot([i-w i+w], [median(g) median(g)], 'k-', 'LineWidth', 2)
        jitter = (rand(size(g)) - 0.5) * w;
        scatter(i + jitter, g, 10, 'filled');
    end
    
    xlim([0.5 numel(G)+0.5])
    set(gca, 'XTick', 1:numel(G), 'XTickLabel', labels)
    ylabel('P(Buy)'); box on; grid on; hold off
    
    r = range(allvals); if r == 0, r = 1; end
    ylim([0, 1])    
    title(sprintf('Set %s: P(Buy) Live vs All Others vs Non-social', setNames{s}));
end

%% Step 5. Choice distribution visualization (number of Buy/Hold/Sell per trial)
if RUN_STEP5
    choiceTypes = {'Buy','Hold','Sell'};
    numChoices = numel(choiceTypes);
    
    % Compute counts per trial (each trial as separate observation)
    ChoiceCount = struct();
    trialLabels = struct(); % store condition index per trial

    % NEW: track file index per row and average totals per file for labeling
    FileIndex = struct();      % maps each row in counts -> file number
    FileNames = struct();      % simple labels per file
    FileAvgTotal = struct();   % mean total count across trials per file

    for c = 1:numel(condNames)
        cname = condNames{c};
        C = Cstruct.(cname);
    
        counts = [];
        condIdxPerTrial = [];
        fileIdxPerTrial = [];
        fileNameList = strings(numel(C.eventTables),1);

        for f = 1:numel(C.eventTables)
            etab = C.eventTables{f};
            if ~ismember('option', etab.Properties.VariableNames)
                continue;
            end
            optsAllTrials = etab.option{1}; % cell array, one cell per trial
            % if string(cname) == "Live"
            %     oppoOptsAllTrials = etab.opponentOption{1};
            %     optsAllTrials = [optsAllTrials; oppoOptsAllTrials];
            % end

            % derive a file label (fallback to index)
            try
                out = regexp(C.files{f}, '([0-9A-Za-z]+)_behaviorOnly\.nex$', 'tokens', 'once');
                out = out{1};
                fileNameList(f) = string(out);
            catch
                fileNameList(f) = string(f);
            end
    
            for t = 1:size(optsAllTrials, 2)
                opts = optsAllTrials(:, t);
                if isnumeric(opts)
                    map = {'Buy','Hold','Sell'};
                    opts = arrayfun(@(x) map{x}, opts, 'UniformOutput', false);
                end
                choices = opts(:);
                valid = ~cellfun(@isempty, choices);
    
                trialCounts = zeros(1,numChoices);
                for k = 1:numChoices
                    trialCounts(k) = sum(strcmp(choices(valid), choiceTypes{k}));
                end
    
                counts = [counts; trialCounts];           % each row = trial
                condIdxPerTrial = [condIdxPerTrial; c];   % condition index
                fileIdxPerTrial = [fileIdxPerTrial; f];   % file index
            end
        end
    
        % compute per-file average TOTAL (Buy) across its ~6 trials
        nFiles_c = numel(C.eventTables);
        avgTot = nan(nFiles_c,1);
        for ff = 1:nFiles_c
            rows_ff = (fileIdxPerTrial==ff);
            if any(rows_ff)
                avgTot(ff) = mean(sum(counts(rows_ff,1),2),'omitnan');
            end
        end

        ChoiceCount.(cname) = counts;
        trialLabels.(cname) = condIdxPerTrial;
        FileIndex.(cname) = fileIdxPerTrial;
        FileNames.(cname) = fileNameList;
        FileAvgTotal.(cname) = avgTot;
    end
    
    % Prepare data for plotting
    allConds = fieldnames(ChoiceCount);
    numConds = numel(allConds);

    figure('Name','Choice distributions by condition', 'Color','w', 'Position',[100, 100, 1200, 600]); hold on;
    typeColors = lines(numChoices); % color per choice type (for boxes only)
    width = 0.2;                % width of each box
    spacing = 1.2;              % space between condition triplets
    jitterAmount = 0.2;         % jitter for scatter points
    
    % Position tracker for x-axis labels
    xTicks = zeros(numConds,1);
    
    for c = 1:numConds
        counts = ChoiceCount.(allConds{c}); % trials × 3
        fIdx   = FileIndex.(allConds{c});
        fileLbls = FileNames.(allConds{c});
        avgTot  = FileAvgTotal.(allConds{c});
        nFiles  = max(max(fIdx,[], 'omitnan'), 0);
        if isnan(nFiles) || nFiles==0, nFiles = 0; end
        fileCols = lines(max(nFiles,1)); % unique color per file

        for k = 1:numChoices
            xCenter = (c-1)*spacing + (k-2)*width*1.5; % Buy/Hold/Sell spacing within condition
    
            % Boxplot for this choice type
            h = boxplot(counts(:,k), 'Positions', xCenter, 'Colors', typeColors(k,:), ...
                'Widths', width, 'Symbol','', 'PlotStyle','traditional');
            set(h, 'LineWidth', 2);
    
            % % Scatter overlay colored by FILE (the 6 dots per file share a color)
            % if nFiles>0
            %     for ff = 1:nFiles
            %         mask = (fIdx==ff);
            %         if ~any(mask), continue; end
            %         scatterX = xCenter + (rand(sum(mask),1)-0.5)*jitterAmount;
            %         scatterY = counts(mask,k)+(rand(sum(mask),1)-0.5)*jitterAmount*3;
            %         % scatter(scatterX, scatterY, 10, fileCols(ff,:), 'filled', 'MarkerFaceAlpha',0.9);
            %         scatter(scatterX, scatterY, 10, 'k', 'filled', 'MarkerFaceAlpha',0.3);
            %     end
            % else
            %     % fallback single color
            %     scatterX = xCenter + (rand(size(counts(:,k)))-0.5)*jitterAmount;
            %     scatterY = counts(mask,k)+(rand(sum(mask),1)-0.5)*jitterAmount*3;
            %     scatter(scatterX, scatterY, 10, typeColors(k,:), 'filled', 'MarkerFaceAlpha',0.6);
            % end

            % % Label files with average TOTAL > BUY_THRES (label once at k==3 to avoid triplicate)
            % if nFiles>0 && k==1
            %     for ff = 1:nFiles
            %         if ff<=numel(avgTot) && ~isnan(avgTot(ff)) && avgTot(ff) > BUY_THRES
            %             mask = (fIdx==ff);
            %             if ~any(mask), continue; end
            %             yLabel = mean(counts(mask,k)) + 0.1;
            %             txt = char(fileLbls(ff));
            %             text(xCenter + 0.2, yLabel + 0.02*ff, txt, 'Color', fileCols(ff,:), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
            %         end
            %     end
            % end
        end
        % Position for condition label
        xTicks(c) = (c-1)*spacing;
    end
    
    % count how many sessions are from each monkey
    numM1 = zeros(length(conds), 1);
    numM2 = zeros(length(conds), 1);
    for c = 1:numConds
        for ii = 1:length(Cstruct.(allConds{c}).eventTables)
            if ~isempty(Cstruct.(allConds{c}).eventTables{ii})
                if unique(Cstruct.(allConds{c}).eventTables{ii}.monkey{1}) == 1
                    numM1(c) = numM1(c) + 1;
                elseif unique(Cstruct.(allConds{c}).eventTables{ii}.monkey{1}) == 2
                    numM2(c) = numM2(c) + 1;
                else
                end
            end
        end
    end

    % X-axis: no label for each box, but condition labels centered above triplets
    set(gca, 'XTick', xTicks);
    numSessionsPerCond = cellfun(@(c) ceil(size(ChoiceCount.(c),1)), allConds);
    set(gca, 'XTickLabel', allConds + " = " + numSessionsPerCond./6 + " sesh (" + numM1 + "+" + numM2 + ")");
    xtickangle(0);
    ylabel('N Trials', 'FontWeight','bold', 'FontSize',12);
    ylim([-1 25]);
    yticks(0:3:15);
    % Count number of sessions per condition
    totalSessions = sum(numSessionsPerCond/6);
    
    % xlabel(sprintf('Total Sesh across conds = %d, BUY THRESH = %d', totalSessions, BUY_THRES));
    xlabel(sprintf('Total Sesh across conds = %d', totalSessions), 'FontWeight','bold', 'FontSize',12);
    title(sprintf('Choice distribution per condition\nCondition set = "%s"', setNames{s}), 'FontWeight','bold', 'FontSize',14);
    
    %% === ANOVA + significance overlay ===
    alpha = 0.05;

    % Compute data max across all conditions/choices
    dataMax = -inf;
    for c = 1:numel(fieldnames(ChoiceCount))
        nm = allConds{c};
        dataMax = max(dataMax, max(ChoiceCount.(nm), [], 'all', 'omitnan'));
    end
    yl = ylim; 
    headroom = max(yl(2) - dataMax, 1); % usable space above data (e.g., ~5)

    % Y-anchors (kept well below the top of the axis to avoid clipping)
    y_buy_text = min(yl(2)-0.5, 24);
    y_buy_pair0 = min(yl(2)-1.0, 20);
    y_within_text = min(yl(2)-0.5, 18);
    y_within_pair0 = min(yl(2)-1.0, 16);
    pair_step = min(0.1*headroom, 0.6); % vertical stacking gap for multiple sig lines

    % ---- 1) Buy across CONDITIONS (one-way ANOVA) ----
    y = []; g = [];
    for c = 1:numConds
        cnt = ChoiceCount.(allConds{c});
        y = [y; cnt(:,1)];
        g = [g; repmat(c,size(cnt,1),1)];
    end
    [p_buy,~,stats_buy] = anova1(y, g, 'off');
    [c_buy,~,~,~] = multcompare(stats_buy, 'Display','off');
    % X locations of the Buy boxes per condition
    buyX = arrayfun(@(cc) (cc-1)*spacing + (1-2)*width*1.5, 1:numConds);
    % Global p-value text
    text(mean(buyX), y_buy_text, sprintf('Buy across conds: global p=%.3g', p_buy), ...
    'HorizontalAlignment','center','FontWeight','bold','HandleVisibility','off');
    % Pairwise lines (Tukey), only if few conditions to avoid clutter
    if numConds <= 5
        sigRows = find(c_buy(:,end) < 1);
        y0 = y_buy_pair0; lvl = 0;
        for rr = 1:numel(sigRows)
            i = c_buy(sigRows(rr),1); j = c_buy(sigRows(rr),2);
            plot_sigline([buyX(i), buyX(j)], y0 + lvl*pair_step, p2stars(c_buy(sigRows(rr),end)));
            lvl = lvl + 1;
        end
    end

    % ---- 2) Within each CONDITION: Buy vs Hold vs Sell ----
    for c = 1:numConds
        cnt = ChoiceCount.(allConds{c}); % trials × 3
        y = [cnt(:,1); cnt(:,2); cnt(:,3)];
        gg = [ones(size(cnt,1),1); 2*ones(size(cnt,1),1); 3*ones(size(cnt,1),1)];
        [p_within,~,stats_within] = anova1(y, gg, 'off');
        cc = multcompare(stats_within, 'Display','off');
        % Per-condition ANOVA p-value
        text(xTicks(c), y_within_text, sprintf('Within cond p=%.3g', p_within), ...
            'HorizontalAlignment','center','FontWeight','normal','HandleVisibility','off');
        
        % Pairwise lines within the triplet
        xBuy  = (c-1)*spacing + (1-2)*width*1.5;
        xHold = (c-1)*spacing + (2-2)*width*1.5;
        xSell = (c-1)*spacing + (3-2)*width*1.5;
        xTrip = [xBuy, xHold, xSell];
        
        sigRows = find(cc(:,end) < 1);
        ybase = y_within_pair0; lvl = 0;
        for rr = 1:numel(sigRows)
            i = cc(sigRows(rr),1); j = cc(sigRows(rr),2);
            plot_sigline([xTrip(i), xTrip(j)], ybase + lvl*pair_step, p2stars(cc(sigRows(rr),end)));
            lvl = lvl + 1;
        end
    end
    
    
    
    % Legend for choice types (boxes)
    for k = 1:numChoices
        h(k) = plot(NaN,NaN,'Color',typeColors(k,:),'LineWidth',2);
    end
    legend(h, choiceTypes, 'Location','northeastoutside');
    grid on;
    box on;    
    saveas(gca, [OUTDIR '/' sprintf('choice_distribution_%s_set.png', setNames{s})])
end

end


% ===== helper functions (place at end of script) =====

function plot_sigline(xpair, y, starstr)
hold on;
line(xpair, [y y], 'Color','k','LineWidth',1.2);
dy = 0.02*diff(ylim);
% line([xpair(1) xpair(1)], [y y-dy], 'Color','k','LineWidth',1.0, 'HandleVisibility','off');
% line([xpair(2) xpair(2)], [y y-dy], 'Color','k','LineWidth',1.0, 'HandleVisibility','off');
text(mean(xpair), y + 0.01*diff(ylim), starstr, 'HorizontalAlignment','center','FontWeight','bold');
end

function s = p2stars(p)
if p < 1e-3, s = '***';
elseif p < 1e-2, s = '**';
elseif p < 0.05, s = '*';
elseif p < 0.01, s = '.';
else, s = 'n.s.';
end
end
