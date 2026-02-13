%% Buy Bias Analysis: Live vs Other Contexts (session-exclude filter version)
% ================================================================
% THREE 3-panel figures per condition set:
%   1) Summary-only (BOXPLOT or BAR+ERROR; see SUMMARY_PLOT_TYPE)
%   2) Scatter only (rainbow color = session, ordered high->low)
%   3) Combined summary + scatter
%
% Split into:
%   A) Both monkeys combined
%   B) Monkey 1 only
%   C) Monkey 2 only
%
% Session filtering:
%   - You specify explicit sessions to EXCLUDE, as strings:
%       'yyyy-mm-dd-sod'  (sod = session-of-day, 2 digits)
%   - Any file whose parsed (date, session-of-day) matches the exclude list
%     is skipped (across all condition sets).
%
% Aggregation:
%   - USE_SESSION_MEAN = 0 : each TRIAL is one observation
%   - USE_SESSION_MEAN = 1 : ONE mean vector [Buy,Hold,Sell] per SESSION
%
% NEW:
%   - SUMMARY_PLOT_TYPE = 'box' or 'bar'
%     If 'bar', plots mean ± error bars (SEM by default; see BAR_ERROR_TYPE).
% ================================================================

clear; clc;

% ---- Global toggles ----
SHOW_SIG         = 0;   % 1 = show significance markers; 0 = hide

DO_BOX_FIG       = 1;   % 1 = make summary-only figure (boxplot OR bar+error; see SUMMARY_PLOT_TYPE)
DO_SCATTER_FIG   = 0;   % 1 = make scatter-only figure
DO_COMBINED_FIG  = 0;   % 1 = make combined summary+scatter figure

USE_SESSION_MEAN = 0;   % 0 = per trial, 1 = one mean value per session

% ---- Summary plot toggle (box vs bar+error) ----
SUMMARY_PLOT_TYPE = 'box';   % 'box' or 'bar'
BAR_ERROR_TYPE    = 'sd';   % 'sem' or 'sd' (only used when SUMMARY_PLOT_TYPE='bar')

% ---- Session EXCLUDE list filter ----
USE_EXCLUDE_LIST = 1;   % 1 = exclude sessions in EXCLUDE_SESSIONS; 0 = keep all

% >>> sessions to exclude (yyyy-mm-dd-sod) <<<
EXCLUDE_SESSIONS = { ...
    '2018-07-07-01', ... % M1, baseline AI
    '2018-07-20-02', ... % M1, baseline Replay
    '2018-07-21-01', ... % M1, baseline Decoy
    '2018-07-06-02', ... % M2, baseline AI
    '2018-07-06-01', ... % M2, baseline AI
    '2018-08-01-03', ... % M2, baseline Decoy
    '2018-08-22-02', ... % M2, OT Decoy
};

% Path to condition packs
OUTDIR = 'C:\Users\plattlab\Tim\Stock_market_experiment\Tim\condition_packs_behavior_only_from_dataTabFULL';

% Condition sets (each "set" is a group of corresponding condition packs)
% Comment out any rows you do NOT want to run.
cond_sets = { ...
    {'AI','Replay','Decoy','Live'}; ...
    {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
    {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
    };
setNames = {'baseline','OT','Saline'};

allResults = struct(); %#ok<NASGU>

for s = 1:numel(cond_sets)
    conds = cond_sets{s};
    fprintf('\n=== Analyzing Set: %s ===\n', setNames{s});

    % ---- Load condition packs for this set ----
    Cstruct = struct();
    for c = 1:numel(conds)
        cname   = conds{c};
        matFile = fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cname));
        if ~exist(matFile,'file')
            warning('Missing file: %s', matFile);
            continue;
        end
        tmp = load(matFile, 'C');
        Cstruct.(matlab.lang.makeValidName(cname)) = tmp.C;
    end
    condNames = fieldnames(Cstruct);  % condition names actually loaded

    %% Step 5. Choice distribution visualization (number of Buy/Hold/Sell per trial or per session)
    choiceTypes = {'Buy','Hold','Sell'};
    numChoices  = numel(choiceTypes);

    % Compute counts (either per trial or per session)
    ChoiceCount  = struct();   % ChoiceCount.(condName): [nObs x 3]
    trialLabels  = struct();   % kept for completeness (per-trial mode)
    FileIndex    = struct();   % FileIndex.(condName): file index per observation
    FileNames    = struct();   % FileNames.(condName): human-readable file label per session
    FileAvgTotal = struct();   % FileAvgTotal.(condName): mean Buy count per session
    FileMonkey   = struct();   % FileMonkey.(condName)(fileIndex,1) = monkey ID (1 or 2)

    for c = 1:numel(condNames)
        cname = condNames{c};
        C     = Cstruct.(cname);

        countsTrial      = [];    % trials x 3
        condIdxPerTrial  = [];
        fileIdxPerTrial  = [];
        fileNameList     = strings(numel(C.eventTables),1);
        monkeyPerFile    = nan(numel(C.eventTables),1); % per-file monkey ID

        for f = 1:numel(C.eventTables)

            % ---- SESSION EXCLUDE FILTER ----
            if USE_EXCLUDE_LIST
                thisFile = C.files{f};
                [~, base, ~] = fileparts(thisFile);

                % Expect base like: mmddyyL1_behaviorOnly  (L/B/X = monkey, 1+ = session of day)
                tok = regexp(base, '(\d{6})([BLX])(\d+)', 'tokens', 'once');

                if ~isempty(tok)
                    mmddyy = tok{1};
                    mm     = str2double(mmddyy(1:2));
                    dd     = str2double(mmddyy(3:4));
                    yy     = str2double(mmddyy(5:6));   % last two digits
                    sod    = str2double(tok{3});        % session number in that day

                    % convert yy -> yyyy (safe for 19xx/20xx)
                    if yy < 70
                        yyyy = 2000 + yy;
                    else
                        yyyy = 1900 + yy;
                    end

                    sessKey = sprintf('%04d-%02d-%02d-%02d', yyyy, mm, dd, sod);

                    % skip this file if in exclude list
                    if any(strcmp(sessKey, EXCLUDE_SESSIONS))
                        continue;
                    end
                else
                    warning('Could not parse date/session from filename (no exclude applied): %s', base);
                end
            end
            % -------------------------------

            etab = C.eventTables{f};

            if ~ismember('option', etab.Properties.VariableNames)
                continue;
            end
            optsAllTrials = etab.option{1};   % cell array: one cell per trial

            % Per-file monkey ID
            mID = NaN;
            if ismember('monkey', etab.Properties.VariableNames) && ~isempty(etab.monkey{1})
                mvals = etab.monkey{1};
                mu    = unique(mvals);
                if numel(mu) == 1
                    mID = mu;                % 1 or 2
                end
            end
            monkeyPerFile(f,1) = mID;

            % Simple file label
            try
                out = regexp(C.files{f}, '([0-9A-Za-z]+)_behaviorOnly\.nex$', 'tokens', 'once');
                out = out{1};
                fileNameList(f) = string(out);
            catch
                fileNameList(f) = string(f);
            end

            % For each trial within this file
            for t = 1:size(optsAllTrials, 2)
                opts = optsAllTrials(:, t);

                % Convert numeric codes to strings if needed
                if isnumeric(opts)
                    map  = {'Buy','Hold','Sell'};
                    opts = arrayfun(@(x) map{x}, opts, 'UniformOutput', false);
                end
                choices = opts(:);
                valid   = ~cellfun(@isempty, choices);

                % Count Buy/Hold/Sell occurrences for this trial
                trialCounts = zeros(1,numChoices);
                for k = 1:numChoices
                    trialCounts(k) = sum(strcmp(choices(valid), choiceTypes{k}));
                end

                countsTrial      = [countsTrial; trialCounts];   %#ok<AGROW>
                condIdxPerTrial  = [condIdxPerTrial; c];         %#ok<AGROW>
                fileIdxPerTrial  = [fileIdxPerTrial; f];         %#ok<AGROW>
            end
        end

        % Per-session average TOTAL (Buy) across its trials (always per-session here)
        nFiles_c = numel(C.eventTables);
        avgTot   = nan(nFiles_c,1);
        for ff = 1:nFiles_c
            rows_ff = (fileIdxPerTrial == ff);
            if any(rows_ff)
                avgTot(ff) = mean(sum(countsTrial(rows_ff,1),2), 'omitnan');
            end
        end

        % Optionally collapse to ONE mean vector per session
        if USE_SESSION_MEAN
            countsSession  = [];
            fileIdxSession = [];
            for ff = 1:nFiles_c
                rows_ff = (fileIdxPerTrial == ff);
                if any(rows_ff)
                    countsSession   = [countsSession; mean(countsTrial(rows_ff,:),1,'omitnan')]; %#ok<AGROW>
                    fileIdxSession  = [fileIdxSession; ff];                                    %#ok<AGROW>
                end
            end
            ChoiceCount.(cname) = countsSession;     % each row = session
            FileIndex.(cname)   = fileIdxSession;    % one entry per session-row
            trialLabels.(cname) = [];                % not used in this mode
        else
            ChoiceCount.(cname) = countsTrial;       % each row = trial
            FileIndex.(cname)   = fileIdxPerTrial;   % one file index per trial
            trialLabels.(cname) = condIdxPerTrial;   %#ok<NASGU>
        end

        FileNames.(cname)    = fileNameList;
        FileAvgTotal.(cname) = avgTot;
        FileMonkey.(cname)   = monkeyPerFile;
    end

    % Prepare data for plotting
    allConds  = fieldnames(ChoiceCount);
    keepConds = false(numel(allConds),1);
    for c = 1:numel(allConds)
        if ~isempty(ChoiceCount.(allConds{c}))
            keepConds(c) = true;
        end
    end
    allConds = allConds(keepConds);

    % ---- Make the three figures (summary / scatter / combined) ----
    if DO_BOX_FIG
        make_choice_figure(ChoiceCount, FileIndex, FileNames, FileAvgTotal, FileMonkey, ...
                           allConds, setNames{s}, 'boxOnly', ...
                           true, false, SHOW_SIG, OUTDIR, USE_SESSION_MEAN, ...
                           SUMMARY_PLOT_TYPE, BAR_ERROR_TYPE);
    end

    if DO_SCATTER_FIG
        make_choice_figure(ChoiceCount, FileIndex, FileNames, FileAvgTotal, FileMonkey, ...
                           allConds, setNames{s}, 'scatterOnly', ...
                           false, true, SHOW_SIG, OUTDIR, USE_SESSION_MEAN, ...
                           SUMMARY_PLOT_TYPE, BAR_ERROR_TYPE);
    end

    if DO_COMBINED_FIG
        make_choice_figure(ChoiceCount, FileIndex, FileNames, FileAvgTotal, FileMonkey, ...
                           allConds, setNames{s}, 'combined', ...
                           true, true, SHOW_SIG, OUTDIR, USE_SESSION_MEAN, ...
                           SUMMARY_PLOT_TYPE, BAR_ERROR_TYPE);
    end

end  % for each condition set


% ===== helper functions ================================================

function make_choice_figure(ChoiceCount, FileIndex, FileNames, FileAvgTotal, FileMonkey, ...
                            allConds, setName, figTag, ...
                            SHOW_BOXPLOT, SHOW_SCATTER, SHOW_SIG, OUTDIR, USE_SESSION_MEAN, ...
                            SUMMARY_PLOT_TYPE, BAR_ERROR_TYPE)

choiceTypes = {'Buy','Hold','Sell'};
numChoices  = numel(choiceTypes);
numConds    = numel(allConds);

typeColors   = lines(numChoices);  % Buy/Hold/Sell colors
width        = 0.2;
spacing      = 1.2;
jitterAmount = 0.05;

figName = sprintf('Choice distributions (%s) - %s', figTag, setName);
if SHOW_BOXPLOT && strcmpi(SUMMARY_PLOT_TYPE,'bar')
    figName = sprintf('%s [bar+err]', figName);
elseif SHOW_BOXPLOT && strcmpi(SUMMARY_PLOT_TYPE,'box')
    figName = sprintf('%s [box]', figName);
end

figure('Name', figName, ...
       'Color','w', 'Position',[100, 100, 1600, 600]);

monkeyGroups = {'both','M1','M2'};
groupTitles  = {'Both monkeys combined', 'Monkey 1 only', 'Monkey 2 only'};

for mg = 1:3
    subplot(1,3,mg); hold on;

    ChoiceCount_sub = struct();
    FileIndex_sub   = struct();

    dataMax     = -inf;
    condHasData = false(numConds,1);

    % ---- Subset observations for this monkey group ----
    for c = 1:numConds
        cname      = allConds{c};
        countsAll  = ChoiceCount.(cname);
        fIdxAll    = FileIndex.(cname);
        fileMonkey = FileMonkey.(cname);

        if isempty(countsAll)
            ChoiceCount_sub.(cname) = [];
            FileIndex_sub.(cname)   = [];
            continue;
        end

        monkeyPerObs = nan(size(fIdxAll));
        validMask    = fIdxAll > 0 & fIdxAll <= numel(fileMonkey);
        monkeyPerObs(validMask) = fileMonkey(fIdxAll(validMask));

        switch monkeyGroups{mg}
            case 'both'
                mask_obs = ~isnan(monkeyPerObs);
            case 'M1'
                mask_obs = (monkeyPerObs == 1);
            case 'M2'
                mask_obs = (monkeyPerObs == 2);
        end

        countsSub = countsAll(mask_obs,:);
        fIdxSub   = fIdxAll(mask_obs);

        ChoiceCount_sub.(cname) = countsSub;
        FileIndex_sub.(cname)   = fIdxSub;

        if ~isempty(countsSub)
            dataMax = max(dataMax, max(countsSub, [], 'all', 'omitnan'));
            condHasData(c) = true;
        end
    end

    usedConds = allConds(condHasData);
    numUsed   = numel(usedConds);

    xTicks = zeros(numUsed,1);

    % --- legend accumulators (one legend per subplot) ---
    sessLegendHandles = gobjects(0);
    sessLegendLabels  = strings(0,1);

    for ic = 1:numUsed
        cname    = usedConds{ic};
        counts   = ChoiceCount_sub.(cname);  % obs x 3 (obs = trials or sessions)
        fIdx     = FileIndex_sub.(cname);
        fileLbls = FileNames.(cname);
        avgTot   = FileAvgTotal.(cname);     %#ok<NASGU>

        if isempty(counts)
            continue;
        end

        % --- SESSION COLORS (RAINBOW) ORDERED HIGH->LOW BY MEAN BUY ---
        validFiles  = fIdx(fIdx > 0);
        uniqueFiles = unique(validFiles);
        nFiles      = numel(uniqueFiles);

        if nFiles > 0
            sessMetric = nan(nFiles,1);
            for jj = 1:nFiles
                fID  = uniqueFiles(jj);
                rows = (fIdx == fID);
                sessMetric(jj) = mean(counts(rows,1), 'omitnan');  % col 1 = Buy
            end

            [~, sortIdx]       = sort(sessMetric, 'descend');
            uniqueFilesSorted  = uniqueFiles(sortIdx);

            cmapBase = jet(nFiles);
            cmap     = cmapBase(sortIdx,:);

            [~, loc] = ismember(fIdx, uniqueFilesSorted);
            loc(loc == 0) = 1;
            sessionColors = cmap(loc, :);

            for jj = 1:nFiles
                fID = uniqueFilesSorted(jj);
                if fID <= 0 || fID > numel(fileLbls)
                    continue;
                end
                if USE_SESSION_MEAN
                    thisLabel = sprintf('%s - %s (sess mean)', cname, fileLbls(fID));
                else
                    thisLabel = sprintf('%s - %s', cname, fileLbls(fID));
                end

                if any(sessLegendLabels == thisLabel)
                    continue;
                end

                hLeg = scatter(NaN, NaN, 30, cmap(jj,:), 'filled');
                sessLegendHandles(end+1,1) = hLeg; %#ok<AGROW>
                sessLegendLabels(end+1,1)  = thisLabel; %#ok<AGROW>
            end
        else
            sessionColors = repmat([0 0 0], size(fIdx,1), 1);
        end
        % ----------------------------------------------------------------

        for k = 1:numChoices
            xCenter = (ic-1)*spacing + (k-2)*width*1.5;

            if SHOW_BOXPLOT
                if strcmpi(SUMMARY_PLOT_TYPE,'box')
                    hBox = boxplot(counts(:,k), 'Positions', xCenter, ...
                                   'Colors', typeColors(k,:), ...
                                   'Widths', width, 'Symbol','', ...
                                   'PlotStyle','traditional');
                    set(hBox, 'LineWidth', 2);

                elseif strcmpi(SUMMARY_PLOT_TYPE,'bar')
                    yk = counts(:,k);
                    mu = mean(yk, 'omitnan');

                    nn = sum(~isnan(yk));
                    if nn <= 1
                        err = 0;
                    else
                        switch lower(BAR_ERROR_TYPE)
                            case 'sd'
                                err = std(yk, 'omitnan');
                            otherwise % 'sem'
                                err = std(yk, 'omitnan') ./ sqrt(nn);
                        end
                    end

                    hb = bar(xCenter, mu, 'BarWidth', width, ...
                             'FaceColor', typeColors(k,:), ...
                             'EdgeColor','k', 'LineWidth', 1.2); %#ok<NASGU>

                    errorbar(xCenter, mu, err, 'k', 'LineStyle','none', ...
                             'LineWidth', 1.2, 'CapSize', 8, ...
                             'HandleVisibility','off');
                else
                    error('SUMMARY_PLOT_TYPE must be ''box'' or ''bar''.');
                end
            end

            if SHOW_SCATTER
                scatterX = xCenter + (rand(size(counts(:,k))) - 0.5)*jitterAmount;
                scatterY = counts(:,k) + (rand(length(counts(:,k)),1) - 0.5)*jitterAmount*3;

                scatter(scatterX, scatterY, 15, sessionColors, ...
                        'filled', 'HandleVisibility','off');
            end
        end

        xTicks(ic) = (ic-1)*spacing;
    end

    % ---- Buy/Hold/Sell + session legend on this subplot ----
    hType = gobjects(numChoices,1);
    for k = 1:numChoices
        hType(k) = scatter(NaN, NaN, 40, typeColors(k,:), 'filled', ...
                           'MarkerEdgeColor','k');
    end

    if ~isempty(sessLegendHandles) && SHOW_SCATTER
        allHandles = [hType(:); sessLegendHandles(:)];
        allLabels  = [choiceTypes, cellstr(sessLegendLabels')];
    else
        allHandles = hType(:);
        allLabels  = choiceTypes;
    end

    legend(allHandles, allLabels, ...
           'Location','northwest', 'Interpreter','none');

    set(gca, 'XTick', xTicks);

    labels = cell(numUsed,1);
    for ic = 1:numUsed
        cname   = usedConds{ic};
        fIdxSub = FileIndex_sub.(cname);
        if isempty(fIdxSub)
            nSess = 0;
        else
            nSess = numel(unique(fIdxSub(fIdxSub > 0)));
        end
        if USE_SESSION_MEAN
            labels{ic} = sprintf('%s\\newline%d sess (means)', cname, nSess);
        else
            labels{ic} = sprintf('%s\\newline%d sess', cname, nSess);
        end
    end
    set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex');
    xtickangle(0);

    if USE_SESSION_MEAN
        if SHOW_BOXPLOT && strcmpi(SUMMARY_PLOT_TYPE,'bar')
            ylabel(sprintf('Mean N per session (%s)', upper(BAR_ERROR_TYPE)), ...
                   'FontWeight','bold', 'FontSize', 12);
        else
            ylabel('Mean N per session', 'FontWeight','bold', 'FontSize', 12);
        end
    else
        if SHOW_BOXPLOT && strcmpi(SUMMARY_PLOT_TYPE,'bar')
            ylabel(sprintf('Mean N per trial (%s)', upper(BAR_ERROR_TYPE)), ...
                   'FontWeight','bold', 'FontSize', 12);
        else
            ylabel('N Trials', 'FontWeight','bold', 'FontSize', 12);
        end
    end

    if isfinite(dataMax)
        ylMax = max(ceil(dataMax + 2), 3);
        ylim([-0.5, ylMax]);
        yticks(0:3:ylMax);
    else
        ylim([0 1]);
        yticks(0:1:1);
    end

    totalSessions = 0;
    for ic = 1:numUsed
        cname   = usedConds{ic};
        fIdxSub = FileIndex_sub.(cname);
        if ~isempty(fIdxSub)
            totalSessions = totalSessions + numel(unique(fIdxSub(fIdxSub > 0)));
        end
    end
    xlabel(sprintf('Total sessions across conditions = %d', totalSessions), ...
           'FontWeight','bold', 'FontSize', 12);

    title(groupTitles{mg}, 'FontWeight','bold', 'FontSize', 14);

    box off;

    if SHOW_SIG
        alpha = 0.05; %#ok<NASGU>

        if isfinite(dataMax)
            yl       = ylim;
            headroom = max(yl(2) - dataMax, 1);
        else
            yl       = ylim;
            headroom = yl(2) - yl(1);
        end

        y_buy_pair0    = min(yl(2) - 1.0, yl(2) - 1.0);
        y_within_pair0 = min(yl(2) - 1.0, yl(2) - 1.0);
        pair_step      = min(0.1*headroom, 0.6);

        % ---- 1) Buy across CONDITIONS (one-way ANOVA) ----
        yANOVA = [];
        gANOVA = [];
        numUsedLocal = numUsed;
        for ic = 1:numUsedLocal
            nm  = usedConds{ic};
            cnt = ChoiceCount_sub.(nm);
            if isempty(cnt)
                continue;
            end
            yANOVA = [yANOVA; cnt(:,1)]; %#ok<AGROW>
            gANOVA = [gANOVA; repmat(ic, size(cnt,1), 1)]; %#ok<AGROW>
        end

        if numel(unique(gANOVA)) > 1 && numel(yANOVA) > 1
            [p_buy,~,stats_buy] = anova1(yANOVA, gANOVA, 'off'); %#ok<ASGLU>
            [c_buy,~,~,~]       = multcompare(stats_buy, 'Display', 'off');

            buyX = arrayfun(@(cc) (cc-1)*spacing + (1-2)*width*1.5, 1:numUsedLocal);

            if numUsedLocal <= 5
                sigRows = find(c_buy(:,end) < 1);
                y0      = y_buy_pair0;
                lvl     = 0;
                for rr = 1:numel(sigRows)
                    i = c_buy(sigRows(rr),1);
                    j = c_buy(sigRows(rr),2);

                    if any(gANOVA == i) && any(gANOVA == j)
                        plot_sigline([buyX(i), buyX(j)], ...
                                     y0 + lvl*pair_step, ...
                                     p2stars(c_buy(sigRows(rr),end)));
                        lvl = lvl + 1;
                    end
                end
            end
        end

        % ---- 2) Within each CONDITION: Buy vs Hold vs Sell ----
        for ic = 1:numUsedLocal
            nm  = usedConds{ic};
            cnt = ChoiceCount_sub.(nm);
            if isempty(cnt) || size(cnt,1) < 2
                continue;
            end

            yWithin = [cnt(:,1); cnt(:,2); cnt(:,3)];
            gg      = [ones(size(cnt,1),1); ...
                       2*ones(size(cnt,1),1); ...
                       3*ones(size(cnt,1),1)];
            [p_within,~,stats_within] = anova1(yWithin, gg, 'off'); %#ok<ASGLU>
            cc_within = multcompare(stats_within, 'Display','off');

            xBuy  = (ic-1)*spacing + (1-2)*width*1.5;
            xHold = (ic-1)*spacing + (2-2)*width*1.5;
            xSell = (ic-1)*spacing + (3-2)*width*1.5;
            xTrip = [xBuy, xHold, xSell];

            sigRows = find(cc_within(:,end) < 1);
            ybase   = y_within_pair0;
            lvl     = 0;
            for rr = 1:numel(sigRows)
                i = cc_within(sigRows(rr),1);
                j = cc_within(sigRows(rr),2);
                plot_sigline([xTrip(i), xTrip(j)], ...
                             ybase + lvl*pair_step, ...
                             p2stars(cc_within(sigRows(rr),end)));
                lvl     = lvl + 1;
            end
        end
    end  % SHOW_SIG

end  % monkey-group loop

ax = findall(gcf,'Type','axes');
for k = 1:numel(ax)
    pos = get(ax(k),'Position');
    pos(2) = pos(2) - 0.05;
    set(ax(k),'Position',pos);
end

sgTitleStr = sprintf('%s: Choice distributions by condition (%s)', setName, figTag);
if SHOW_BOXPLOT
    if strcmpi(SUMMARY_PLOT_TYPE,'bar')
        sgTitleStr = sprintf('%s [bar+%s]', sgTitleStr, upper(BAR_ERROR_TYPE));
    elseif strcmpi(SUMMARY_PLOT_TYPE,'box')
        sgTitleStr = sprintf('%s [box]', sgTitleStr);
    end
end
sgtitle(sgTitleStr);

saveTag = figTag;
if SHOW_BOXPLOT && strcmpi(SUMMARY_PLOT_TYPE,'bar')
    saveTag = sprintf('%s_bar', figTag);
end

saveas(gcf, fullfile(OUTDIR, ...
    sprintf('choice_distribution_%s_%s_set_by_monkey.png', setName, saveTag)));

end  % make_choice_figure

function plot_sigline(xpair, y, starstr)
if ~contains(starstr, '*')
    return;
end
hold on;
line(xpair, [y y], 'Color','k', 'LineWidth', 1.2, 'HandleVisibility','off');
text(mean(xpair), y + 0.01*diff(ylim), starstr, ...
    'HorizontalAlignment','center', 'FontWeight','bold', 'Color','r');
end

function s = p2stars(p)
if p < 1e-3
    s = '***';
elseif p < 1e-2
    s = '**';
elseif p < 0.05
    s = '*';
elseif p < 0.01
    s = '.';
else
    s = 'n.s.';
end
end
