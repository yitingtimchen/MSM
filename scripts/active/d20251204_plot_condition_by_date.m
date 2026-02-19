% File: d20251204_plot_condition_by_date.m
% plot_condition_by_date_per_monkey.m
% Load per-condition packs created by behavior_only_from_dataTabFULL
% and plot condition vs. date separately for each monkey (L and B),
% with dot size indicating how many sessions are on that date/condition.

clear; clc;

% ---- Paths ----
outdir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_behavior_only_from_dataTabFULL';
% outdir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_251118'; % neural branch

if ~isfolder(outdir)
    error('Output directory not found: %s', outdir);
end

files = dir(fullfile(outdir, '*_condition_pack.mat'));
if isempty(files)
    error('No *_condition_pack.mat files found in: %s', outdir);
end

% ---- Accumulators for session-level info ----
allDates   = datetime.empty(0,1);
allMonkeys = strings(0,1);
allConds   = strings(0,1);

for k = 1:numel(files)
    packPath = fullfile(outdir, files(k).name);
    S = load(packPath, 'C');
    if ~isfield(S, 'C')
        warning('File %s does not contain struct C; skipping.', files(k).name);
        continue;
    end
    C = S.C;

    if ~isfield(C, 'files') || ~isfield(C, 'cond')
        warning('Struct C in %s missing required fields; skipping.', files(k).name);
        continue;
    end

    condLabel = string(C.cond);
    fileList  = C.files(:);

    for f = 1:numel(fileList)
        thisFile = fileList{f};
        [~, base, ~] = fileparts(thisFile);

        % Expect base like: mmddyyL1_behaviorOnly
        tok = regexp(base, '(\d{6})([BLX])(\d+)', 'tokens', 'once');
        if isempty(tok)
            warning('Could not parse placeholder filename: %s', base);
            continue;
        end

        mmddyy = tok{1};
        mChar  = tok{2};
        % sessOfDay = str2double(tok{3}); %#ok<NASGU>

        % Only keep L/B; ignore X or anything else
        if ~(mChar == 'L' || mChar == 'B')
            continue;
        end

        mm = str2double(mmddyy(1:2));
        dd = str2double(mmddyy(3:4));
        yy = str2double(mmddyy(5:6));

        if yy > 69
            yyyy = 1900 + yy;
        else
            yyyy = 2000 + yy;
        end

        try
            d = datetime(yyyy, mm, dd);
        catch
            warning('Invalid date parsed from %s; skipping.', base);
            continue;
        end

        allDates(end+1,1)   = d;              %#ok<SAGROW>
        allMonkeys(end+1,1) = string(mChar);  %#ok<SAGROW>
        allConds(end+1,1)   = condLabel;      %#ok<SAGROW>
    end
end

if isempty(allDates)
    error('No valid sessions found from condition pack files.');
end

% ---- Build ordered condition categories ----
allCondLevels = { ...
    'AI','Replay','Decoy','Live', ...
    'Saline AI','Saline Replay','Saline Decoy','Saline Live', ...
    'OT AI','OT Replay','OT Decoy','OT Live'};

present = ismember(allCondLevels, unique(allConds));
condLevels = allCondLevels(present);

condCat = categorical(allConds, condLevels);
condIdx = double(condCat); %#ok<NASGU> % kept for reference if needed

% ---- Plot per monkey ----
monkeys = unique(allMonkeys);
monkeys = sort(monkeys);  % e.g. 'B','L' -> 'B','L'

figure('Name','Condition vs Date per Monkey','Color','w');
nM = numel(monkeys);

for iM = 1:nM
    mID  = monkeys(iM);
    mask = allMonkeys == mID;

    dM     = allDates(mask);
    cM     = condCat(mask);
    cIdxM  = double(cM);

    if isempty(dM)
        continue;
    end

    % Group by (date, condition) and count sessions
    X = [datenum(dM), cIdxM];
    [uniqRows, ~, ic] = unique(X, 'rows');
    counts = accumarray(ic, 1);

    uDates   = datetime(uniqRows(:,1), 'ConvertFrom', 'datenum');
    uCondIdx = uniqRows(:,2);

    % Marker size encodes number of sessions (1,2,3,...) on that date/condition
    baseSize = 20;   % area for 1 session
    scale    = 20;   % additional area per extra session
    sizes    = baseSize + scale*(counts - 1);

    subplot(nM, 1, iM);
    scatter(uDates, uCondIdx, sizes, 'filled');
    grid on;

    set(gca, ...
        'YTick', 1:numel(condLevels), ...
        'YTickLabel', condLevels, ...
        'YLim', [0.5, numel(condLevels)+0.5]);
    datetick(gca, 'x', 'yyyy-mm-dd', 'keeplimits', 'keepticks');

    title(sprintf('Monkey %s', mID));
end

% Link x-axes across subplots for easier comparison
ax = findall(gcf, 'Type', 'axes');
linkaxes(ax, 'x');
