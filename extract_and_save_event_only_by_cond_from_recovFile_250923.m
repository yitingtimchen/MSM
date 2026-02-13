
% =========================================================================
% MATLAB PIPELINE — Behavior-Only (Ignores .nex Completely)
%                   Build per-condition packs with behavior matrices only.
%
% Version: behavior_only_from_recoveryFile
% Author: Platt Lab (Tim Chen) — minimal .nex-free variant
%
% WHAT THIS DOES (high level):
%   1) Scans condition sheets (recoveryFile) to list sessions.
%   2) Skips .nex entirely — no reading, no neurons, no events.
%   3) For each session, builds an empty placeholder file-output struct but
%      ATTACHES behavior variables from dataTabFULL.xlsx.
%   4) Keeps per-condition pack structure identical to original script:
%        C.cond, C.dt, C.files, C.rawSpikes, C.S, C.neuronNames, C.Tbins,
%        C.trialStartEvent, C.eventTables, C.eventNamesPerFile, C.commonEventNames,
%        C.events3D  (numeric behavior vars only; text vars skipped).
%   5) Writes <COND>_condition_pack.mat + diary + event_name_variations.tsv
%
% IMPORTANT CHANGES:
%   - Entire Part (B) "PER-FILE PROCESSING: SPIKES & EVENTS" is replaced with
%     a behavior-only block that creates empty placeholders and attaches behavior.
%   - No neuron matching; FILTER_TO_SHEET_NEURONS is ignored.
%   - No .nex path existence check; placeholder filenames are synthesized.
%
% OUTPUTS:
%   OUTDIR/
%     - <COND>_condition_pack.mat
%     - pipeline_diary_<timestamp>.txt
%     - event_name_variations.tsv  (empty/headers-only in this variant)
%
% DEPENDENCIES:
%   - dataTabFULL.xlsx (behavior variables extraction)
%   - recoveryFile spreadsheet: MSM_neuronSortMASTER_recoveryFile.xlsx
% =========================================================================

%% ----------------------------- CONFIG ---------------------------------- %%
clear; clc;

is_mac = 0;

% Conditions and inputs (char cell, not string array)
conditions       = {'AI','Replay','Decoy','Live', ...
                    'OT AI', 'OT Replay', 'OT Decoy', 'OT Live', ...
                    'Saline AI', 'Saline Replay', 'Saline Decoy', 'Saline Live'};
% conditions       = {'AI','Replay','Decoy','Live'}; % test

spreadsheet_path = 'C:\Users\plattlab\Tim\Stock_market_experiment\MSM_neuronSortMASTER_recoveryFile_modified_filtered.xlsx';
behavior_path    = 'C:\Users\plattlab\Tim\Stock_market_experiment\MATLAB_Copy\dataTabFULL.xlsx';
root_path        = 'C:\Users\plattlab\Tim\Stock_market_experiment\MSM_Sorted\';
OUTDIR = fullfile('C:\Users\plattlab\Tim\Stock_market_experiment\Tim', 'condition_packs_behavior_only_v2');

% Binning and trial structure (kept for compatibility with original packs)
BIN_DT        = 0.001 ;   % 1 ms
TR            = 15;       % trials per market
MKT           = 6;        % markets per session
TRIALS_TOTAL  = TR*MKT;   % 90

% Output directory
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

% Diary (log) file
diaryFile = fullfile(OUTDIR, sprintf('pipeline_diary_%s.txt', datestr(now,'yyyymmdd_HHMMSS')));
try
    diary(diaryFile);
    cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>
    fprintf('Diary started: %s\n', diaryFile);
catch ME
    warning('Failed to start diary: %s', ME.message);
end

%% -------------------- OUTPUT CONTAINERS (IN-MEMORY) -------------------- %%
cond2files          = containers.Map('KeyType','char','ValueType','any');
fileOutputs         = containers.Map('KeyType','char','ValueType','any');
eventVariantTally   = containers.Map('KeyType','char','ValueType','any'); % will remain empty

% Load behavior table once (required)
behaviorTable = [];
if exist(behavior_path, 'file')
    try
        behaviorTable = readtable(behavior_path, 'Sheet','Sheet1');
    catch ME
        warning('Could not read behavior table: %s', ME.message);
    end
else
    warning('Behavior table path does not exist; behavior variables will be missing.');
end

%% --------------------- (A) BUILD CONDITION→FILES MAP ------------------- %%
for c = 1:numel(conditions)
    cond = conditions{c};   % char

    % Read this condition's sheet (drop rows missing key fields)
    Tsheet = readtable(spreadsheet_path, 'Sheet', cond);
    Tsheet = Tsheet(~any(ismissing(Tsheet(:, {'Date','Session','Monkey'})), 2), :);

    % get 15 sessions/90 blocks per condition
    % Tsheet = Tsheet(Tsheet.Monkey == 2, :);
    % Tsheet = sortrows(Tsheet, {'Monkey', 'Date'}, {'ascend', 'descend'});
    % Tsheet = Tsheet(1:min(15, size(Tsheet, 1)), :);

    % Build unique per-session keys
    sessionKeys = {};  % char cell
    for i = 1:height(Tsheet)
        dt = normalizeDate(Tsheet.Date(i));
        date_str = datestr(dt, 'mmddyy');
        monkey_map = {'B','L'};
        if Tsheet.Monkey(i) < 1 || Tsheet.Monkey(i) > 2
            warning('Unknown Monkey code at row %d; skipping row.', i);
            continue
        end
        mchar = monkey_map{Tsheet.Monkey(i)};
        sess  = num2str(Tsheet.Session(i));
        key   = sprintf('%s%s%s', upper(date_str), mchar, sess);
        sessionKeys{end+1,1} = key; %#ok<AGROW>
    end
    sessionKeys = unique(sessionKeys);

    % Synthesize placeholder ".nex" paths (they need only to encode the key)
    selected = cell(numel(sessionKeys),1);
    for i = 1:numel(sessionKeys)
        selected{i} = fullfile(root_path, sprintf('%s_behaviorOnly.nex', sessionKeys{i})); %#ok<AGROW>
    end

    cond2files(cond) = selected;
end

%% -------------- (B) BEHAVIOR-ONLY PER-FILE PROCESSING ------------------ %%
allConds = cond2files.keys;
for c = 1:numel(allConds)
    cond = allConds{c};   % char key
    fileList = cond2files(cond);

    for f = 1:numel(fileList)
        nexPath = fileList{f};
        if isKey(fileOutputs, nexPath), continue; end % already processed

        fprintf('Processing (behavior-only): %s\n', nexPath);

        % Build behavior-only event table
        evTable = table();
        if ~isempty(behaviorTable)
            evTable = addBehaviorVarsToEventTable(evTable, nexPath, behaviorTable, TR, MKT);
        end

        % Cache outputs (empty placeholders for spikes/events)
        out = struct();
        out.neuronNames     = {};          % empty
        out.rawSpikes       = {};          % empty
        out.spikeBinMatrix  = [];          % empty
        out.binWidthSec     = BIN_DT;      % keep for compatibility
        out.T_bins          = 1;           % minimal placeholder
        out.trialStartEvent = '';          % none
        out.eventIndexTable = evTable;     % behavior variables only
        out.alignedInfo     = struct();    % none

        fileOutputs(nexPath) = out;
    end
end

%% ---------------- (C) SAVE PER-CONDITION PACKS TO DISK ----------------- %%
condNames = cond2files.keys;
for c = 1:numel(condNames)
    cond = condNames{c};   % char

    files = cond2files(cond);
    Fcount = numel(files);

    S_cells = cell(Fcount,1);
    rawSpikes = cell(Fcount,1);
    names   = cell(Fcount,1);
    Tbins   = zeros(Fcount,1);
    trialStartEvent = cell(Fcount,1);
    evTables = cell(Fcount,1);
    evVarNames = cell(Fcount,1);
    dt = [];

    for f = 1:Fcount
        out = fileOutputs(files{f});
        S_cells{f}         = out.spikeBinMatrix;
        rawSpikes{f}       = out.rawSpikes;
        names{f}           = out.neuronNames;
        Tbins(f)           = out.T_bins;
        trialStartEvent{f} = out.trialStartEvent;
        evTables{f}        = out.eventIndexTable;
        if ~isempty(evTables{f})
            evVarNames{f} = evTables{f}.Properties.VariableNames;
        else
            evVarNames{f} = {};
        end
        if isempty(dt), dt = out.binWidthSec; end
    end

    if Fcount > 0
        common = evVarNames{1};
        for f = 2:Fcount, common = intersect(common, evVarNames{f}, 'stable'); end
    else
        common = {};
    end

    % Build events3D only for NUMERIC behavior variables (text vars skipped)
    events3D = struct();
    if ~isempty(common) && Fcount > 0 && ~isempty(evTables{1})
        for e = 1:numel(common)
            X0 = evTables{1}.(common{e}){1};
            if isnumeric(X0)
                [TR_ok, MKT_ok] = size(X0);
                arr = nan(TR_ok, MKT_ok, Fcount);
                for f = 1:Fcount
                    arr(:,:,f) = evTables{f}.(common{e}){1};
                end
                events3D.(common{e}) = arr;
            else
                % skip non-numeric (likely text/categorical) variables
            end
        end
    end

    C = struct();
    C.cond               = cond;
    C.dt                 = dt;
    C.files              = files(:);
    C.rawSpikes          = rawSpikes;
    C.S                  = S_cells;
    C.neuronNames        = names;
    C.Tbins              = Tbins(:);
    C.trialStartEvent    = trialStartEvent(:);
    C.eventTables        = evTables;
    C.eventNamesPerFile  = evVarNames;
    C.commonEventNames   = common;
    C.events3D           = events3D;

    save(fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cond)), 'C', '-v7.3');
end

fprintf('Done. Files written to: %s\n', OUTDIR);
variants = finalizeEventVariantReport(eventVariantTally, OUTDIR); %#ok<NASGU>
try, diary('off'); catch, end

%% filtered behavior tables
filteredBehaviorPath = writeFilteredBehaviorTable(spreadsheet_path, behavior_path, conditions, OUTDIR);

%% ============================== HELPERS ================================ %%
function dt = normalizeDate(raw)
    if isdatetime(raw)
        dt = raw;
    elseif isnumeric(raw) && isscalar(raw)
        dt = datetime(raw, 'ConvertFrom','excel');
    else
        try
            dt = datetime(raw, 'InputFormat','MM/dd/yy');
        catch
            error('Unable to parse date: %s', string(raw));
        end
    end
end

function evTableOut = addBehaviorVarsToEventTable(evTableIn, nexPath, Tfull, TR, MKT)
% Attach behavior matrices only. Assumes per-session identification from
% placeholder filename: <mmddyy><B|L><session>_behaviorOnly.nex
    evTableOut = evTableIn;
    if isempty(Tfull), return, end

    % reqVars = {'year','month','day','SessionOfDay','session', ...
    %            'option','optimal','marketOrig','tResp','bubbleMarket', ...
    %            'monkey','postPortfolio','divPerShare','sizeB'};
    % missing = setdiff(reqVars, Tfull.Properties.VariableNames);
    % if ~isempty(missing)
    %     warning('Behavior table missing variables: %s', strjoin(missing,', '));
    %     return
    % end

    [~, base, ~] = fileparts(nexPath);
    tok = regexp(base,'(\d{6})([BL])(\d+)','tokens','once');
    if isempty(tok), warning('Cannot parse date/session from placeholder: %s', base); return, end
    mm = str2double(tok{1}(1:2));
    dd = str2double(tok{1}(3:4));
    yy = str2double(tok{1}(5:6));
    sessOfDay = str2double(tok{3});
    yyyy = 2000 + yy; if yy > 69, yyyy = 1900 + yy; end

    rows = Tfull.year==yyyy & Tfull.month==mm & Tfull.day==dd & Tfull.SessionOfDay==sessOfDay;
    Tsel = Tfull(rows,:);
    nrows = height(Tsel);
    if nrows ~= 180, warning('Expected 180 behavior rows for %s, got %d.', base, nrows); end
    if nrows == 0, return, end

    sessions = unique(Tsel.session, 'sorted');
    if numel(sessions) ~= MKT, sessions = unique(Tsel.session, 'stable'); end

    % varsToExtract = {'session','option','optimal','marketOrig','tResp','bubbleMarket','monkey','postPortfolio','divPerShare','sizeB'};
    varsToExtract = Tsel.Properties.VariableNames;
    for v = 1:numel(varsToExtract)
        vn = varsToExtract{v};
        col = Tsel.(vn);
        isTextLike = iscellstr(col) || isstring(col) || iscategorical(col);
        if isTextLike 
            M = cell(TR, MKT); 
            Mopp = cell(TR, MKT);
        else
            M = nan(TR, MKT); 
            Mopp = nan(TR, MKT);
        end
        for si = 1:min(MKT, numel(sessions))
            s = sessions(si);
            blk = Tsel(Tsel.session == s, :);
            vals = blk.(vn);
            if iscategorical(vals), vals = cellstr(vals); end
            if isstring(vals),     vals = cellstr(vals);  end
            if iscell(vals)
                vodd = vals(1:2:end);
                M(1:min(TR,numel(vodd)), si) = vodd(1:min(TR,numel(vodd)));
                veven = vals(2:2:end);
                Mopp(1:min(TR,numel(veven)), si) = veven(1:min(TR,numel(veven)));
            else
                vodd = double(vals(1:2:end));
                M(1:min(TR,numel(vodd)), si) = vodd(1:min(TR,numel(vodd)));
                veven = double(vals(2:2:end));
                Mopp(1:min(TR,numel(veven)), si) = veven(1:min(TR,numel(veven)));
            end
        end
        evTableOut.(vn) = {M};
        evTableOut.([vn 'Opp']) = {Mopp};
    end
end

function variants = finalizeEventVariantReport(tally, OUTDIR)
% Minimal no-op report for compatibility
    bases = tally.keys; %#ok<NASGU>
    try
        T = cell2table(cell(0,4), 'VariableNames', {'base','nVariants','variant','count'});
        outf = fullfile(OUTDIR, 'event_name_variations.tsv');
        writetable(T, outf, 'FileType','text','Delimiter','\t');
        fprintf('\nWrote empty variation table: %s\n', outf);
    catch ME
        fprintf('Could not write TSV: %s\n', ME.message);
    end
    variants = struct('base',{},'nVariants',{},'variants',{},'counts',{});
end

function outPath = writeFilteredBehaviorTable(spreadsheet_path, behavior_path, conditions, OUTDIR)
if ~exist(behavior_path,'file')
error('Behavior table not found: %s', behavior_path);
end
Tfull = readtable(behavior_path,'Sheet','Sheet1');

needCols = {'year','month','day','SessionOfDay'};
miss = setdiff(needCols, Tfull.Properties.VariableNames);
if ~isempty(miss)
    error('Behavior table missing required columns: %s', strjoin(miss,', '));
end

% Collect session keys [YYYY-MM-DD_SessionOfDay] from recoveryFile sheets
sessKeys = strings(0,1);
for c = 1:numel(conditions)
    cond = conditions{c};
    Tsheet = readtable(spreadsheet_path, 'Sheet', cond);
    Tsheet = Tsheet(~any(ismissing(Tsheet(:, {'Date','Session','Monkey'})), 2), :);
    for i = 1:height(Tsheet)
        dt = normalizeDateRF(Tsheet.Date(i));
        sess = Tsheet.Session(i);
        sessKeys(end+1,1) = sprintf('%04d-%02d-%02d_%d', year(dt), month(dt), day(dt), sess); %#ok<AGROW>
    end
end
sessKeys = unique(sessKeys);

% Build keys for behavior rows
behKeys = strings(height(Tfull),1);
for i = 1:height(Tfull)
    behKeys(i) = sprintf('%04d-%02d-%02d_%d', Tfull.year(i), Tfull.month(i), Tfull.day(i), Tfull.SessionOfDay(i));
end

mask = ismember(behKeys, sessKeys);
Tfiltered = Tfull(mask, :);

% Write out
outPath = fullfile(OUTDIR, 'dataTabFULL_filtered_from_recoveryFile.xlsx');
writetable(Tfiltered, outPath, 'Sheet','Sheet1');

% Console summary
nSess = numel(sessKeys);
expected = nSess * 180; % expected rows per session (90 trials * 2 rows)
fprintf('Filtered %d/%d rows across %d sessions (expected ~%d).\n', height(Tfiltered), height(Tfull), nSess, expected);
fprintf('Wrote: %s\n', outPath);


end

function dt = normalizeDateRF(raw)
if isdatetime(raw)
dt = raw;
elseif isnumeric(raw) && isscalar(raw)
dt = datetime(raw, 'ConvertFrom','excel');
else
try
dt = datetime(raw, 'InputFormat','MM/dd/yy');
catch
error('Unable to parse date: %s', string(raw));
end
end
end