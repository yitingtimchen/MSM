% Tim 251201:
% While this does include more sessions than the ones listed in AM's
% recoveryFile, it doesn't produce meaningful analysis results. Don't use
% it...

% =========================================================================
% MATLAB PIPELINE — Behavior-Only from dataTabFULL Only
%                   Build per-condition packs with behavior matrices only.
%
% Version: behavior_only_from_dataTabFULL
% Author: Platt Lab (Tim Chen) — dataTabFULL-only variant
%
% WHAT THIS DOES (high level):
%   1) Reads dataTabFULL.xlsx (behavior table) once.
%   2) Uses only columns in dataTabFULL to:
%        - Define sessions (year, month, day, SessionOfDay).
%        - Map condition codes to labels (AI/Replay/Decoy/Live).
%        - Map OT codes to sub-conditions (none / Saline / OT).
%   3) For each session, builds a placeholder file-output entry and
%      attaches behavior variables from dataTabFULL.xlsx.
%   4) Produces a compact per-condition pack struct with:
%        C.cond, C.dt, C.files, C.eventTables,
%        C.eventNamesPerFile, C.commonEventNames, C.events3D.
%   5) Writes <COND>_condition_pack.mat + diary.
%
% CONDITION / OT MAPPING:
%   - Column "condition" is formatted like #.#, but only the digit before
%     the decimal matters:
%         1 = AI
%         2 = Replay
%         3 = Decoy
%         4 = Live
%
%   - Column "OT":
%         0 = (not OT or Saline)
%         1 = Saline
%         2 = OT
%
%   - Final condition labels:
%         OT=0:  'AI','Replay','Decoy','Live'
%         OT=1:  'Saline AI','Saline Replay','Saline Decoy','Saline Live'
%         OT=2:  'OT AI','OT Replay','OT Decoy','OT Live'
%
% MONKEY COLUMN:
%   - Column "monkey": 1=L, 2=B, discard 3 (rows dropped).
%   - Monkey letter is encoded in placeholder filenames but not used in
%     behavior selection.
%
% DEPENDENCIES:
%   - dataTabFULL.xlsx (behavior variables extraction; all info comes here)
% =========================================================================

%% ----------------------------- CONFIG ---------------------------------- %%
clear; clc;

behavior_path    = 'C:\Users\plattlab\Tim\Stock_market_experiment\MATLAB_Copy\dataTabFULL.xlsx';
root_path        = 'C:\Users\plattlab\Tim\Stock_market_experiment\MSM_Sorted\';
% OUTDIR = fullfile('C:\Users\plattlab\Tim\Stock_market_experiment\Tim', 'condition_packs_behavior_only_from_dataTabFULL'); % all sessions included
OUTDIR = fullfile('C:\Users\plattlab\Tim\Stock_market_experiment\Tim', 'condition_packs_behavior_only_from_dataTabFULL_v2'); % exclude odd dates

% Binning and trial structure (kept for compatibility with original packs)
BIN_DT        = 0.001 ;   % 1 ms
TR            = 15;       % trials per market
MKT           = 6;        % markets per session

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
cond2files    = containers.Map('KeyType','char','ValueType','any');
fileOutputs   = containers.Map('KeyType','char','ValueType','any');

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

%% -------------- (A) BUILD CONDITION→FILES MAP FROM DATATABFULL -------- %%
if ~isempty(behaviorTable)

    % Required columns
    reqVars = {'year','month','day','SessionOfDay', ...
               'condition','OT','monkey','session'};
    miss = setdiff(reqVars, behaviorTable.Properties.VariableNames);
    if ~isempty(miss)
        error('Behavior table missing required columns: %s', strjoin(miss,', '));
    end

    % Drop monkey==3 rows
    if any(behaviorTable.monkey == 3)
        badMask = behaviorTable.monkey == 3;
        warning('Dropping %d rows with monkey==3.', nnz(badMask));
        behaviorTable(badMask,:) = [];
    end

    nRows = height(behaviorTable);

    % Build per-row session keys: YYYY-MM-DD_SessionOfDay
    sessKey = strings(nRows,1);
    for i = 1:nRows
        sessKey(i) = sprintf('%04d-%02d-%02d_%d', ...
            behaviorTable.year(i), behaviorTable.month(i), ...
            behaviorTable.day(i), behaviorTable.SessionOfDay(i));
    end
    behaviorTable.sessionKey = sessKey;

    % Parse condition codes (#.# -> floor integer before decimal)
    condRaw = behaviorTable.condition;
    condCode = nan(nRows,1);
    if iscell(condRaw)
        for i = 1:nRows
            v = condRaw{i};
            if isempty(v)
                condCode(i) = NaN;
            else
                condCode(i) = floor(str2double(char(v)));
            end
        end
    elseif isstring(condRaw) || iscategorical(condRaw)
        for i = 1:nRows
            v = char(condRaw(i));
            condCode(i) = floor(str2double(v));
        end
    else
        condCode = floor(double(condRaw));
    end
    behaviorTable.condCode = condCode;

    % Parse OT codes (0/1/2)
    otRaw = behaviorTable.OT;
    otCode = nan(nRows,1);
    if iscell(otRaw)
        for i = 1:nRows
            v = otRaw{i};
            if isempty(v)
                otCode(i) = NaN;
            else
                otCode(i) = str2double(char(v));
            end
        end
    elseif isstring(otRaw) || iscategorical(otRaw)
        for i = 1:nRows
            v = char(otRaw(i));
            otCode(i) = str2double(v);
        end
    else
        otCode = double(otRaw);
    end
    behaviorTable.OTcode = otCode;

    % Build map: condition label -> list of placeholder "nex" paths
    [sessKeys,~,sessIdx] = unique(behaviorTable.sessionKey);

    for s = 1:numel(sessKeys)
        rowsS = (sessIdx == s);
        if ~any(rowsS), continue; end

        % Condition code per session
        condCodes = unique(behaviorTable.condCode(rowsS));
        condCodes(isnan(condCodes)) = [];
        if isempty(condCodes)
            warning('Session %s has no valid condition code; skipping.', sessKeys(s));
            continue;
        end
        if numel(condCodes) > 1
            warning('Session %s has multiple condition codes: %s; using first.', ...
                sessKeys(s), mat2str(condCodes));
        end
        condCodeS = condCodes(1);

        % OT code per session
        otCodes = unique(behaviorTable.OTcode(rowsS));
        otCodes(isnan(otCodes)) = [];
        if isempty(otCodes)
            warning('Session %s has no OT code; skipping.', sessKeys(s));
            continue;
        end
        if numel(otCodes) > 1
            warning('Session %s has multiple OT codes: %s; using first.', ...
                sessKeys(s), mat2str(otCodes));
        end
        otCodeS = otCodes(1);

        % Map to condition label
        condLabel = mapConditionLabel(condCodeS, otCodeS);
        if isempty(condLabel)
            warning('Session %s has unmapped condition/OT (cond=%g, OT=%g); skipping.', ...
                sessKeys(s), condCodeS, otCodeS);
            continue;
        end

        % Representative row for this session
        idxFirst = find(rowsS, 1, 'first');
        yearK      = behaviorTable.year(idxFirst);
        monthK     = behaviorTable.month(idxFirst);
        dayK       = behaviorTable.day(idxFirst);
        sessOfDayK = behaviorTable.SessionOfDay(idxFirst);
        monkeyK    = behaviorTable.monkey(idxFirst);

        switch monkeyK
            case 1
                mChar = 'B';
            case 2
                mChar = 'L';
            otherwise
                mChar = 'X';
        end

        yy2 = mod(yearK, 100);
        baseName = sprintf('%02d%02d%02d%s%d', monthK, dayK, yy2, mChar, sessOfDayK);
        placeholderPath = fullfile(root_path, sprintf('%s_behaviorOnly.nex', baseName));

        if ~isKey(cond2files, condLabel)
            cond2files(condLabel) = {placeholderPath};
        else
            existing = cond2files(condLabel);
            existing{end+1,1} = placeholderPath;
            cond2files(condLabel) = existing;
        end
    end

else
    warning('Behavior table is empty; no condition→files mapping built.');
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

        % Cache outputs (behavior-only)
        out = struct();
        out.binWidthSec     = BIN_DT;      % keep for compatibility
        out.eventIndexTable = evTable;     % behavior variables only

        fileOutputs(nexPath) = out;
    end
end

%% ---------------- (C) SAVE PER-CONDITION PACKS TO DISK ----------------- %%
condNames = cond2files.keys;
for c = 1:numel(condNames)
    cond = condNames{c};   % char

    files = cond2files(cond);
    Fcount = numel(files);

    evTables   = cell(Fcount,1);
    evVarNames = cell(Fcount,1);
    dt = [];

    for f = 1:Fcount
        out = fileOutputs(files{f});
        evTables{f} = out.eventIndexTable;
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
    C.cond              = cond;
    C.dt                = dt;
    C.files             = files(:);
    C.eventTables       = evTables;
    C.eventNamesPerFile = evVarNames;
    C.commonEventNames  = common;
    C.events3D          = events3D;

    save(fullfile(OUTDIR, sprintf('%s_condition_pack.mat', cond)), 'C', '-v7.3');
end

fprintf('Done. Files written to: %s\n', OUTDIR);
try, diary('off'); catch, end

%% ============================== HELPERS ================================ %%
function evTableOut = addBehaviorVarsToEventTable(evTableIn, nexPath, Tfull, TR, MKT)
% Attach behavior matrices only. Assumes per-session identification from
% placeholder filename: <mmddyy><B|L><SessionOfDay>_behaviorOnly.nex
    evTableOut = evTableIn;
    if isempty(Tfull), return, end

    [~, base, ~] = fileparts(nexPath);
    tok = regexp(base,'(\d{6})([BLX])(\d+)','tokens','once');
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

function label = mapConditionLabel(condCode, otCode)
% Map numeric condition/OT codes to string labels.
%   condition (integer part of "condition" column):
%       1 = AI
%       2 = Replay
%       3 = Decoy
%       4 = Live
%
%   OT:
%       0 = none
%       1 = Saline
%       2 = OT

    base = '';
    switch condCode
        case 1
            base = 'AI';
        case 2
            base = 'Replay';
        case 3
            base = 'Decoy';
        case 4
            base = 'Live';
        otherwise
            label = '';
            return;
    end

    switch otCode
        case 0
            label = base;
        case 1
            label = ['Saline ' base];
        case 2
            label = ['OT ' base];
        otherwise
            label = '';
    end
end
