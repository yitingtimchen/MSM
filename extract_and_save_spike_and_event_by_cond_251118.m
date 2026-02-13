% =========================================================================
% MATLAB PIPELINE — Build Condition→File Map, Create 1 ms Spike Matrices,
%                   Global Tru ncation @ m1trialStart(1), ON/OFF Pair Validation
%                   (bidirectional OFF shifting), Build Choice from aligned ON,
%                   Size Consistency Report, Diary Logs, Save Per-Condition Packs
%
% Version: aligned_pairs_from_off_on + global_truncation_m1_trialstart_first
% Author: Platt Lab (Tim Chen) — updated per clarified requirements
%
% WHAT THIS DOES (high level):
%   1) Scans condition sheets to select the best .nex file per session.
%   2) Reads neurons & events from each .nex and builds 1 ms spike
%   matrices. NO FILTERING NEURONS.
%   3) Normalizes event names (keep only m1*/m2*, unify case and aliases).
%   4) *** GLOBAL TRUNCATION *** Drop ALL event timestamps occurring earlier
%      than the FIRST entry of m1trialStart. (Your clarified rule.)
%   5) For each reward family (rew1, rew2) and monkey (m1, m2):
%      - Validate ON/OFF with the rule 0 < (OFF - ON) < 12 seconds.
%      - If lengths differ or pairs invalid, attempt to re-align by
%        shifting OFF indices forward/backward (relative to pair index)
%        to find valid OFF for each ON. Keep only valid (ON,OFF) pairs.
%      - Truncate both ON and OFF afterward so the two vectors have
%        the SAME length (dangling stamps dropped).
%   6) Build derived variables:
%      - m[12]choice := the aligned m[12]rew1on timestamps - 250 ms.
%        (Uses ONLY ONs that have a valid OFF pair after alignment.)
%      - m[12]rew1time := aligned (OFF - ON) durations for rew1 (in seconds).
%      - m[12]rew2time := aligned (OFF - ON) durations for rew2 (in seconds).
%   7) Convert timestamp events to TR×MKT (15×6) bin-index matrices
%      (1-based bins using BIN_DT) per trial window (this monkey's trialStart(i)
%      to the next same-monkey trialStart(i+1)). Durations are stored as raw
%      seconds in TR×MKT matrices.
%   8) After building the table, check that all variables have the same size
%      and report inconsistent finite counts.
%   9) Aggregates per-condition packs and writes .mat files + logs.
%
% OUTPUTS:
%   OUTDIR/
%     - <COND>_condition_pack.mat
%     - pipeline_diary_<timestamp>.txt
%     - event_name_variations.tsv
%
% DEPENDENCIES:
%   - readNexFile: returns struct F with fields:
%       F.neurons{i}.name, F.neurons{i}.timestamps
%       F.events{i}.name,  F.events{i}.timestamps
%   - dataTabFULL.xlsx (optional; behavior variables extraction)
% =========================================================================

%% ----------------------------- CONFIG ---------------------------------- %%
clear; clc;

is_mac = 0;

% Conditions and inputs (char cell, not string array)
conditions       = {'AI','Replay','Decoy','Live', ...
                    'OT AI', 'OT Replay', 'OT Decoy', 'OT Live', ...
                    'Saline AI', 'Saline Replay', 'Saline Decoy', 'Saline Live'};
% conditions       = {'Replay'}; % test

spreadsheet_path = 'C:\Users\plattlab\Tim\Stock_market_experiment\MSM_neuronSortMASTER_recoveryFile.xlsx';
behavior_path    = 'C:\Users\plattlab\Tim\Stock_market_experiment\MATLAB_Copy\dataTabFULL.xlsx';
root_path        = 'C:\Users\plattlab\Tim\Stock_market_experiment\MSM_Sorted\';
OUTDIR = fullfile('C:\Users\plattlab\Tim\Stock_market_experiment\Tim', 'condition_packs_251118');

% Binning and trial structure
BIN_DT        = 0.001 ;   % 1 ms
TR            = 15;       % trials per market
MKT           = 6;        % markets per session
TRIALS_TOTAL  = TR*MKT;   % 90

% Options
USE_SPARSE                = false;  % set true to save memory on long sessions
FILTER_TO_SHEET_NEURONS   = true;   % include only neurons listed in the spreadsheet

% Trial-start candidates (first present will be used) AFTER canonicalization
TRIAL_START_CANDIDATES = ["m1trialStart","m1start","m1trialON"];

% Output directory
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

% Diary (log) file
diaryFile = fullfile(OUTDIR, sprintf('pipeline_diary_%s.txt', datestr(now,'yyyymmdd_HHMMSS')));
try
    diary(diaryFile);
    cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU> ensure diary closes on exit
    fprintf('Diary started: %s\n', diaryFile);
catch ME
    warning('Failed to start diary: %s', ME.message);
end

%% -------------------- OUTPUT CONTAINERS (IN-MEMORY) -------------------- %%
cond2files          = containers.Map('KeyType','char','ValueType','any');
fileOutputs         = containers.Map('KeyType','char','ValueType','any');
eventVariantTally   = containers.Map('KeyType','char','ValueType','any');

% Load behavior table once (optional if you don't have it)
behaviorTable = [];
if exist(behavior_path, 'file')
    try
        behaviorTable = readtable(behavior_path, 'Sheet','Sheet1');
    catch ME
        warning('Could not read behavior table: %s', ME.message);
    end
else
    warning('Behavior table path does not exist; skipping behavior attachment.');
end

%% --------------------- (A) BUILD CONDITION→FILES MAP ------------------- %%
for c = 1:numel(conditions)
    cond = conditions{c};   % char

    % Read this condition's sheet (drop rows missing key fields)
    Tsheet = readtable(spreadsheet_path, 'Sheet', cond);
    Tsheet = Tsheet(~any(ismissing(Tsheet(:, {'Date','Session','Monkey'})), 2), :);

    % Build per-session neuron allowlist
    perSessionNeuronSet = containers.Map('KeyType','char','ValueType','any');

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

        if FILTER_TO_SHEET_NEURONS && ismember('NeuronID', Tsheet.Properties.VariableNames)
            try
                nid = Tsheet.NeuronID{i};
                if ~isempty(nid)
                    nm = normalizeNeuronID(nid);
                    if isKey(perSessionNeuronSet, key)
                        perSessionNeuronSet(key) = unique([perSessionNeuronSet(key); {nm}]);
                    else
                        perSessionNeuronSet(key) = unique({nm});
                    end
                end
            catch
                % ignore malformed NeuronID; continue
            end
        end
        sessionKeys{end+1,1} = key; %#ok<AGROW>
    end
    sessionKeys = unique(sessionKeys);

    % Group .nex files by session key
    d = dir(fullfile(root_path, '*.nex'));
    sessMap = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(d)
        fname = d(k).name;
        tok = regexp(fname,'(\d{6})([BL])(\d+)','tokens','once');
        if isempty(tok), continue; end
        thisKey = sprintf('%s%s%s', tok{1}, tok{2}, tok{3});
        if ~ismember(thisKey, sessionKeys), continue; end
        if ~isKey(sessMap, thisKey), sessMap(thisKey) = {fname};
        else, sessMap(thisKey) = [sessMap(thisKey), {fname}];
        end
    end

    % Select best file per session
    selected = {};
    for i = 1:numel(sessionKeys)
        key = sessionKeys{i};
        if isKey(sessMap, key)
            bestName = selectBestFile(sessMap(key));
            selected{end+1} = fullfile(root_path, bestName); %#ok<SAGROW>
        end
    end

    cond2files(cond) = selected;
    if FILTER_TO_SHEET_NEURONS, cond2files([cond '_NeuronSets']) = perSessionNeuronSet; end
end

%% -------------- (B) PER-FILE PROCESSING: SPIKES & EVENTS --------------- %%
allConds = cond2files.keys;
for c = 1:numel(allConds)
    cond = allConds{c};   % char key
    if endsWith(cond,'_NeuronSets'), continue; end

    fileList = cond2files(cond);
    if FILTER_TO_SHEET_NEURONS
        perSessionNeuronSet = cond2files([cond '_NeuronSets']);
    else
        perSessionNeuronSet = [];
    end

    for f = 1:numel(fileList)
        nexPath = fileList{f};
        if isKey(fileOutputs, nexPath), continue; end % already processed

        fprintf('Processing: %s\n', nexPath);

        % Read .nex content
        F = readNexFile(nexPath);

        % Normalize containers
        [neuronNamesRaw, neuronSpikesRaw] = normalizeNeurons(F);
        [eventNames, eventTimes]          = normalizeEvents(F);
        updateEventVariantTally(eventVariantTally, eventNames);

        % Optional neuron filtering to sheet allowlist
        if FILTER_TO_SHEET_NEURONS && ~isempty(perSessionNeuronSet)
            sessKey = parseSessionKeyFromFilename(nexPath);
            if isKey(perSessionNeuronSet, sessKey)
                allowed = perSessionNeuronSet(sessKey);
                keep = ismember(neuronNamesRaw, allowed);
                neuronNames  = neuronNamesRaw(keep);
                neuronSpikes = neuronSpikesRaw(keep);
            else
                neuronNames  = neuronNamesRaw;
                neuronSpikes = neuronSpikesRaw;
            end
        else
            neuronNames  = neuronNamesRaw;
            neuronSpikes = neuronSpikesRaw;
        end

        % Determine total duration from spikes/events
        vals_spk = cellfun(@(v) ifemptymax(v), neuronSpikes);
        tmax_spk = max([0; vals_spk(:)]);
        vals_ev  = cellfun(@(v) ifemptymax(v), eventTimes);
        tmax_ev  = max([0; vals_ev(:)]);
        tmax = max(tmax_spk, tmax_ev);
        if ~isfinite(tmax) || tmax <= 0, tmax = 1; end

        % Build N×T binary spike matrix
        [Sbin, Tbins] = buildSpikeBinMatrix(neuronSpikes, BIN_DT, tmax, USE_SPARSE);

        % === GLOBAL TRUNCATION: drop ALL events earlier than FIRST m1trialStart ===
        [eventTimes] = truncateAllEventsBeforeFirstM1Start(eventNames, eventTimes, nexPath);

        % Identify trial-starts (we'll still need per-window bins)
        [trialStartName, trialStarts] = findTrialStart(eventNames, eventTimes, TRIAL_START_CANDIDATES);

        % === ALIGN ON/OFF (after truncation) & BUILD DERIVED ===
        [eventNames, eventTimes, alignedInfo] = alignAllRewardPairs(eventNames, eventTimes, nexPath);

        % Build event table (timestamps->bins, durations raw seconds)
        if isempty(trialStartName)
            warning('No trial-start event found in %s; creating empty event table.', nexPath);
            evTable = table();
        else
            evTable = buildEventIndexTable_fromAligned(eventNames, eventTimes, trialStarts, BIN_DT, TR, MKT, Tbins, nexPath);
            checkEventTableConsistency(evTable, TR, MKT, nexPath);
        end

        % Attach behavior variables (optional)
        if ~isempty(behaviorTable)
            evTable = addBehaviorVarsToEventTable(evTable, nexPath, behaviorTable, TR, MKT);
        end

        % Cache outputs
        out = struct();
        out.neuronNames     = neuronNames;
        out.rawSpikes       = neuronSpikes;
        out.spikeBinMatrix  = Sbin;
        out.binWidthSec     = BIN_DT;
        out.T_bins          = Tbins;
        out.trialStartEvent = trialStartName;
        out.eventIndexTable = evTable;
        out.alignedInfo     = alignedInfo;  % diagnostics: counts kept/dropped

        fileOutputs(nexPath) = out;
    end
end

%% ---------------- (C) SAVE PER-CONDITION PACKS TO DISK ----------------- %%
condNames = cond2files.keys;
for c = 1:numel(condNames)
    cond = condNames{c};   % char
    if endsWith(cond,'_NeuronSets'), continue; end

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

    events3D = struct();
    if ~isempty(common) && Fcount > 0 && ~isempty(evTables{1})
        X0 = evTables{1}.(common{1}){1};
        [TR_ok, MKT_ok] = size(X0);
        for e = 1:numel(common)
            arr = nan(TR_ok, MKT_ok, Fcount);
            for f = 1:Fcount
                arr(:,:,f) = evTables{f}.(common{e}){1};
            end
            events3D.(common{e}) = arr;
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
variants = finalizeEventVariantReport(eventVariantTally, OUTDIR);
try, diary('off'); catch, end

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

function best = selectBestFile(files)
    idx = find(~cellfun(@isempty, regexp(files,'-sorted-02','once')),1,'first');
    if ~isempty(idx), best = files{idx}; return; end
    idx = find(~cellfun(@isempty, regexp(files,'-sorted-01','once')),1,'first');
    if ~isempty(idx), best = files{idx}; return; end
    idx = find(~cellfun(@isempty, regexp(files,'-sorted','once')),1,'first');
    if ~isempty(idx), best = files{idx}; return; end
    idx = find(~cellfun(@isempty, regexp(files,'-autosorted','once')),1,'first');
    if ~isempty(idx), best = files{idx}; return; end
    best = files{1};
end

function name = normalizeNeuronID(nid)
    if isstring(nid) || ischar(nid), nid = char(nid); else, nid = char(string(nid)); end
    num_part    = regexp(nid,'\d+','match','once');
    letter_part = regexp(nid,'[a-zA-Z]+','match','once');
    if isempty(num_part)
        error('NeuronID missing numeric part: "%s"', nid);
    end
    if isempty(letter_part), letter_part = ''; end
    name = ['SPK' sprintf('%03d', str2double(num_part)) lower(letter_part(1))];
end

function key = parseSessionKeyFromFilename(nexPath)
    [~, base, ~] = fileparts(nexPath);
    tok = regexp(base,'(\d{6})([BL])(\d+)','tokens','once');
    if isempty(tok)
        error('Filename format not recognized: %s', nexPath);
    end
    key = sprintf('%s%s%s', tok{1}, tok{2}, tok{3});
end

function [names, spikes] = normalizeNeurons(F)
    if isfield(F,'neurons') && ~isempty(F.neurons)
        N = numel(F.neurons);
        names  = cell(N,1);
        spikes = cell(N,1);
        for i = 1:N
            names{i}  = char(F.neurons{i}.name);
            spikes{i} = double(F.neurons{i}.timestamps(:))';
        end
    else
        names = {};
        spikes = {};
    end
end

function [names, times] = normalizeEvents(F)
    % Keep ONLY the requested events (others from dataTabFULL instead):
    keepList = { ...
        'm1trialStart','m1rew1on','m1rew1off','m1rew2on','m1rew2off','m1trialStop', ...
        'm2trialStart','m2rew1on','m2rew1off','m2rew2on','m2rew2off','m2trialStop'};

    names = {};
    times = {};
    if ~isfield(F,'events') || isempty(F.events), return; end

    E = numel(F.events);
    raw_names = cell(E,1);
    raw_times = cell(E,1);
    for i = 1:E
        raw_names{i} = char(F.events{i}.name);
        raw_times{i} = double(F.events{i}.timestamps(:))';
    end

    canon2times = containers.Map('KeyType','char','ValueType','any');
    for i = 1:E
        nm_canon = canonicalizeEventName(raw_names{i});
        if ismember(nm_canon, keepList)
            ts = raw_times{i};
            if ~isKey(canon2times, nm_canon)
                canon2times(nm_canon) = ts;
            else
                canon2times(nm_canon) = [canon2times(nm_canon), ts]; %#ok<AGROW>
            end
        end
    end

    keys = canon2times.keys;
    names = keys(:);
    times = cell(numel(keys),1);
    for k = 1:numel(keys)
        tk = canon2times(keys{k});
        tk = sort(tk(:)');
        dt = [true, diff(tk) > 1e-6];
        times{k} = tk(dt);
    end
end

function out = canonicalizeEventName(nm)
    nml = lower(nm);
    if any(strcmp(nml, {'m1start','m1trialon','m1trialstart'})), out = 'm1trialStart'; return; end
    if any(strcmp(nml, {'m2start','m2trialon','m2trialstart'})), out = 'm2trialStart'; return; end
    if ~isempty(regexp(nml, '^m[12]rew[12](on|off)$', 'once')), out = nml; return; end
    if any(strcmp(nml, {'m1trialstart','m2trialstart'})), out = [nml(1:2) 'trialStart']; return; end
    out = nm;
end

function m = ifemptymax(v), if isempty(v), m = 0; else, m = max(v); end, end

function [Sbin, Tbins] = buildSpikeBinMatrix(neuronSpikes, bin_dt, tmax, use_sparse)
    N = numel(neuronSpikes);
    Tbins = max(1, ceil(tmax/bin_dt));
    if use_sparse, Sbin = sparse(N, Tbins); else, Sbin = false(N, Tbins); end
    invdt = 1/bin_dt;
    for n = 1:N
        ts = neuronSpikes{n};
        if isempty(ts), continue; end
        idx = floor(ts*invdt) + 1;
        idx = idx(idx>=1 & idx<=Tbins);
        if isempty(idx), continue; end
        if use_sparse, Sbin(n, unique(idx)) = 1; else, Sbin(n, unique(idx)) = true; end
    end
end

function [trialName, trialStarts] = findTrialStart(eventNames, eventTimes, candidates)
    trialName = '';
    trialStarts = [];
    for k = 1:numel(candidates)
        cand = char(candidates(k));
        m = find(strcmpi(eventNames, cand));
        if ~isempty(m)
            trialName   = eventNames{m(1)};
            trialStarts = eventTimes{m(1)};
            trialStarts = trialStarts(:)';
            return
        end
    end
end

function [eventTimesOut] = truncateAllEventsBeforeFirstM1Start(eventNames, eventTimesIn, nexPath)
% Drop ALL event timestamps earlier than the FIRST entry of m1trialStart.
    eventTimesOut = eventTimesIn;
    idxM1 = find(strcmpi(eventNames,'m1trialStart'),1,'first');
    if isempty(idxM1)
        fprintf('GLOBAL TRUNCATION: No m1trialStart found in %s; skipping.\n', nexPath);
        return;
    end
    tsM1 = eventTimesOut{idxM1};
    if isempty(tsM1)
        fprintf('GLOBAL TRUNCATION: m1trialStart is empty in %s; skipping.\n', nexPath);
        return;
    end
    cut = tsM1(1);
    changed = {};
    for e = 1:numel(eventTimesOut)
        ts = eventTimesOut{e};
        if isempty(ts), continue; end
        before = numel(ts);
        ts = ts(ts >= cut);
        if numel(ts) < before
            eventTimesOut{e} = ts;
            changed{end+1} = eventNames{e}; %#ok<AGROW>
        end
    end
    if ~isempty(changed)
        fprintf('GLOBAL TRUNCATION @ first m1trialStart=%.6f s in %s. Affected: %s\n', ...
            cut, nexPath, strjoin(changed, ', '));
    end
end

function [eventNames, eventTimes, alignedInfo] = alignAllRewardPairs(eventNames, eventTimes, nexPath)
% For each (monkey, rew#), align ON/OFF using 0<dt<12. Try shifting OFF
% index forward/backward around each ON to find the nearest valid OFF.
% After alignment, truncate both ON and OFF so lengths match.
% Also build: m[12]choice (aligned rew1on), m[12]rew1time, m[12]rew2time.

    alignedInfo = struct();
    name2idx = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(eventNames), name2idx(lower(eventNames{i})) = i; end

    monkeys = {'m1','m2'};
    rewards = {'rew1','rew2'};

    % Ensure outputs exist for times
    for mi = 1:2
        for ri = 1:2
            base = [monkeys{mi} rewards{ri}];
            if ~isKey(name2idx, [base 'on']),  eventNames{end+1} = [base 'on'];  eventTimes{end+1} = []; name2idx([base 'on'])  = numel(eventNames); end %#ok<AGROW>
            if ~isKey(name2idx, [base 'off']), eventNames{end+1} = [base 'off']; eventTimes{end+1} = []; name2idx([base 'off']) = numel(eventNames); end %#ok<AGROW>
        end
    end

    for mi = 1:2
        mkey = monkeys{mi};
        for ri = 1:2
            base = [mkey rewards{ri}];
            ion  = name2idx([base 'on']);
            ioff = name2idx([base 'off']);
            on   = eventTimes{ion};
            off  = eventTimes{ioff};

            [on_al, off_al, kept_idx_on, kept_idx_off] = pair_on_off_bidirectional(on, off);

            % Truncate to same length (safety; pairer already returns equal lengths)
            L = min(numel(on_al), numel(off_al));
            on_al  = on_al(1:L);  off_al = off_al(1:L);

            % Write back aligned vectors
            eventTimes{ion}  = on_al;
            eventTimes{ioff} = off_al;

            % Save debug info
            alignedInfo.([base '_n_on_raw'])  = numel(on);
            alignedInfo.([base '_n_off_raw']) = numel(off);
            alignedInfo.([base '_n_on_kept']) = numel(on_al);
            alignedInfo.([base '_n_off_kept'])= numel(off_al);

            % Build *time durations from aligned pairs
            name_time = [base 'time'];
            dur = off_al - on_al;  % already guaranteed 0<dt<12 by pairer
            % attach or create variable
            idx_time = find(strcmpi(eventNames, name_time), 1, 'first');
            if isempty(idx_time)
                eventNames{end+1} = name_time; eventTimes{end+1} = dur; %#ok<AGROW>
            else
                eventTimes{idx_time} = dur;
            end

            % Build m[12]choice from aligned rew1on (only for rew1)
            if strcmp(rewards{ri}, 'rew1')
                cho_name = [mkey 'choice'];
                cho_ts   = on_al - 0.250;   % 250 ms shift
                idx_cho = find(strcmpi(eventNames, cho_name), 1, 'first');
                if isempty(idx_cho)
                    eventNames{end+1} = cho_name; eventTimes{end+1} = cho_ts; %#ok<AGROW>
                else
                    eventTimes{idx_cho} = cho_ts;
                end
            end
        end
    end

    fprintf('ALIGNMENT SUMMARY [%s]:\n', nexPath);
    disp(alignedInfo);
end

function [on_out, off_out, kept_on_idx, kept_off_idx] = pair_on_off_bidirectional(on_in, off_in)
% Greedy, per-sequence alignment:
%   - Inputs are already globally truncated.
%   - For each ON, search OFF indices around the current j pointer to find the
%     nearest OFF that satisfies 0 < (OFF-ON) < 12. We allow moving OFF pointer
%     forward OR backward relative to its current position to minimize |Δt|, but
%     Δt must remain positive and < 12.
%   - If multiple OFFs satisfy the rule, pick the one with the smallest Δt.
%   - OFF and ON used once each (1-to-1 pairing). Unpaired stamps are dropped.
% Returns paired vectors of the same length and the indices kept from inputs.

    on  = on_in(:)';  off = off_in(:)';
    on_out = []; off_out = [];
    kept_on_idx  = []; kept_off_idx = [];

    if isempty(on) || isempty(off), return; end

    i = 1; j_used = false(1, numel(off));
    while i <= numel(on)
        on_t = on(i);

        % Candidate OFFs that can pair with this ON (OFF > ON and OFF-ON < 12)
        mask_valid = (off > on_t) & ((off - on_t) < 12) & ~j_used;
        cand_idx   = find(mask_valid);
        if isempty(cand_idx)
            % No valid OFF for this ON → drop ON
            i = i + 1;
            continue;
        end
        % Choose the OFF with minimal (OFF-ON)
        [~, k] = min(off(cand_idx) - on_t);
        j = cand_idx(k);

        % Commit the pair
        on_out(end+1)   = on_t;    %#ok<AGROW>
        off_out(end+1)  = off(j);  %#ok<AGROW>
        kept_on_idx(end+1)  = i;   %#ok<AGROW>
        kept_off_idx(end+1) = j;   %#ok<AGROW>
        j_used(j) = true;
        i = i + 1;
    end

    % Ensure same length (should be by construction)
    L = min(numel(on_out), numel(off_out));
    on_out  = on_out(1:L);
    off_out = off_out(1:L);
end

function evTable = buildEventIndexTable_fromAligned(eventNames, eventTimes, trialStarts, bin_dt, TR, MKT, Tbins, nexPath)
% Convert:
%   - Timestamp events → per-window FIRST occurrence index (bin-based)
%   - Duration events (*time, numeric scalars per valid pair) → per-window FIRST duration
% Windows are defined by this monkey's trialStart: [t0(i), t1(i)), i=1..90

    totalNeeded = TR*MKT;  % usually 90
    t0 = trialStarts(length(trialStarts)-90+1:end);
    % if numel(t0) < totalNeeded
    %     t0(end+1:totalNeeded) = NaN; %#ok<AGROW>
    % else
    %     t0 = t0(1:totalNeeded);
    % end
    t1 = [t0(2:end), Inf];

    % Map name → idx
    name2idx = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(eventNames), name2idx(lower(eventNames{i})) = i; end

    varNames = matlab.lang.makeValidName(eventNames);
    evCells  = cell(1, numel(eventNames));

    for e = 1:numel(eventNames)
        ename = eventNames{e};
        ename_l = lower(ename);
        isDuration = ~isempty(regexp(ename_l, '^m[12]rew[12]time$', 'once'));

        if isDuration
            % Raw durations; pick first per window (these are per-pair values, not stamps)
            vdur = eventTimes{e}(:)';
            % We need to assign per trial window -> use aligned ONs as anchors:
            % Find the corresponding aligned ON vector to know which trials they fall in.
            on_name = regexprep(ename_l, 'time$', 'on');
            if ~isKey(name2idx, on_name)
                V = nan(length(t0),1);
            else
                on_ts = eventTimes{name2idx(on_name)}(:)';
                V = windowFirstScalar(on_ts, vdur, t0, t1);
            end
            evCells{e} = reshape(V, [TR, MKT]);
            continue;
        else
            % Timestamp → first time in window → bin index
            ts = eventTimes{e};
            if isempty(ts)
                v = nan(length(t0),1);
                counts = zeros(length(t0),1);
            else
                [v, counts] = pickFirstInWindowsCount(ts, t0, t1);
            end

            bad = find(counts > 1);
            for bi = 1:numel(bad)
                i = bad(bi);
                a = t0(i); b = t1(i); if ~isfinite(b), b = inf; end
                fprintf('WARNING [%s] event "%s": trial %d has %d occurrence(s) in [%.6f, %.6f)\n', ...
                    nexPath, ename, i, counts(i), a, b);
            end

            idx = round(v./bin_dt) + 1;
            idx(~isfinite(idx) | idx<1 | idx>Tbins) = NaN;
            evCells{e} = reshape(idx, [TR, MKT]);
        end
    end

    evTable = table();
    for e = 1:numel(eventNames), evTable.(varNames{e}) = {evCells{e}}; end
end

function [v, counts] = pickFirstInWindowsCount(ts, t0, t1)
    n = numel(t0);
    v = nan(n,1); counts = zeros(n,1);
    if isempty(ts), return; end
    ts = ts(:)';
    for i = 1:n
        a = t0(i); b = t1(i);
        if ~isfinite(a), counts(i) = 0; continue; end
        in = (ts >= a) & (ts < b);
        counts(i) = sum(in);
        if counts(i) >= 1, v(i) = ts(find(in, 1, 'first')); end
    end
end

function V = windowFirstScalar(on_ts, dur_vec, t0, t1)
% Assign the FIRST duration to the window that contains the matching ON.
% If multiple ONs fall in the same window, take the first one.
    n = numel(t0);
    V = nan(n,1);
    if isempty(on_ts) || isempty(dur_vec), return; end
    on_ts = on_ts(:)';
    K = min(numel(on_ts), numel(dur_vec));
    on_ts = on_ts(1:K); dur_vec = dur_vec(1:K);
    for i = 1:n
        a = t0(i); b = t1(i);
        if ~isfinite(a), continue; end
        in = (on_ts >= a) & (on_ts < b);
        idx = find(in, 1, 'first');
        if ~isempty(idx)
            V(i) = dur_vec(idx);
        end
    end
end

function evTableOut = addBehaviorVarsToEventTable(evTableIn, nexPath, Tfull, TR, MKT)
% (unchanged from your earlier version; attaches behavior matrices)
    evTableOut = evTableIn;
    if isempty(Tfull), return, end

    [~, base, ~] = fileparts(nexPath);
    tok = regexp(base,'(\d{6})([BL])(\d+)','tokens','once');
    if isempty(tok), warning('Cannot parse date/session from filename: %s', base); return, end
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

    % varsToExtract = {'session','option','optimal','marketOrig','tResp','bubbleMarket','monkey'};
    varsToExtract = Tfull.Properties.VariableNames;
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

function updateEventVariantTally(tally, eventNames)
    if isempty(eventNames), return; end
    keep = cellfun(@(s) startsWith(lower(s),'m1') || startsWith(lower(s),'m2'), eventNames);
    names = eventNames(keep);
    for i = 1:numel(names)
        raw = names{i};
        base = lower(raw);
        if ~isKey(tally, base)
            sub = containers.Map('KeyType','char','ValueType','double');
            sub(raw) = 1; tally(base) = sub;
        else
            sub = tally(base);
            if isKey(sub, raw), sub(raw) = sub(raw) + 1; else, sub(raw) = 1; end
            tally(base) = sub;
        end
    end
end

function variants = finalizeEventVariantReport(tally, OUTDIR)
    bases = tally.keys;
    variants = struct('base',{},'nVariants',{},'variants',{},'counts',{});
    fprintf('\n=== Event Name Variation Report (m1*/m2*) ===\n');
    fprintf('Base form = lowercase exact name (case-insensitive grouping)\n\n');
    bases = sort(bases);
    for i = 1:numel(bases)
        base = bases{i};
        sub  = tally(base);
        vraw = sub.keys;
        vcnt = cellfun(@(k) sub(k), vraw);
        [vcnt, order] = sort(vcnt, 'descend');
        vraw = vraw(order);
        variants(end+1).base     = base; %#ok<AGROW>
        variants(end).nVariants  = numel(vraw);
        variants(end).variants   = vraw;
        variants(end).counts     = vcnt;
        fprintf('%s  | variants: %d\n', base, numel(vraw));
        for j = 1:numel(vraw)
            fprintf('   - %-20s  count: %d\n', vraw{j}, vcnt(j));
        end
    end
    try
        outRows = {};
        for i = 1:numel(variants)
            base  = variants(i).base;
            for j = 1:numel(variants(i).variants)
                outRows(end+1, :) = {base, variants(i).nVariants, variants(i).variants{j}, variants(i).counts(j)}; %#ok<AGROW>
            end
        end
        T = cell2table(outRows, 'VariableNames', {'base','nVariants','variant','count'});
        outf = fullfile(OUTDIR, 'event_name_variations.tsv');
        writetable(T, outf, 'FileType','text','Delimiter','\t');
        fprintf('\nWrote variation table: %s\n', outf);
    catch ME
        fprintf('Could not write TSV: %s\n', ME.message);
    end
end

function checkEventTableConsistency(evTable, TR, MKT, nexPath)
    if isempty(evTable)
        fprintf('CONSISTENCY CHECK: empty event table for %s\n', nexPath);
        return
    end
    vnames = evTable.Properties.VariableNames;
    shapes_ok = true;
    nnz_counts = nan(numel(vnames),1);
    for i = 1:numel(vnames)
        X = evTable.(vnames{i}){1};
        sz_ok = isequal(size(X), [TR, MKT]);
        if ~sz_ok
            shapes_ok = false;
            fprintf('CONSISTENCY CHECK: %s shape %s != [%d %d] in %s\n', ...
                vnames{i}, mat2str(size(X)), TR, MKT, nexPath);
        end
        nnz_counts(i) = nnz(~isnan(X));
    end
    if shapes_ok
        fprintf('CONSISTENCY CHECK: all variables have shape [%d %d] in %s\n', TR, MKT, nexPath);
    end
    u = unique(nnz_counts(~isnan(nnz_counts)));
    if numel(u) == 1
        fprintf('CONSISTENCY CHECK: all variables have the same finite-count = %d in %s\n', u, nexPath);
    else
        fprintf('CONSISTENCY CHECK: finite-counts differ in %s:\n', nexPath);
        for i = 1:numel(vnames)
            fprintf('  %s: %d\n', vnames{i}, nnz_counts(i));
        end
    end
end
 