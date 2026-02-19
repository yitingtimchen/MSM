% File: d20250804_review_neuron_trial_pileup.m
% This code was made to directly view each neuron in the .nex files. Tim,
% 251104

clear; close all; clc;

%% Load Data
root_path = "C:\Users\plattlab\MSM\data\neural\raw_external\MSM_Sorted\nex\";
sheet = "C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx";

animal = "L";

interactiveNeuronBrowser(root_path, sheet, animal)


function interactiveNeuronBrowser(rootPath, sheetPath, animal)
% INTERACTIVENEURONBROWSER Browse neurons with dropdown menus and smoothed PSTH
%   rootPath: directory containing .nex files
%   sheetPath: full path to dataTabFULL.xlsx
%   animal: character for file filtering, e.g. 'B'

% Parameters
alignedEventLabels = ["m1rew1On","m1rew1on","m1rew1ON"];  % extendable list of acceptable labels
beforeEvent      = 1;
afterEvent       = 2;
temp_res         = 0.01;
sigma            = 0;

% Load condition table (Sheet1)
opts = detectImportOptions(sheetPath, 'Sheet','Sheet1');
opts.VariableNamingRule    = 'preserve';
opts.SelectedVariableNames = opts.VariableNames;
condTable = readtable(sheetPath, opts);

% Persistent cache file on disk
cacheFile = fullfile(rootPath, sprintf('%s_validFiles.mat', animal));
if exist(cacheFile, 'file')
    s = load(cacheFile, 'fileList');
    fileList = s.fileList;
else
    % Validate files once and save
    rawFiles = dir(fullfile(rootPath, ['*', char(animal), '*.nex']));
    validList = {};
    for i = 1:numel(rawFiles)
        fname = rawFiles(i).name;
        data  = d20260213_read_nex_file(fullfile(rootPath,fname));
        evNames = string(cellfun(@(e)e.name, data.events, 'UniformOutput',false));
        evTimes = cellfun(@(e)e.timestamps, data.events, 'UniformOutput',false);
        % check membership
        idxA = find(ismember(evNames, alignedEventLabels), 1);
        if ~isempty(idxA) && ~isempty(evTimes{idxA})
            validList{end+1} = fname; %#ok<AGROW>
        else
            fprintf('Excluding %s: no valid alignment event.\n', fname);
        end
    end
    if isempty(validList)
        error('No files with valid alignment labels found.');
    end
    fileList = validList;
    save(cacheFile, 'fileList');
end

numFiles = numel(fileList);
fileIdx   = 1;
neuronIdx = 1;

% Create figure and axes
hFig = figure('Name','Neuron Browser','NumberTitle','off','Position',[200 200 900 550]);
hAx  = axes('Parent',hFig,'Units','normalized','Position',[0.1 0.25 0.85 0.7]);

% File dropdown
hFileDropdown = uicontrol('Style','popupmenu','String',fileList, ...
    'Units','normalized','Position',[0.1 0.15 0.3 0.05],'Callback',@fileDropdownCallback);
% Neuron dropdown
hNeuronDropdown = uicontrol('Style','popupmenu','String',{''}, ...
    'Units','normalized','Position',[0.5 0.15 0.3 0.05],'Callback',@neuronDropdownCallback);

% Navigation buttons
uicontrol('Style','pushbutton','String','Prev','Units','normalized', ...
    'Position',[0.05 0.05 0.1 0.05],'Callback',@prevCallback);
uicontrol('Style','pushbutton','String','Next','Units','normalized', ...
    'Position',[0.85 0.05 0.1 0.05],'Callback',@nextCallback);

% Initialize dropdowns and plot
enumUpdate();
updatePlot();

    function fileDropdownCallback(src,~)
        fileIdx = src.Value;
        neuronIdx = 1;
        enumUpdate();
        updatePlot();
    end

    function neuronDropdownCallback(src,~)
        neuronIdx = src.Value;
        updatePlot();
    end

    function prevCallback(~,~)
        neuronIdx = neuronIdx - 1;
        if neuronIdx < 1
            fileIdx = fileIdx - 1;
            if fileIdx < 1, fileIdx = numFiles; end
            enumUpdate();
            neuronIdx = numel(get(hNeuronDropdown,'String'));
        end
        set(hFileDropdown,'Value',fileIdx);
        set(hNeuronDropdown,'Value',neuronIdx);
        updatePlot();
    end

    function nextCallback(~,~)
        neuronIdx = neuronIdx + 1;
        if neuronIdx > numel(get(hNeuronDropdown,'String'))
            fileIdx = mod(fileIdx,numFiles) + 1;
            enumUpdate();
            neuronIdx = 1;
        end
        set(hFileDropdown,'Value',fileIdx);
        set(hNeuronDropdown,'Value',neuronIdx);
        updatePlot();
    end

    function enumUpdate()
        data = d20260213_read_nex_file(fullfile(rootPath,fileList{fileIdx}));
        names = cellfun(@(n)n.name, data.neurons, 'UniformOutput',false);
        set(hNeuronDropdown,'String',names,'Value',neuronIdx);
    end

    function updatePlot()
        cla(hAx); hold(hAx,'off');
        fname = fileList{fileIdx};
        data  = d20260213_read_nex_file(fullfile(rootPath,fname));
        % Condition lookup
        try
            condStr = getConditionFromFileName(fname, condTable, animal);
        catch
            condStr = 'Unknown';
        end
        % Alignment times using explicit acceptable labels list
        evNames = string(cellfun(@(e)e.name,data.events,'UniformOutput',false));
        evTimes = cellfun(@(e)e.timestamps,data.events,'UniformOutput',false);
        idxA = find(ismember(evNames, alignedEventLabels), 1);
        if isempty(idxA)
            error('Alignment event not found. Acceptable labels: %s', strjoin(alignedEventLabels, ', '));
        end
        tAlign = evTimes{idxA};
        % Bin edges
        edges   = -beforeEvent:temp_res:afterEvent;
        centers = edges(1:end-1)+temp_res/2;
        % Spikes and binning
        neuron = data.neurons{neuronIdx};
        spikes = neuron.timestamps;
        binned = zeros(numel(tAlign),numel(centers));
        for t = 1:numel(tAlign)
            dts = spikes - tAlign(t);
            dts = dts(dts>=-beforeEvent & dts<afterEvent);
            binned(t,:) = histcounts(dts,edges);
        end
        % Smooth PSTH
        sm = gaussianSmooth1DPerNeuron(binned,sigma);
        mu = mean(sm,1) ./ temp_res; % Hz
        sdv = std(sm,[],1);
        % Plot
        hold on
        fill(hAx,[centers fliplr(centers)],[mu+sdv fliplr(mu-sdv)], ...
             [0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5);
        plot(hAx,centers,mu,'k','LineWidth',1.5);
        hold off
        xline(hAx,0,'r','LineWidth',1.5);
        xlim(hAx,[-beforeEvent afterEvent]);
        ylim([0 200])
        xlabel(hAx,'Time (s)'); 
        ylabel(hAx,'Firing rate (Hz)');
        title(hAx,sprintf('%s | %s | %s (%d/%d)',fname,unique(condStr),neuron.name,neuronIdx,numel(data.neurons)),'Interpreter','none');
    end
end

function condition = getConditionFromFileName(fileName,T,animal)
    tokens = regexp(fileName,['(\d{2})(\d{2})(\d{2})', char(animal), '(\d+)'],'tokens');
    yy = str2double(tokens{1}{3}); mm = str2double(tokens{1}{1}); dd = str2double(tokens{1}{2}); sess = str2double(tokens{1}{4}); year = 2000+yy;
    idx = T.year==year & T.month==mm & T.day==dd & T.SessionOfDay==sess;
    condition = string(T.condition(idx));
end

function smoothed = gaussianSmooth1DPerNeuron(data,sigma)
    if sigma==0, smoothed=data; return; end
    half = ceil(3*sigma); x = -half:half;
    ker = exp(-x.^2/(2*sigma^2)); ker = ker/sum(ker);
    pad = padarray(data,[0 half],'replicate','both');
    c = conv2(pad,ker(:)','same');
    smoothed = c(:,half+1:end-half);
end
