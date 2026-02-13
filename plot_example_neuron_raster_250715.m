%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike Raster Plot with Event Alignment
%
% Description:
%   Loads a .nex neural recording file, aligns trials to a specific event, 
%   bins spike data, annotates other task events, and plots trial-by-trial 
%   rasters for a specified neuron with visual overlays for key intervals.
%
% Inputs:
%   - root_path: Path to the folder containing the .nex file
%   - file_name: Name of the .nex file to load
%   - aligned_events: Event to align all trials to (e.g., 'm1trialStart')
%   - other_events: Events to annotate in the raster (e.g., reward on/off)
%   - beforeEvent: Time (s) before the aligned event to include per trial
%   - afterEvent: Time (s) after the aligned event to include per trial
%   - temp_res: Temporal resolution (bin width in seconds)
%
% Output:
%   - Trial-aligned spike raster plot for the specified neuron
%   - Event marker overlays and interval highlights on each trial
%
% Dependencies:
%   - Requires NexTools or similar implementation for readNexFile()
%
% Author: Tim Chen
% Lab: Platt Lab
% Date: 7/16/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear; close all; clc;

%% Parameters
root_path = "C:\Users\plattlab\Tim\Stock_market_experiment\MSM_Sorted\";
file_name = "082118B2-sorted.nex";

aligned_events = ["m1trialStart"]; % this is what event to align each trial to
other_events = ["m1rew1On", "m1rew1Off", "m1rew2On", "m1rew2Off", "m1trialStop", "m2trialStart"]; % this is what other events to annotate during each trial
beforeEvent = 1;  % seconds, the window before the aligned event
afterEvent = 10;  % seconds, the window after the aligned event
temp_res = 0.001; % seconds, bin resolution

%% Load Data
file = readNexFile(root_path + file_name);

% Extract events and neurons
events = cell(length(file.events), 2);
for i = 1:length(file.events)
    events{i, 1} = string(file.events{i}.name);
    events{i, 2} = file.events{i}.timestamps;
end

neurons = cell(length(file.neurons), 2);
for i = 1:length(file.neurons)
    neurons{i, 1} = string(file.neurons{i}.name);
    neurons{i, 2} = file.neurons{i}.timestamps;
end

%% Define Trial Intervals
aligned_events_idx = find(ismember(cell2mat({events{:,1}}'), aligned_events));
other_events_idx = find(ismember(cell2mat({events{:,1}}'), other_events));

intervals = cell(length(aligned_events_idx), 1);
for i = 1:length(aligned_events_idx)
    ts = events{aligned_events_idx(i), 2};
    intervals{i} = [ts - beforeEvent, ts + afterEvent];
end

% Get aligned reference times (per trial)
aligned_times = events{aligned_events_idx, 2};

% For each "other" event type
aligned_other_event_times = containers.Map();

for i = 1:length(other_events)
    event_name = other_events(i);
    event_ts = events{other_events_idx(i), 2};

    % Map each timestamp to the nearest aligned_time (1-to-1 assumption)
    aligned_rel_times = zeros(size(aligned_times));

    for t = 1:length(aligned_times)
        % Find the closest event timestamp after the aligned time
        idx = find(event_ts >= intervals{1}(t,1) & event_ts <= intervals{1}(t,2), 1, 'first');
        if ~isempty(idx)
            aligned_rel_times(t) = event_ts(idx) - aligned_times(t);
        else
            aligned_rel_times(t) = NaN;
        end
    end

    % Store relative times
    aligned_other_event_times(event_name) = aligned_rel_times;
end

%% Bin Spikes
neuronData = cell(length(aligned_events_idx), 1);

for evt = 1:length(aligned_events_idx)
    num_trials = size(intervals{evt}, 1);
    neuronData{evt} = cell(num_trials, 1);

    for trial = 1:num_trials
        neuronData{evt}{trial} = cell(length(neurons), 2);
        t0 = intervals{evt}(trial, 1);
        t1 = intervals{evt}(trial, 2);
        edges = t0:temp_res:t1;

        for n = 1:length(neurons)
            spikes = neurons{n, 2};
            trial_spikes = spikes(spikes >= t0 & spikes < t1);
            neuronData{evt}{trial}{n, 1} = trial_spikes;
            neuronData{evt}{trial}{n, 2} = histcounts(trial_spikes, edges);
        end
    end
end


%% Plot All Trials for One Neuron
neuron_idx = 21;    % SPK011c for 082118B2-sorted.nex
evt = 1;
sigma = 1;

num_trials = length(neuronData{evt});
edges = intervals{evt}(1,1):temp_res:intervals{evt}(1,2);
centers = edges(1:end-1) + temp_res/2;

% Collect raster + binned
raster_spikes = cell(num_trials, 1);
binned_matrix = zeros(num_trials, length(centers));

for trial = 1:num_trials
    spikes = neuronData{evt}{trial}{neuron_idx, 1};
    event_time = events{aligned_events_idx(evt), 2}(trial);
    raster_spikes{trial} = spikes - event_time;
    binned_matrix(trial, :) = neuronData{evt}{trial}{neuron_idx, 2};
end

% Get sorting durations
t_start_list = zeros(num_trials, 1); % all 0 since aligned
t_stop_sort = aligned_other_event_times(char("m2trialStart"));
[~, sort_idx] = sort(t_stop_sort);  % sort by duration

% For highlight use m1rew2Off
t_stop_highlight = aligned_other_event_times(char("m1trialStop"));

figure; hold on;

keys = {'m1rew1On', 'm1rew1Off', 'm1rew2On', 'm1rew2Off', 'm2trialStart'};
markers = {'v', '^', 'v', '^', 'o'};
colors = {[0 0.6 0], [0 0.6 0], [0.8 0 0], [0.8 0 0], [0 0.2 1]};
event_marker_map = containers.Map(keys, markers);
event_color_map = containers.Map(keys, colors);

for i = 1:num_trials
    trial = sort_idx(i);
    spikes = raster_spikes{trial};

    if ~isempty(spikes)
        x = repmat(spikes(:)', 3, 1); x(3,:) = NaN;
        y = [repmat(i - 0.3, 1, length(spikes));
             repmat(i + 0.3, 1, length(spikes));
             NaN(1, length(spikes))];
        plot(x(:), y(:), 'k', 'HandleVisibility', 'off');
    end

    % Overlay events
    for j = 1:length(other_events)
        evt_name = other_events(j);
        if evt_name == "m1trialStop"
            continue;
        end
        rel_time = aligned_other_event_times(char(evt_name));
        t_rel = rel_time(trial);
        if ~isnan(t_rel)
            plot(t_rel, i, ...
                'Marker', event_marker_map(char(evt_name)), ...
                'MarkerEdgeColor', event_color_map(char(evt_name)), ...
                'MarkerFaceColor', event_color_map(char(evt_name)), ...
                'MarkerSize', 4, 'HandleVisibility', 'off');
        end
    end

    % Highlight m1trialStart to m1rew2Off
    t_high = t_stop_highlight(trial);
    if ~isnan(t_high)
        x_patch = [0, t_high, t_high, 0];
        y_patch = [i - 0.5, i - 0.5, i + 0.5, i + 0.5];
        patch(x_patch, y_patch, [1 1 0], ...
              'EdgeColor', 'none', ...
              'FaceAlpha', 0.2, ...
              'HandleVisibility', 'off');
    end
end

% xline(0, '-r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlim([-beforeEvent, afterEvent]);
ylim([0, num_trials + 1]);
xlabel('Time (s)'); ylabel('Trial');

% Manual legend
plot(NaN, NaN, '^', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], ...
     'DisplayName', 'm1rew1On', 'LineStyle', 'none');
plot(NaN, NaN, 'v', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], ...
     'DisplayName', 'm1rew1Off', 'LineStyle', 'none');
plot(NaN, NaN, '^', 'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], ...
     'DisplayName', 'm1rew2On', 'LineStyle', 'none');
plot(NaN, NaN, 'v', 'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], ...
     'DisplayName', 'm1rew2Off', 'LineStyle', 'none');
plot(NaN, NaN, 'o', 'Color', [0 0.2 1], 'MarkerFaceColor', [0 0.2 1], ...
     'DisplayName', 'm2trialStart', 'LineStyle', 'none');
legend('Location', 'bestoutside');

title(neurons{neuron_idx, 1})
