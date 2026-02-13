clc
close all

%% ==================== CONFIG ====================
r_thresh = 0.9;         % correlation threshold
variant  = 'raw';       % 'raw' | 'base' | 'z' (flexible)

% Condition sets
groupsets = {
    {'AI','Replay','Decoy','Live'}, ...
    {'OT_AI','OT_Replay','OT_Decoy','OT_Live'}, ...
    {'Saline_AI','Saline_Replay','Saline_Decoy','Saline_Live'}
};
setnames = {'AIReplayDecoyLive','OT_AIOT_ReplayOT_DecoyOT_Live','Saline_AISaline_ReplaySaline_DecoySaline_Live'};

root_dir = 'C:\Users\plattlab\Tim\Stock_market_experiment\Tim\condition_packs\';

%% ==================== LOAD UNITS ====================
units_list = load_units();   % Loads and trims units per group

%% ==================== MAIN PROCESS ====================
keep_indices = struct();
keep_assign  = struct();
corr_data    = struct();

for ss = 1:numel(setnames)
    setname = setnames{ss};
    groups  = groupsets{ss};
    
    for g = 1:numel(groups)
        grp = groups{g};
        
        % --- Load & orient Units ---
        U = orient_neurons_by_bins(units_list.(grp), 400, grp, 'Units');
        
        % --- Load & orient AllSpikes ---
        Sraw = AllSpikes.(setname).m1rew1on.(grp).(variant);
        S = orient_neurons_by_bins(Sraw, 400, grp, 'AllSpikes');
        
        % --- Standardize for robust matching ---
        U = zscore_rows(U);
        S = zscore_rows(S);
        
        % --- Correlation matrix ---
        C = corr(U.', S.', 'Rows','pairwise');   % [Units x AllSpikes]
        C(~isfinite(C)) = -Inf;
        corr_data.(setname).(grp).(variant) = C; % Store correlations
        
        % --- Best match per AllSpikes neuron ---
        [bestR, bestUnitIdx] = max(C, [], 1);
        bestR       = bestR(:);
        bestUnitIdx = bestUnitIdx(:);
        
        % --- Threshold & deduplicate ---
        N = size(S,1);
        pass = find(bestR >= r_thresh);
        keep_mask  = false(N,1);
        used_units = false(size(U,1),1);
        
        if ~isempty(pass)
            [~, order] = sort(bestR(pass), 'descend');
            for j = pass(order).'
                i = bestUnitIdx(j);
                if ~used_units(i)
                    keep_mask(j) = true;
                    used_units(i) = true;
                end
            end
        end
        
        % --- Collect outputs ---
        keep_j = find(keep_mask);
        keep_indices.(setname).(grp).(variant) = keep_j;
        
        % Audit table
        keep_assign.(setname).(grp).(variant) = table( ...
            keep_j, bestUnitIdx(keep_j), bestR(keep_j), ...
            'VariableNames', {'AllSpikesNeuron','MatchedUnitsNeuron','Correlation'} ...
        );
        
        fprintf('%s: kept %d / %d (r >= %.2f)\n', grp, numel(keep_j), N, r_thresh);
    end
end

%% ==================== VISUALIZATION ====================
for ss = 1:numel(setnames)
    setname = setnames{ss};
    groups  = groupsets{ss};
    
    figure('Name', setname, 'Position', [100 100 1700 600]);
    
    for g = 1:numel(groups)
        grp = groups{g};
        C = corr_data.(setname).(grp).(variant);
        keep_j = keep_indices.(setname).(grp).(variant);

        % --- Mask correlations below threshold ---
        C_masked = C;
        C_masked(C_masked < r_thresh) = NaN;   % or = 0 if you prefer

        % --- Heatmap subplot ---
        subplot(2,numel(groups),g);
        [~, bestUnitIdx] = max(C, [], 1);        % 1 x N
        [~, sortOrder] = sort(bestUnitIdx);      % sort AllSpikes neurons by matched Units
        C_sorted = C_masked(:, sortOrder);
        imagesc(C_sorted);
        caxis([r_thresh 1]);                      % scale from threshold to 1
        colorbar;
        
        % Axis labels with neuron counts
        numUnits = size(C_sorted, 1);
        numAllSpikes = size(C_sorted, 2);
        xlabel(sprintf('AllSpikes neurons (sorted, N=%d)', numAllSpikes), 'Interpreter','none');
        ylabel(sprintf('Units neurons (N=%d)', numUnits), 'Interpreter','none');
        
        title(sprintf('%s - Heatmap', grp), 'Interpreter','none');

        % --- Histogram subplot ---
        subplot(2,numel(groups),g + numel(groups));
        bestR = max(C,[],1);
        if ~isempty(keep_j)
            bestR_keep = bestR(keep_j);
            histogram(bestR_keep, r_thresh:0.005:1, 'FaceColor', [0.2 0.6 0.5]);
        else
            histogram([], 20);
        end
        xline(r_thresh,'r--','Threshold','LineWidth',1.5);
        xlabel('Correlation', 'Interpreter','none');
        ylabel('Count', 'Interpreter','none');
        xlim([r_thresh 1]);
        title(sprintf('%s - Kept neurons', grp), 'Interpreter','none');
    end
    
    sgtitle(sprintf('%s - %s', setname, variant), 'Interpreter','none');

    % --- Save after plotting each figure ---
    outdir = fullfile(pwd, 'neuron_validation');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    saveas(gcf, fullfile(outdir, sprintf('%s_%s.png', setname, variant)));
    savefig(gcf, fullfile(outdir, sprintf('%s_%s.fig', setname, variant)));
end


save([outdir, '/keep_indices.mat'], 'keep_indices', 'variant', '-v7.3');  % persist for plotting

%% ==================== LOCAL FUNCTIONS ====================

function units = load_units()
    base = 'C:\Users\plattlab\Tim\Stock_market_experiment\MATLAB_Copy\MSM_SpikeData';
    files = struct( ...
        'AI',     fullfile(base,'AI','m1rew1on.mat'), ...
        'Replay', fullfile(base,'Replay','m1rew1on.mat'), ...
        'Decoy',  fullfile(base,'Decoy','m1rew1on.mat'), ...
        'Live',   fullfile(base,'Live','m1rew1on.mat'), ...
        'OT_AI',     fullfile(base,'OT','AI','m1rew1on.mat'), ...
        'OT_Replay', fullfile(base,'OT','Replay','m1rew1on.mat'), ...
        'OT_Decoy',  fullfile(base,'OT','Decoy','m1rew1on.mat'), ...
        'OT_Live',   fullfile(base,'OT','Live','m1rew1on.mat'), ...
        'Saline_AI',     fullfile(base,'Saline','AI','m1rew1on.mat'), ...
        'Saline_Replay', fullfile(base,'Saline','Replay','m1rew1on.mat'), ...
        'Saline_Decoy',  fullfile(base,'Saline','Decoy','m1rew1on.mat'), ...
        'Saline_Live',   fullfile(base,'Saline','Live','m1rew1on.mat') ...
        );
    units = struct();
    for f = fieldnames(files).'
        grp = f{1};
        S = load(files.(grp),'spikes');
        units.(grp) = S.spikes;
    end
    % % Optional: do not delete! Apply trimming rules
    % units.AI     = units.AI(:,1:160);
    % units.Decoy  = [units.Decoy(:,1:114), units.Decoy(:,129:139), units.Decoy(:,end-34:end)];
    % units.Live   = [units.Live(:,1:62), units.Live(:,end-97:end)];
end

function Xo = orient_neurons_by_bins(X, bins, grp, tag)
    if size(X,2) == bins
        Xo = X;
    elseif size(X,1) == bins
        Xo = X.';
    else
        error('"%s" %s has no %d-bin dimension.', grp, tag, bins);
    end
end

function Xz = zscore_rows(X)
    mu = mean(X, 2, 'omitnan');
    Xc = X - mu;
    sd = std(Xc, 0, 2, 'omitnan');
    sd(~isfinite(sd) | sd < 1e-12) = 1;
    Xz = Xc ./ sd;
end
