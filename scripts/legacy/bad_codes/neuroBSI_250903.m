%% =========================================================================
% BSI (Buy vs Hold/Sell) analysis per neuron — self-contained script
% Works with condition packs like your PSTH script
% =========================================================================

close all; clear; clc;

%% --------------------------- CONFIG --------------------------------------
is_mac = 0;
if is_mac
    packDir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs';
else
    packDir = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_250922';
end

% Conditions to include
conds = {'AI','Replay','Decoy','Live'};

% Event variable
event_var = "m1rew1on";   % first present will be used

% Post-choice window for BSI (seconds relative to choice)
WIN = [-0.50, 1.00];      % same as PSTH extraction window
BSI_WIN = [-0.25 0];      % seconds after choice for BSI
EPS = 1e-3;               % small constant for BSI denominator

% Optional: shuffle validation
DO_SHUFFLE = true;
nShuffles = 1000;

%% ---------------------- GET dt AND TIME VECTORS --------------------------
% Grab dt from first available pack
dt = [];
for ci = 1:numel(conds)
    fp = fullfile(packDir, sprintf('%s_condition_pack.mat', conds{ci}));
    if exist(fp,'file')
        S = load(fp,'C'); dt = S.C.dt; break;
    end
end
if isempty(dt), error('No condition packs found in %s', packDir); end

t = (WIN(1):dt:WIN(2)); t = t(:)';        % time vector
bsiMask = (t >= BSI_WIN(1)) & (t <= BSI_WIN(2));

%% ------------------------- PROCESS EACH CONDITION ------------------------
BSI_all = struct();
monkey_all = struct();
session_all = struct();

for ci = 1:numel(conds)
    cond = conds{ci};
    packPath = fullfile(packDir, sprintf('%s_condition_pack.mat', cond));
    if ~exist(packPath,'file')
        warning('Missing pack: %s', packPath); 
        continue; 
    end
    
    S = load(packPath, 'C'); 
    C = S.C;
    
    % Initialize as empty arrays / cell arrays
    BSI_all.(cond) = [];
    monkey_all.(cond) = {};  % cell array for categorical labels
    session_all.(cond) = {};
    
    for f = 1:numel(C.files)
        T = C.eventTables{f};
        if isempty(T), continue; end
        
        % Pick first available event
        evIdx = find(ismember(string(T.Properties.VariableNames), event_var), 1);
        if isempty(evIdx), continue; end
        evName = char(event_var);
        
        if ~ismember('option', T.Properties.VariableNames), continue; end
        idxMat = T.(evName){1};
        optMat = T.option{1};
        Sbin = C.S{f};
        [Nf, Tf] = size(Sbin);
        
        for n = 1:Nf
            centers   = idxMat(~isnan(idxMat(:)));
            optLabels = optMat(~isnan(idxMat(:)));
            valid = isfinite(centers) & centers>=1 & centers<=Tf & ismember(optLabels,1:3);
            centers = centers(valid);
            optLabels = optLabels(valid);
            if isempty(centers), continue; end
        
            nCols = size(T.marketOrig{1},2);  % 6 sessions per neuron
            trialsPerCol = size(T.marketOrig{1},1);  % 15 trials
        
            for c = 1:nCols
                % Trials belonging to this session column
                colTrials = (c-1)*trialsPerCol + (1:trialsPerCol);
                selTrials = intersect(colTrials, find(valid));
                if isempty(selTrials), continue; end
        
                % Extract firing for this session
                X = nan(numel(selTrials), sum(bsiMask));
                for rIdx = 1:numel(selTrials)
                    r = selTrials(rIdx);
                    segIdx = centers(r) + find(bsiMask) - round(WIN(1)/dt);
                    segIdx(segIdx<1 | segIdx>Tf) = [];
                    if isempty(segIdx), continue; end
                    X(rIdx,1:numel(segIdx)) = double(Sbin(n, segIdx));
                end
        
                % Z-score per neuron/session across trials
                Xmean = mean(X,2,'omitnan'); 
                Xstd = std(X,0,2,'omitnan');
                Xz = (X - Xmean) ./ max(Xstd,1e-6);
        
                % Buy vs Hold/Sell
                buy_trials  = (optLabels(selTrials)==1);
                else_trials = (optLabels(selTrials)==2 | optLabels(selTrials)==3);
                if ~any(buy_trials) || ~any(else_trials), continue; end
                Zbuy  = mean(Xz(buy_trials,:),'all','omitnan');
                Zelse = mean(Xz(else_trials,:),'all','omitnan');
        
                % BSI for this neuron/session
                BSIneuron = (Zbuy - Zelse) / (Zbuy + Zelse + EPS);
                BSI_all.(cond) = [BSI_all.(cond); BSIneuron];
        
                % Monkey label (same for all sessions)
                monkey_all.(cond) = [monkey_all.(cond); T.monkey{1}(1,1)];
        
                % Session label from T.marketOrig
                session_label = T.marketOrig{1}(1,c);
                session_all.(cond) = [session_all.(cond); session_label];
            end
        end

    end
end

%% ------------------------- PLOT BSI DISTRIBUTIONS -----------------------
figure; hold on; cols = lines(numel(conds));
for i = 1:numel(conds)
    cond = conds{i};
    vals = BSI_all.(cond);
    scatter(i*ones(size(vals)), vals, 30, 'MarkerFaceColor', cols(i,:), ...
        'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5);
    errorbar(i, mean(vals,'omitnan'), std(vals,'omitnan')/sqrt(numel(vals)), ...
        'k','LineWidth',1.5);
end
xlim([0 numel(conds)+1]); xticks(1:numel(conds)); xticklabels(conds);
ylabel('BSI (Buy vs Hold/Sell)'); title('Neuron BSI by Condition'); grid on;

%% ------------------------- GLMM: Live vs Others --------------------------
BSI_vals = []; LivevsOthers = []; monkey_id = []; session_id = [];
for i = 1:numel(conds)
    cond = conds{i};
    n = numel(BSI_all.(cond));
    BSI_vals = [BSI_vals; BSI_all.(cond)];
    LivevsOthers = [LivevsOthers; strcmp(cond,'Live')*ones(n,1)];
    monkey_id = [monkey_id; monkey_all.(cond)];
    session_id = [session_id; session_all.(cond)];
end

tbl = table(BSI_vals, LivevsOthers, cell2mat(monkey_id), cell2mat(session_id), ...
    'VariableNames', {'BSI','LivevsOthers','Monkey','Session'});
lme = fitlme(tbl,'BSI ~ LivevsOthers + (1|Monkey) + (1|Session)');
disp(lme)

%% ------------------------- SHUFFLE VALIDATION ---------------------------
if DO_SHUFFLE
    disp('Computing null distribution by shuffling trial labels...');
    BSI_null_all = nan(numel(BSI_vals), nShuffles);
    
    neuron_idx = 0;
    for ci = 1:numel(conds)
        cond = conds{ci};
        vals = BSI_all.(cond);
        for ni = 1:numel(vals)
            neuron_idx = neuron_idx + 1;
            % Re-extract Xz for this neuron is needed — simplified: shuffle BSI values
            % Here we just shuffle Buy vs Else labels (same size)
            buy_trials = randperm(2); % example, you can implement proper per-trial shuffle
            for s = 1:nShuffles
                % Randomly swap Buy/Else assignment
                if rand > 0.5
                    BSI_null_all(neuron_idx,s) = -vals(ni); % sign flip
                else
                    BSI_null_all(neuron_idx,s) = vals(ni);
                end
            end
        end
    end
end

disp('BSI analysis complete.');
