% Combined correlation analysis: Z_fr vs Opponent Portfolio AND Subject Portfolio
% Fully annotated version that runs both analyses side-by-side and saves figures.
% -----------------------------------------------------------------------------
clear; clc; close all;

% ========================== USER SETTINGS =====================================
% Path to condition packs
OUTDIR = 'C:\Users\plattlab\MSM\outputs_local\condition_packs_250922';

% If true, collapse all non-Live trials into a single "not_Live" bin
LIVE_NON_LIVE_ONLY = 0;

% Events and windows of interest (ms)
TARGET_EVTS     = {'m1rew1on', 'm1rew2on'};
TARGET_EVT_WINS = {[-250, 250], [-250, 250]}; % ms

% Condition sets to analyze (edit as needed)
% cond_sets = { ...
%     {'AI','Replay','Decoy','Live'}; ...
%     {'OT AI','OT Replay','OT Decoy','OT Live'}; ...
%     {'Saline AI','Saline Replay','Saline Decoy','Saline Live'} ...
% };
% setNames = {'baseline','OT','Saline'};

cond_sets = { {'AI','Replay','Decoy','Live'} };
setNames  = {'baseline'};

% Where to save figures
figdir = fullfile(pwd, 'figs', 'corr');
if ~exist(figdir, 'dir'), mkdir(figdir); end

% =============================================================================
allResults = struct();

for pair = 1:numel(TARGET_EVTS)
    TARGET_EVT     = TARGET_EVTS{pair};
    TARGET_EVT_WIN = TARGET_EVT_WINS{pair};

    for s = 1:numel(cond_sets)
        conds = cond_sets{s};
        fprintf('\n=== Analyzing Set: %s | Event: %s ===\n', setNames{s}, TARGET_EVT);
        
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
    
        %% ---- Build trial-level table ----
        % Each row = one trial; carries covariates + aligned neural features
        trialData = table(); 
        PBuy = struct(); %#ok<NASGU>
        for c = 1:numel(condNames)
            cname = condNames{c};
            C = Cstruct.(cname);
    
            for f = 1:numel(C.eventTables)
                etab = C.eventTables{f};
                if ~ismember('option', etab.Properties.VariableNames)
                    continue; % skip if option column not present
                end
                if isempty(etab.option{1})
                    continue; % skip if option column exists but is empty
                end

                % --- Neural firing feature (Z_fr): event-locked and baseline z-scored
                Z_fr = get_Z_fr(C, f, TARGET_EVT, [1 2], TARGET_EVT_WIN, 'm1rew1on', [-1000, -700]);
    
                % --- Subject choices (option) as UPPER strings: BUY/HOLD/SELL
                subjChoice = mapChoicesToLabels(etab.option{1});  % TR × 1 string
    
                % Compute p(BUY) as a quick sanity check (stored but unused here)
                valid = ~ismissing(subjChoice);
                pBuy = mean(subjChoice(valid) == "BUY");
                PBuy.(cname)(f,1) = pBuy; %#ok<STRNU>
    
                subjBuy = subjChoice(valid) == "BUY";
    
                % Condition / flags
                condLabel = repmat(string(cname), sum(valid), 1);
                if LIVE_NON_LIVE_ONLY
                    condLabel(condLabel ~= "Live") = "not_Live";
                end
    
                % Covariates (carry through as-is)
                monkey    = string(etab.monkey{1}(:));
                session   = double(etab.session{1}(:));
                trial     = double(etab.trialNum{1}(:));
                market    = string(etab.marketOrig{1}(:));
                bubble    = string(etab.bubbleMarket{1}(:));
                portfolio = double(etab.prePortfolio{1}(:));           % SUBJECT portfolio
    
                % Opponent choice mapped to UPPER strings and time-aligned
                oppRaw            = mapChoicesToLabels(etab.optionOpp{1});
                opponentChoice    = [missing; oppRaw(1:end-1)]; % opponent responds to the subject
                opponentBuy       = opponentChoice == "BUY";
                opponentPortfolio = etab.postPortfolioOpp{1}(:);
                opponentPortfolio = [missing; opponentPortfolio(1:end-1)]; % opponent responds to subject
    
                % Previous subject choice and reward (shifted)
                prevChoice        = [missing; subjChoice(1:end-1)];
                prevBuy           = prevChoice == "BUY";
                choiceReward      = etab.sizeB{1}(:) .* etab.divPerShare{1}(:);
                prevChoiceReward  = [missing; choiceReward(1:end-1)];
    
                % Assemble row-block (strings now; convert to categorical later)
                tmpT = table( ...
                    subjBuy, prevBuy, opponentBuy, ...
                    condLabel(valid), string(monkey(valid)), session(valid), trial(valid), ...
                    string(market(valid)), string(bubble(valid)), ...
                    string(subjChoice(valid)), string(opponentChoice(valid)), string(prevChoice(valid)), ...
                    choiceReward(valid), prevChoiceReward(valid), ...
                    portfolio(valid), opponentPortfolio(valid), ...
                    Z_fr(valid), ...
                    'VariableNames', {'subjBuy', 'prevBuy', 'opponentBuy', ...
                        'Condition','Monkey','Session','Trial','Market','Bubble', ...
                        'subjChoice','opponentChoice','prevChoice', ...
                        'choiceReward', 'prevChoiceReward', ...
                        'portfolio', 'opponentPortfolio', ...
                        'Z_fr'});
    
                trialData = [trialData; tmpT]; %#ok<AGROW>
            end
        end
        
        % ---- Categorical encodings ----
        trialData.Monkey = categorical(trialData.Monkey);
        trialData.Market = categorical(trialData.Market);
        trialData.Bubble = categorical(trialData.Bubble);
        
        % Ensure categorical predictors use HOLD as baseline (first level)
        if LIVE_NON_LIVE_ONLY
            trialData.Condition = setcats(categorical(trialData.Condition), {'not_Live','Live'});        
        else
            trialData.Condition = setcats(categorical(trialData.Condition), {'AI','Replay','Decoy','Live'});
        end
        trialData.subjChoice     = setcats(categorical(trialData.subjChoice),     {'HOLD','BUY','SELL'});
        trialData.opponentChoice = setcats(categorical(trialData.opponentChoice), {'HOLD','BUY','SELL'});
        trialData.prevChoice     = setcats(categorical(trialData.prevChoice),     {'HOLD','BUY','SELL'});
    
        % ---- Keep trials that have neural data (and opponent portfolio for first analysis)
        trialData = trialData(~isnan(trialData.Z_fr) & ~isnan(trialData.opponentPortfolio), :);

        % ---- Partition into condition groups ----
        trialDataAI       = trialData(trialData.Condition == "AI", :);
        trialDataReplay   = trialData(trialData.Condition == "Replay", :);
        trialDataDecoy    = trialData(trialData.Condition == "Decoy", :);
        trialDataLive     = trialData(trialData.Condition == "Live", :);

        trialDataNonLive  = trialData(trialData.Condition ~= "Live", :);
        trialDataSocial   = trialData(trialData.Condition == "Live" | trialData.Condition == "Decoy", :);
        trialDataNonSocial= trialData(trialData.Condition == "AI"   | trialData.Condition == "Replay", :);

        condition_grouped = {
            % trialDataAI;
            % trialDataReplay;
            % trialDataDecoy;
            trialDataLive;
            trialDataNonLive;
            trialDataSocial;
            trialDataNonSocial;
        };
        condition_grouped_names = [
            % "AI";
            % "Replay";
            % "Decoy";
            "Live";
            "Non-Live";
            "Social";
            "Non-social";
        ];

        % ==================== ANALYSIS A: Opponent Portfolio ====================
        % Correlations of Z_fr with opponentPortfolio, by choice class.
        choice_labels = ["Buy","Hold","Sell","Non-buy"];
        rhos_opp      = zeros(numel(condition_grouped), 4); % BUY, HOLD, SELL, NON-BUY
        pvals_rho_opp = zeros(numel(condition_grouped), 4);
        zs_opp        = zeros(numel(condition_grouped), 1); % Fisher z for BUY vs NON-BUY
        pvals_z_opp   = zeros(numel(condition_grouped), 1);

        for gp = 1:numel(condition_grouped)
            curr_gp      = condition_grouped{gp};
            buy_trials   = curr_gp(curr_gp.subjChoice == "BUY", :);
            hold_trials  = curr_gp(curr_gp.subjChoice == "HOLD", :);
            sell_trials  = curr_gp(curr_gp.subjChoice == "SELL", :);
            nonbuy_trials= curr_gp(~curr_gp.subjBuy, :);

            % Use 'Rows','complete' to safely handle any residual NaNs per pair
            [rhos_opp(gp,1), pvals_rho_opp(gp,1)] = corr(buy_trials.Z_fr,   buy_trials.opponentPortfolio, "Type","Pearson", "Rows","complete");
            [rhos_opp(gp,2), pvals_rho_opp(gp,2)] = corr(hold_trials.Z_fr,  hold_trials.opponentPortfolio,"Type","Pearson", "Rows","complete");
            [rhos_opp(gp,3), pvals_rho_opp(gp,3)] = corr(sell_trials.Z_fr,  sell_trials.opponentPortfolio,"Type","Pearson", "Rows","complete");
            [rhos_opp(gp,4), pvals_rho_opp(gp,4)] = corr(nonbuy_trials.Z_fr,nonbuy_trials.opponentPortfolio,"Type","Pearson","Rows","complete");

            n_buy    = sum(curr_gp.subjChoice=="BUY" & ~isnan(curr_gp.Z_fr) & ~isnan(curr_gp.opponentPortfolio));
            n_nonbuy = sum(~curr_gp.subjBuy           & ~isnan(curr_gp.Z_fr) & ~isnan(curr_gp.opponentPortfolio));
            [zs_opp(gp), pvals_z_opp(gp)] = corr_diff_indep(rhos_opp(gp,1), n_buy, rhos_opp(gp,4), n_nonbuy, 'two');
        end

        % ---- Plot: Opponent Portfolio ----
        fh1 = figure('Color','w');        
        ax1 = subplot(1, 2, 1);
        imagesc(rhos_opp);
        colorbar;
        n = 256; % number of colors
        rwb = [linspace(0,1,n/2)' linspace(0,1,n/2)' ones(n/2,1);  % blue→white
               ones(n/2,1) linspace(1,0,n/2)' linspace(1,0,n/2)']; % white→red
        colormap(ax1, rwb);        
        clim([-1 1]);
        set(ax1,'XTick',1:numel(choice_labels),'XTickLabel',choice_labels, ...
            'YTick',1:numel(condition_grouped_names),'YTickLabel',condition_grouped_names, ...
            'TickDir','out','LineWidth',1);
        xlabel('Choice'); ylabel('Condition');
        title(['Correlation: Z_fr @ ' TARGET_EVT ' vs Opponent Portfolio'], 'Interpreter','none');

        % Annotate cells with r and significance stars
        [nRows,nCols] = size(rhos_opp);
        for i = 1:nRows
            for j = 1:nCols
                v = rhos_opp(i,j);
                p = pvals_rho_opp(i,j);
                star = pstars(p);
                lbl = sprintf('%.3f (%s)', v, star);
                tc = 'k';
                if ~isnan(v) && p < 0.05, tc = 'r'; end
                text(j, i, lbl, 'HorizontalAlignment','center', 'FontWeight','bold', 'Color', tc, 'FontSize', 11);
            end
        end

        subplot(1, 2, 2); hold on
        bar(1:numel(condition_grouped_names), zs_opp);
        for x = 1:numel(condition_grouped_names)
            tc = 'k';
            if pvals_z_opp(x) < 0.05, tc = 'r'; end
            text(x, zs_opp(x)+signOrOne(zs_opp(x))*0.3, pstars(pvals_z_opp(x)), 'HorizontalAlignment','center', 'FontWeight','bold', 'Color', tc, 'FontSize', 14);
        end
        xticks(1:numel(condition_grouped_names))
        xticklabels(condition_grouped_names)
        ylabel("z")
        ylim([-5 5])
        title([string(['Correlation: Z_fr @ ' TARGET_EVT ' vs Opponent Portfolio']); ...
               "Fisher r-to-z test for BUY - NON-BUY"], 'Interpreter',"none")
        hold off
        
        % Save Opponent Portfolio figure
        exportgraphics(fh1, fullfile(figdir, [TARGET_EVT '_fr_vs_opponPort_corr.svg']), 'Resolution', 300);
        savefig(fh1,      fullfile(figdir, [TARGET_EVT '_fr_vs_opponPort_corr.fig']));

        % ==================== ANALYSIS B: Subject Portfolio =====================
        % Correlations of Z_fr with SUBJECT's own portfolio, by choice class.
        rhos_subj      = zeros(numel(condition_grouped), 4); % BUY, HOLD, SELL, NON-BUY
        pvals_rho_subj = zeros(numel(condition_grouped), 4);
        zs_subj        = zeros(numel(condition_grouped), 1); % Fisher z for BUY vs NON-BUY
        pvals_z_subj   = zeros(numel(condition_grouped), 1);

        for gp = 1:numel(condition_grouped)
            curr_gp       = condition_grouped{gp};
            buy_trials    = curr_gp(curr_gp.subjChoice == "BUY", :);
            hold_trials   = curr_gp(curr_gp.subjChoice == "HOLD", :);
            sell_trials   = curr_gp(curr_gp.subjChoice == "SELL", :);
            nonbuy_trials = curr_gp(~curr_gp.subjBuy, :);

            % Pairwise-complete handling of NaNs for subject portfolio
            [rhos_subj(gp,1), pvals_rho_subj(gp,1)] = corr(buy_trials.Z_fr,    buy_trials.portfolio, "Type","Pearson", "Rows","complete");
            [rhos_subj(gp,2), pvals_rho_subj(gp,2)] = corr(hold_trials.Z_fr,   hold_trials.portfolio,"Type","Pearson", "Rows","complete");
            [rhos_subj(gp,3), pvals_rho_subj(gp,3)] = corr(sell_trials.Z_fr,   sell_trials.portfolio,"Type","Pearson", "Rows","complete");
            [rhos_subj(gp,4), pvals_rho_subj(gp,4)] = corr(nonbuy_trials.Z_fr, nonbuy_trials.portfolio,"Type","Pearson","Rows","complete");

            n_buy_subj    = sum(curr_gp.subjChoice=="BUY" & ~isnan(curr_gp.Z_fr) & ~isnan(curr_gp.portfolio));
            n_nonbuy_subj = sum(~curr_gp.subjBuy           & ~isnan(curr_gp.Z_fr) & ~isnan(curr_gp.portfolio));
            [zs_subj(gp), pvals_z_subj(gp)] = corr_diff_indep(rhos_subj(gp,1), n_buy_subj, rhos_subj(gp,4), n_nonbuy_subj, 'two');
        end

        % ---- Plot: Subject Portfolio ----
        fh2 = figure('Color','w');        
        ax2 = subplot(1, 2, 1);
        imagesc(rhos_subj);
        colorbar;
        colormap(ax2, rwb);
        clim([-1 1]);
        set(ax2,'XTick',1:numel(choice_labels),'XTickLabel',choice_labels, ...
            'YTick',1:numel(condition_grouped_names),'YTickLabel',condition_grouped_names, ...
            'TickDir','out','LineWidth',1);
        xlabel('Choice'); ylabel('Condition');
        title(['Correlation: Z_fr @ ' TARGET_EVT ' vs Subject Portfolio'], 'Interpreter','none');

        % Annotate cells
        [nRows2,nCols2] = size(rhos_subj);
        for i = 1:nRows2
            for j = 1:nCols2
                v = rhos_subj(i,j);
                p = pvals_rho_subj(i,j);
                star = pstars(p);
                lbl = sprintf('%.3f (%s)', v, star);
                tc = 'k';
                if ~isnan(v) && p < 0.05, tc = 'r'; end
                text(j, i, lbl, 'HorizontalAlignment','center', 'FontWeight','bold', 'Color', tc, 'FontSize', 11);
            end
        end

        subplot(1, 2, 2); hold on
        bar(1:numel(condition_grouped_names), zs_subj);
        for x = 1:numel(condition_grouped_names)
            tc = 'k';
            if pvals_z_subj(x) < 0.05, tc = 'r'; end
            text(x, zs_subj(x)+signOrOne(zs_subj(x))*0.3, pstars(pvals_z_subj(x)), 'HorizontalAlignment','center', 'FontWeight','bold', 'Color', tc, 'FontSize', 14);
        end
        xticks(1:numel(condition_grouped_names))
        xticklabels(condition_grouped_names)
        ylabel("z")
        ylim([-8 8])
        title([string(['Correlation: Z_fr @ ' TARGET_EVT ' vs Subject Portfolio']); ...
               "Fisher r-to-z test for BUY - NON-BUY"], 'Interpreter',"none")
        hold off
        
        % Save Subject Portfolio figure
        exportgraphics(fh2, fullfile(figdir, [TARGET_EVT '_fr_vs_subjPort_corr.svg']), 'Resolution', 300);
        savefig(fh2,      fullfile(figdir, [TARGET_EVT '_fr_vs_subjPort_corr.fig']));

        % ---- (Optional) store outputs into allResults ----
        resKey = matlab.lang.makeValidName([setNames{s} '_' TARGET_EVT]);
        allResults.(resKey).rhos_opp    = rhos_opp;
        allResults.(resKey).pvals_opp   = pvals_rho_opp;
        allResults.(resKey).zs_opp      = zs_opp;
        allResults.(resKey).pvals_z_opp = pvals_z_opp;

        allResults.(resKey).rhos_subj    = rhos_subj;
        allResults.(resKey).pvals_subj   = pvals_rho_subj;
        allResults.(resKey).zs_subj      = zs_subj;
        allResults.(resKey).pvals_z_subj = pvals_z_subj;
    end
end

%% ---------- Helpers ----------

% Fisher r-to-z test for difference between two independent correlations
function [z,p] = corr_diff_indep(r1,n1,r2,n2,tail)
    if nargin<5, tail = 'two'; end
    if any([n1 n2] < 4) || any(abs([r1 r2])>=1)
        z = NaN; p = NaN; return;
    end
    z1 = atanh(r1);
    z2 = atanh(r2);
    se = sqrt(1/(n1-3) + 1/(n2-3));
    z = (z1 - z2) / se;
    Phi = @(x) 0.5*erfc(-x/sqrt(2)); % standard normal CDF
    switch lower(tail)
        case 'two',   p = 2*Phi(-abs(z));
        case 'right', p = 1 - Phi(z); % H1: r1 > r2
        case 'left',  p = Phi(z);     % H1: r1 < r2
        otherwise,    p = 2*Phi(-abs(z));
    end
end

function s = pstars(p)
    if isnan(p), s = ''; return; end
    if p < 1e-4,      s = '****';
    elseif p < 1e-3,  s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'ns';
    end
end

function out = mapChoicesToLabels(x)
% Returns a string array with values "BUY","HOLD","SELL" (or missing)
    if iscell(x); x = x(:); end
    if isnumeric(x)
        out = strings(numel(x),1);
        out(:) = missing;
        out(x==1) = "BUY";
        out(x==2) = "HOLD";
        out(x==3) = "SELL";
    elseif isstring(x)
        out = upper(x(:));
    elseif iscellstr(x)
        out = upper(string(x(:)));
    else
        out = strings(numel(x),1); out(:) = missing;
    end
end

function fr = get_Z_fr(C, f, target_event, neuron_type, target_win, baseline_event, baseline_win)
% Extract per-trial firing rate (FR) from S, z-scored against a baseline window.
% - target_event: string (field name in eventTables), e.g., 'm1rew1on'
% - neuron_type: vector of int codes (e.g., [1 2]) to include
% - target_win/baseline_win: [t_start, t_end] in ms
% Returns: fr (TR×1) mean z-scored FR across selected neuron types.
    target_interval   = (target_win   ./ 1000) ./ C.dt;
    baseline_interval = (baseline_win ./ 1000) ./ C.dt;
    fr = [];
    try
        target_anchor   = C.eventTables{f}.(target_event){1,1}(:);
        target_ts       = target_anchor + target_interval(1); 
        target_te       = target_anchor + target_interval(2);
        baseline_anchor = C.eventTables{f}.(baseline_event){1,1}(:);
        baseline_ts     = baseline_anchor + baseline_interval(1); 
        baseline_te     = baseline_anchor + baseline_interval(2);
    catch
        fr = nan(90, 1);
        return
    end
    for i = 1:numel(neuron_type)
        mask  = C.neuronType{f} == neuron_type(i);
        Ssel  = C.S{f}(mask,:);
        ntr   = numel(target_ts); 
        S_interval = nan(sum(mask), ntr);

        % baseline mean and std
        baseline_values = nan(sum(mask), ntr);
        for t = 1:ntr
            a = baseline_ts(t); b = baseline_te(t);
            if ~isnan(a) && ~isnan(b)
                baseline_values(:,t) = sum(Ssel(:,a:b),2) ./ ((baseline_win(2) - baseline_win(1))/1000/C.dt);
            end
        end
        baseline_mean = mean(baseline_values, 2, 'omitnan');
        baseline_std  = std(baseline_values, 0,  2, 'omitnan');

        for t = 1:ntr
            a = target_ts(t); b = target_te(t);
            if ~isnan(a) && ~isnan(b) 
                S_interval(:,t) = sum(Ssel(:,a:b),2) ./ ((target_win(2) - target_win(1))/1000/C.dt);
                S_interval(:,t) = (S_interval(:,t) - baseline_mean) ./ baseline_std;
            end
        end
        fr = [fr; S_interval]; %#ok<AGROW>
    end
    fr = mean(fr, 1, 'omitnan');
    fr = fr(:);
end

function s = signOrOne(x)
% Utility: returns sign(x) but coerces 0 to +1 (for text offset purposes)
    if x == 0 || isnan(x)
        s = 1;
    else
        s = sign(x);
    end
end
