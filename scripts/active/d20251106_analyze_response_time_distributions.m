% File: d20251106_analyze_response_time_distributions.m
clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tResp_sweep_cSup_heatmap_251106.txt  (MATLAB script)
% -----------------------------------------------------------------------------
% Sweep both SHIFT and FLOOR (with FLOOR > SHIFT) and fit:
%   log(tResp_corr) ~ cSup + (1|monkey) + (1|session)
% Produce heat maps of cSup fixed-effect estimates (Replay, Decoy, Live)
% across the (shift, floor) grid, with text overlays (stars by significance).
%
% Save this file as .m if you want to run directly (rename .txt -> .m).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ================================ USER SETTINGS ================================
dataFile  = 'C:\Users\plattlab\MSM\data\tooling\MATLAB_Copy\dataTabFULL.xlsx';
sheetName = 'Sheet1';

% Sweep values (seconds)
shiftsSec  = 0:0.05:0.3;   % global shifts added to tResp
floorsSec  = 0:0.05:0.5;   % hard floors applied after shift (must be > shift)

% Core analysis toggles
useOddTrialBlocksOnly = true;     % Player1 only (odd trialBlock)
excludeMonkey3        = true;     % drop monkey==3 if present
doSessionTrim         = true;     % session-aware 1–99% trimming before LME
maxRT_cap             = 5.00;     % upper cap (s) for trimming

makeFigs              = 1;        % make heat maps
figTitlePrefix        = 'cSup effects: log RT LME';  % interpreter=none enforced below

%% ================================ LOAD DATA ====================================
fprintf('\n[LOAD] Reading table: %s (sheet %s)\n', dataFile, sheetName);
T = readtable(dataFile, 'Sheet', sheetName);

needCols = {'tResp','condition','trialBlock','monkey','session','time'};
assert(all(ismember(needCols, T.Properties.VariableNames)), ...
    'Input table must contain columns: %s', strjoin(needCols, ', '));

% Super-condition (integer 1..4): AI=1, Replay=2, Decoy=3, Live=4
if ~ismember('supCond', T.Properties.VariableNames)
    T.supCond = floor(T.condition);
end

% Odd blocks only (subject / player 1)
if useOddTrialBlocksOnly
    T = T(mod(T.trialBlock,2)==1, :);
end

% Exclude monkey 3
if excludeMonkey3
    T = T(T.monkey ~= 3, :);
end

% Types / convenience
T.tResp   = safe_to_seconds(T.tResp);
T.sessKey = build_session_key(T.time, T.session);
T.cSup    = categorical(T.supCond, 1:4, {'AI','Replay','Decoy','Live'});
T.cMky    = categorical(T.monkey);

fprintf('[LOAD] N trials after prep: %d\n', height(T));
if height(T)==0, error('No trials after prep.'); end

% ---- Estimate hardware latency shift from UNFLOORED/UNSHIFTED tResp ----
% Adds back the helper to estimate s* that best normalizes log(tResp + s).
capForShift = max(shiftsSec);  % bound the search to the sweep's max shift
sHat = estimate_shift_lognormal(T.tResp, capForShift);
fprintf('[SHIFT*] estimated by lognormal fitting: %.0f ms (cap %.0f ms)\n', 1000*sHat, 1000*capForShift);

% LME availability
hasLME = (exist('fitlme','file')==2);

%% ============================ PREP SWEEP CONTAINERS ============================
S = numel(shiftsSec); F = numel(floorsSec);
levs = {'Replay','Decoy','Live'}; L = numel(levs);

betaMat   = nan(S, F, L);   % fixed-effect estimate for each cSup level vs AI baseline
pvalMat   = nan(S, F, L);   % p-values for each effect
nUsedMat  = nan(S, F);      % #obs used in the fit
failMat   = false(S, F);    % fit failure flags

%% ================================== SWEEP =====================================
for si = 1:S
    shift = shiftsSec(si);
    for fi = 1:F
        floorV = floorsSec(fi);
        if floorV < shift
            % invalid combo (floor has no effect if <= shift)
            continue;
        end

        % Apply transform + hard floor
        Twork = T;
        Twork.tResp_corr = Twork.tResp + shift;
        Twork = Twork(Twork.tResp >= floorV, :);
        if isempty(Twork)
            fprintf('[SKIP] shift=%g floor=%g -> no trials after floor.\n', shift, floorV);
            continue;
        end

        % Session-aware trimming (1..99 pct) to reduce extreme leverage
        if doSessionTrim
            [G,~] = findgroups(Twork.sessKey);
            Q = splitapply(@(x) prctile(x(isfinite(x)), [1 99]), Twork.tResp_corr, G);
            qL = Q(:,1); qH = Q(:,2);
            Lb = max(qL(G), floorV);
            Ub = min(qH(G), maxRT_cap);
            keep = Twork.tResp_corr >= Lb & Twork.tResp_corr <= Ub;
            Twork = Twork(keep,:);
        end

        if height(Twork) < 5
            fprintf('[SKIP] shift=%g floor=%g -> too few trials after trim.\n', shift, floorV);
            continue;
        end

        % Fit LME: logRT ~ cSup + (1|monkey) + (1|session)
        if hasLME
            try
                Twork.cSess = categorical(Twork.sessKey);
                mdl = fitlme(Twork, 'log(tResp_corr) ~ cSup + (1|cMky) + (1|cSess)');
                nUsedMat(si,fi) = height(mdl.Variables);

                % Pull fixed effects for cSup_*; baseline is AI
                coef = mdl.Coefficients;
                for li = 1:L
                    nm = ['cSup_' levs{li}];
                    row = strcmp(coef.Name, nm);
                    if any(row)
                        betaMat(si,fi,li) = coef.Estimate(row);
                        pvalMat(si,fi,li) = coef.pValue(row);
                    else
                        betaMat(si,fi,li) = NaN;
                        pvalMat(si,fi,li) = NaN;
                    end
                end
            catch ME
                failMat(si,fi) = true;
                fprintf('[FAIL] shift=%g floor=%g : %s\n', shift, floorV, ME.message);
                continue;
            end
        else
            warning('fitlme not available. Aborting sweep.');
            break;
        end
    end
end

%% ================================ HEAT MAPS ====================================
if makeFigs
    % Build labels
    xLabels = compose('%d', round(1000*floorsSec)); % floors in ms
    yLabels = compose('%d', round(1000*shiftsSec)); % shifts in ms

    % Determine shared CLim across panels (ignore NaNs)
    allVals = betaMat(:);
    allVals = allVals(isfinite(allVals));
    if isempty(allVals)
        climVals = [-1 1];
    else
        q = quantile(abs(allVals), 0.95);
        q = max(q, 0.05);
        climVals = [-q q];
    end

    figure('Color','w','Name','cSup Estimates Heatmaps');
    tlo = tiledlayout(1, L, 'TileSpacing','compact', 'Padding','compact');

    % Compute y-position (index space) for the estimated shift line
    yLine = interp1(shiftsSec, 1:S, sHat, 'linear', 'extrap');

    for li = 1:L
        nexttile;
        M = squeeze(betaMat(:,:,li));  % S x F
        imagesc(1:F, 1:S, M, 'AlphaData', ~isnan(M));
        axis equal tight;
        % Custom diverging colormap centered at 0
        nC = 100;
        cmap = interp1([0.3 0 -0.3], [1 0.3 0.3; 1 1 1; 0.3 0.3 1], linspace(-0.3,0.3,nC), 'linear');
        colormap(cmap);
        clim(climVals);  % ensure symmetric scaling around 0
        cb = colorbar; set(cb.Label,'String','Estimate (log-RT vs AI)','Interpreter','none');
        set(gca,'XTick',1:F,'XTickLabel',xLabels);
        set(gca,'YTick',1:S,'YTickLabel',yLabels);
        xlabel('FLOOR (ms)','Interpreter','none');
        ylabel('SHIFT (ms)','Interpreter','none');
        ttl = title(sprintf('%s — %s', figTitlePrefix, levs{li})); set(ttl,'Interpreter','none');

        % Overlay significance stars instead of numeric values
        hold on;
        for si = 1:S
            for fi = 1:F
                if fi<=numel(floorsSec) && si<=numel(shiftsSec)
                    p = pvalMat(si,fi,li);
                    if isfinite(p)
                        if     p < 0.001, stars = '***';
                        elseif p < 0.01,  stars = '**';
                        elseif p < 0.05,  stars = '*';
                        else,  stars = '';
                        end
                        if ~isempty(stars)
                            text(fi, si, stars, 'HorizontalAlignment','center', ...
                                 'VerticalAlignment','middle', 'Color','k', ...
                                 'FontSize',8, 'FontWeight','bold', 'Interpreter','none');
                        end
                    end
                end
            end
        end

        % --- Dashed horizontal line at the estimated shift (in index coords) ---
        if isfinite(yLine) && yLine >= 1 && yLine <= S
            plot([1 F], [yLine yLine], 'k--', 'LineWidth', 1.2);
            % annotate on the left margin
            text(1, yLine, sprintf('  s*≈%d ms', round(1000*sHat)), ...
                 'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                 'Color','k','FontSize',8,'Interpreter','none');
        end
        hold off;

        set(gca,'YDir','normal');    % data pre-flipped -> keep y increasing upward
        set(gca,'XDir','reverse');   % visually right-to-left floors
    end
    t = title(tlo, 'cSup fixed-effect estimates across SHIFT×FLOOR');
    set(t,'Interpreter','none');
end

fprintf('\n[NOTE] Estimates are on log scale; exp(beta) gives multiplicative change vs AI.\n');
fprintf('[DONE] Sweep complete.\n');

%% ================================ HELPERS ======================================
function x = safe_to_seconds(x)
    if isduration(x)
        x = seconds(x);
    elseif iscell(x)
        try
            x = seconds(duration(x));
        catch
            x = str2double(x);
        end
    elseif isdatetime(x)
        x = seconds(x - x(1));
    end
    x = double(x);
end

function key = build_session_key(t, s)
    if isdatetime(t)
        d = dateshift(t,'start','day');
        key = string(datestr(d,'yyyy-mm-dd')) + "_" + string(s);
    else
        key = string(t) + "_" + string(s);
    end
    key = categorical(key);
end

% Significance color map for text overlays (by p-value)
% p<0.001 red, p<0.01 magenta, p<0.05 blue, else gray
function c = p2color(p) %#ok<DEFNU>
    if isnan(p), c = [0.5 0.5 0.5];
    elseif p < 1e-3, c = [0.85 0.1 0.1];
    elseif p < 1e-2, c = [0.7 0.1 0.7];
    elseif p < 5e-2, c = [0.1 0.2 0.8];
    else, c = [0.5 0.5 0.5];
    end
end

function s = estimate_shift_lognormal(y, capSec)
%ESTIMATE_SHIFT_LOGNORMAL  Nonnegative latency shift s to stabilize low tail.
%   s = argmin_{0<=s<=capSec} NLL( log(y+s) ) under a Gaussian model.
    if nargin<2 || isempty(capSec), capSec = 0.25; end
    y = double(y(:));
    y = y(isfinite(y) & y>0);
    if isempty(y), s = 0; return; end

    function val = obj(s0)
        if s0 < 0 || s0 > capSec, val = Inf; return; end
        z = log(y + s0);
        mu = mean(z);
        sg = std(z) + eps;
        val = numel(z)*log(sg) + 0.5*sum(((z-mu)./sg).^2);
    end

    try
        s = fminbnd(@obj, 0, capSec);
    catch
        gs = linspace(0,capSec,161);
        vals = arrayfun(@obj, gs);
        [~,ix] = min(vals); s = gs(ix);
    end
end
