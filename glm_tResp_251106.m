%% tResp LME: customizable model formulas (no block/trialInBlock in models)
% Self-contained pipeline to load data, build previous-choice regressors,
% transform RT, fit a LIST of user-defined formulas, and compare fits.
% Customize the model list at the TOP as strings. No random slopes unless you put them.
%
% Expected columns in input table T:
%   session (categorical/char/string)  : session identifier (used as grouping if included in formula)
%   trialBlock (int, 1..30, resets)    : within-block index; odd=Player1, even=Player2
%   tResp (numeric, ms)                : response time
%   option (numeric or categorical)    : chosen option on the trial
% Optional covariates (detected/aliased):
%   cSup (categorical)   -> aliased to 'condition' if 'condition' is not present
%   cMky (categorical)   -> aliased to 'monkey'    if 'monkey'    is not present
%
% Notes:
% - Previous choices computed WITHIN (session, blockId), where blockId resets when trialBlock==1.
% - Player 1 trials = odd trialBlock. On odd rows:
%       p2_prev = lag1 (opponent even from previous trial)
%       p1_prev = lag2 (self odd two-back)
% - Interpreter is set to 'none' for all labels/titles to avoid TeX parsing.
%
% -------------------------------------------------------------------------
%                           USER CONFIGURATION
% -------------------------------------------------------------------------
SHIFT   = 0;      % global shift added to tResp before flooring
FLOOR   = 0;      % hard floor threshold (s)
FLOOR_MODE = 'exclude';% 'clamp' or 'exclude'
MIN_CAT_CT = 2;        % drop categorical levels with < MIN_CAT_CT rows
DEFAULT_PROMPT = false;% prompt for file if T not in workspace

% >>>>>>>>>>>>>>>>>>>>>  CUSTOMIZE YOUR MODEL FORMULAS HERE  <<<<<<<<<<<<<<<<<<<<<<
% Put any valid 'fitlme' Wilkinson formula strings below.
% You can use variables: p1_prev, p2_prev, condition, monkey, session, tResp_t.
% Examples include the user's request: no slopes; random intercepts only.
FORMULAE = [ ...
    "log(tResp_t) ~ cSup + (1|monkey) + (1|session)";             % Condition only + REs
    "log(tResp_t) ~ p1_prev + p2_prev + (1|monkey) + (1|session)";                                           % History only
    "log(tResp_t) ~ p1_prev + p2_prev + (1|cSup) + (1|session)";                                           % History only
    "log(tResp_t) ~ p1_prev + p2_prev + (1|cSup:monkey) + (1|session)";                                           % History only
];  %#ok<NASGU>
% ==========================================================================

%% Global: turn off TeX interpreter for labeling
set(groot, 'defaultTextInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter', 'none'); 
set(groot, 'defaultLegendInterpreter', 'none');

%% ================================ LOAD DATA ====================================
if ~exist('T','var') || ~istable(T) || DEFAULT_PROMPT
    [fn, fp] = uigetfile({'*.mat;*.csv','MAT or CSV files (*.mat, *.csv)'}, 'Select data file containing table T');
    if isequal(fn,0), error('No file selected.'); end
    fpath = fullfile(fp, fn);
    [~,~,ext] = fileparts(fpath);
    switch lower(ext)
        case '.mat'
            S = load(fpath);
            if isfield(S,'T') && istable(S.T)
                T = S.T;
            else
                fns = fieldnames(S); got = false;
                for i=1:numel(fns)
                    if istable(S.(fns{i})), T = S.(fns{i}); got = true; break; end
                end
                if ~got, error('No table found in MAT file.'); end
            end
        case '.csv'
            T = readtable(fpath);
        otherwise
            error('Unsupported file type: %s', ext);
    end
end

% Required columns
reqCols = {'trialBlock','tResp','option'};
for k=1:numel(reqCols)
    assert(ismember(reqCols{k}, T.Properties.VariableNames), 'Missing required column: %s', reqCols{k});
end

% session (for grouping if used)
if ~ismember('session', T.Properties.VariableNames)
    warning('Column "session" not found. Creating a dummy session.');
    T.session = repmat("S1", height(T), 1);
end
if ~iscategorical(T.session), T.session = categorical(string(T.session)); end

% Optional: cSup/cMky; create aliases 'condition' and 'monkey' if not present
if ismember('cSup', T.Properties.VariableNames) && ~ismember('condition', T.Properties.VariableNames)
    T.condition = T.cSup;
end
if ismember('cMky', T.Properties.VariableNames) && ~ismember('monkey', T.Properties.VariableNames)
    T.monkey = T.cMky;
end
if ismember('condition', T.Properties.VariableNames) && ~iscategorical(T.condition)
    T.condition = categorical(string(T.condition));
end
if ismember('monkey', T.Properties.VariableNames) && ~iscategorical(T.monkey)
    T.monkey = categorical(string(T.monkey));
end

% Ensure types
assert(isnumeric(T.tResp), 'tResp must be numeric (ms).');
assert(isnumeric(T.trialBlock), 'trialBlock must be numeric integers in 1..30.');
T.row = (1:height(T))';  % preserve original order

%% --------------- BUILD blockId (for lagging only) & LAG REGRESSORS ---------------
Gs = findgroups(T.session);
blockId = zeros(height(T),1);
uGs = unique(Gs);
for gi = reshape(uGs,1,[])
    idx = find(Gs==gi);
    [~,ord] = sortrows([T.row(idx) T.trialBlock(idx)], [1 2]);
    idx = idx(ord);
    b = cumsum(T.trialBlock(idx)==1);
    if any(b==0)
        if T.trialBlock(idx(1)) ~= 1
            b = 1 + cumsum([0; T.trialBlock(idx(2:end))==1]);
        end
    end
    blockId(idx) = b;
end

GB = findgroups(T.session, blockId);
lag1 = nan(height(T),1);
lag2 = nan(height(T),1);
for gb = unique(GB(:))'
    idx = find(GB==gb);
    [~,ord] = sort(T.trialBlock(idx));
    idx = idx(ord);
    oi = T.option(idx);
    if iscategorical(oi), oi = removecats(oi); end
    if numel(idx) >= 2
        lag1(idx(2:end)) = double(oi(1:end-1));
    end
    if numel(idx) >= 3
        lag2(idx(3:end)) = double(oi(1:end-2));
    end
end
T.lag1 = lag1;  % prev trial within block
T.lag2 = lag2;  % two-back within block

% Player 1 regressors on odd rows only
isOdd = mod(T.trialBlock,2)==1;
p2_prev = nan(height(T),1);   % opponent previous (even) -> lag1 on odd rows
p1_prev = nan(height(T),1);   % self two-back (odd)     -> lag2 on odd rows
p2_prev(isOdd) = T.lag1(isOdd);
p1_prev(isOdd) = T.lag2(isOdd);
T.p1_prev = categorical(p1_prev);
T.p2_prev = categorical(p2_prev);

%% ------------------------------ RT TRANSFORM ------------------------------------
t = T.tResp + SHIFT;
switch lower(FLOOR_MODE)
    case 'clamp'
        t = max(t, FLOOR);
        keepFloor = true(height(T),1);
    case 'exclude'
        keepFloor = t >= FLOOR;
    otherwise
        error('Unknown FLOOR_MODE: %s', FLOOR_MODE);
end
T.tResp_t = t;

% Final analysis set: odd trials with both prevs defined and finite RT
keep = isOdd & keepFloor & ~isundefined(T.p1_prev) & ~isundefined(T.p2_prev) & isfinite(T.tResp_t);
Tf = T(keep,:);

% Drop sparse categorical levels
Tf = drop_sparse_levels(Tf, MIN_CAT_CT);

%% ------------------------------ FIT & COMPARE -----------------------------------
% If FORMULAE is not in workspace (e.g., cleared), rebuild defaults:
if ~exist('FORMULAE','var') || isempty(FORMULAE)
    FORMULAE = [ ...
        "tResp_t ~ 1 + (1|session)";
        "tResp_t ~ p1_prev + p2_prev + (1|session)";
        "tResp_t ~ p1_prev + p2_prev + (1|monkey) + (1|condition:monkey) + (1|condition)";
        "tResp_t ~ p1_prev + p2_prev + condition + (1|monkey) + (1|condition:monkey) + (1|condition)";
        "tResp_t ~ condition + (1|monkey) + (1|condition:monkey) + (1|condition)";
    ];
end

nM = numel(FORMULAE);
results = struct([]);

for i = 1:nM
    form = FORMULAE(i);
    fprintf('\n=== Model %d ===\nFormula: %s\n', i, form);
    [lme, ok, msg] = fitModelSafe(Tf, form);
    results(i).form = form;
    results(i).ok   = ok;
    results(i).msg  = msg;
    if ok
        results(i).lme   = lme;
        results(i).anova = anova(lme,'DFMethod','Residual');
        crit = lme.ModelCriterion;
        results(i).AIC   = crit.AIC;
        results(i).BIC   = crit.BIC;
        results(i).logLik= crit.LogLikelihood;
        results(i).N     = height(lme.Variables);

        % Print a concise summary
        disp(lme);
        disp('ANOVA (Residual DF):'); disp(results(i).anova);

        % Plots: fixed effects forest & diagnostics
        try, plotFixedEffects(lme, sprintf('Model %d',i)); catch, end
        try
            figure('Name',sprintf('Diagnostics: Model %d',i),'Color','w');
            subplot(1,2,1); plotResiduals(lme, 'histogram'); title(sprintf('Model %d Residuals',i));
            subplot(1,2,2); plotResiduals(lme, 'fitted');    title(sprintf('Model %d Residuals vs Fitted',i));
        catch
        end
    else
        results(i).lme   = [];
        results(i).anova = [];
        results(i).AIC   = Inf;
        results(i).BIC   = Inf;
        results(i).logLik= -Inf;
        results(i).N     = sum(~isnan(Tf.tResp_t));
        fprintf('[WARN] Model %d failed: %s\n', i, msg);
    end
end

% Comparison table & AIC/BIC bar plot
cmp = table( (1:nM)', string({results.form})', [results.AIC]', [results.BIC]', [results.logLik]', [results.N]', ...
    'VariableNames', {'idx','formula','AIC','BIC','logLik','N'});
disp('=== Model Comparison ==='); disp(cmp);

figure('Name','Model Comparison AIC/BIC','Color','w');
subplot(1,2,1);
bar(cmp.idx, cmp.AIC); title('AIC'); xlabel('Model'); ylabel('AIC'); xticks(cmp.idx); xticklabels(cmp.idx); grid on;
subplot(1,2,2);
bar(cmp.idx, cmp.BIC); title('BIC'); xlabel('Model'); ylabel('BIC'); xticks(cmp.idx); xticklabels(cmp.idx); grid on;

[~,bestIdx] = min(cmp.BIC);
disp('Best (by BIC):'); disp(cmp(bestIdx,:));

if isfinite(cmp.BIC(bestIdx)) && results(bestIdx).ok
    fprintf('Best model formula: %s\n', results(bestIdx).form);
end

%% ============================== FUNCTIONS ==============================
function Tf = drop_sparse_levels(Tf, minCount)
    cats = {'session','p1_prev','p2_prev','condition','monkey'};
    for i=1:numel(cats)
        v = cats{i};
        if ismember(v, Tf.Properties.VariableNames) && iscategorical(Tf.(v))
            Tf.(v) = removecats(Tf.(v));
            C = categories(Tf.(v));
            toDrop = false(size(C));
            for j=1:numel(C)
                toDrop(j) = sum(Tf.(v)==C{j}) < minCount;
            end
            if any(toDrop)
                keepMask = true(height(Tf),1);
                for j=1:numel(C)
                    if toDrop(j)
                        keepMask = keepMask & ~(Tf.(v)==C{j});
                    end
                end
                Tf = Tf(keepMask,:);
                Tf.(v) = removecats(Tf.(v));
            end
        end
    end
end

function [lme, ok, msg] = fitModelSafe(Tf, form)
    ok = true; msg = ''; lme = [];
    try
        lme = fitlme(Tf, form, 'FitMethod','REML');
    catch
        try
            vars = {'session','p1_prev','p2_prev','condition','monkey'};
            for v = vars
                vn = v{1};
                if ismember(vn, Tf.Properties.VariableNames) && iscategorical(Tf.(vn))
                    Tf.(vn) = removecats(Tf.(vn));
                end
            end
            lme = fitlme(Tf, form, 'FitMethod','REML');
        catch ME2
            ok = false;
            msg = ME2.message;
            % Escape percent signs for safe warning printf
            msgSafe  = strrep(msg,  '%', '%%');
            formSafe = strrep(form, '%', '%%');
            warning('Model failed: %s | Form: %s', msgSafe, formSafe);
            lme = [];
        end
    end
end

function plotFixedEffects(lme, tag)
    tbl = lme.Coefficients;  % table
    names = string(tbl.Name);
    b = tbl.Estimate;
    se = tbl.SE;
    ciLo = b - 1.96*se;
    ciHi = b + 1.96*se;
    figure('Name',['Fixed Effects: ' char(tag)],'Color','w');
    hold on;
    y = numel(b):-1:1;
    for i=1:numel(b)
        plot([ciLo(i) ciHi(i)], [y(i) y(i)], '-','LineWidth',2);
        plot(b(i), y(i), 'o','MarkerSize',6,'MarkerFaceColor',[0.3 0.3 0.3]);
    end
    yticklabels(names(end:-1:1));
    yticks(1:numel(b));
    xlabel('Estimate (ms)'); title(['Fixed Effects (95% CI): ' char(tag)]);
    grid on; box on; hold off;
end
