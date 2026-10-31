%% photocurrent_vs_power.m
% Photocurrent of an optically pumped sample against pump power:
% - Parses the power from each file name (mw / w / no-laser)
% - Uses 0W files as baseline (mean + std), propagates uncertainty
% - Keeps the last 70% of the points for the steady-state statistics
% - Log-log fits + regime detection with slope & R^2 per regime
%
% Outputs:
%   Figures:
%     1) Photocurrent vs power (no error bars)
%     2) Photocurrent vs power (error bars)
%     3) Log-log (all)
%     4) Log-log (P <= 4W)
%     5) Local slope vs P (regime map)
%     6) Log-log + regime fits (each regime slope & R2)
%   CSV:
%     powerSweep_summary.csv
%     regime_fits.csv

clc; clear; close all;

%% ================= USER SETTINGS =================
dataDir = uigetdir(pwd,'Select folder with CSV files');
if isequal(dataDir,0), error('No folder selected.'); end

% Steady-state definition: keep last (1-transientFrac) of points
transientFrac = 0.30;   % same as old code (remove first 30%, keep last 70%)

% Sign handling for log-log:
% If your photocurrent is negative by convention, set signFlip = true (recommended)
signFlip = false;       % if true, fits use Y = -(Iphoto)
useAbsForLoglog = false; % if true, fits use abs(Iphoto). Not as defensible as signFlip.

% Log-log fit cutoff
PmaxFit = 4.0;          % W

% Regime detection knobs
winPts          = 9;     % sliding window points for local slope (odd)
maxRegimes      = 3;     % up to 3 regimes
minPtsPerRegime = 8;     % minimum points per regime fit
slopeJumpThresh = 0.20;  % fallback split sensitivity if findchangepts unavailable

% Exports
exportSummaryCSV = true;
summaryCSVName   = 'powerSweep_summary.csv';
exportRegimeCSV  = true;
regimeCSVName    = 'regime_fits.csv';

%% ================= LOAD FILES =================
files = dir(fullfile(dataDir,'*.csv'));
if isempty(files), error('No CSV files found in selected folder.'); end

% Preallocate vectors
laserPower_W = [];
avgCurrent_A = [];
stdCurrent_A = [];
fileNameList = strings(0,1);

%% ================= LOOP THROUGH FILES =================
for i = 1:length(files)
    filename = files(i).name;
    filepath = fullfile(dataDir, filename);
    fname = lower(filename);

    % -------- PARSE LASER POWER FROM NAME (old behavior) --------
    if contains(fname,'no-laser') || contains(fname,'no laser') || contains(fname,'nolaser')
        powerW = 0;

    elseif contains(fname,'mw')
        token = regexp(fname,'([\d_\.]+)mw','tokens','once');
        if isempty(token), continue; end
        powerStr = strrep(token{1},'_','.');
        powerW = str2double(powerStr)/1000;   % mW → W

    elseif contains(fname,'w')
        token = regexp(fname,'([\d_\.]+)w','tokens','once');
        if isempty(token), continue; end
        powerStr = strrep(token{1},'_','.');
        powerW = str2double(powerStr);        % W

    else
        % if power isn't encoded in name, skip
        continue
    end

    % -------- READ CURRENT TRACE ROBUSTLY --------
    S = readKeithleyCSV_preserve(filepath);  % uses "Reading," header if present
    I = S.I;

    % Remove Keithley overflow (old behavior)
    I(I > 9e36) = NaN;

    % Remove transient: keep last 70%
    N = numel(I);
    iStart = max(1, round(transientFrac*N));
    I = I(iStart:end);

    avgI = mean(I,'omitnan');
    stdI = std(I,'omitnan');

    % -------- STORE --------
    laserPower_W(end+1,1) = powerW; %#ok<SAGROW>
    avgCurrent_A(end+1,1) = avgI;   %#ok<SAGROW>
    stdCurrent_A(end+1,1) = stdI;   %#ok<SAGROW>
    fileNameList(end+1,1) = string(filename); %#ok<SAGROW>
end

if isempty(laserPower_W)
    error('No usable files after parsing power from filenames.');
end

%% ================= SORT BY LASER POWER =================
[laserPower_W, idx] = sort(laserPower_W);
avgCurrent_A = avgCurrent_A(idx);
stdCurrent_A = stdCurrent_A(idx);
fileNameList = fileNameList(idx);

%% ================= BASELINE (0W) SUBTRACTION (IMPROVED) =================
darkIdx = (laserPower_W == 0);
if ~any(darkIdx)
    warning('No 0W / no-laser files found. Baseline subtraction will be skipped (Iphoto = I).');
    baseMean = 0;
    baseStd  = 0;
else
    % baseline = mean of 0W means (like old), but also compute baseline std
    baseMean = mean(avgCurrent_A(darkIdx), 'omitnan');

    % For baseline uncertainty: combine measurement stds + scatter between baseline files
    % (This is more realistic than using just one number.)
    scatterStd = std(avgCurrent_A(darkIdx), 'omitnan');
    measStd    = mean(stdCurrent_A(darkIdx), 'omitnan'); % typical within-file noise
    baseStd    = sqrt(scatterStd^2 + measStd^2);
end

photoCurrent_A = avgCurrent_A - baseMean;
photoStd_A     = sqrt(stdCurrent_A.^2 + baseStd.^2);  % propagate baseline uncertainty

% convert to microamp for linear plots
photoCurrent_uA = photoCurrent_A * 1e6;
photoStd_uA     = photoStd_A     * 1e6;

%% ================= 1) LINEAR PLOT (NO ERROR BARS) =================
figure('Color','w','Name','Photocurrent vs Power (No Error Bars)');
plot(laserPower_W, photoCurrent_uA, 'o-', 'LineWidth',1.5,'MarkerSize',8);
xlabel('Laser Power (W)');
ylabel('Photocurrent (\muA)');
title('Photocurrent vs Laser Power at 300 V (baseline-subtracted)');
grid on;

%% ================= 2) LINEAR PLOT (WITH ERROR BARS) =================
figure('Color','w','Name','Photocurrent vs Power (Error Bars)');
errorbar(laserPower_W, photoCurrent_uA, photoStd_uA, 'o-', 'LineWidth',1.5,'MarkerSize',8);
xlabel('Laser Power (W)');
ylabel('Photocurrent (\muA)');
title('Photocurrent vs Laser Power at 300 V (baseline-subtracted + propagated std)');
grid on;

%% ================= PREP LOG-LOG VECTORS =================
P = laserPower_W(:);
Iph = photoCurrent_A(:);

% Apply sign handling for fitting
YfitBase = Iph;
yLabel = 'I_{photo} (A)';

if signFlip
    YfitBase = -YfitBase;
    yLabel = '-I_{photo} (A)';
end
if useAbsForLoglog
    YfitBase = abs(YfitBase);
    yLabel = '|I_{photo}| (A)';
end

% Positive-only for log space
goodAll = (P > 0) & (YfitBase > 0) & isfinite(P) & isfinite(YfitBase);
xAll = P(goodAll);
yAll = YfitBase(goodAll);

%% ================= 3) LOG-LOG FIT (ALL POSITIVE POINTS) =================
figure('Color','w','Name','Log-Log (All)');
loglog(xAll, yAll, 'o','MarkerSize',8,'LineWidth',1.5); hold on; grid on;

[alphaAll, Aall, R2all] = loglogFitWithR2(xAll, yAll);
fitYall = Aall * xAll.^alphaAll;
loglog(xAll, fitYall, '--k','LineWidth',1.5);

xlabel('Laser Power (W)');
ylabel(yLabel);
title(sprintf('Log–Log (ALL positive points): \\alpha = %.2f, R^2 = %.3f', alphaAll, R2all));
legend('Data','Power-law fit','Location','best');

fprintf('\nALL fit: alpha = %.6f, A = %.6g, R^2 = %.6f, N=%d\n', alphaAll, Aall, R2all, numel(xAll));

%% ================= 4) LOG-LOG FIT (P <= 4 W) =================
good4 = goodAll & (P <= PmaxFit);
x4 = P(good4);
y4 = YfitBase(good4);

figure('Color','w','Name','Log-Log (P <= 4W)');
loglog(x4, y4, 'o','MarkerSize',8,'LineWidth',1.5); hold on; grid on;

[alpha4, A4, R24] = loglogFitWithR2(x4, y4);
fitY4 = A4 * x4.^alpha4;
loglog(x4, fitY4, '--k','LineWidth',1.5);

xlabel('Laser Power (W)');
ylabel(yLabel);
title(sprintf('Log–Log (P \\le %.1f W): \\alpha = %.2f, R^2 = %.3f', PmaxFit, alpha4, R24));
legend('Data','Power-law fit','Location','best');

fprintf('P<=%.1f fit: alpha = %.6f, A = %.6g, R^2 = %.6f, N=%d\n', PmaxFit, alpha4, A4, R24, numel(x4));

%% ================= 5) LOCAL SLOPE vs POWER (REGIME MAP) =================
% For regime detection we use (xAll, yAll) sorted by power already.
if numel(xAll) >= max(10, winPts)
    if mod(winPts,2)==0, winPts = winPts+1; end
    halfW = floor(winPts/2);

    lx = log10(xAll);
    ly = log10(yAll);

    slopeLocal = nan(size(lx));
    r2Local    = nan(size(lx));

    for i = 1:numel(lx)
        i1 = max(1, i-halfW);
        i2 = min(numel(lx), i+halfW);
        if (i2-i1+1) < 5, continue; end

        xi = lx(i1:i2);
        yi = ly(i1:i2);

        c = polyfit(xi, yi, 1);
        yhat = polyval(c, xi);

        SSres = sum((yi-yhat).^2);
        SStot = sum((yi-mean(yi)).^2);
        R2 = 1 - SSres/SStot;

        slopeLocal(i) = c(1);
        r2Local(i) = R2;
    end

    figure('Color','w','Name','Local slope vs Power');
    semilogx(xAll, slopeLocal, 'o-','LineWidth',1.2,'MarkerSize',6);
    grid on;
    xlabel('Laser Power (W)');
    ylabel('Local slope dlog(I)/dlog(P)');
    title(sprintf('Local slope vs P (winPts=%d). Flat = single regime.', winPts));

    %% ================= 6) REGIME BOUNDARIES + REGIME FITS =================
    idxValid = find(isfinite(slopeLocal));
    Pm = xAll(idxValid);
    sm = slopeLocal(idxValid);
    sm_s = movmedian(sm, max(3, round(winPts/3)));

    cp = [];
    try
        cp = findchangepts(sm_s, 'MaxNumChanges', maxRegimes-1, 'Statistic', 'rms');
    catch
        cp = [];
    end

    bounds = [1; cp(:)+1; numel(Pm)+1];

    % Enforce min points per regime by merging tiny regimes
    boundsClean = bounds(1);
    for b = 2:numel(bounds)
        if (bounds(b) - boundsClean(end)) >= minPtsPerRegime
            boundsClean(end+1,1) = bounds(b); %#ok<SAGROW>
        end
    end
    if boundsClean(end) ~= numel(Pm)+1
        boundsClean(end+1,1) = numel(Pm)+1;
    end
    bounds = boundsClean;

    % Fallback split if only 1 regime survived
    if numel(bounds) <= 2
        ds = abs(diff(sm_s));
        j = find(ds > slopeJumpThresh, 1, 'first');
        if ~isempty(j) && j >= minPtsPerRegime && (numel(Pm)-j) >= minPtsPerRegime
            bounds = [1; j+1; numel(Pm)+1];
        else
            bounds = [1; numel(Pm)+1];
        end
    end

    % Draw boundaries on slope plot
    figure(findobj('Name','Local slope vs Power')); hold on;
    for k = 2:(numel(bounds)-1)
        xline(Pm(bounds(k)), '--', 'LineWidth', 1.2);
    end
    hold off;

    % Log-log + regime fits
    figure('Color','w','Name','Log-Log + Regime Fits');
    loglog(xAll, yAll, 'o', 'MarkerSize', 6); grid on; hold on;

    regID = []; PminR=[]; PmaxR=[]; NptsR=[]; slopeR=[]; intR=[]; R2R=[];

    fprintf('\n================ REGIME FIT RESULTS ================\n');
    fprintf('Regime |   Pmin(W)   Pmax(W) |    slope       intercept      R^2      N\n');
    fprintf('----------------------------------------------------------------------\n');

    for r = 1:(numel(bounds)-1)
        a = bounds(r);
        b = bounds(r+1)-1;

        Pr = Pm(a:b);
        if numel(Pr) < minPtsPerRegime, continue; end

        xr = log10(Pr);

        % get corresponding y-values by matching Pr into xAll
        % since xAll is sorted and unique-ish, use ismembertol
        tf = ismembertol(xAll, Pr, 1e-12);
        Yr = yAll(tf);
        % but tf might include extra points if duplicates; safer:
        Yr = interp1(xAll, yAll, Pr, 'linear', 'extrap');
        yr = log10(Yr);

        c = polyfit(xr, yr, 1);
        m = c(1); b0 = c(2);

        yhat = polyval(c, xr);
        SSres = sum((yr - yhat).^2);
        SStot = sum((yr - mean(yr)).^2);
        R2 = 1 - SSres/SStot;

        regID(end+1,1)  = r; %#ok<SAGROW>
        PminR(end+1,1)  = min(Pr);
        PmaxR(end+1,1)  = max(Pr);
        NptsR(end+1,1)  = numel(Pr);
        slopeR(end+1,1) = m;
        intR(end+1,1)   = b0;
        R2R(end+1,1)    = R2;

        Pline = logspace(log10(min(Pr)), log10(max(Pr)), 200);
        Yline = 10^(b0) .* (Pline.^m);
        loglog(Pline, Yline, '-', 'LineWidth', 2);

        Pmid = 10^(mean(log10([min(Pr), max(Pr)])));
        Ymid = 10^(b0) * (Pmid^m);
        text(Pmid, Ymid, sprintf('Reg %d: m=%.2f, R^2=%.2f', r, m, R2), ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','left');

        fprintf('%6d | %8.3g  %8.3g | %10.4f  %12.4f  %8.4f  %5d\n', ...
            r, min(Pr), max(Pr), m, b0, R2, numel(Pr));
    end

    xlabel('Laser Power (W)');
    ylabel(yLabel);
    title('Log–log with regime fits');
    legend('Data','Regime fits','Location','best');
    hold off;

    % Export regime table
    RegimeTable = table(regID, PminR, PmaxR, NptsR, slopeR, intR, R2R, ...
        'VariableNames', {'Regime','Pmin_W','Pmax_W','N','Slope','Intercept','R2'});
    disp(' ');
    disp(RegimeTable);

    if exportRegimeCSV
        outRegPath = fullfile(dataDir, regimeCSVName);
        writetable(RegimeTable, outRegPath);
        fprintf('Saved regime fit table: %s\n', outRegPath);
    end
else
    warning('Not enough positive points for regime detection (need at least ~winPts).');
end

%% ================= EXPORT SUMMARY TABLE =================
if exportSummaryCSV
    Tout = table(fileNameList, P, avgCurrent_A, stdCurrent_A, ...
                 repmat(baseMean,size(P)), repmat(baseStd,size(P)), ...
                 photoCurrent_A, photoStd_A, ...
        'VariableNames', {'File','Power_W','I_mean_A','I_std_A','Baseline_mean_A','Baseline_std_A','Iphoto_mean_A','Iphoto_std_A'});

    outPath = fullfile(dataDir, summaryCSVName);
    writetable(Tout, outPath);
    fprintf('\nSaved summary: %s\n', outPath);
end

%% ===================== LOCAL FUNCTIONS =====================

function S = readKeithleyCSV_preserve(filePath)
% Robust reader for Keithley-style CSV:
% If file contains a header line that starts with "Reading,", use it.
% Otherwise fall back to normal readtable.

    L = readlines(filePath);
    hdrIdx = find(startsWith(L, "Reading,"), 1, 'first');

    if isempty(hdrIdx)
        % fallback (simple readtable)
        T = readtable(filePath,'VariableNamingRule','preserve');
    else
        opts = detectImportOptions(filePath, 'NumHeaderLines', hdrIdx-1);
        opts.VariableNamesLine = hdrIdx;
        opts.DataLines = [hdrIdx+1, Inf];
        opts.VariableNamingRule = 'preserve';
        T = readtable(filePath, opts);
    end

    % Current column: "Reading"
    if any(strcmp(T.Properties.VariableNames, "Reading"))
        I = T.("Reading");
    else
        % sometimes variable naming differs; try to find a column containing "Reading"
        vn = T.Properties.VariableNames;
        k = find(contains(vn, "Reading", 'IgnoreCase', true), 1, 'first');
        if isempty(k)
            error('Cannot find "Reading" column in %s', filePath);
        end
        I = T.(vn{k});
    end

    S.I = double(I);
    S.T = T;
end

function [alpha, A, R2] = loglogFitWithR2(x, y)
% Fit y = A*x^alpha using linear fit in log10 space.
% Returns alpha, A, and R^2 in log-space.

    x = x(:); y = y(:);
    good = isfinite(x) & isfinite(y) & (x>0) & (y>0);
    x = x(good); y = y(good);

    lx = log10(x);
    ly = log10(y);

    c = polyfit(lx, ly, 1);
    alpha = c(1);
    intercept = c(2);
    A = 10^(intercept);

    yhat = polyval(c, lx);
    SSres = sum((ly - yhat).^2);
    SStot = sum((ly - mean(ly)).^2);
    R2 = 1 - SSres/SStot;
end
