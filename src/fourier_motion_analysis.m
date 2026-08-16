%% Fourier "best-fit" for 2D motion of 8 body points using FFT
clear; clc; close all;

%% ----------------- 1. Load data -----------------

% Find project folder from this script's location
scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);

% Data file inside project/data
filename = fullfile(projectFolder, 'data', 'Oliver Running Data.csv');

% Preserve original column names
T = readtable(filename, "VariableNamingRule","preserve");

% Time vector
time = T.("VideoAnalysis: Time (s)");
time = time - time(1);
N_samp = numel(time);

% Estimate sampling frequency
dt = mean(diff(time));
Fs = 1/dt;

fprintf('Estimated frame rate: %.3f Hz\n', Fs);

%% ----------------- 2. FFT frequency axis -----------------

Nfft  = N_samp;
f     = (0:Nfft-1) * (Fs/Nfft);
k_max = floor(Nfft/2);
f_pos = f(1:k_max);

%% ----------------- 3. Settings -----------------

numHarm = 6;
maxFreq = Inf;

% Outlier handling
doOutlierClean      = true;
outlierWindow       = 15;
outlierThreshFactor = 3;

% Body points
points = { ...
    "Wrist", ...
    "Elbow", ...
    "Shoulder", ...
    "Ear", ...
    "Hip", ...
    "Knee", ...
    "Ankle", ...
    "Toe"};

%% ----------------- 4. Global 3D figure -----------------

figAll = figure('Name','All points 3D trajectory','Color','w');
hold on;
grid on;

xlabel('Time (s)');
ylabel('X (m)');
zlabel('Y (m)');
title('All body points – 3D trajectory (Time, X, Y)');
view(135,25);

colors = lines(numel(points));

%% ----------------- 5. Loop over each body point -----------------

for p = 1:numel(points)

    pt = points{p};

    xName = "VideoAnalysis: X " + pt + " (m)";
    yName = "VideoAnalysis: Y " + pt + " (m)";

    if ~ismember(xName, T.Properties.VariableNames) || ...
       ~ismember(yName, T.Properties.VariableNames)

        warning('Skipping %s (columns not found).', pt);
        continue;
    end

    X_raw = T.(xName);
    Y_raw = T.(yName);

    %% ----- 5a. Outlier cleaning -----

    if doOutlierClean

        idxOutX = isoutlier( ...
            X_raw, ...
            'movmedian', ...
            outlierWindow, ...
            'ThresholdFactor', ...
            outlierThreshFactor);

        idxOutY = isoutlier( ...
            Y_raw, ...
            'movmedian', ...
            outlierWindow, ...
            'ThresholdFactor', ...
            outlierThreshFactor);

        X_clean = X_raw;
        Y_clean = Y_raw;

        X_clean(idxOutX) = NaN;
        Y_clean(idxOutY) = NaN;

        X_clean = fillmissing(X_clean, 'pchip');
        Y_clean = fillmissing(Y_clean, 'pchip');

    else

        X_clean = X_raw;
        Y_clean = Y_raw;

    end

    %% ----- 5b. FFT + reconstruction for X(t) -----

    X_fft = fft(X_clean);
    X_fft_pos = X_fft(1:k_max);

    ampX   = (2/Nfft) * abs(X_fft_pos);
    phaseX = angle(X_fft_pos);

    validIdx = 2:k_max;

    if ~isinf(maxFreq)
        validIdx = validIdx(f_pos(validIdx) <= maxFreq);
    end

    [~, orderX] = sort(ampX(validIdx), 'descend');

    keepCountX = min(numHarm, numel(orderX));
    keepIdxX   = validIdx(orderX(1:keepCountX));

    X_recon = mean(X_clean) * ones(size(time));

    for k = 1:keepCountX

        n   = keepIdxX(k);
        Ak  = ampX(n);
        fk  = f_pos(n);
        phk = phaseX(n);

        X_recon = X_recon + ...
            Ak * cos(2*pi*fk*time + phk);

    end

    %% ----- 5c. FFT + reconstruction for Y(t) -----

    Y_fft = fft(Y_clean);
    Y_fft_pos = Y_fft(1:k_max);

    ampY   = (2/Nfft) * abs(Y_fft_pos);
    phaseY = angle(Y_fft_pos);

    validIdx = 2:k_max;

    if ~isinf(maxFreq)
        validIdx = validIdx(f_pos(validIdx) <= maxFreq);
    end

    [~, orderY] = sort(ampY(validIdx), 'descend');

    keepCountY = min(numHarm, numel(orderY));
    keepIdxY   = validIdx(orderY(1:keepCountY));

    Y_recon = mean(Y_clean) * ones(size(time));

    for k = 1:keepCountY

        n   = keepIdxY(k);
        Ak  = ampY(n);
        fk  = f_pos(n);
        phk = phaseY(n);

        Y_recon = Y_recon + ...
            Ak * cos(2*pi*fk*time + phk);

    end

    %% ----- 5d. Per-point figure -----

    figure('Name', pt, 'Color','w');

    sgtitle( ...
        sprintf('%s – Fourier reconstruction (%d harmonics)', ...
        pt, numHarm), ...
        'Interpreter','none');

    % X(t)
    subplot(3,1,1);

    plot(time, X_raw, 'b.-', ...
        'DisplayName','X raw');

    hold on;

    plot(time, X_recon, 'r', ...
        'LineWidth',1.5, ...
        'DisplayName','X Fourier fit');

    ylabel('X (m)');
    title(pt + " : X(t)");
    legend('Location','best');
    grid on;

    % Y(t)
    subplot(3,1,2);

    plot(time, Y_raw, 'b.-', ...
        'DisplayName','Y raw');

    hold on;

    plot(time, Y_recon, 'r', ...
        'LineWidth',1.5, ...
        'DisplayName','Y Fourier fit');

    ylabel('Y (m)');
    title(pt + " : Y(t)");
    legend('Location','best');
    grid on;

    % 3D trajectory
    subplot(3,1,3);

    plot3( ...
        time, X_raw, Y_raw, ...
        'c.-', ...
        'DisplayName','Original');

    hold on;

    plot3( ...
        time, X_recon, Y_recon, ...
        'm', ...
        'LineWidth',1.5, ...
        'DisplayName','Fourier fit');

    xlabel('Time (s)');
    ylabel('X (m)');
    zlabel('Y (m)');
    title(pt + " : 3D trajectory (Time, X, Y)");
    legend('Location','best');
    grid on;
    view(135,25);

    %% ----- 5e. Add to global figure -----

    figure(figAll);

    plot3( ...
        time, X_recon, Y_recon, ...
        'LineWidth',1.5, ...
        'Color',colors(p,:), ...
        'DisplayName',char(pt));

end

figure(figAll);
legend('Location','best');
grid on;