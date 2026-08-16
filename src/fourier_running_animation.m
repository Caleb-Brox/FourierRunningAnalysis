%% Fourier running model + mirrored left side + continuous animation
% Reconstructs right-side joint motion using Fourier components,
% estimates stride period, generates the left side using a half-stride
% phase shift, and continuously animates the reconstructed runner.

clear; clc; close all;

%% ----------------- 1. Load data -----------------

% Find project folder based on this script's location
scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);

% Locate running data
filename = fullfile( ...
    projectFolder, ...
    'data', ...
    'Oliver Running Data.csv');

T = readtable(filename, ...
    "VariableNamingRule","preserve");

% Experimental timestamps
time = T.("VideoAnalysis: Time (s)");
time = time - time(1);

N_samp = numel(time);

% Estimate average sampling interval and frame rate
dt = mean(diff(time));
Fs = 1/dt;

fprintf("Estimated frame rate: %.3f Hz\n", Fs);

%% ----------------- 2. FFT frequency axis -----------------

Nfft = N_samp;

f = (0:Nfft-1) * (Fs/Nfft);

k_max = floor(Nfft/2);

f_pos = f(1:k_max);

%% ----------------- 3. Fourier settings -----------------

% Number of strongest Fourier components retained
numHarm = 2;

% Optional maximum allowed frequency
maxFreq = Inf;

% Right-side body points
pointsR = [ ...
    "Wrist", ...
    "Elbow", ...
    "Shoulder", ...
    "Ear", ...
    "Hip", ...
    "Knee", ...
    "Ankle", ...
    "Toe"];

nR = numel(pointsR);

% Reconstructed joint positions
XR = zeros(N_samp,nR);
YR = zeros(N_samp,nR);

% Store Fourier information
freqX  = cell(1,nR);
ampX   = cell(1,nR);
phaseX = cell(1,nR);
meanX  = zeros(1,nR);

freqY  = cell(1,nR);
ampY   = cell(1,nR);
phaseY = cell(1,nR);
meanY  = zeros(1,nR);

%% ----------------- 4. Fourier reconstruction -----------------

for p = 1:nR

    pt = pointsR(p);

    xName = "VideoAnalysis: X " + pt + " (m)";
    yName = "VideoAnalysis: Y " + pt + " (m)";

    if ~ismember(xName,T.Properties.VariableNames) || ...
       ~ismember(yName,T.Properties.VariableNames)

        warning("Missing data for %s",pt);
        continue;

    end

    % X coordinate
    [XR(:,p), ...
     freqX{p}, ...
     ampX{p}, ...
     phaseX{p}, ...
     meanX(p)] = recon_fft( ...
        T.(xName), ...
        time, ...
        Nfft, ...
        k_max, ...
        f_pos, ...
        numHarm, ...
        maxFreq);

    % Y coordinate
    [YR(:,p), ...
     freqY{p}, ...
     ampY{p}, ...
     phaseY{p}, ...
     meanY(p)] = recon_fft( ...
        T.(yName), ...
        time, ...
        Nfft, ...
        k_max, ...
        f_pos, ...
        numHarm, ...
        maxFreq);

end

%% ----------------- 5. Print Fourier equations -----------------

fprintf('\n');
fprintf('================ FOURIER EQUATIONS =================\n');

for p = 1:nR

    name = char(pointsR(p));

    % X equation
    fprintf('\nX_%s(t) = %.5g', ...
        name, ...
        meanX(p));

    for k = 1:numel(freqX{p})

        fprintf( ...
            ' + %.5g*cos(2*pi*%.5g*t + %.5g)', ...
            ampX{p}(k), ...
            freqX{p}(k), ...
            phaseX{p}(k));

    end

    fprintf('\n');

    % Y equation
    fprintf('Y_%s(t) = %.5g', ...
        name, ...
        meanY(p));

    for k = 1:numel(freqY{p})

        fprintf( ...
            ' + %.5g*cos(2*pi*%.5g*t + %.5g)', ...
            ampY{p}(k), ...
            freqY{p}(k), ...
            phaseY{p}(k));

    end

    fprintf('\n');

end

fprintf('\n');
fprintf(['Units: position in meters, time in seconds, ' ...
         'frequency in Hz, phase in radians.\n\n']);

%% ----------------- 6. Estimate stride period -----------------

idxAnkleR = find(pointsR=="Ankle");

YA = YR(:,idxAnkleR);

YA_fft = fft(YA);
YA_fft_pos = YA_fft(1:k_max);

ampYA = (2/Nfft) * abs(YA_fft_pos);

% Ignore DC component
[~,idxMax] = max(ampYA(2:end));

idxMax = idxMax + 1;

f0 = f_pos(idxMax);

Tstride = 1/f0;

fprintf( ...
    "Estimated stride frequency = %.3f Hz\n", ...
    f0);

fprintf( ...
    "Estimated stride period = %.3f s\n", ...
    Tstride);

%% ----------------- 7. Generate left side -----------------

% Opposite limbs should be roughly half a stride out of phase
halfStrideTime = Tstride/2;

halfShiftSamples = round(halfStrideTime/dt);

fprintf( ...
    "Half-period shift = %d samples\n", ...
    halfShiftSamples);

symR = [ ...
    "Wrist", ...
    "Elbow", ...
    "Knee", ...
    "Ankle", ...
    "Toe"];

idxSymR = zeros(size(symR));

for k = 1:numel(symR)

    idxSymR(k) = find( ...
        pointsR==symR(k), ...
        1);

end

%% ----------------- 8. Assemble skeleton -----------------

% Joint indices
idxW   = find(pointsR=="Wrist");
idxE   = find(pointsR=="Elbow");
idxS   = find(pointsR=="Shoulder");
idxEar = find(pointsR=="Ear");
idxH   = find(pointsR=="Hip");
idxK   = find(pointsR=="Knee");
idxA   = find(pointsR=="Ankle");
idxT   = find(pointsR=="Toe");

% 13 joints total
XJ = zeros(N_samp,13);
YJ = zeros(N_samp,13);

% Right side
XJ(:,1) = XR(:,idxW);
YJ(:,1) = YR(:,idxW);

XJ(:,2) = XR(:,idxE);
YJ(:,2) = YR(:,idxE);

XJ(:,3) = XR(:,idxS);
YJ(:,3) = YR(:,idxS);

XJ(:,4) = XR(:,idxEar);
YJ(:,4) = YR(:,idxEar);

XJ(:,5) = XR(:,idxH);
YJ(:,5) = YR(:,idxH);

XJ(:,6) = XR(:,idxK);
YJ(:,6) = YR(:,idxK);

XJ(:,7) = XR(:,idxA);
YJ(:,7) = YR(:,idxA);

XJ(:,8) = XR(:,idxT);
YJ(:,8) = YR(:,idxT);

%% Generate left limbs using half-stride shift

XR_L = zeros(N_samp,numel(symR));
YR_L = zeros(N_samp,numel(symR));

for k = 1:numel(symR)

    idxR = idxSymR(k);

    XR_L(:,k) = circshift( ...
        XR(:,idxR), ...
        halfShiftSamples);

    YR_L(:,k) = circshift( ...
        YR(:,idxR), ...
        halfShiftSamples);

end

% Left arm
XJ(:,9)  = XR_L(:,symR=="Wrist");
YJ(:,9)  = YR_L(:,symR=="Wrist");

XJ(:,10) = XR_L(:,symR=="Elbow");
YJ(:,10) = YR_L(:,symR=="Elbow");

% Left leg
XJ(:,11) = XR_L(:,symR=="Knee");
YJ(:,11) = YR_L(:,symR=="Knee");

XJ(:,12) = XR_L(:,symR=="Ankle");
YJ(:,12) = YR_L(:,symR=="Ankle");

XJ(:,13) = XR_L(:,symR=="Toe");
YJ(:,13) = YR_L(:,symR=="Toe");

%% ----------------- 9. Skeleton connections -----------------

segments = [ ...

    1 2;      % Right wrist → elbow
    2 3;      % Right elbow → shoulder

    3 5;      % Shoulder → hip

    5 6;      % Right hip → knee
    6 7;      % Right knee → ankle
    7 8;      % Right ankle → toe

    4 3;      % Ear → shoulder

    9 10;     % Left wrist → elbow
    10 3;     % Left elbow → shoulder

    5 11;     % Hip → left knee
    11 12;    % Left knee → ankle
    12 13];   % Left ankle → toe

%% ----------------- 10. Create animation -----------------

fig = figure( ...
    "Name","Fourier Running Animation", ...
    "Color","w");

axis equal;
grid on;
hold on;

xlabel("X (m)");
ylabel("Y (m)");

title("Fourier-Reconstructed Running Motion");

% Set viewing limits
margin = 0.05;

xMin = min(XJ(:));
xMax = max(XJ(:));

yMin = min(YJ(:));
yMax = max(YJ(:));

xRange = xMax-xMin;
yRange = yMax-yMin;

xlim([ ...
    xMin-margin*xRange, ...
    xMax+margin*xRange]);

ylim([ ...
    yMin-margin*yRange, ...
    yMax+margin*yRange]);

%% Create joints

i0 = 1;

hJoints = plot( ...
    XJ(i0,:), ...
    YJ(i0,:), ...
    "ko", ...
    "MarkerFaceColor","k", ...
    "MarkerSize",6);

%% Create body segments

hSeg = gobjects(size(segments,1),1);

for s = 1:size(segments,1)

    idxA = segments(s,1);
    idxB = segments(s,2);

    hSeg(s) = line( ...
        XJ(i0,[idxA idxB]), ...
        YJ(i0,[idxA idxB]), ...
        "LineWidth",2);

end

%% ----------------- 11. Continuous animation -----------------

% Display every second recorded sample
frameStep = 2;

frameIndices = 1:frameStep:N_samp;

% Continue until the animation window is closed
while isgraphics(fig)

    for j = 1:numel(frameIndices)

        if ~isgraphics(fig)
            break;
        end

        i = frameIndices(j);

        % Update joint locations
        set( ...
            hJoints, ...
            "XData",XJ(i,:), ...
            "YData",YJ(i,:));

        % Update body segments
        for s = 1:size(segments,1)

            idxA = segments(s,1);
            idxB = segments(s,2);

            set( ...
                hSeg(s), ...
                "XData",XJ(i,[idxA idxB]), ...
                "YData",YJ(i,[idxA idxB]));

        end

        drawnow;

        % -----------------------------------------
        % Follow the ACTUAL experimental timestamps
        % -----------------------------------------

        if j < numel(frameIndices)

            currentIndex = frameIndices(j);
            nextIndex = frameIndices(j+1);

            waitTime = ...
                time(nextIndex) - time(currentIndex);

        else

            % At the end of the recording,
            % use the average sampling timing before
            % starting the next loop.
            waitTime = frameStep * dt;

        end

        pause(waitTime);

    end

end


%% ----------------- Local Fourier function -----------------

function [sigRecon,freqs,amps,phases,meanVal] = recon_fft( ...
    sigRaw,time,Nfft,k_max,f_pos,numHarm,maxFreq)

    %% FFT

    Sig_fft = fft(sigRaw);

    Sig_fft_pos = Sig_fft(1:k_max);

    %% Amplitude and phase

    amp = ...
        (2/Nfft) * abs(Sig_fft_pos);

    phase = ...
        angle(Sig_fft_pos);

    %% Ignore DC component

    validIdx = 2:k_max;

    if ~isinf(maxFreq)

        validIdx = validIdx( ...
            f_pos(validIdx)<=maxFreq);

    end

    %% Find strongest frequencies

    [~,order] = sort( ...
        amp(validIdx), ...
        "descend");

    keepCount = min( ...
        numHarm, ...
        numel(order));

    keepIdx = validIdx( ...
        order(1:keepCount));

    %% Save Fourier parameters

    freqs = f_pos(keepIdx);

    amps = amp(keepIdx);

    phases = phase(keepIdx);

    meanVal = mean(sigRaw);

    %% Reconstruct signal

    sigRecon = ...
        meanVal * ones(size(time));

    for kk = 1:keepCount

        Ak = amps(kk);
        fk = freqs(kk);
        ph = phases(kk);

        sigRecon = ...
            sigRecon + ...
            Ak*cos( ...
                2*pi*fk*time + ph);

    end

end