%% Continuous hand-drawn shape -> Fourier epicycles
% Toolbox-free version using mouse callbacks in base MATLAB

clear; clc; close all;

%% ----------------- 1. Draw shape with mouse -----------------

fig = figure( ...
    'Color','k', ...
    'Name','Draw a Shape');

ax = axes(fig);

set(ax, ...
    'Color','k', ...
    'XColor','w', ...
    'YColor','w');

axis(ax,[-1 1 -1 1]);
axis(ax,'equal');
axis(ax,'manual');

grid(ax,'on');
hold(ax,'on');

xlabel(ax,'x','Color','w');
ylabel(ax,'y','Color','w');

title(ax, ...
    {'Click and drag to draw a shape', ...
     'Release the mouse to close it automatically'}, ...
    'Color','w');

% Line used to display drawing
hDraw = plot(ax,NaN,NaN, ...
    'LineWidth',1.5);

% Store drawing information inside figure
fig.UserData.isDrawing = false;
fig.UserData.x = [];
fig.UserData.y = [];
fig.UserData.hDraw = hDraw;
fig.UserData.ax = ax;

% Mouse callbacks
fig.WindowButtonDownFcn = @startDrawing;
fig.WindowButtonMotionFcn = @continueDrawing;
fig.WindowButtonUpFcn = @stopDrawing;

% Wait until drawing is finished
uiwait(fig);

% Retrieve drawing
xClicks = fig.UserData.x;
yClicks = fig.UserData.y;

if numel(xClicks) < 3
    error('Not enough points were drawn to define a shape.');
end

%% ----------------- 2. Resample curve uniformly -----------------

N = 400;

dx = diff(xClicks);
dy = diff(yClicks);

ds = hypot(dx,dy);
ds = ds(:);

% Cumulative distance along drawn path
s = [0; cumsum(ds)];

% Remove duplicate path positions
% These can occur if the mouse reports the same location twice
[s, uniqueIdx] = unique(s,'stable');

xClicks = xClicks(uniqueIdx);
yClicks = yClicks(uniqueIdx);

Ltotal = s(end);

if Ltotal == 0
    error('The drawn shape has zero total length.');
end

% Uniform spacing along shape
s_uniform = linspace(0,Ltotal,N).';

x = interp1( ...
    s, ...
    xClicks, ...
    s_uniform, ...
    'linear');

y = interp1( ...
    s, ...
    yClicks, ...
    s_uniform, ...
    'linear');

%% ----------------- 3. Convert to complex signal -----------------

z = x + 1i*y;

% Recenter around origin
z = z - mean(z);

%% ----------------- 4. Fourier coefficients -----------------

Z = fft(z)/N;

% Number of positive and negative harmonics
M = 12;

% Frequency order:
% 0, +1, -1, +2, -2, ...
n_list = 0;

for m = 1:M
    n_list = [n_list, m, -m]; %#ok<AGROW>
end

modes = struct('n',{},'C',{});

for j = 1:numel(n_list)

    n = n_list(j);

    % Convert signed Fourier frequency to MATLAB FFT index
    if n >= 0
        k_idx = n + 1;
    else
        k_idx = N + n + 1;
    end

    modes(j).n = n;
    modes(j).C = Z(k_idx);

end

n_vec = [modes.n];
C_vec = [modes.C];

K = numel(modes);

%% ----------------- 5. Reconstruct Fourier curve -----------------

t_curve = linspace(0,1,800);

z_curve = zeros(size(t_curve));

for j = 1:K

    n = n_vec(j);
    C = C_vec(j);

    z_curve = z_curve + ...
        C .* exp(1i*2*pi*n*t_curve);

end

%% ----------------- 6. Create animation window -----------------

fig2 = figure( ...
    'Color','k', ...
    'Name','Fourier Epicycles');

%% Main shape plot

ax1 = subplot(1,2,1);

set(ax1, ...
    'Color','k', ...
    'XColor','w', ...
    'YColor','w');

axis(ax1,'equal');
grid(ax1,'on');
hold(ax1,'on');

xlabel(ax1,'x','Color','w');
ylabel(ax1,'y','Color','w');

title(ax1, ...
    'Original Shape and Fourier Reconstruction', ...
    'Color','w');

% Original resampled shape
plot(ax1, ...
    real(z), ...
    imag(z), ...
    'LineWidth',1);

% Fourier reconstruction
plot(ax1, ...
    real(z_curve), ...
    imag(z_curve), ...
    'LineWidth',1.5);

%% Epicycle segments

segLines = gobjects(1,K);

for j = 1:K

    segLines(j) = plot(ax1, ...
        [0 0], ...
        [0 0], ...
        'w-');

end

endPoint = plot(ax1, ...
    0, ...
    0, ...
    'o', ...
    'MarkerFaceColor','w');

%% Viewing limits

minx = min(real(z));
maxx = max(real(z));

miny = min(imag(z));
maxy = max(imag(z));

rangeX = maxx-minx;
rangeY = maxy-miny;

marginX = max(0.2,0.2*rangeX);
marginY = max(0.2,0.2*rangeY);

axis(ax1,[ ...
    minx-marginX ...
    maxx+marginX ...
    miny-marginY ...
    maxy+marginY]);

%% ----------------- 7. X coordinate plot -----------------

ax2 = subplot(2,2,2);

set(ax2, ...
    'Color','k', ...
    'XColor','w', ...
    'YColor','w');

hold(ax2,'on');
grid(ax2,'on');

xlabel(ax2,'t','Color','w');
ylabel(ax2,'x(t)','Color','w');

title(ax2, ...
    'Endpoint x(t)', ...
    'Color','w');

xTimeLine = plot(ax2,NaN,NaN);

%% ----------------- 8. Y coordinate plot -----------------

ax3 = subplot(2,2,4);

set(ax3, ...
    'Color','k', ...
    'XColor','w', ...
    'YColor','w');

hold(ax3,'on');
grid(ax3,'on');

xlabel(ax3,'t','Color','w');
ylabel(ax3,'y(t)','Color','w');

title(ax3, ...
    'Endpoint y(t)', ...
    'Color','w');

yTimeLine = plot(ax3,NaN,NaN);

%% ----------------- 9. Animate epicycles -----------------

Fs = 60;

t_anim = t_curve;

x_hist = zeros(size(t_anim));
y_hist = zeros(size(t_anim));

for ii = 1:length(t_anim)

    if ~isgraphics(fig2)
        break;
    end

    tau = t_anim(ii);

    xCurr = 0;
    yCurr = 0;

    xs = zeros(1,K+1);
    ys = zeros(1,K+1);

    xs(1) = xCurr;
    ys(1) = yCurr;

    for j = 1:K

        n = n_vec(j);
        C = C_vec(j);

        z_j = C * ...
            exp(1i*2*pi*n*tau);

        dx = real(z_j);
        dy = imag(z_j);

        xCurr = xCurr + dx;
        yCurr = yCurr + dy;

        xs(j+1) = xCurr;
        ys(j+1) = yCurr;

        set(segLines(j), ...
            'XData',[xs(j) xs(j+1)], ...
            'YData',[ys(j) ys(j+1)]);

    end

    x_hist(ii) = xCurr;
    y_hist(ii) = yCurr;

    set(endPoint, ...
        'XData',xCurr, ...
        'YData',yCurr);

    set(xTimeLine, ...
        'XData',t_anim(1:ii), ...
        'YData',x_hist(1:ii));

    set(yTimeLine, ...
        'XData',t_anim(1:ii), ...
        'YData',y_hist(1:ii));

    drawnow;

    pause(1/Fs);

end


%% ----------------- Mouse callback functions -----------------

function startDrawing(src,~)

    data = src.UserData;

    data.isDrawing = true;

    point = data.ax.CurrentPoint;

    x = point(1,1);
    y = point(1,2);

    data.x = x;
    data.y = y;

    set(data.hDraw, ...
        'XData',data.x, ...
        'YData',data.y);

    src.UserData = data;

end


function continueDrawing(src,~)

    data = src.UserData;

    if ~data.isDrawing
        return;
    end

    point = data.ax.CurrentPoint;

    x = point(1,1);
    y = point(1,2);

    % Only record points inside drawing window
    if x >= -1 && x <= 1 && ...
       y >= -1 && y <= 1

        data.x(end+1) = x;
        data.y(end+1) = y;

        set(data.hDraw, ...
            'XData',data.x, ...
            'YData',data.y);

    end

    src.UserData = data;

end


function stopDrawing(src,~)

    data = src.UserData;

    data.isDrawing = false;

    % Snap final point back to starting point
    if numel(data.x) >= 2

        data.x(end+1) = data.x(1);
        data.y(end+1) = data.y(1);

        set(data.hDraw, ...
            'XData',data.x, ...
            'YData',data.y);

    end

    src.UserData = data;

    % Show closed shape briefly
    drawnow;
    pause(0.3);

    % Allow main script to continue
    uiresume(src);

end