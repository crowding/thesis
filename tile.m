function tile(m,n,monitor, figs)
%tile([m],[n],[monitor], [figs])
%
%Tiles matlab figures in the monitor in rows starting at the top-right.
%
%Inputs:
%   m  #of rows of figures (default is 3)
%   n  #of columns (default is 3)
%   monitor monitor number (1 is desktop, 2 is extended, default is
%       highest available)
%   handles array of figure handles to tile (defaults to all open figures)

%Written 3/30/12 gmb (after getting frustrated by Matlab's positioning of
%figures outside the monitor by default)

% Updated 8/2012 PBM because matlab's idea of "screen coordinates" is
% more horribly broken than sensible people like GMB suspect.

% Did you know, MonitorPositions gives you coordinates that grow DOWN
% and to the right, while FIgure placement uses coordinates that grow
% UP and to the right? So, at least, we want the screen coordinates
% going the same way that we're going to place the figures.

monitorPos = get(0,'MonitorPositions');
boxsize = max(monitorPos(:,[1 2]) + monitorPos(:, [3 4]), [], 1);
monitorPos(:,2) = boxsize(2) - (monitorPos(:,2) + monitorPos(:,4));

% Furthermore, on the Mac the horizontal origin of the Figure window
% coordinate system is the left edge of screen 1 -- not the left edge
% of the bounding box that MonitorPositions reported. But the vertical
% origin is usually the bottom of the bouhnding box. Sometimes. And
% sometimes it's in some other damned place.

if strncmp(computer, 'MAC', 3)
    monitorPos(:,1) = monitorPos(:,1) - monitorPos(1,1);
end

%This situation might be totally different on different
%platforms. Caveat utilitor matlabacus.

nMonitors = size(monitorPos,1);

if ~exist('monitor','var') || isempty(monitor)
    monitor= nMonitors;
end

if monitor>nMonitors
    error(sprintf('No monitor #%d available',monitor));
end

sz = monitorPos(monitor,[3,4]);

if ~exist('figs', 'var') || isempty(figs)
    figs = sort(get(0,'Children'))';
end

if ~exist('n','var') || isempty(n)
    if ~exist('m', 'var') || isempty(m)
        n = ceil(sqrt(length(figs)));
    else
        n = ceil(length(figs) / m);
    end
end

if ~exist('m','var') || isempty(m)
    m = ceil(length(figs) / n);
end



x = 0;

if monitor==1
    x0=-2;
    y0 = sz(2)+2;
    dx = round((sz(1)+2)/n);
    dy = round((sz(2)-38)/m);
else
    x0 = monitorPos(monitor,1)-2;
    y0 = monitorPos(monitor,2)+2;
    dx = round((sz(1)+2)/n);
    dy = round(sz(2)/m);
end
x =x0;
y =y0;


set(0,'Units','pixels');

[xi, yi] = ind2sub([n m], 1+mod(1:numel(figs), m*n));
x = x0 + (xi-1)*dx;
y = y0 + (yi-1)*dy;

arrayfun(@place, figs(:), x(:), y(:))
function place(fig,x,y)
    set(fig,'units','pixels');
    set(fig,'OuterPosition',[x,y,dx,dy]);
end

end
