function theAxis = tightsubplot(nrows, ncols, thisPlot, varargin)
%TIGHTSUBPLOT Create axes in tiled positions.
%Copy of SUBPLOT with no internal insets.

narg = nargin;
kill_siblings = 0;
create_axis = true;
move_axis = false;
delay_destroy = false;
useAutoLayout = true;
tol = sqrt(eps);
parent = get(0,'CurrentFigure');
parentfigure = parent;
if ~isempty(parent) && ~isempty(get(parent,'CurrentAxes'))
  parent = get(get(parent,'CurrentAxes'),'Parent');
  parentfigure = parent;
  if ~strcmp(get(parentfigure,'Type'),'figure')
    parentfigure = ancestor(parent,'figure');
  end
end
pvpairs = {};
preventMove = false;
% This is the percent offset from the subplot grid of the plotbox.
%inset = [.2 .18 .04 .1]; % [left bottom right top]
inset = [0 0 0 0]; % [left bottom right top]


if narg == 0 % make compatible with 3.5, i.e. subplot == subplot(111)
  nrows = 111;
  narg = 1;
end

%check for encoded format
handle = '';
position = '';
explicitParent = false;

if narg == 1
  % The argument could be one of 3 things:
  % 1) a 3-digit number 100 < num < 1000, of the format mnp
  % 2) a 3-character string containing a number as above
  % 3) an axis handle
  code = nrows;

  % turn string into a number:
  if(ischar(code)) 
    code = str2double(code);
  end

  % number with a fractional part can only be an identifier:
  if(rem(code,1) > 0)
    handle = code;
    if ~strcmp(get(handle,'type'),'axes')
      error(id('InvalidAxesHandle'),'Requires valid axes handle for input.')
    end
    create_axis = false;
    % all other numbers will be converted to mnp format:
  else
    thisPlot = rem(code, 10);
    ncols = rem( fix(code-thisPlot)/10,10);
    nrows = fix(code/100);
    if nrows*ncols < thisPlot
      error(id('SubplotIndexTooLarge'),'Index exceeds number of subplots.');
    end
    kill_siblings = 1;
    if(code == 111)
      create_axis   = false;
      delay_destroy = true;
    else
      create_axis   = true;
      delay_destroy = false;
    end
  end
    
elseif narg == 2
  % The arguments MUST be the string 'position' and a 4-element vector:
  if(strcmpi(nrows, 'position'))
    pos_size = size(ncols);
    if(pos_size(1) * pos_size(2) == 4)
      position = ncols;
    else
      error(id('InvalidPositionParameter'),...
            'Position must be of the form [left bottom width height].')
    end
  else
    error(id('UnknownOption'),'Unknown command option.')
  end
  kill_siblings = 1; % Kill overlaps here also.
  useAutoLayout = false;
  
elseif narg == 3
  % passed in subplot(m,n,p) -- we should kill overlaps
  % here too:
  kill_siblings = 1;
  
elseif narg >= 4
  arg = varargin{1};
  if ~ischar(arg)
    % passed in subplot(m,n,p,H,...)
    handle = arg;
    parent = get(handle,'Parent');
    parentfigure = ancestor(handle,'figure');
    % If the parent is passed in explicitly, don't create a new figure
    % when the "NextPlot" property is set to "new" in the figure.
    explicitParent = true;
    set(parentfigure,'CurrentAxes',handle);
    move_axis = true;
    create_axis = false;
    if narg >= 5 && strcmpi(varargin{2},'PreventMove')
        preventMove = true;
        pvpairs = varargin(3:end);
    else
        pvpairs = varargin(2:end);
    end
  elseif strncmpi(arg,'replace',1)
    % passed in subplot(m,n,p,'replace')
    kill_siblings = 2; % kill nomatter what
  elseif strcmpi(arg,'align') || strcmpi(arg,'v6')
    % passed in subplot(m,n,p,'align') or 'v6'
    % since obeying position will remove the axes from the grid just set
    % useAutoLayout to false to skip adding it to the grid to start with
    useAutoLayout = false; 
    kill_siblings = 1; % kill if it overlaps stuff
  else
    % passed in prop-value pairs
    kill_siblings = 1;
    pvpairs = varargin;
    par = find(strncmpi('Parent',pvpairs(1:2:end),6));
    if any(par)
        % If the parent is passed in explicitly, don't create a new figure
        % when the "NextPlot" property is set to "new" in the figure.
        explicitParent = true;
        parent = varargin{2*par(1)};
        parentfigure = ancestor(parent,'figure');
    end
  end
end

% if we recovered an identifier earlier, use it:
if ~isempty(handle) && ~move_axis
  parent = get(handle,'Parent');
  parentfigure = ancestor(handle,'figure');
  set(parentfigure,'CurrentAxes',handle);
else  % if we haven't recovered position yet, generate it from mnp info:
  if isempty(parent), 
    parent = gcf; 
    parentfigure = parent; 
  end
  if isempty(position)
    if min(thisPlot) < 1
      error(id('SubplotIndexTooSmall'),'Illegal plot number.')
    elseif max(thisPlot) > ncols*nrows
      error(id('SubplotIndexTooLarge'),'Index exceeds number of subplots.');
    else

      row = (nrows-1) -fix((thisPlot-1)/ncols);
      col = rem (thisPlot-1, ncols);

      % get default axis position in normalized units
      % If we have checked this quanitity once, cache it.
      if ~isappdata(parentfigure,'SubplotDefaultAxesLocation')
          if ~strcmp(get(parentfigure,'defaultaxesunits'),'normalized')
              tmp = axes;
              set(tmp,'units','normalized')
              def_pos = get(tmp,'position');
              delete(tmp)
          else
              def_pos = get(parentfigure,'DefaultAxesPosition');
          end
          setappdata(parentfigure,'SubplotDefaultAxesLocation',def_pos);
      else
          def_pos = getappdata(parentfigure,'SubplotDefaultAxesLocation');
      end

      % compute outerposition and insets relative to figure bounds
      rw = max(row)-min(row)+1;
      cw = max(col)-min(col)+1;
      width = def_pos(3)/(ncols - inset(1) - inset(3));
      height = def_pos(4)/(nrows - inset(2) - inset(4));
      inset = inset.*[width height width height];
      outerpos = [def_pos(1) + min(col)*width - inset(1), ...
                  def_pos(2) + min(row)*height - inset(2), ...
                  width*cw height*rw];

      % adjust outerpos and insets for axes around the outside edges
      if min(col) == 0
        inset(1) = def_pos(1); 
        outerpos(3) = outerpos(1)+outerpos(3); 
        outerpos(1) = 0; 
      end
      if min(row) == 0, 
        inset(2) = def_pos(2);
        outerpos(4) = outerpos(2)+outerpos(4); 
        outerpos(2) = 0;
      end
      if max(col) == ncols-1, 
        inset(3) = max(0,1-def_pos(1)-def_pos(3));
        outerpos(3) = 1-outerpos(1);
      end
      if max(row) == nrows-1, 
        inset(4) = max(0,1-def_pos(2)-def_pos(4));
        outerpos(4) = 1-outerpos(2);
      end

      % compute inner position
      position = [outerpos(1:2) + inset(1:2),...
                  outerpos(3:4) - inset(1:2) - inset(3:4)];

    end
  end
end

% kill overlapping siblings if mnp specifier was used:
nextstate = get(parentfigure,'nextplot');

if strncmp(nextstate,'replace',7)
  nextstate = 'add'; 
elseif strncmp(nextstate,'new',3)
  kill_siblings = 0;
end

if(kill_siblings)
  if delay_destroy
    if nargout
      error(id('TooManyOutputs'),...
            'Function called with too many output arguments.')
    else
      set(parentfigure,'NextPlot','replace'); 
      return
    end
  end
  sibs = datachildren(parent);
  newcurrent = [];
  for i = 1:length(sibs)
    % Be aware that handles in this list might be destroyed before
    % we get to them, because of other objects' DeleteFcn callbacks...
    if(ishandle(sibs(i)) && strcmp(get(sibs(i),'Type'),'axes'))
      units = get(sibs(i),'Units');
      sibpos = get(sibs(i),'Position');
      if ~strcmp(units,'normalized')
          sibpos = hgconvertunits(parentfigure,sibpos,units,'normalized',parent);
      end
      intersect = 1;
      if((position(1) >= sibpos(1) + sibpos(3)-tol) || ...
         (sibpos(1) >= position(1) + position(3)-tol) || ...
         (position(2) >= sibpos(2) + sibpos(4)-tol) || ...
         (sibpos(2) >= position(2) + position(4)-tol))
        intersect = 0;
      end
      if intersect
        % if already found the position match axis (got_one)
        % or if this intersecting axis doesn't have matching pos
        % or if we know we want to make a new axes anyway (kill_siblings==2)
        % delete it
        if any(abs(sibpos - position) > tol) || kill_siblings==2
            % if this axis is in the grid the don't delete it
            grid = getappdata(parent,'SubplotGrid');
            if isempty(grid) ||  ~any(grid(:) == sibs(i)) || ...
                size(grid,1) ~= nrows || size(grid,2) ~= ncols || ...
                ~isscalar(row) || ~isscalar(col) || kill_siblings==2
                delete(sibs(i));
            end
        end
        if ishandle(sibs(i))
            % if this axes overlaps the other one exactly then
            if ~isempty(newcurrent) && ishandle(newcurrent)
                delete(newcurrent);
            end
            newcurrent = sibs(i);
        end
      end
    end
  end
  if ~isempty(newcurrent) && ishandle(newcurrent)
      set(parentfigure,'CurrentAxes',newcurrent);
      create_axis = false;
  end
  set(parentfigure,'NextPlot',nextstate);
end

% create the axis:
if create_axis
  if strcmp(nextstate,'new') && ~explicitParent
    parent = figure;
    parentfigure = parent;
  end
  ax = axes('units','normal','Position',position,...
            'LooseInset',inset,'Parent',parent,pvpairs{:});
  set(ax,'units',get(parentfigure,'defaultaxesunits'))
  if useAutoLayout
    addAxesToGrid(ax,nrows,ncols,row,col,position);
  end
elseif move_axis && ~preventMove
  ax = handle;
  units = get(handle,'units');
  set(handle,'units','normal','Position',position,...
             'LooseInset',inset,'Parent',parent,pvpairs{:});
  set(handle,'units',units);
  if useAutoLayout
    addAxesToGrid(ax,nrows,ncols,row,col,position);
  end
else
  % this should only happen with subplot(H)
  ax = get(parentfigure,'CurrentAxes'); 
end

% return identifier, if requested:
if(nargout > 0)
    theAxis = ax;
end

% Create subplot listeners to align plot boxes automatically
function createListeners(p,axlisth)
setappdata(p,'SubplotListeners',[])
fig = p;
if ~isequal(get(fig,'Type'),'figure')
  fig = ancestor(fig,'figure');
end
list = [...
    handle.listener(axlisth,findprop(axlisth(1),'Units'), ... % must be first
                    'PropertyPostSet',@axesUnitsPostSet);
    handle.listener(axlisth,findprop(axlisth(1),'Units'), ...
                    'PropertyPreSet',@axesUnitsPreSet);
    handle.listener(axlisth,findprop(axlisth(1),'Position'), ...
                    'PropertyPostSet',@axesMoved);
    handle.listener(axlisth,findprop(axlisth(1),'ActivePositionProperty'), ...
                    'PropertyPreSet',@axesMoved);
    handle.listener(axlisth,findprop(axlisth(1),'Parent'), ...
                    'PropertyPreSet',@axesMoved);
    handle.listener(axlisth,'AxisInvalidEvent',{@subplotlayoutInvalid,p});
    handle.listener(handle(fig),'FigureUpdateEvent',{@subplotlayout,p})];
for k=1:length(axlisth)
  ax = axlisth(k);
  if ~isappdata(double(ax),'SubplotDeleteListener')
    setappdata(double(ax),'SubplotDeleteListener',...
                          handle.listener(ax,'ObjectBeingDestroyed', ...
                                          @axesDestroyed));
  end
end
setappdata(p,'SubplotListeners',list)

% Add ax to a matrix of handles in the specified location.
% The grid is stored on the parent appdata.
% Also store the insets in ax appdata.
% Only stores the axes if it is in a 1-by-1 cell and
% the grid size matches any existing grid.
function addAxesToGrid(ax,nrows,ncols,row,col,position)
p = get(ax,'parent');
grid = getappdata(p,'SubplotGrid');
if isempty(grid)
  grid = zeros(nrows,ncols);
end
if any(size(grid) ~= [nrows ncols]), return; end
if length(row) ~= 1 || length(col) ~= 1, return; end
if round(row) ~= row || round(col) ~= col, return; end
if grid(row+1,col+1) == ax, return, end
grid(row+1,col+1) = ax;
list = grid(:);
list(list == 0) = []; % remove root
list(~ishandle(list)) = []; % remove invalid handles
createListeners(p,handle(list));
setappdata(p,'SubplotGrid',grid)
setappdata(ax,'SubplotPosition',position); % normalized
subplotlayoutInvalid(handle(ax),[],p);

% Remove ax from grid of subplots in p
function removeAxesFromGrid(p,ax)
grid = getappdata(p,'SubplotGrid');
if ~isempty(grid)
  n = grid == ax;
  if any(n(:))
    grid(n) = 0;
    list = grid(:);
    list(list == 0) = []; % remove root
    list(~ishandle(list)) = [];
    if isempty(list)
      rmappdata(p,'SubplotListeners');
      rmappdata(p,'SubplotGrid');
    else
      setappdata(p,'SubplotGrid',grid);
    end
  end
end

% Callback when axis moves to remove it from subplot layout grid
function axesMoved(hSrc,evdata) %#ok
ax = double(evdata.affectedObject);
% If the legend or colorbar is causing the move, do not remove the axes
% from the subplot grid. Do, however, update it's cached position:
if isappdata(ax,'inLayout') && ~isempty(getappdata(ax,'inLayout'))
    setappdata(ax,'SubplotPosition',get(ax,'Position'));
else
    removeAxesFromGrid(get(ax,'Parent'),ax);
end

% Callback when axis changes units
function axesUnitsPreSet(hSrc,evdata) %#ok
ax = double(evdata.affectedObject);
p = get(ax,'Parent');
list = getappdata(p,'SubplotListeners');
if ~isempty(list)
    set(list(2:end),'enable','off');
end

% Callback when axis done changing units
function axesUnitsPostSet(hSrc,evdata) %#ok
ax = double(evdata.affectedObject);
p = get(ax,'Parent');
list = getappdata(p,'SubplotListeners');
if ~isempty(list)
    set(list(2:end),'enable','on');
end

% Callback when axis destroyed 
function axesDestroyed(hSrc,evdata) %#ok
ax = double(hSrc);
p = get(ax,'Parent');
if strcmp(get(p,'BeingDeleted'),'off')
  removeAxesFromGrid(p,ax);
elseif isappdata(p,'SubplotListeners')
  rmappdata(p,'SubplotListeners');
  rmappdata(p,'SubplotGrid');
end

function str = id(str)
str = ['MATLAB:subplot:' str];
