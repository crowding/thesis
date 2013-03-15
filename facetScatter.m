function handles = ...
    facetScatter(data, varargin)
    %function handles = ...
    %    facetScatter(data, figVar, xVar, yVar, rowVar, ...
    %                 colVar, colorVar, sizeVar, morePlotting);
    %
    %Build up a big lattice-type scatterplot.
    %figures, etc.
    %
    %use 'const_' in a var to not marginalize over that var.
    %
    %'morePlotting', if present is a function handle that will take
    %the data chunk and list of graphics handles and be called once
    %per subplot.

    %Does not yet DTRT with making a legend and titles and axis labels in
    %all the right places.
    pars = inputParser();
    pars.addParamValue('fig', 'const_');
    pars.addParamValue('x', 'const_');
    pars.addParamValue('y', 'const_');
    pars.addParamValue('row', 'const_');
    pars.addParamValue('col', 'const_');
    pars.addParamValue('color', 'const_');
    pars.addParamValue('size', 'const_');
    pars.addParamValue('title', 'insert title here')
    pars.addParamValue('morePlotting', []);
    pars.addParamValue('newFigures', false);
    pars.parse(varargin{:});

    figVar = pars.Results.fig;
    xVar = pars.Results.x;
    yVar = pars.Results.y;
    rowVar = pars.Results.row;
    colVar = pars.Results.col;
    colorVar = pars.Results.color;
    sizeVar = pars.Results.size;
    titleLabel = pars.Results.title;
    morePlotting = pars.Results.morePlotting;
    newFigures = pars.Results.newFigures;

    data.const_ = ones(size(data, 1),1);

    %need to set the xlim globally.

    [~, figno, figix] = unique(data.(figVar), 'first');
    if newFigures
        fignumbers = max([1; get(0, 'Children')]) + (1:numel(figno));
    else
        fignumbers = 1:numel(figno);
    end
    fighandles = arrayfun(@openFigure, fignumbers(:));
    function no = openFigure(no)
        figure(no);
        clf;
    end
    data.fighandle_ = fighandles(figix);

    %I want the rows and columns to be assigned per figure. xlim and ylim
    %should vary per figure also.
    data = groupfun(data, figVar, @setupMappings);
    function figData = setupMappings(figData)
        [rowVals, ~, figData.figRow_] = unique(figData.(rowVar));
        [colVals, ~, figData.figCol_] = unique(figData.(colVar));
        [colorVals, ~, figData.color_] = unique(figData.(colorVar));
        figData.x_ = figData.(xVar);
        figData.y_ = figData.(yVar);
        figData.size_ = figData.(sizeVar);
        figData.nRows_ = zeros(size(figData, 1),1) + length(rowVals);
        figData.nCols_ = zeros(size(figData, 1),1) + length(colVals);

        if length(colorVals) > 1
            colormap = cool(length(colorVals));
        else
            colormap = [0 0 0];
        end
        figData.color_ = colormap(figData.color_, :);

        figData.xmin_ = zeros(size(figData, 1),1) + min(figData.(xVar));
        figData.xmax_ = zeros(size(figData, 1),1) + max(figData.(xVar));
        figData.ymin_ = zeros(size(figData, 1),1) + min(figData.(yVar));
        figData.ymax_ = zeros(size(figData, 1),1) + max(figData.(yVar));
    end

    %compute lowest-most and leftest-most graphs
    lowlefts = unique(data(:,{'fighandle_' 'figRow_' 'figCol_'}));
    %lowlefts = groupfun(lowlefts, {'fighandle_', 'figCol_'}, @addLowest);
    %lowlefts = groupfun(lowlefts, {'fighandle_', 'figRow_'}, @addLeftest);
    lowlefts = groupfun(lowlefts, {'fighandle_'}, @addLowest);
    lowlefts = groupfun(lowlefts, {'fighandle_'}, @addLeftest);
    function chunk = addLowest(chunk)
        chunk.lowest = (chunk.figRow_ == max(chunk.figRow_));
    end
    function chunk = addLeftest(chunk)
        chunk.leftest = (chunk.figCol_ == min(chunk.figCol_));
    end

    handles = groupfun(data, {figVar rowVar colVar}, @subplot);
    function handles = subplot(chunk)
        handles.fig = chunk.fighandle_(1);
        figure(handles.fig);
        handles.ax = subplotix(chunk.nRows_(1), chunk.nCols_(1), ...
                               chunk.figRow_(1), chunk.figCol_(1));
        handles.handle = {...
            gscatter(chunk.x_, chunk.y_, ...
                     1:size(chunk,1), chunk.color_, '.', ...
                     sqrt(chunk.size_) * 5, 0)};
        set(handles.ax, 'XLim', [chunk.xmin_(1) chunk.xmax_(1)], ...
                        'YLim', [chunk.ymin_(1) chunk.ymax_(1)]);
        %enable axes, legend, h/v based on which subplot.
        lowleft = join(chunk(1, {'fighandle_' 'figRow_' 'figCol_'}), ...
                       lowlefts, 'Type', 'Inner', 'MergeKeys', true);

        annotation = '';
        xl = xVar;
        if ~strcmp(colVar, 'const_')
            val = chunk.(colVar)(1);
            if iscell(val)
                val = val{1};
            end
            %                xl = sprintf('%s\n(%s = %s)', ...
            %             xl, colVar, num2str(val,3));
            annotation = sprintf('%s = %s\n%s', ...
                                 colVar, num2str(val,3), annotation);
        end

        yl = yVar;
        if ~strcmp(rowVar, 'const_')
            val = chunk.(rowVar)(1);
            if iscell(val)
                val = val{1};
            end
            %yl = sprintf('%s\n(%s = %s)', ...
            %             yl, rowVar, num2str(val,3));
            annotation = sprintf('%s = %s\n%s', ...
                                 rowVar, num2str(val,3), annotation);
        end

        if lowleft.lowest
            xlabel(handles.ax, xl);
        else
            set(handles.ax, 'XTickLabel', []);
        end

        if lowleft.leftest
            ylabel(handles.ax, yl);
        else
            set(handles.ax, 'YTickLabel', []);
        end

        text(chunk.xmin_(1), chunk.ymax_(1), annotation, ...
             'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
             'interpreter', 'none');
        %'other' is for other functions you may want to plot while in this axis.
        if exist('morePlotting', 'var') && ~isempty(morePlotting)
            handles = morePlotting(handles, chunk);
        end
    end

    % create a central plot title
    groupfun(handles, {'fig'}, @addTitle);
    function addTitle(chunk)
        f = chunk.fig(1);

        if isa(titleLabel, 'function_handle')
            tit = titleLabel(handles, chunk);
        else
            tit = titleLabel;
        end
        if ~strcmp(figVar, 'const_')
            d = join(data(:,{'fighandle_', figVar}), dataset({f, 'fighandle_'}), ...
                     'Type', 'inner', 'mergeKeys', true);
            val = d.(figVar)(1);
            if iscell(val)
                val = val{1};
            end
            tit = sprintf('%s\n%s = %s', tit, figVar, num2str(val,3));
        end
        set(0, 'CurrentFigure', f);
        set(f,'NextPlot','add');
        ax = axes();
        h = title(tit);
        set(gca,'Visible','off');
        set(h,'Visible','on');

        %now while we're at it we'll make a global legend.
        colorset = unique(join(data(:, {'fighandle_', 'color_', colorVar}), ...
                        dataset({f, 'fighandle_'}), ...
                        'Type', 'inner', 'mergeKeys', true));

        %How to make a global legend? Plot invisible lines using the needed
        %colors into the invisible axes we just made.
        hold(ax, 'on');
        labels = maybe_num2str(colorset.(colorVar));
        for i = 1:size(colorset, 1)
            handle(i) = plot(ax, [1 2], [1 2], '.-', ...
                             'Color', colorset{i, 'color_'}, ...
                             'DisplayName', labels{i});
            set(handle(i), 'Visible', 'off');
        end
        l = legend(ax, 'Location', 'West');
        pos = get(l, 'Position');
        pos(1) = [0.05];
        set(l, 'Position', pos);
        set(get(l, 'Title'), 'string', colorVar);
        han = struct('fig', f, 'ax', ax, 'handle', handle, 'type', {'title'});
    end
end