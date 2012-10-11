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
    pars.addParamValue('morePlotting', []);
    pars.parse(varargin{:});

    figVar = pars.Results.fig;
    xVar = pars.Results.x;
    yVar = pars.Results.y;
    rowVar = pars.Results.row;
    colVar = pars.Results.col;
    colorVar = pars.Results.color;
    sizeVar = pars.Results.size;
    morePlotting = pars.Results.morePlotting;

    data.const_ = ones(size(data, 1),1);

    %need to set the xlim globally.

    [~, figno, figix] = unique(data.(figVar), 'first');
    fighandles = arrayfun(@figure, (1:numel(figno))');
    data.fighandle_ = fighandles(figix);

    %I want the rows and columns to be assigned per figure. xlim and ylim
    %should vary per figure also.
    data = groupfun(data, figVar, @setupMappings);
    function figData = setupMappings(figData)
        [rowVals, ~, figData.figRow_] = unique(figData.(rowVar));
        [colVals, ~, figData.figCol_] = unique(figData.(colVar));
        [colorVals, ~, figData.color_] = unique(figData.(colorVar));
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

    handles = groupfun(data, {figVar rowVar colVar}, @subplot);
    function handles = subplot(chunk)
        handles.fig = chunk.fighandle_(1);
        figure(handles.fig);
        handles.ax = subplotix(chunk.nRows_(1), chunk.nCols_(1), ...
                               chunk.figRow_(1), chunk.figCol_(1));
        handles.handle = {...
            gscatter(chunk.(xVar), chunk.(yVar), ...
                     1:size(chunk,1), chunk.color_, '.', ...
                     sqrt(chunk.(sizeVar)) * 5, 0)};
        set(handles.ax, 'XLim', [chunk.xmin_(1) chunk.xmax_(1)], ...
                        'YLim', [chunk.ymin_(1) chunk.ymax_(1)]);
        %enable axes, legend, h/v based on which subplot.
        if chunk.figRow_(1) == chunk.nRows_(1);
            xl = xVar;
            if ~strcmp(colVar, 'const_')
                val = chunk.(colVar)(1);
                if iscell(val)
                    val = val{1};
                end
                xl = sprintf('%s\n(%s = %s)', ...
                             xl, colVar, num2str(val,3));
            end
            xlabel(handles.ax, xl);
        else
            set(handles.ax, 'XTickLabel', []);
        end

        if chunk.figCol_(1) == 1
            yl = yVar;
            if ~strcmp(rowVar, 'const_')
                val = chunk.(rowVar)(1);
                if iscell(val)
                    val = val{1};
                end
                yl = sprintf('%s\n(%s = %s)', ...
                             yl, rowVar, num2str(val,3));
            end
            ylabel(handles.ax, yl);
        else
            set(handles.ax, 'YTickLabel', []);
        end

        if chunk.figRow_(1) == 1
            tit = 'insert title here';
            if ~strcmp(figVar, 'const_')
                val = chunk.(figVar)(1);
                if iscell(val)
                    val = val{1};
                end
                tit = sprintf('%s\n%s = %s', tit, figVar, num2str(val,3));
            end
            if ~strcmp(colorVar, 'const_')
                %in lieu of placing a legend just yet
                tit = sprintf('%s\n(colors indicate %s)', tit, colorVar);
            end
            title(handles.ax, tit);
        end

        %'other' is for other functions you may want to plot while in this axis.
        if exist('morePlotting', 'var') && ~isempty(morePlotting)
            morePlotting(handles, chunk);
        end
        %TODO; indicate what size means.
    end
end
