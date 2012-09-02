function handles = ...
    plotResiduals(data, figVar, xVar, yVar, rowVar, colVar, colorVar)
    %function handles = ...
    %    plotResiduals(data, figVar, xVar, yVar, rowVar, colVar, colorVar)
    %Build up a big lattice-type scatterplot.
    %and here's a function for plotting the residuals over some
    %variables. Returns a dataset of handles to axes, linesseries,
    %figures, etc.
    %use 'const_' to not marginalize over that var.
    data.const_ = ones(size(data, 1),1);

    %need to set the xlim globally.
    xlim = [min(data.(xVar)) max(data.(xVar))];
    ylim = [min(data.(yVar)) max(data.(yVar))];

    [~, figno, figix] = unique(data.(figVar), 'first');
    fighandles = arrayfun(@figure, 1:numel(figno));
    data.fighandle_ = fighandles(figix);

    [rowVals, ~, data.figRow_] = unique(data.(rowVar));
    [colVals, ~, data.figCol_] = unique(data.(colVar));
    [colorVals, ~, data.color_] = unique(data.(colorVar));
    colormap = cool(length(colorVals));
    data.color_ = colormap(data.color_, :);
    nrow = length(rowVals);
    ncol = length(colVals);

    handles = groupfun(data, {figVar rowVar colVar}, @subplot);
    function handles = subplot(chunk)
        handles.fig = chunk.fighandle_(1);
        figure(handles.fig);
        handles.ax = subplotix(length(rowVals), length(colVals), ...
                               chunk.figRow_(1), chunk.figCol_(1));
        handles.handle = {...
            gscatter(chunk.(xVar), chunk.(yVar), ...
                     1:size(chunk,1), chunk.color_, '.', ...
                     sqrt(chunk.n_obs) * 5, 0)};
        set(handles.ax, 'XLim', xlim, 'YLim', ylim);
        %enable axes, legend, h/v based on which subplot.
        if chunk.figRow_(1) == numel(rowVals)
            xlabel(handles.ax, ...
                   sprintf('%s\n(%s = %s)', ...
                           xVar, colVar, num2str(chunk.(colVar)(1),3)));
        else
            set(handles.ax, 'XTickLabel', []);
        end

        if chunk.figCol_(1) == 1
            ylabel(handles.ax, ...
                   sprintf('(%s = %s)\n%s', ...
                           rowVar, num2str(chunk.(rowVar)(1),3), yVar));
        else
            set(handles.ax, 'YTickLabel', []);
        end

        if chunk.figRow_(1) == 1
            tit = 'Residual values for model fit';
            if ~strcmp(colorVar, 'const_')
                tit = sprintf('%s\n%s = %s', tit, figVar, chunk.(figVar)(1));
            end
            if ~strcmp(colorVar, 'const_')
                %in lieu of placing a legend just yet
                tit = sprintf('%s\n(colors indicate %s)', tit, colorVar);
            end
            title(handles.ax, tit);
        end
    end
end
