function f = plotFits(fitfile, datafile, outfile)

    if ~exist('fitfile', 'var') || ~isempty(fitfile)
       fitfile = fullfile(fileparts(mfilename('fullpath')), 'fits.mat');
    end

    if ~exist('datafile', 'var') || ~isempty(datafile)
       datafile = fullfile(fileparts(mfilename('fullpath')), 'data.mat');
    end

    if ~exist('outfile', 'var')
        outfile = fullfile(fileparts(mfilename('fullpath')), 'fit_plots.list')
    end

    load(fitfile, 'fits');
    load(datafile, 'data');
    data = doRename(data);
    split = {'subject', 'dx', 'content', 'spacing', 'exp_type'};


    %plot only data that matches one of our fits.
    [~, ia, ib] = join(data, fits, 'type', 'inner');
    data = data(unique(ia),:);
    fits = fits(unique(ib),:);

    %break down the trials in terms of number yes and number no per unique
    %condition
    rates = countdata(data, split);

    %make a lattice plot for all the data. It says plotResiduals but it's
    %just a scatter plot.
    %                              figure,    x     y    row
    handles = plotResiduals(rates, 'subject', 'dx', 'p', 'exp_type', ...
    ...%                    col        color      size
                            'content', 'spacing', 'n', @addPredictions);


    %for each subplot, plot over dx and add handles.
    function [handle] = addPredictions(handle, chunk)
        %make prediction lines on this subplot based on the data
        chunk.dx = [];
        chunk.n = [];
        chunk.p = [];
        chunk = unique(chunk);
        chunk.line = (1:size(chunk, 1))';

        dx_eval = dataset(...
            {linspace(chunk.xmin_(1), chunk.xmax_(1), 100)', 'dx'},...
            {ones(100, 1), 'const_'});
        pred = groupfun(chunk, 'line', @predictLine);
        function line = predictLine(line)
            param = join(line, fits, 'type', 'inner', 'MergeKeys', true);
            line = join(line, dx_eval, 'type', 'inner', 'MergeKeys', true);
            line.p = MotionModel(param, line);
        end

        %and add these predictions to the graph. They are keyed by color.
        hold(handle.ax(1), 'on');
        groupfun(pred, 'line', @predLine);
        %this doesn't return the actual graphics handles corresponding to the
        %predicted line though.
        function predLine(line)
            [x, o] = sort(line.dx);
            y = line.p(o);
            plot(handle.ax(1), x, y, 'color', line.color_(1,:));
        end
        hold(handle.ax(1), 'off');
    end

    %finally, print the figures
    fh = fopen(outfile, 'w');
    C = onCleanup(@() close(fh));
    groupfun(handles, 'fig', @writeFigure);
    function writeFigure(handles)
        fig = handles.fig(1);
        sub = handles.subject{1};
        set(fig, 'PaperOrientation', 'landscape', ...
                 'PaperPosition', [0.25 0.25 10.5 8] );
        outfile = fullfile(fileparts(outfile), ...
                           sprintf('fit_%s.pdf', sub));
        print(fig,'-dpdf', outfile);
        fprintf(fh, '%s\n', outfile);
    end

end