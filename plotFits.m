function f = plotFits(fitfile, outfile)

    if ~exist('fitfile', 'var') || ~isempty(fitfile)
       fitfile = fullfile(fileparts(mfilename('fullpath')), 'fits.mat');
    end

    if ~exist('outfile', 'var')
        outfile = fullfile(fileparts(mfilename('fullpath')), 'fit_plots.list');
    end

    load(fitfile, 'models', 'modelNames');

    split = {'subject', 'dx', 'content', 'spacing', 'exp_type'};

    fh = fopen(outfile, 'w');
    C = onCleanup(@() close(fh));

    cellfun(@plotModel, models, modelNames)
    function plotModel(model, modelName)
        %break down the trials in terms of number yes and number no per unique
        %condition
        rates = countdata(model.data, split);

        %make a lattice plot for all the data. It says plotResiduals but it's
        %just a scatter plot.
        %                              figure,    x     y    row
        handles = plotResiduals(rates, 'subject', 'dx', 'p', 'exp_type', ...
                                'content', 'spacing', 'n', @addPredictions);
        %                       col        color      size

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
                line = join(line, dx_eval, 'type', 'inner', 'MergeKeys', true);
                line.p = model.fullPredict([], line);
            end

            %and add these predictions to the graph. They are keyed by color.
            hold(handle.ax(1), 'on');
            groupfun(pred, 'line', @predLine);
            %this doesn't return the actual graphics handles corresponding to
            %the predicted line though.
            function predLine(line)
                [x, o] = sort(line.dx);
                y = line.p(o);
                plot(handle.ax(1), x, y, 'color', line.color_(1,:));
            end
            hold(handle.ax(1), 'off');
        end

        %finally, print the figures
        groupfun(handles, 'fig', @writeFigure);
        function writeFigure(handles)
            fig = handles.fig(1);
            sub = handles.subject{1};
            set(fig, 'PaperOrientation', 'landscape', ...
                     'PaperPosition', [0.25 0.25 10.5 8] );
            outfile = fullfile(fileparts(outfile), ...
                               sprintf('fit_%s_%s.pdf', modelName, sub));
            print(fig,'-dpdf', outfile);
            fprintf(fh, '%s\n', outfile);
        end
    end
end