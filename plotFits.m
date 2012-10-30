function f = plotFits(fitfile, outfile)

    if ~exist('fitfile', 'var') || ~isempty(fitfile)
       fitfile = fullfile(fileparts(mfilename('fullpath')), 'fits.mat');
    end

    if ~exist('outfile', 'var')
        outfile = fullfile(fileparts(mfilename('fullpath')), 'fit_plots.list');
    end

    load(fitfile, 'models', 'modelNames');
    modelNames = modelNames(:);

    fh = fopen(outfile, 'w');
    C = onCleanup(@() fclose(fh));

    %plot all models, one window per subject
    handles = cellfun(@(x,y)plotModel(x,'plotBy', ...
                                      {'fig', 'subject', 'title', y, 'col', 'exp_type', ...
                                       'x', 'dx', 'y', 'p', 'row', 'spacing', ...
                                       'color', 'content', 'size', 'n'}), ...
                      models, modelNames, 'UniformOutput', 0);

    tile(1,1);

    for i = 1:numel(modelNames)
        ln = size(handles{1},1);
        handles{i}.modelName = modelNames(zeros(ln,1) + i);
    end

    handles = cat(1, handles{:});

    %finally, print the figures
    groupfun(handles, 'fig', @writeFigure);
    function writeFigure(handles)
        fig = handles.fig(1);
        sub = handles.subject{1};
        set(fig, 'PaperOrientation', 'landscape', ...
                 'PaperPosition', [0.25 0.25 10.5 8] );
        outfile = fullfile(fileparts(outfile), ...
                           sprintf('fit_%s_%s.pdf', handles.modelName{1}, sub))
        print(fig,'-dpdf', outfile);
        fprintf(fh, '%s\n', outfile);
    end

end