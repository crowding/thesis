function makeDiagnostics(fitfile, outfiles)

% There are some aspects of the data that I think the model does not
% capture well, and I want to show it.

% What I'd like are residual or deviance plots.

% Side note: you know, this is going to take renamed data because the
% model takes renamed data. Perhaps I should make "folding" or
% "unfolding" a model paramter.

% these are the variables the model operates over.

%now plot this data over these marginals

vars = {'subject', 'content', 'spacing', 'dx', 'response', 'exp_type'};

if ~exist('fitfile', 'var') || isempty(fitfile)
    fitfile = fullfile(fileparts(mfilename('fullpath')), 'fits.mat');
end

if ~exist('outfile', 'var')
    outfile = fullfile(fileparts(mfilename('fullpath')), 'diagnostics.list');
end

fh = fopen(outfile, 'w');
onCleanup(@(x)close(fh));

load(fitfile, 'models', 'modelNames');

model = models{1};
modelName = modelNames{1};

split = {'subject', 'spacing', 'exp_type'};
binsize = 25;
binvar = 'dx';

resid = model.residuals(split, binvar, binsize);
handles = plotResiduals(resid, 'subject', 'dx', 'pearson_resid', ...
                        'exp_type', 'const_', 'spacing', 'n_obs');


tile([],[],[],unique(handles.fig));

groupfun(handles, 'fig', @writeFigure);
function writeFigure(handles)
    fig = handles.fig(1);
    sub = handles.subject{1};
    set(fig, 'PaperOrientation', 'landscape', ...
             'PaperPosition', [0.25 0.25 10.5 8] );
    outfile = fullfile(fileparts(outfile), ...
                       sprintf('fit_%s_%s.pdf', ...
                               modelName, sub));
    print(fig,'-dpdf', outfile);
    fprintf(fh, '%s\n', outfile);
end

end
