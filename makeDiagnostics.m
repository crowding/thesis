function makeDiagnostics(outfiles)

% There are some aspects of the data that I think the model does not
% capture well, and I want to show it.

% What I'd like are residual or deviance plots.

% Side note: you know, this is going to take renamed data because the
% model takes renamed data. Perhaps I should make "folding" or
% "unfolding" a model paramter.

% these are the variables the model operates over.


%now plot this data over these marginals

vars = {'subject', 'content', 'spacing', 'dx', 'response', 'exp_type'};

if ~exist('data', 'var') || isempty(data)
    load('data.mat', 'data')
    data = rename(data, ...
                  'target_spacing', 'spacing', ...
                  'folded_direction_content', 'content', ...
                  'folded_displacement', 'dx', ...
                  'folded_response_with_carrier', 'response');
    data = data( strcmp(data.exp_type, 'spacing') ...
                & strcmp(data.subject, 'pbm'), vars);
end

if ~exist('params', 'var') || isempty(params)
    load('modelResults.mat', 'p')
    params = p;
end

if ~exist('model', 'var') || isempty(model)
    model = @MotionModel;
end

if ~isa(data, 'dataset')
    data = dataset(data);
end

if ~exist('binsize', 'var')
    binsize = 25;
end

if ~exist('binvar', 'var')
    binvar = 'dx';
end

if ~exist('split', 'var') || split == 0
    split = {'subject', 'content', 'spacing'};
end

resid = residuals(data, model, params, split, binvar, binsize);
handles = plotResiduals(resid, 'const_', 'dx', 'pearson_resid', ...
                        'content', 'const_', 'spacing');
nfigs = unique(unique(handles.fig));
ncol = ceil(sqrt(nfigs));
nrow = ceil(nfigs/ncol);
tile(nrow,ncol,[],unique(handles.fig));

%list the outputs produced
if exist('outfile', 'var')
    fh = open(outfile, 'w');
    for i = unique(handles.fig)
        outfile = fullfile(fileparts(outfile), sprintf('residuals_%d', i));
        fprintf(fh, '%s\n', outfile);
        print('-depsc', outfile);
    end
end

end
 