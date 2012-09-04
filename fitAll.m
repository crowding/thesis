function fitAll(infile, outfile, splits, subset)

if ~exist('filename', 'var')
    infile = fullfile(fileparts(mfilename), 'data.mat');
end

if ~exist('outfile', 'var')
    outfile = fullfile(fileparts(mfilename), 'fits.mat');
end

if ~exist('splits', 'var')
    %By default do independent fits for each subject and experiment type.
    splits = {'subject', 'exp_type'};
end

if ~exist('subset', 'var')
    %By default restrict to only content and spacing experiments.
    subset = dataset({{'content';'spacing'},'exp_type'});
end

load(infile, 'data');
data = doRename(data);

%Select the subset.
data = join(data, subset, 'type', 'inner', 'MergeKeys', true);

%Look up some initial parameters. Fit the model independently per
%subject and experiment type.
params = initialParams(data, splits);

freeParams = {'sig0','sigk','siga','mua','mukc','muks'};

fits = groupfun(data, splits, @makeFit)
function f = makeFit(chunk)
    split = chunk(1,splits)
    chunkParams = join(params, split, 'type', 'inner', 'MergeKeys', true);
    f = fit(@fitMotionModel, chunkParams, freeParams, chunk);
end

save(outfile, 'fits');

end