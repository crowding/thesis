function fitAll(infile, outfile, splits, subset)

if ~exist('infile', 'var')
    infile = fullfile(fileparts(mfilename('fullpath')), 'data.mat');
end

if ~exist('outfile', 'var')
    outfile = fullfile(fileparts(mfilename('fullpath')), 'fits.mat');
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

models = {BoyntonModel()};
models = cellfun(@(x) x.fit(data), models, 'UniformOutput', 0);
modelNames = {'BoyntonModel'};

%for backward compatibility right now.
fits = models{1}.parameters;

save(outfile, 'models', 'modelNames');

end