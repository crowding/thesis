function fitAll(infile, outfile, splits, subset)

if ~exist('infile', 'var')
    infile = fullfile(fileparts(mfilename('fullpath')), 'data.mat');
end

if ~exist('outfile', 'var')
    outfile = fullfile(fileparts(mfilename('fullpath')), 'fits.mat');
end

if ~exist('splits', 'var')
    %By default do independent fits for each subject.
    splits = {'subject'};
end

if ~exist('subset', 'var')
    %By default restrict to only content and spacing experiments.
    subset = dataset({{'content';'spacing'},'exp_type'});
end

trial_threshold = 2000;
%Only bother with subjects that have this many trials.

%model fits are without folding.
load(infile, 'data');
data = doRename(data, false);

%Select the subset.
data = join(data, subset, 'type', 'inner', 'MergeKeys', true);

%Filter on trial count.
trial_counts = grpstats(data, 'subject', @numel, 'DataVars', 'response');
good_subjects = trial_counts.subject(trial_counts.GroupCount > trial_threshold);
data=data(ismember(data.subject, good_subjects),:);

models = {SlopeModel('splits', splits, ...
                     'freeParams', {'mu_0', 'beta_0', 'cs', 'beta_summation', ...
                                    'saturating_induced', 'wiggle_induced'}, ...
                     'data', data)};

models = cellfun(@(x) x.fit(), models, 'UniformOutput', 0);
modelNames = {'SlopeModel'};

%A separate variable with the parameters.
fits = models{1}.parameters;

save(outfile, 'models', 'modelNames', 'fits');

end