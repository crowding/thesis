load('data.mat');
data = doRename(data, false); %no folding

%only take data from spacing and direction-content experiments
data = data(ismember(data.exp_type, {'spacing', 'content'}), :);

%and only use subjects for whon we have >2000 trials
trial_counts = grpstats(data, 'subject', @numel, 'DataVars', 'response')
good_subjects = trial_counts.subject(trial_counts.GroupCount > 2000);
data=data(ismember(data.subject, good_subjects),:);

%fit the model(s)
M = fit(SlopeModel('splits', {'subject'}, ...
                   'freeParams', {'mu_0', 'beta_0', 'cs', 'beta_summation', ...
                                  'saturating_induced', 'wiggle_induced'}, ...
                   'data', data));

%make plots of the models
plotModel(M, 'plotBy', {'fig', 'subject', 'col', 'exp_type', ...
                       'x', 'dx', 'y', 'p', 'row', 'spacing', ...
                       'color', 'content', 'size', 'n'});
tile(1,1);

%let's show contours of direction content vs. spacing (fixed
%at delta-x=0) for each subject.

sMin = min(data.spacing)
sMax = max(data.spacing)
cMin = min(data.content)
cMax = max(data.content)

[content, spacing, const_, dx] = ndgrid(linspace(cMin, cMax, 50)', ...
                                        linspace(sMin, sMax, 50)', ...
                                        1, 1);

mesh = dataset({content(:), 'content'}, {spacing(:), 'spacing'}, ...
               {dx(:), 'dx'}, {const_(:), 'const_'});
subjs = unique(data(:,'subject'));
subjs.const_ = repmat(1, size(subjs, 1), 1);
mesh = join(mesh, subjs, 'Type', 'Inner', 'MergeKeys', true);
mesh.pred = M.fullPredict(M.parameters,mesh);


for s in unique(pred.subject)'
    figure()
    subset = mesh(strcmp(mesh.subject, s), :);
    pred = reshape(subset.pred, size(content));
    contour(spacing, content, pred)
    xlabel('Direction Content')
    ylabel('Spacing');
    title('Carrier motion induced bias (response rates where delta-x = 0)')
    colorbar;
end

    tile()
%let's show contours of Dx versis spacing for each subject.