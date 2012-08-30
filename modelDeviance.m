function residuals = modelDeviance(model, params, data)

% There are some aspects of the data that I think the model does not
% capture well, and I want to show it.

% What I'd like are residual or deviance plots.

% Side note: you know, this is going to take renamed data because the
% model takes renamed data. Perhaps I should make "folding" or
% "unfolding" a model paramter.

% these are the variables the model operates over.
vars = {'subject', 'content', 'spacing', 'dx', 'exp_type'};

if ~exist('data', 'var')
    load('data.mat', 'data')
    data = rename(data, ...
                  'target_spacing', 'spacing', ...
                  'folded_direction_content', 'content', ...
                  'folded_displacement', 'dx', ...
                  'folded_response_with_carrier', 'response');
    data = data( strcmp(data.exp_type, 'spacing') ...
                & strcmp(data.subject, 'pbm'), vars);
end

if ~exist('params', 'var')
    load('modelResults.mat', 'p')
    params = p;
end

if ~exist('model', 'var')
    model = @MotionModel;
end

if ~isa(data, 'dataset')
    data = dataset(data);
end

params = dataset(params);

% make the predictions (the model is fit already)
data.p_pred = model(params, data);

% the residual is the number of subject responses, minus the number of
% target responses.
% We'll try different margins along which to collapse the data.

% I have a suspision that there is a motion-energy effect that is
% independent of spacing and not captured.

margin = {'spacing', 'subject'};

% could also make standard error of prediction (neglecting error of
% model fit) using arcsine transform.
residuals = grpstats(data, setdiff(vars, margin), ...
                     @(ds) dataset({length(ds.response), 'n_obs', },...
                                   {sum(ds.response),    'ncw_obs'}, ...
                                   {sum(ds.pred),        'n_pred'} ));




end