function resid = modelDeviance(data, margin, model, params)

% There are some aspects of the data that I think the model does not
% capture well, and I want to show it.

% What I'd like are residual or deviance plots.

% Side note: you know, this is going to take renamed data because the
% model takes renamed data. Perhaps I should make "folding" or
% "unfolding" a model paramter.

% these are the variables the model operates over.
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

%if ~isa(params, 'dataset')
%    params = dataset(params);
%end

% make the predictions (the model is fit already)
data.p_pred = model(params, data);

% The residual is the number of subject responses, minus the number of
% target responses.  For a binomial response variable that's almost
% useless, so, we'll marginalize over different parameters and plot
% the residuals of the sums.

if ~exist('margin', 'var') || margin == 0
    margin = {'spacing', 'subject', 'response', 'exp_type'};
end
% WIth predictions in hand, we compute expected number of yesses,
% and residual deviance, groupwise, over some marginal variables.
function resid = residuals(margin)
    resid = groupfun(data, setdiff(vars, margin), ...
         @(x) struct('n_pred', sum(x.p_pred), ...
                     'n_yes', sum(x.response), ...
                     'n_total', size(x,1), ...
                     'n_pred_sd', sqrt(sum(x.p_pred.*(1-x.p_pred))), ...
                     'deviance', 2*sum(log(x.response.*x.p_pred ...
                                           + (1-x.response).*(1-x.p_pred)))));
    %TODO: take the deviance and n_samples and compute a chi-squared test
    %statistic.
end

%and here's a function for plotting the residuals over some variables.
    function handles = plotResid(data, figVar, xVar, rowVar, colVar, colorVar)
    %use 'const_' to not marginalize over that var.
    data.const_ = repmat(1, size(data, 1));

    %need to set the xlim globally.
    xlim = [min(data.(xVar)) max(data.(xVar))];
    ylim = [min(data.(yVar)) max(data.(yVar))];

    [figs, ~, data.fighandle_] = arrayfun(@figure, unique(data.(figVar))) + 97;
    [rows, ~, data.figRow_] = unique(data.(rowVar));
    [color, ~, data.figCol_] = unique(data.(colVar));
    [~, ~, data.color_] = unique(data.(colorVar));
    colormap = cool(length(color));
    data.color_ = colormap(data.color_, :);
    nrow = length(rows);
    ncol = length(cols);

    handles = groupfun(data, {figVar rowVar colVar}, @subplot);
    function handles = subplot(chunk)
        figure(data.fighandle_);
        a = subplotix(length(rows), length(cols), data.figRow_, data.figCol_);
        gscatter(data.(xVar), data.n_yes - data.n_pred, ...
                 1:size(data,1), data.color_, '.', sqrt(data.n_obs) * 5);
        h = data.(xVar), data.n_yes - data.n_pred
        handles = struct('figure', {fh}, 'axis', {a}, 'handles', {h});

        %TODO: enable axes, legend, h/v based on which subplot.
    end
end

%now plot this data over these marginals
resid = residuals(margin);

% plotResid(resid, 'subject', 'dx', 'const_', 'const_', 'const_');

end