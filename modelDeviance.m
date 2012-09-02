function [resid, handles] = ...
        modelDeviance(data, split, model, params, binvar, binsize)

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

if ~exist('binsize', 'var')
    binsize = 25;
end
if ~exist('binvar', 'var')
    binvar = 'dx';
end

% make the predictions (the model is fit already)
data.p_pred = model(params, data);

% The residual is the number of subject responses, minus the number of
% target responses.  For a binomial response variable that's almost
% useless, so, we'll marginalize over different parameters and plot
% the residuals of the sums.


if ~exist('split', 'var') || split == 0
    split = {'subject', 'content', 'spacing'};
end

resid = residuals(data, split, binvar, binsize);

handles = plotResid(resid, 'const_', 'dx', 'pearson_resid', ...
                    'content', 'const_', 'spacing');
tile(1,1,[],unique(handles.fig));

%now plot this data over these marginals

% Computing residuals.
% WIth predictions in hand, we compute expected number of yesses,
% and residual deviance, groupwise, over some marginal variables.
% THe pearson residu
function r = residuals(d, split, binvar, bin_size)

    %per point residuals and variance of the predicted mean
    d.resid = logical(d.response) - d.p_pred;
    d.pred_var = (d.p_pred).*(1-d.p_pred);

    %Most logistic regression examples use toy problems where "the mean of
    %observations as some x value)" makes sense.  Gelman and
    %Hill(2007) is a great book because it uses real data instead,
    %like binary observations where the abscissa is not binned (the
    %water well dataset.) They suggest the expedient of binning data
    %across. There is also a library written in R to do this sort of
    %thing for you, and all this code here is for if you don't use R.
    %
    %Most statisticians use R today, so you should learn it. In my
    %experience, working with datasets and making plots, a line of R
    %is worth, on average, two or three lines of MATLAB. Remember
    %Wichmann and Hill, who wrote about out how best to fit
    %psychometric functions, and produced the psignifit library? They
    %don't maintain psignifit for MATLAB any more and they publish
    %code in R now. As an example.

    d = groupfun(d, split, @bin);
    function x = bin(x)
        [~, ord] = sort(x.(binvar));
        x = x(ord,:);
        x.bin__ = ceil((1:size(x,1))'./binsize);
    end

    %Now sum and compute summary residuals over bins and splits
    r = groupfun(d, union(split, 'bin__'), @group_residual);
    function d = group_residual(x)
        d.(binvar) = mean(x.(binvar));
        d.n_pred = sum(x.p_pred);
        d.n_yes = sum(x.response);
        d.n_obs = size(x,1);
        d.n_pred_sd = sqrt(sum(x.pred_var));
        d.total_resid = sum(x.resid);
        d.pearson_resid = sum(x.resid) / sqrt(sum(x.pred_var));
        %something about the log-likelihood residual deviance, but I don't
        %understand that stat yet enough to say what it is when you
        %bin over observations.

        %d.deviance = -2.*(sum(log(x.response.*x.p_pred ...
        %                          + ~x.response.*(1-x.p_pred))) ...
        %                  - sum(log(x.response*mean(x.response) ...
        %                            + ~x.response*(1-mean(x.response)))));
    end
end

%and here's a function for plotting the residuals over some
%variables. Returns a dataset of handles to axes, linesseries,
%figures, etc.
function handles = plotResid(data, figVar, xVar, yVar, rowVar, colVar, colorVar)
    %use 'const_' to not marginalize over that var.
    data.const_ = ones(size(data, 1),1);

    %need to set the xlim globally.
    xlim = [min(data.(xVar)) max(data.(xVar))];
    ylim = [min(data.(yVar)) max(data.(yVar))];

    [~, figno, figix] = unique(data.(figVar), 'first');
    fighandles = arrayfun(@figure, 1:numel(figno));
    data.fighandle_ = fighandles(figix);

    [rowVals, ~, data.figRow_] = unique(data.(rowVar));
    [colVals, ~, data.figCol_] = unique(data.(colVar));
    [colorVals, ~, data.color_] = unique(data.(colorVar));
    colormap = cool(length(colorVals));
    data.color_ = colormap(data.color_, :);
    nrow = length(rowVals);
    ncol = length(colVals);

    handles = groupfun(data, {figVar rowVar colVar}, @subplot);
    function handles = subplot(chunk)
        handles.fig = chunk.fighandle_(1);
        figure(handles.fig);
        handles.ax = subplotix(length(rowVals), length(colVals), ...
                               chunk.figRow_(1), chunk.figCol_(1));
        handles.handle = {...
            gscatter(chunk.(xVar), chunk.(yVar), ...
                     1:size(chunk,1), chunk.color_, '.', ...
                     sqrt(chunk.n_obs) * 5, 0)};
        set(handles.ax, 'XLim', xlim, 'YLim', ylim);
        %enable axes, legend, h/v based on which subplot.
        if chunk.figRow_(1) == numel(rowVals)
            xlabel(handles.ax, ...
                   sprintf('%s\n(%s = %s)', ...
                           xVar, colVar, num2str(chunk.(colVar)(1),3)));
        else
            set(handles.ax, 'XTickLabel', []);
        end

        if chunk.figCol_(1) == 1
            ylabel(handles.ax, ...
                   sprintf('(%s = %s)\n%s', ...
                           rowVar, num2str(chunk.(rowVar)(1),3), yVar));
        else
            set(handles.ax, 'YTickLabel', []);
        end

        if chunk.figRow_(1) == 1 && ~strcmp(colorVar, 'const_')
            %in lieu of placing a legend just yet
            title(handles.ax, sprintf('(colors indicate %s)', colorVar));
        end
    end
end



end