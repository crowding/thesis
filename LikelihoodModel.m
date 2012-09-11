classdef LikelihoodModel
    properties
        parameters = dataset();
        data = dataset();

        freeParams = {};

        splits = {'subject', 'exp_type'}

        initialParamsTable = dataset();
        initialParamsDefaults = dataset();
    end

    methods
        function M = LikelihoodModel()
        end

        function params = initialParams(M, data)
            conditions = unique(data(:,M.splits));

            if ~isempty(M.initialParamsTable)
                [params, ~, ib] = join(M.initialParamsTable, conditions, ...
                                       'type', 'inner', 'MergeKeys', true);
            else
                params = dataset();
                ib = [];
            end

            %any conditions that weren't matched get default...
            %Ugh, doing it thisway was more trouble than it was worth.
            unmatched = conditions(~ismember(1:size(conditions, 1), ib), :);
            unmatched.dsfargeg = ones(size(unmatched, 1), 1);
            param_defaults = M.initialParamsDefaults;
            param_defaults.dsfargeg = ones(size(param_defaults, 1), 1);

            defaulted = join(param_defaults, unmatched, ...
                             'type', 'inner', 'MergeKeys', true);
            defaulted(:,'dsfargeg') = [];

            params = [params; defaulted];
        end

        function M = fit(M, data)
            M.data = data;
            if isempty(M.parameters)
                params = M.initialParams(data);
                M.parameters = params;
            end

            params = groupfun(data, M.splits, @makeFit);

            function [f, err] = makeFit(chunk)
                split = chunk(1,M.splits)
                chunkParams = join(params, split, ...
                                   'type', 'inner', 'MergeKeys', true);
                [f, f.err] = fit(@M.negLogLikelihood, chunkParams, ...
                                 M.freeParams, chunk);
                f.n_obs = size(chunk, 1);
            end
            M.parameters = params;
        end

        function err = negLogLikelihood(M, params, data)
            if ~exist('data', 'var') || isempty(data)
                data = M.data;
            end
            if ~exist('params', 'var') || isempty(params)
                params=M.params;
            end
            prob = M.predict(params, data);

            prob = prob*.99+.005;
            err = -sum(  data.response.*log(prob) ...
                       + (1-data.response).*log(1-prob));
        end

        function p = fullPredict(M, params, data)
            if ~exist('data', 'var') || isempty(data)
                data = M.data;
            end
            if ~exist('params', 'var') || isempty(params)
                params=M.parameters;
            end
            %use the model's predict method but first join the parameter table
            %to account for varying parameters.
            param = join(data, params, 'type', 'inner', 'MergeKeys', true);
            p = M.predict(param, data);
        end

        function r = residuals(M, split, binvar, binsize)
        %function r = residuals(model, splitvars, binvar, binsize)
        % Having already fit a model, compute the resudials, marginalizing
        % over some variables and binning over another. Right now I compute
        % Pearson residuals, should add deviance values also.


        % The residual is the number of subject responses, minus the number of
        % target responses.  For a binomial response variable that's almost
        % useless, so, we'll marginalize over different parameters and plot
        % the residuals of the sums, sort of group wise residuals.

        % make the predictions (the model is fit already)
            d = M.data;
            d.p_pred = M.fullPredict();

            %per point residuals and variance of the predicted mean
            d.resid = logical(d.response) - d.p_pred;
            d.pred_var = (d.p_pred).*(1-d.p_pred);

            %Most logistic regression examples use toy problems where
            %"the mean of observations at some x value" makes sense.
            %Gelman and Hill(2007) is a great book because it uses
            %real data instead, like binary observations where the
            %abscissa is not binned (the water well dataset.) They
            %suggest the expedient of binning data across a
            %variable. There is also a library written in R to do this
            %sort of thing for you, and all this code here is for if
            %you don't use R.

            %this bin-splitter function is pretty damn slow. Hmm.
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
                d.pearson_resid = d.total_resid ./ sqrt(sum(x.pred_var));

                % the deviance residual is basically the difference between "expected"
                % log likelihood and observed log likelihood. It should follow
                % a chi-square distribution, so we can convert it to p-values.
                %
                % d.expected_log_likelihood = d.  d.likelihood_deviance =
                %d.deviance = -2.*(sum(log(x.response.*x.p_pred ...
                %                          + ~x.response.*(1-x.p_pred))) ...
                %                  - sum(log(x.response*mean(x.response) ...
                %                            + ~x.response*(1-mean(x.response)))));
            end
        end
    end
end
