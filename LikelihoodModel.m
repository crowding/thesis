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

            [params, ~, ib] = join(M.initialParamsTable, conditions, ...
                                   'type', 'inner', 'MergeKeys', true);

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
                chunkParams = join(params, split, 'type', 'inner', 'MergeKeys', true);
                [f, f.err] = fit(@M.negLogLikelihood, chunkParams, M.freeParams, chunk);
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
            err = -sum(data.response.*log(prob) + (1-data.response).*log(1-prob));
        end
    end
end
