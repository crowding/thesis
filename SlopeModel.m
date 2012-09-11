classdef SlopeModel < LikelihoodModel
% this is PBM's first attempt at a model.

    methods
        function SM = SlopeModel()
            SM.freeParams = {'ke', 'ks', 'kdx', 'cs'};

            args = cellfun(@(x) {0 x}, SM.freeParams, 'UniformOutput', 0);
            SM.initialParamsDefaults = dataset(args{:});
            SM.initialParamsDefaults.cs = 4;
        end

        function prob = predict(SM, p, data)
            if ~exist('data', 'var') || isempty(data)
                data = BM.data;
            end
            if ~exist('p', 'var') || isempty(p)
                if isempty(SM.parameters)
                    error('bModel:needParams', ...
                          'Need parameters to predict responses');
                else
                    p = SM.parameters;
                end
            end

            [s, c, dx] = deal(data.spacing, data.content, data.dx);

            logit = @(x)1./(1+exp(-x));

            %more motion energy ("e") affects bias
            r = p.ke .* c ./ s;

            %more motion inside "critical distance" (ks for summation) affects bias
            r = r + p.ks .* c .* (1 + p.cs ./ s);

            %more motion inside "critical spacing" affects sensitivity to dx
            r = r + dx .* p.kdx ./ (1 + p.cs ./ s);

            prob = logit(r);
        end
    end
end