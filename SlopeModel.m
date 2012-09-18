classdef SlopeModel < LikelihoodModel
% this is PBM's first attempt at a model.

    methods
        function SM = SlopeModel()
            SM.initialParamsDefaults = dataset(struct(...
                'mu_0', 0, ...
                'beta_0', 0, ...
                'beta_small', 0, ...
                'cs', 0, ...
                'beta_summation', 0, ...  %or BOTH these parameters
                'beta_induced', 0, ...    %(not all three, not another combo)
                'sbeta_summation', 0, ...
                'sbeta_induced', 0, ...
                'content_compression', 5 ...
                ));

            SM.freeParams = {'mu_0', 'beta_0', 'cs', 'sbeta_summation', 'sbeta_induced'};
        end

        function prob = predict(SM, p, data)
            if ~exist('data', 'var') || isempty(data)
                data = SM.data;
            end
            if ~exist('p', 'var') || isempty(p)
                if isempty(SM.parameters)
                    error('bModel:needParams', ...
                          'Need parameters to predict responses');
                else
                    p = SM.parameters;
                end
            end

            %the sensitivity to is controlled by the critical spacing.
            sens = (2 - 2./(1+exp(-p.cs./data.spacing)));

            %alternately there is the degree of "summation"
            summation = 1./data.spacing;

            scaled_content = tanh(data.content .* p.content_compression);

            bias = p.mu_0 ...
                   + p.beta_induced.*data.content ...
                   + p.beta_summation.*summation.*data.content ...
                   + p.sbeta_induced.*scaled_content ...
                   + p.sbeta_summation.*summation.*scaled_content;


            prelink = bias + p.beta_0.*data.dx.*sens + p.beta_small.*data.dx.*(1-sens);

            logit = @(x)0.98./(1+exp(-x)) + 0.01;
            prob = logit(prelink);
        end
    end
end