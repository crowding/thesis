classdef SlopeModel < LikelihoodModel
    methods
        function SM = SlopeModel(varargin)
            SM = SM@LikelihoodModel(varargin{:});

            SM.initialParamsDefaults = dataset(struct(...
                'mu_0', 0, ...
                'beta_0', 8, ...
                'beta_small', 0, ...
                'cs', 4, ...
                'beta_summation', 4, ...
                'beta_induced', 0, ...
                'saturating_induced', 0, ...
                'wiggle_induced', 0, ...
                'induced_scale', 6 ...
                ));

            if isempty(SM.freeParams)
                SM.freeParams = {'mu_0', 'beta_0', 'cs', 'saturating_induced', 'wiggle_induced'};
            end
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

            %the "induced motion* comes in two types, which I'll fit with a
            %logistic plus the third derivative of a logistic (no real
            %justification here other than hte combination looks loke
            %the data.) The scale parameter of the logit is chosen
            %arbitrarily, not fit.
            logit = @(x) 1./(exp(x)+1);
            logit3 = @(x) exp(x)*(exp(x)-1)./(exp(x)+1).^3;
            induced = (p.saturating_induced .* logit(p.induced_scale.*data.content) ...
                       + p.wiggle_induced .* logit(p.induced_scale.*data.content));

            bias = p.mu_0 ...
                   + p.beta_induced.*data.content ...
                   + p.beta_summation.*summation.*data.content ...
                   + induced;

            prelink = bias + p.beta_0.*data.dx.*sens + p.beta_small.*data.dx.*(1-sens);

            logit = @(x)0.98./(1+exp(-x)) + 0.01;
            prob = logit(prelink);
        end
    end
end