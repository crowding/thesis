function r = residuals(d, model, params, split, binvar, binsize)
    %function r = residuals(data, model, params, splitvars, binvar, binsize)
    % Having already fit a model, compute the resudials, marginalizing
    % over some variables and binning over another. Right now I compute
    % Pearson residuals, should add deviance values also.


    % The residual is the number of subject responses, minus the number of
    % target responses.  For a binomial response variable that's almost
    % useless, so, we'll marginalize over different parameters and plot
    % the residuals of the sums, sort of group wise residuals.

    % make the predictions (the model is fit already)
    d.p_pred = model(params, d);

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
        d.pearson_resid = d. / sqrt(sum(x.pred_var));

        % the deviance residual is basically the difference between "expected"
        % log likelihood and observed log likelihood. It should follow
        % a chi-quoare distribution, so we can convert it to p-values.
        %
        % d.expected_log_likelihood = d.  d.likelihood_deviance =
        % d.deviance_x2 =
        %something about the log-likelihood residual deviance, but I don't
        %understand that stat yet enough to say what it is when you
        %bin over observations.

        %d.deviance = -2.*(sum(log(x.response.*x.p_pred ...
        %                          + ~x.response.*(1-x.p_pred))) ...
        %                  - sum(log(x.response*mean(x.response) ...
        %                            + ~x.response*(1-mean(x.response)))));
    end
end