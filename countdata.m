function count = countdata(data, split)

%This is the matlab translation of this R code.
%
%count = ddply(data, split, function(x) c(n=nrow(x), n.yes = sum(x.response), p = mean(x.response)

    count = grpstats(data, split, @(x)mean(logical(x)), 'DataVars', 'response');
    count = rename(count, 'GroupCount', 'n', 'Fun1_response', 'p');
end
