function d = groupfun(data, groupvars, fun, varargin)
%Split a dataset up into groups, apply a function to each chunk, and
%concatenate the results.

%You'd think 'grpstats' would do this kind of thing. But it only wants
%to compute column-wise statistics! Who ever computes column-wise
%statistics?! Different columns mean different things!
if ~all(ismember(groupvars, get(data, 'VarNames')))
    error('groupfun:noSuchVars',...
          ['Can''t split on variables that aren''t in the dataset']);
end

[~, ~, ix] = unique(data(:,groupvars));
splits = separate(ix);

mapped_splits = cellfun(@map_split, splits, 'UniformOutput', 0);
function V = map_split(ix)
    V = fun(data(ix,:), varargin{:});
    if ~isa(V, 'dataset')
        V = dataset(V);
    end
    fill = setdiff(groupvars, get(V, 'VarNames'));
    if ~isempty(fill)
        V(:,fill) = data(repmat(ix(1), size(V,1), 1),fill);
    end
end

d = cat(1, mapped_splits{:});

end
