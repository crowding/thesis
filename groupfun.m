function varargout = groupfun(data, groupvars, fun, varargin)
%function varargout = groupfun(data, groupvars, fun, varargin)
%
%Split a dataset up into groups, apply a function to each group, and
%concatenate the results into a dataset.
%
%Example:
%
% >> load fisheriris
% >> iris = dataset({nominal(species),'species'},...
%                {meas,'SepalLength','SepalWidth','PetalLength','PetalWidth'}...
%                );
%
% >> MeanAspects = ...
%     groupfun(iris, 'species',  @(x)struct( ...
%                       'SepalAspect', mean(x.SepalLength./x.SepalWidth), ...
%                       'PetalAspect', mean(x.PetalLength./x.PetalWidth)));
%
% MeanAspects =
%
%     SepalAspect    PetalAspect    species
%     1.4702          6.908         setosa
%     2.1604         3.2428         versicolor
%     2.2305         2.7807         virginica

%You'd think 'grpstats' would do this kind of thing. But it only wants
%to compute column-wise statistics! Who ever computes column-wise
%statistics?! Different columns mean different things and somtimes you
%want statistics that pull together several columns.

if ~all(ismember(groupvars, get(data, 'VarNames')))
    error('groupfun:noSuchVars',...
          ['Can''t split on variables that aren''t in the dataset']);
end

[~, firstixs, ix] = unique(data(:,groupvars));
splits = separate(ix);
[varargout{1:nargout}] = deal([]);

[varargout{1:nargout}] = cellfun(@map_split, splits, 'UniformOutput', 0);
function varargout = map_split(ix)
    [varargout{1:nargout}] = fun(data(ix,:), varargin{:});
end
varargout = cellfun( @mkdatasets, varargout, 'UniformOutput', 0);
function Vs=mkdatasets(Vs)
    Vs = cellfun(@mkdataset, Vs, 'UniformOutput', 0);
    function V = mkdataset(V)
        if ~isa(V, 'dataset')
            V = dataset(V);
        end
    end
end

varargout = cellfun(@fillout, varargout, 'UniformOutput', 0);
function Vs=fillout(Vs)
    Vs = cellfun(@fillgrp,Vs,num2cell(firstixs),'UniformOutput', 0);
    function V = fillgrp(V,firstix)
        fill = setdiff(groupvars, get(V, 'VarNames'));
        if ~isempty(fill)
            V(:,fill) = data(repmat(firstix, size(V,1), 1),fill);
        end
    end
end

varargout = cellfun(@rbind, varargout, 'UniformOutput', 0);
function V = rbind(V)
    V = cat(1, V{:});
end

end
