function s = rename(s, varargin)
% function s = rename(s, oldName, newName, oldName, newName...)
% function s = rename(s, {oldName, newName, oldName, newName...})
% function s = rename(s, struct(oldName, newName, oldname, newname...))
%
% Replace a struct's or dataset's field names with other names.

if isstruct(varargin{1})
    names = fieldnames(varargin{1});
    values = struct2cell(varargin{1}, 1);
elseif iscell(varargin{1})
    names = varargin{1}(1:2:end);
    newnames = varargin{1}(2:2:end)
else
    names = varargin(1:2:end);
    newnames = varargin(2:2:end);
end

if (numel(names) ~= numel(newnames))
    error('rename:badArguments', 'number of names and new names must match.');
end
    
if isa(s, 'dataset')
    vals = s(:,names);
    s(:,names) = [];
    s(:,newnames) = vals(:,names);
else
    for i = [names(:) newnames(:)]';
        [name, newname] = i{:};
        [value{1:numel(s)}] = s.(name);
        s = rmfield(s, name);
        [s.(newname)] = value{:};
    end
end