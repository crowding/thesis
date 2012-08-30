function UnpackPeter(varargin)
% loads all the data into structs.
% output variables are S, a struct containing all the experiment
% types, and 'data', in which experiment type is another column.

if isempty(varargin)
    d = dir('*_trials.csv');
    filenames = {d.name};
else
    filenames = varargin;
end

S = struct();

cellfun(@loadFile, filenames)
function loadFile(filename)
    groupname = regexp(filename, '([a-zA-Z])*_series_trials.*', 'tokens');
    S.([groupname{1}{1} 'Data']) = read_csv(filename);
end

fields = fieldnames(S.contentData);
S.allData = cellfun(@(x) struct2cell(orderfields(x, fields)), ...
                  {S.contentData, S.spacingData}, 'UniformOutput', 0);
S.allData = cell2struct(cellfun(@vertcat, S.allData{:}, 'UniformOutput', 0)...
                        , fields, 1);



%or, with datasets
data = cellfun(@load_data_set, filenames, 'UniformOutput', 0);
function d = load_data_set(filename)
    groupname = regexp(filename, '([a-zA-Z])*_series_trials.*', 'tokens');
    d = dataset(read_csv(filename));
    d.exp_type = repmat(groupname{1}, size(d,1),1);
end

S.data = cat(1, data{:});

save('data.mat', '-struct', 'S')

end
