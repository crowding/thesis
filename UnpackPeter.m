function UnpackPeter(varargin)

if isempty(varargin)
    d = dir('*_trials.csv');
    varargin = {d.name};
end

S = struct();

cellfun(@loadFile, varargin)
function loadFile(filename)
    groupname = regexp(filename, '([a-zA-Z])*_series_trials.*', 'tokens');
    S.([groupname{1}{1} 'Data']) = read_csv(filename);
end

fields = fieldnames(S.contentData);
S.allData = cellfun(@(x) struct2cell(orderfields(x, fields)), ...
                  {S.contentData, S.spacingData}, 'UniformOutput', 0);
S.allData = cell2struct(cellfun(@vertcat, S.allData{:}, 'UniformOutput', 0), fields, 1);

save('data.mat', '-struct', 'S')
end
