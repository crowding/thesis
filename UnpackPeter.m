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
    varname = [groupname{1}{1} 'Data'];
    S.(varname) = read_csv(filename);
    %since we end up wanting to flip the response, might as well do that
    %here for all candidate models. THe actual problem is that
    %displacement is measured clockwise, while carrier is measured
    %anticlockwise. THat's an old consequence of stuff I did early in
    %writing the graphics :/
    S.(varname).folded_displacement = -S.(varname).folded_displacement;
    S.(varname).abs_displacement = -S.(varname).abs_displacement;
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
    d.folded_displacement = -d.folded_displacement;
    d.abs_displacement = -d.abs_displacement;
end

S.data = cat(1, data{:});

save('data.mat', '-struct', 'S')

end
