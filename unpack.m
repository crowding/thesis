function unpack(varargin)
% loads all the data into structs.
% output variables are S, a struct containing all the experiment
% types, and 'data', in which experiment type is another column.

if isempty(varargin)
    d = dir('*_trials.csv');
    filenames = {d.name};
else
    filenames = varargin;
end

data = cellfun(@load_data_set, filenames, 'UniformOutput', 0);
function d = load_data_set(filename)
    groupname = regexp(filename, '([a-zA-Z])*_series_trials.*', 'tokens');
    d = dataset(read_csv(filename));
    d.exp_type = repmat(groupname{1}, size(d,1),1);
    %This flip is because the knob recorded values with a flipped mapping.
    d.folded_response_with_carrier = ~d.folded_response_with_carrier;
    d.abs_response_cw = ~d.abs_response_cw;
end

data = dsvertcat(1, data{:});

save(fullfile(fileparts(mfilename('fullpath')), 'data.mat'), 'data')

end
