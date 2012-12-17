function a = dsvertcat(varargin)

data = varargin;
varnames = cellfun(@(x) get(x, 'VarNames'), data, 'UniformOutput', 0);

cellfun(@check, unique(cat(1, varnames{:})))
function check(varname)
    structful = cellfun(@(x) iscell(x.(varname)), data);
    if ~all(structful) & any(structful)
        for i = 1:numel(data)
            if ~structful(i)
                data{i}.(varname) = num2cell(data{i}.(varname));
            end
        end
    end
end

a = vertcat(data{:});
end
