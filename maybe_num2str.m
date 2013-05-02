function x = maybe_num2str(x)
    if (iscell(x))
        if (numel(x) > 1)
            x = cellfun(@maybe_num2str, x, 'UniformOutput', 0);
        end
    elseif (isnumeric(x))
        if numel(x) > 1
            x = arrayfun(@maybe_num2str, x, 'UniformOutput', 0);
        else
            x = num2str(x);
        end
    end
end
