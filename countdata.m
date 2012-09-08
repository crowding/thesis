function count = countdata(data, split)
%This is the matlab translation of this R code.
%
%count = ddply(data, split, function(x) c(n=nrow(x), n.yes = sum(x.response), p = mean(x.response)

%this is several times faster than the shorter version (countdata2), but way too
%much code.

%Count the responses and probability per condition. Because doing it
%straightforwardly like you would with ddply in plyr/R just takes too
%long with MATLAB's poky dataset arrays.

    [values, ~, ix] = datasetfun(@unique, data(:,split), 'UniformOutput', 0);
    [cases, ~, caseix] = unique(data(:,split));
    sz = cellfun('prodofsize', values);

    coverage = size(cases, 1) / prod(sz);

    %accumulating in sparse arrays will be the order of the day. So
    %collapse everything except the row index.
    if coverage < 0.1
        sparse = true;
    else
        sparse = false;
    end

    if length(values) > 2
        %OH FOR FUCK'S SAKE, sparse doesn't do >2-d arrays. So we
        %have to pack this into a 2-d array.
        pack_ix = [ix{1} sub2ind(sz(2:end), ix{2:end})];
        pack_sz = [sz(1) prod(sz(2:end))];
    else
        pack_ix = [ix{:}];
        pack_sz = sz;
    end

    n_yes = accumarray(pack_ix, double(logical(data.response)), pack_sz,...
                       [], 0, sparse);
    n_no = accumarray(pack_ix, double(~data.response), pack_sz,...
                      [], 0, sparse);
    n = n_yes + n_no;

    %now unpack n and n_yes back into a dataset.

    [i,j, n] = find(n);
    n_yes = n_yes(sub2ind(size(n_yes), i,j));

    if length(values) > 2
        [jj{1:numel(sz)-1}] = ind2sub(sz(2:end), j);
        sub = {i jj{:}};
    else
        sub = {i j};
        sub = sub(1:length(sz));
    end

    groupdata = cellfun(@(V, I)V(I), values, sub, 'UniformOutput', 0);
    ds_args = cellfun(@(data, label){data label}, groupdata, split, 'UniformOutput', 0);
    ds_args = [ds_args {{full(n) 'n'}, {full(n_yes ./ n) 'p'}}];

    count = dataset(ds_args{:});
end