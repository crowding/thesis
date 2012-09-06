function f = plotFits(fitfile, datafile, outfile)

    if ~exist('fitfile', 'var') || ~isempty(fitfile)
       fitfile = fullfile(fileparts(mfilename), 'fits.mat');
    end

    if ~exist('datafile', 'var') || ~isempty(datafile)
       datafile = fullfile(fileparts(mfilename), 'data.mat');
    end

    load(fitfile, 'fits');
    load(datafile, 'data');
    data = doRename(data);
    split = {'subject', 'dx', 'content', 'spacing', 'exp_type'};

    %break down the trials in terms of number yes and number no
    rates = countdata(data, split);

    %plot only data that matches one of our fits.
    [~, ia, ib] = join(data, fits, 'type', 'inner');
    data = data(ia,:);
    fits = fits(ib,:);

    %make a lattice plot for all the data. It says plotResiduals but it's
    %just a scatter plot.
    handles = plotResiduals(rates, 'subject', 'dx', 'p', 'exp_type', ...
                            'content', 'spacing', 'n');

    %evaluate model predictions
    %dx_eval = linspace(min(data.dx), max(data.dx), 100);
    %dx_eval = dataset({linspace(min(data.dx), max(data.dx), 100), 'dx'} ...
    %dx_eval = {ones(100, 1), 'dx'};
    %preds = groupfun(rates, setdiff(split, 'dx'), @(chunk)addPredictions)
    function pred = addPredictions(chunk)
        chunk = chunk{1};
        %
        %and add these predictions to the graph.
        fit = join(params, chunk(1,:), 'type', 'inner', 'MergeKeys', 'true');
        graphics = join(handles, chunk(1,:), 'type', 'inner', 'MergeKeys', 'true');

        chunk.dummy_ = ones(size(chunk, 1),1);

        pred = MotionModel(chunk, params)
    end
    
    %Plot on top of that the model fits. which involves making
    %predictions. I wish I had ggplot2 for this.

end