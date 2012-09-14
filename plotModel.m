function handles = plotModel(model, modelName, split)

if ~exist('split', 'var')
    split = {'subject', 'dx', 'content', 'spacing', 'exp_type'};
end
%break down the trials in terms of number yes and number no per unique
%condition
rates = countdata(model.data, split);

%make a lattice plot for all the data. It says plotResiduals but it's
%just a scatter plot.
%                              figure,    x     y    row
handles = plotResiduals(rates, 'subject', 'dx', 'p', 'exp_type', ...
                        'content', 'spacing', 'n', @addPredictions);
%                       col        color      size

%for each subplot, plot over dx and add handles.
function [handle] = addPredictions(handle, chunk)
%make prediction lines on this subplot based on the data
    chunk.dx = [];
    chunk.n = [];
    chunk.p = [];
    chunk = unique(chunk);
    chunk.line = (1:size(chunk, 1))';

    dx_eval = dataset(...
        {linspace(chunk.xmin_(1), chunk.xmax_(1), 100)', 'dx'},...
        {ones(100, 1), 'const_'});

    pred = groupfun(chunk, 'line', @predictLine);
    function line = predictLine(line)
        line = join(line, dx_eval, 'type', 'inner', 'MergeKeys', true);
        line.p = model.fullPredict([], line);
    end
    %and add these predictions to the graph. They are keyed by color.
    hold(handle.ax(1), 'on');
    groupfun(pred, 'line', @plotLine);
    %this doesn't return the actual graphics handles corresponding to
    %the predicted line though.
    function plotLine(line)
        [x, o] = sort(line.dx);
        y = line.p(o);
        plot(handle.ax(1), x, y, 'color', line.color_(1,:));
    end
    hold(handle.ax(1), 'off');
end

end
