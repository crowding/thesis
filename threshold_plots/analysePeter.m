function analysePeter(varargin)
    
clf;

if nargin == 0
    filenames = {
          'collections/spacing_series_calculations.csv' ...
        , 'collections/content_series_calculations.csv', ...
        };
else
    filenames = varargin;
end

datasets = cellfun(@read_csv, filenames, 'UniformOutput', 0);
datasets = cellfun(@dataset, datasets, 'UniformOutput', 0);

all_data = cat(1, datasets{:});

for i = 1:numel(datasets)
    subplot(1, numel(datasets), i);
    plotData(datasets{i});
end

print('ione/analysePeter.fig')

end

function plotData(data)

sublist = unique(data.subject);
dirlist = unique(data.folded_direction_content);
spacinglist = unique(data.target_spacing);

colorlist = hot(length(dirlist)+3);

symlist={'o', 's' '*' 'x' 'd' '^' 'v' '>' '<' 'p' 'h' 'v' '>' '<' 'p' 'h'};

for s=1:length(sublist)
    for d=1:length(dirlist)
        for sp=1:length(spacinglist)
            ind=find(   data.folded_direction_content == dirlist(d) ...
                      & data.target_spacing == spacinglist(sp) ...
                      & strcmp(data.subject, sublist{s}));
            if ~isempty(ind)
                plot3(  dirlist(d), spacinglist(sp), data.bias(ind), [symlist{s}]...
                      , 'MarkerSize', 5, 'Color', colorlist(d, :) ...
                      , 'MarkerFaceColor', colorlist(d, :));
                hold on
            end
        end
    end
end

line([0 0],[0 25], [ 0 0 ])
xlabel('carrier direction content')
ylabel('envelope spacing')
zlabel('bias (log-odds)')
set(gca, 'XLim', [0 1])
set(gca, 'YLim', [0 25])
set(gca, 'ZLim',[-4 4])
%     text(.5, 20, -2, 'local motion direction biases percept')
%     text(.5, 20, 2, 'induced observer motion')
%     text(.5, 20 ,0, 'no effect of local motion on percept')
end