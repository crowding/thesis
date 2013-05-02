function plotter = connectLines(orderBy, groupBy, varargin)
%function plotter = connectLines(orderBy, groupBy)
%Connect the lines in teh current subplot according to group.
%(think about supporting linetypes, etc...?

    plotter = @doConnectLines;
    function handles = doConnectLines(handles, chunk)
    %this is getting called once per subplot
        hold(handles.ax, 'on')
        linehandles = groupfun(chunk, groupBy, @addLine);
        function lh = addLine(lineChunk)
            [~,order] = sortrows(lineChunk(:,orderBy));
            lh = plot(handles.ax, lineChunk.x_(order), lineChunk.y_(order), ...
                      'LineStyle', '-', ...
                      'Color', lineChunk.color_(order(1),:), ...
                      varargin{:});
        end
        handles.linehandles = {linehandles.V}
        hold(handles.ax, 'off')
    end
end