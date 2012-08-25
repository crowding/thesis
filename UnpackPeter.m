% UnpackPeter.m

% Here's CSV files with individual trial data.
% 
% I'm keeping it in two files, 
% "content_series_trials": 
% more values of direction content and fewer values of spacing were used in
% each session

% "spacing_series_trials" 
%  many values of spacing were used but only one value of content

% Subjects seem to adapt to the strength of direction content, so you will
% see a bigger effect of the magnitude of direction content with the first
% dataset than the second.

% "eccentricity" the eccentricity of the spots 

% "abs_direction_content", The strength of the carrier motion, ranges from
% -1 (counterclockwise) to 1 (clockwise).

% "abs_displacement", the spatial step size between each motion pulse. This
% variable was controlled by a staircase procedure.

% "abs_response_cw" TRUE if the subject responded clockwise.

% "folded_direction_content", 
% "folded_displacement"
% "folded_response_with_carrier"
% same data folded over so that carrier direction content is always
% positive. The response variable is TRUE if it agrees with the sign of the
% carrier.

% "target_spacing" The distance between targets along the circumference of
% circle.

% "target_number" The number of targets. target_spacing =
% 2*pi*eccentricity/target_number

% "subject" Subject initials.

% Another note, this only includes trials where the subject gave a response
% within a certain window from stimulus onset. They had feedback on whether
% their response latency were in the window. Some subjects' responses
% varied with the response latency, which isn't included in this data.

contentData = read_csv('content_series_trials.csv');
spacingData = read_csv('spacing_series_trials.csv');

fields = fieldnames(contentData);
allData = cellfun(@(x) struct2cell(orderfields(x, fields)), ...
                  {contentData, spacingData}, 'UniformOutput', 0)
allData = cell2struct(cellfun(@vertcat, allData{:}, 'UniformOutput', 0), fields, 1) 

save data allData contentData spacingData

