## The `*_trials.csv` files contain individual trial data.
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

## The `*_calculations.csv` files contain psychometric function fits.

These are fits for the same data in the `*_trials.csv` files.

Here's a couple of CSV files. There are two, `spacing_series` and
`contrast_series.` The difference is that in `spacing_series` the
directional content was held constant for each entire session and the
spacing varied, while in "contrast_series" I used only 2 values of
spacing and four values of directional content.

The two situations are not directly comparable because subjects appear
to adapt to the average amount of directional content used in any
session. So varying directional content within a session is much more
effective than varying it across sessions, which is why I've separated
out the data.

`subject`, `folded_direction_content`, and `target_spacing` should be
self explanatory.

The psychometric functions are fit by logistic CDFs with a 5%
allowance for lapses in either direction.

`bias` is just the intercept coefficient you get from logistic
regression; that is, it measures (in log-odds) how often subjects
respond `clockwise` to a stimulus with clockwise carrier and no
envelope motion. `bias_sem` gives you a standard error of the bias
measurement.

`sensitivity` is the slope parameter for the logistic CDF and
`sensitivity_sem` its standard error.

`pse` is the envelope-displacement-per-step at which the subject
equivocates between directions. Since steps have intervals of 100 ms,
multiply by 10 to get envelope speed in degrees/sec.

`pse_25` and `pse_75` gives a confidence intervals on the PSE
measurement (these are calculated by parametric bootstrap; I draw
random psychometric functions using the fitted likelihood function,
and find the 25% and 75% quantiles of the PSE)

`threshold` measures how much envelope motion it takes to move between
50% and 75% on the psychometric function. `threshold_25` and
`threshold_75` give a confidence interval.

Note that when sensitivity is low enough that you aren't entirely sure
if sensitivity is positive, the confidence intervals for PSE and
threshold will blow up (It's basically the problem of estimating the
x-intercept of a linear regression)
