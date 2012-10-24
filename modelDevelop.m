%here is what we will do. Let's work on a subset of the data, to begin with.

%I'll work with unfolded data, here, since it allows me to separate an
%absolute "bias" (ML ans SK show an abolute bias) from the nonmonotonic
%effects of direction content (similar to those found in Murakami & Shomojo,
%1993, Bex and Dakin 2010).

%%
load('data.mat');
data = doRename(data, false); %no folding
%start with my data, to begin with.
subset = data(strcmp(data.subject, 'pbm') & ...
              strcmp(data.exp_type, 'spacing'), :);

%%
%To start off, we'll just fit a psychometric function: constant slope,
%constant bias. If the subjects perceive the envelope motion
%viridically this should work.
M = SlopeModel();
M.splits = {};
M.freeParams = {'beta_0', 'mu_0'};
M = M.fit(subset);
plotModel(M);
%plot those fits

%%
%Now, it's clear that the slope is generally too shallow, but it might
%be more clear what's going on if I plot residuals binned over dx,
%conditional on spacing.
resid_base = M.residuals({'spacing'}, 'dx', 25);
facetScatter(resid_base, 'x', 'dx', 'y', 'pearson_resid', 'color', ...
             'spacing', 'size', 'n_obs');

%Here the horizontal axis is dx (or "global speed") and the vertical
%axis is spacing. Colors indicate spacing, with cool colors indicating
%narrow spacings and warm colors indicating wide spacings. Of note, at
%most spacings there is a positive slope to the residuals, but at wide
%spacings there is a negative slope.

%This makes it clear: the model is underestimating the slope of the
%psychometric function on larger spacings and overestimating ths slope
%on smaller spacings. And, there seems to be a "breaking point"
%between the majority of the data, and the two or three smallest
%spacings. So the slope term needs to be dependent on the spacing.

%Models of what are going on with crowding differ, but the general
%descriptive consensus is that "threshold is constant above critical
%spacing and threshold increases when critical spacing decreases below
%critical.

%%
%So we want seneitivity to be asymptotic when spacing is large, and
%close to zero when spacing is small. However, it would not make sense
%to have sensitivity ever go negative. So I'll use a sigmoid: The
%seneitivity will be multiplied by 2 - 2/(1+exp(-cs/s)), where "s" is
%the spacing and "cs" is a model parameter. This has the
%attractive feature of only requiring one additional parameter. Here
%"cs" can be interpreted as a "critical spacing" number; it
%corresponds to the spacing at which threshold multiplies by a
%constant factor, versus the threshold for unflanked elements. (which
%is a frequently used empirical measure of critical spacing in
%crowding literature.) A CS of 0 eliminates dependency on spacing,
%which makes the model nicely nested, and the model does reasonable
%thengs even if you let CS go negative, which makes fminsearch
%happier.

%And the residuals look better:
resid_slopechange = M.residuals({'spacing'}, 'dx', 25);
facetScatter(resid_slopechange, ...
             'x', 'dx', ...
             'y', 'pearson_resid', ...
             'color', 'spacing', ...
             'size', 'n_obs' ...
             );

%% And the plot now can have a slope change:
M.freeParams = {'beta_0', 'mu_0', 'cs'}
M = M.fit();
plotModel(M);

%Note that at this point the model does not care at all about
%direction content; as far as the model is concerned there is none.
%We expect that as small spacings, the responses will be determined by
%direction content, and we see this: the small spacings (cool colors,
%here) have a positive slope. There are residual "ccw" responses when
%direction content is "ccw" and residual "cw" responses when direction
%content is "cw."

%%
%Now what you see is that -- particularly for the cool colors (low
%spacing)--the fits are off. That's reasonable as there is not yet any
%dependence on of direction content in this model.

% So let's add some dependence on direction content.

%Here, I'll show residuals as a function of direction content.
resid_content = M.residuals({'content', 'spacing'}, 'dx', 10000);

facetScatter(resid_content, ...
             'x', 'content', ...
             'y', 'pearson_resid', ...
             'color', 'spacing', ...
             'size', 'n_obs', ...
             'morePlotting', connectLines('content', 'spacing', ...
                                          'LineWidth', 3));
%The horizontal axis of this plot is direction content, while the
%vertical indicates residuals against the model fit (for a model
%without any dependence on direction content.) Colors indicates
%spacing -- warm is wide and cool is narrow.

%%

%Here's another way of looking at this, this time putting spacing on
%the x axis using colors for the direction content:
plotResiduals(resid_content, 'const_', ...
              'spacing', 'pearson_resid', ...
              'const_', 'const_', ...
              'content', 'n_obs');
%In bosth cases, we see curves crossing over each other.The influence
%of direction content plain just plain _reverses_ at some spacing.

%So we can see that the poorness of the model fit is a function of
%spacing. What we don't intuitively expect, but is also made plain by
%these residual plots, is that the effect actually REVERSES; at large
%spacings, there there is a negative slope to the residuals. At large
%spacings, more direction content makes you answer _against_ the
%direction content! And this effect actually looks just as large as
%the effect of direction content.

%Now let's think about the problem in cue combination terms. Suppose
%that both envelope and carrier motion systems give you a noisy
%measurement on a common variable velocity, (I'm not saying there
%necessarily is such a variable but bear with me.) Suppose the
%envelope-motion-detection is viridical but its precision is limited
%by crowding-the variability increases as flankers cross below the
%critical distance.  On the other hand you have a carrier motion
%mechanism that can detect an imbalance in motion energy no matter the
%target spacing. The variability of its signal is fixed but its
%magnitude is proportional to the degree of imbalance (let's say.)

%Then to get the overall response, the visual system multiplies these
%two disagreeing distributions (or, of you prefer not to think of
%explicitly represented distributions, adds the signals in inverse
%proportion to their reliability.)

%If you do this the straighforward way, it just works out to adding a
%term for diretion content; it determines come of the response at wide
%spacing, and when envelope response is attenuated at close spacings,
%it determines more of the response. The new constant is called
%'beta_content' in the model.

M.freeParams = {'mu_0', 'beta_0', 'cs', 'beta_induced'};
M = M.fit();
plotModel(M);

%Now, that didn't do much at all. See why? Look back at the residual
%plots: The dependency on direction content goes in both directions,
%and the model fit tries to split the difference. This is not easy to
%see when plotting curve fits directly (since the steep slope at large
%spacing visually obscures residuals that are just as large is those
%around shallow slopes) but plotting the residuals brings out the
%problem.

%So I'm going to borrow a point from Murakami and Shimojo (1993),
%Mareschal, Morgon and Solomon (1992) and note that repulsion at large
%distances becomes assimilation at small distances. Also, there is in
%the crowding literature an idea of "oblibatory summation;" that is
%you are obligated to sum up the signals within the criical distance
%of your object no matter if they don't come from your object. So
%"beta_induced" controls the repulsion (or "induced motion") portion
%of the response and I need another one to dominate at close
%range. The residual plot looks suggestively like an inverse function
%of spacing, which is consonant with the idea of summing up however
%many targets fall near the critical distance, so I'll use
%that. (Ponder: how does this square with the idea of critical
%distance which I used for the envelope response?)

%So there is just a new constant 'beta_summation' and the mean
%response is going (you could take an alternate tact with cue
%combination and say that the variance decreases with more summation,
%but it turns out the same as far as prediction probabilities of
%response.
M.freeParams = {'mu_0', 'beta_0', 'cs', 'beta_summation', 'beta_induced'};
M = M.fit();
M_content.parameters
plotModel(M);

%Ah, this starts looking pretty good. Let's see those residuals as a
%function spacing, and of direction content:
resid_with_switch = M.residuals({'content', 'spacing'}, 'dx', Inf);

plotResiduals(resid_with_switch, 'const_', ...
              'spacing', 'pearson_resid', ...
              'const_', 'const_', ...
              'content', 'n_obs');

plotResiduals(resid_with_switch, 'const_', ...
              'content', 'pearson_resid', ...
              'const_', 'const_', ...
              'spacing', 'n_obs');

%I don't seem much of a pattern left over either way.

resid_dx = M.residuals({}, 'dx', 50);
facet_scatter(resid_with_switch, ...
              'x', 'dx', 'y', 'pearson_resid', 'size', 'n_obs');
%Not sure if there's anything there with that last one.

%How about if we see how this works for the content-based
%experiment. Note that the direction-content experiment only tests at
%two spacings, so it is pointless to try fitting the critical spacing.
%I will leave 'cs' fixed at the value we already have:

M_content = M;
content_subset = data(strcmp(data.subject, 'pbm') ...
                      & strcmp(data.exp_type, 'content'), :);
M_content.freeParams = {'mu_0', 'beta_0', 'beta_summation', 'beta_induced'}
M_content = fit(M_content, content_subset);
M_content.parameters
%that's interesting, it came up with some different numbers. Much less
%difference between summation and induced motion for one.

%Here's the model fit for PBM, 'content' experiment.
plotModel(M_content);
%that all looks fairly reasonable.

r = M_content.residuals({'content', 'spacing'}, 'dx', Inf);

facetScatter(r, 'x', 'content', 'y', 'pearson_resid', ...
             'color', 'spacing', 'size', 'n_obs');
%That one starts to look like a linear dependence on direction content
%isn't quite right. There's a sort of S-shaped pattern to thsi
%graph. I could believe it's noise, until I see the same thing in
%another subject, where I'll put an asterisk (*)

%How about fitting both experiment types?
M_content = M;
content_subset = data(strcmp(data.subject, 'pbm') ...
                      & strcmp(data.exp_type, 'content'), :);

%Now we should stop playing around with this particularly well behaved
%subject and start looking at some others. Let's try JB since he's
%such a treat

subset = data(strcmp(data.subject, 'jb') ...
              & strcmp(data.exp_type, 'spacing'), :);
%Exceeding the maximum iterations is a bad sign and it looks like a bad fit.
%Let's add parameters one at a dime and see what the data look like.
M = SlopeModel;
M.freeParams = {'mu_0', 'beta_0'}; M = M.fit(subset);
plotModel(M, 'plotBy', {'x', 'dx', 'y', 'p', 'row', 'spacing', ...
                    'color', 'content', 'size', 'n'});

%You'll need a very tall window for that.  And really, look at
%that. Here, cool colors denote CCW direction contents, while warm
%colors denote CW direction contents.  You see that the influence of
%direction content, again, reverses from close to wide spacing -- Same
%thing as with PBM but here the slope change is not much, for some
%reason. 

M.freeParams = {'mu_0', 'beta_0', 'beta_summation', 'beta_induced'}; M = M.fit();
plotModel(M, 'plotBy', {'x', 'dx', 'y', 'p', 'row', 'spacing', ...
                    'color', 'content', 'size', 'n'});

%Hmm. adding direction content was not as nice as I hoped -- still underestimating the slope everywhere.
M.freeParams = {'mu_0', 'beta_0', 'cs', 'beta_summation', 'beta_induced'}; M = M.fit();
plotModel(M, 'plotBy', {'x', 'dx', 'y', 'p', 'row', 'spacing', ...
                    'color', 'content', 'size', 'n'});

%So my puzzle for you is, what is making JB so sensitive at the small
%spacings while so insensitive at the large spacings? I'll think this
%over but in the meantime add an extra parameter ('beta_small') to
%reflect this kind of anomalous sensitivity to envelope motion when
%you should be crowded.
M.parameters.beta_small = 0;
M.parameters.cs = 4;
M.freeParams = {'mu_0', 'beta_0', 'beta_summation', 'beta_induced', 'beta_small'}; M = M.fit();
plotModel(M, 'plotBy', {'x', 'dx', 'y', 'p', 'row', 'spacing', ...
                    'color', 'content', 'size', 'n'});
%Interestingly, that parameter didn't work. Ah, I see a potential
%problem. Look what's happening at spacing = 10.5: There's a lot of
%direction contents being tested, but they all collapse onto these
%'leftward' or 'rightward' functions. So what's the deal there? The
%magnitude of induced motion isn't really sensitive to content?

%So, what the heck is going on there? Notice that the response
%actually seems nonmonotonic, with the blues/reds actually being less
%important than the grays.

%Does that happen for the content experiment as well?
subset = data(strcmp(data.subject, 'jb') ...
              & strcmp(data.exp_type, 'content'), :);
Mjb_c = SlopeModel();
Mjb_c.initialParamsDefaults.cs = 4;
Mjb_c.freeParams = {'mu_0', 'beta_0', 'beta_summation', 'beta_induced'};
Mjb_c = Mjb_c.fit(subset);
%Let's look at the residuals against content here.
Rjb_c = Mjb_c.residuals({'content', 'spacing'}, 'dx', Inf);
facetScatter(Rjb_c, 'x', 'content', 'y', 'pearson_resid', 'color', 'spacing', 'size', 'n_obs');
%Wow, there's that S-shaped thing again (*). and it looks too similar at
%both spacings to be a coincidence.  Comparing the residuals to the
%earlier fit plot, it appears to be saying that middling direction
%contents cause more shift than the model really tolerates -- it's not
%linear in direction content. So the fit is not reacting strongly
%enough to middling direction contents.

%Let's try again with log scaling the direction content.
Mjb_c.parameters.beta_summation = 0;
Mjb_c.parameters.beta_induced = 0;
Mjb_c.freeParams = {'mu_0', 'beta_0', 'sbeta_summation', 'sbeta_induced', 'beta_small'}; Mjb_c = Mjb_c.fit();
Rjb_c = Mjb_c.residuals({'content', 'spacing'}, 'dx', Inf);
facetScatter(Rjb_c, 'x', 'content', 'y', 'pearson_resid', 'color', 'spacing', 'size', 'n_obs');
plotModel(Mjb_c, 'plotBy', {'x', 'dx', 'y', 'p', 'row', 'spacing', 'color', 'content'});

%So I've put a compressive transform on direction content. Does that
%help us with the modeling? Looks like answer is YES.

M.parameters.beta_summation = 0;
M.parameters.beta_induced = 0;
M.freeParams = {'mu_0', 'beta_0', 'cs', 'beta_small', 'sbeta_summation', 'sbeta_induced'}; M = M.fit();
plotModel(M, 'plotBy', {'x', 'dx', 'y', 'p', 'row', 'spacing', 'color', 'content'});

%So we've learned to put a compressive transform on direction content. 
%And it looks like beta_small might not be needed after all.
%(Brainstorm: depending on the interaction between dx and motion
%energy, and people's differing amounts of summation/induced motion
%tradeoff, that may explain why beta_small looks needed sometimes.

%Let's look at a bunch of residuals here, binned over dx, colored by spacing.
R_jb = M.residuals({'content', 'spacing'}, 'dx', Inf);
facetScatter(R_jb, 'x', 'content', 'y', 'pearson_resid', 'color', 'spacing');

%Hmm. looks like direction content response may actually be nonmonotonic...

%Before messing any more with the compressive term, it would be
%instructive to look at actual motion energy.

