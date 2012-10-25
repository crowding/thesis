%here is what we will do. Let's work on a subset of the data, to begin with.

%I'll work with unfolded data, here, since it allows me to separate an
%absolute "bias" (ML ans SK show an abolute bias) from the nonmonotonic
%effects of direction content (similar to those found in Murakami & Shomojo,
%1993, Bex and Dakin 2010).

%%
load('data.mat');
data = doRename(data, false); %no folding
extract_experiment = @(subj, exp) data(strcmp(data.subject, subj) & ...
                                       strcmp(data.exp_type, exp), :);

%%
%To start off, we'll just fit a single psychometric function: constant
%slope, constant bias. If the subjects perceive the envelope motion
%viridically this should work.
%Start with my data, to begin with.
M = SlopeModel('freeParams', {'beta_0', 'mu_0'});
M = M.fit(extract_experiment('pbm', 'spacing'));
plotModel(M);

%%
%Now, it's clear that the slope is generally too shallow,
%be more clear what's going on if I plot residuals binned over dx,
%conditional on spacing.
resid_base = M.residuals({'spacing'}, 'dx', 25);
facetScatter(resid_base, 'x', 'dx', 'y', 'pearson_resid', 'color', ...
             'spacing', 'size', 'n_obs');

%Here the horizontal axis is dx (or "envelope speed") and the vertical
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

%This is already in the model (SlopeModel.m) as 'cs' which was
%previsouly set to 0. Adding it to the list of free parameters, the
%model now fits with a change of slope.
M.freeParams = {'beta_0', 'mu_0', 'cs'}
M = M.fit();
plotModel(M);

%%
%And the residuals binned over displacement look less patterned
resid_slopechange = M.residuals({'spacing'}, 'dx', 25);
facetScatter(resid_slopechange, 'x', 'dx', 'y', 'pearson_resid', 'color', ...
             'spacing', 'size', 'n_obs');

%% after glancing at the residuals, look at the curve fits again.
plotModel(M);
%%
%Note that at this point the model does not care at all about
%direction content We expect that as small spacings, the responses
%will be determined by direction content. Indeed, looking at the small
%spacings, we see that the model is fitting well below the data for
%positive direction content, and above the data for negative.

resid_content = M.residuals({'content', 'spacing'});
%since I'll be plotting variations on this particular residual plot,
%i'll save a shortcut for these arguments.
content_vs_residual = { ...
    'x', 'content', 'y', 'pearson_resid', ...
    'color', 'spacing', 'size', 'n_obs', ...
    'morePlotting', connectLines('content', 'spacing', 'LineWidth', 3)};

facetScatter(resid_content, content_vs_residual{:});
%The horizontal axis of this plot is direction content, while the
%vertical indicates residuals against the model fit (for a model
%without any dependence on direction content.) Colors indicates
%spacing -- warm is wide and cool is narrow.

%So, yes, at narrow spacings (blue lines) there is positive dependence
%on direction content, as we expected. What's perhaps unexpected is
%the opposite pattern ad wide spacings, which appears to be just as
%strong (trick question: why doesn't it look just as strong in the
%graph?)
%%

%So before we add direction content to the model there's another thing to
%look at whcih is how the direction content interacts with
%spacing. Here is a plot of the same where the x=axis is spacing:
spacing_vs_residual = { ...
    'x', 'spacing', 'y', 'pearson_resid', ...
    'color', 'content', 'size', 'n_obs', ...
    'morePlotting', connectLines('spacing', 'content', 'LineWidth', 3)};
facetScatter(resid_content, spacing_vs_residual{:});

%So, yes, there is a strong effect of direction content (color) ad the
%narrow (left) end of the scale. To me these curves look like
%1/spacing, at the left end of the scale, which makes sense to me, and
%that's what I code in the model. (Ponder: how does this square with
%the idea of critical distance which I used for the envelope
%response?)

%%
% One way to think about the model is in cue combination terms. Suppose
% that both envelope and carrier motion systems give you a noisy
% measurement on a common variable velocity, Suppose the
% envelope-motion-detection is viridical but its precision is limited
% by crowding-the variability increases as flankers cross below the
% critical distance.  On the other hand you have a carrier motion
% mechanism that can detect an imbalance in motion energy no matter the
% target spacing. The variability of its signal is fixed but its
% magnitude is proportional to the degree of imbalance (let's say.)

% If you do this the straighforward way, it just works out to adding a
% term for direction content; it determines some of the response at wide
% spacing, and when envelope response is attenuated at close spacings,
% it determines more of the response. The new constant is called
% 'beta_summation' in SlopeModel.m and looks like (beta_summation *
% content / spacing). Adding it to the model,

M.freeParams = {'mu_0', 'beta_0', 'cs', 'beta_summation'};
M = M.fit();
plotModel(M);

% But now you can see that the behavior at wide spacings (warm colors)
% isn't fitting. To me it looks like for positive direction content
% (top row), the data for side spacing (reddish) is above the fit, and
% vice versa at for counterclockwise direction content. So at the wide
% end of the scale the effect of direction content is reversed. (we
% already saw this in the residual plots.)

%%
%Plotting residuals as a function of spacing: for different direction
%contents:
resid_spacing = M.residuals({'spacing', 'content'});
facetScatter(resid_spacing, spacing_vs_residual{:})

%Looks like we really can't ignore that. What we've been looking at,
%in several different ways, at the fact that the effect of direction
%content on decisions actually reverses as a function of spacing. It's
%not a small effect either; in the way that it exerts leverage on the
%model, it's equally large as the local summation effect that we're
%more interested in. So to obtain a better fit for the parameter we
%are interested in, we're going to have to add this large-spacing
%effect to the model.

%So I'm going to borrow a point from Murakami and Shimojo (1993),
%Mareschal, Morgon and Solomon (1992) and note that repulsion at large
%distances becomes assimilation at small distances (induced motion
%becomes assimilation; tilt illusion becomes crowding, etc.)

%%
% So here I'll add another parameter, called "beta_induced". Again,
% this was in SlopeModel.m all along, but was previously set to 0.
M.freeParams = {'mu_0', 'beta_0', 'cs', 'beta_summation', 'beta_induced'};
M = M.fit();
plotModel(M);

%This is starting to look good.

%%
%The residual is less suggestive too.
resid_spacing = M.residuals({'spacing', 'content'});
facetScatter(resid_spacing, content_vs_residual{:});
%I don't see much of a pattern left over.

%%
%Now, that's just one subject, on one experiment. The next thing to
%look for is whether this also works for the direction-content
%experiment (where direction content is varied within session, rather
%than between session. Note that the direction-content experiment only
%tests at two spacings, so it is pointless to try fitting the critical
%spacing.  I will leave 'cs' fixed at the value we already have:

M_content = M;
content_subset = data(strcmp(data.subject, 'pbm') ...
                      & strcmp(data.exp_type, 'content'), :);
M_content.freeParams = {'mu_0', 'beta_0', 'beta_summation', 'beta_induced'}
M_content = fit(M_content, content_subset);
M_content.parameters
%That's interesting, it came up with some different numbers. Summation
%and induced motion both much less strong, which may indicate some
%kind of adaptation to the stronger direction contents used in this
%experiment.
plotModel(M_content);
%%
%that looks reasonable, but the "large spacing" model fits look off
%for some reason, low in some places and high in others. Look, in
%particular, at the second, third, second-from-bottom and
%third-from-bottom plots -- the purple curve looks off. This doesn't
%harm my fit so much but other subjects like JB are more pronounced here.
resid_content = M_content.residuals({'content', 'spacing'}, 'dx', Inf);
facetScatter(resid_content, content_vs_residual{:});
%That one starts to look like a linear dependence on direction content
%isn't quite right. There's a sort of S-shaped pattern to this
%graph. I could believe it's noise, until I see the same thing in
%another subject, which I will.

%%
%Now we should stop playing around with this particularly well behaved
%subject and start looking at some others. Let's try JB since he's
%got a lot of data, and seems markedly different.
Mjb = SlopeModel('freeParams', ...
                 {'mu_0', 'beta_0', 'cs', 'beta_induced', 'beta_summation'});
Mjb = Mjb.fit(extract_experiment('jb', 'spacing'));
%Exceeding the maximum iterations is a bad sign and it looks like not
%a great fit.
plotModel(Mjb);
%%
%Let's interchange the colors and row-facets; so the next plot will
%group subplots by spacing and color by direction content.
facet_row_spacing = {'x', 'dx', 'y', 'p', 'row', 'spacing', ...
                    'color', 'content', 'size', 'n'};
plotModel(Mjb, 'plotBy', facet_row_spacing);

%Here, cool colors denote CCW direction contents, while warm colors
%denote CW direction contents.  You see that the influence of
%direction content, again, completely reverses from close to wide
%spacing as we saw in PBM, but stronger. The "induced" and "summation"
%terms capture some of the switchover, but not all, and the slopes are
%too shallow as a result.
%%
% So what's not being captured? Let's look at residuals as a function
% of spacing and direction content.

resid_jb = Mjb.residuals({'content', 'spacing'});
facetScatter(resid_jb, content_vs_residual{:});
% Hey, look, it's that same wiggle again. Especially at wide spacings.
% So, at wide spacings, the "beta_summation" is not capturing it. It
% does look like hte wiggle disappears at narrow spacings, though --
% so may be the wiggle is tied in to (this would be consistent with
% previous literature on induced motion and the tilt illusion.

% So the "wiggle" is evident in at least two subjects. As we'll come
% back to, it's actually present in all subjects.

%%
%I poked around and decided that replacing the linear dependence a
%logistic to the third derivative of a logistic should be able to fit
%the wiggles. This is implemented in two parameters, 'saturating_induced',
%and 'wiggle_induced', which replace 'beta_induced' in the code. I
%select the arbitrarily.
Mjb = SlopeModel();
Mjb.freeParams = {'mu_0', 'beta_0', 'cs', 'beta_summation', ...
                  'saturating_induced', 'wiggle_induced'};
Mjb = Mjb.fit(subset);
plotModel(Mjb, 'plotBy', facet_row_spacing);
%I think that successfully captured some wiggle. Slopes are still not
%quite steep enough though.

%%
resids = Mjb.residuals({'content', 'spacing'});
facetScatter(resids, content_vs_residual{:})
%I'm still not totally satisfied with that, but it looks better than
%before. I wonder if there needs to be some nonlinearity in the
%summation, or a component that scales with the number of elements on
%screen.

%%
%Now, can we fit both experiment types with the same parameters?
extract_subject = @(subj) data( strcmp(data.subject, subj) ...
                              & ( strcmp(data.exp_type, 'spacing')  ...
                                | strcmp(data.exp_type, 'content') ) ...
                                , :);
M = SlopeModel('splits', {'subject'}, ...
               'freeParams', {'mu_0', 'beta_0', 'cs', 'beta_summation', ...
                              'saturating_induced', 'wiggle_induced'}, ...
               'data', extract_subject('pbm'));
M = M.fit();
plotModel(M, 'plotBy', {facet_row_spacing{:}, 'col', 'exp_type'});
%Heyyyyyy that's actually pretty okay looking. Does it work for other subjects?

%%
%So let's generate figures for everyone. It's the full model, and I'm
%not differentiating between content and spacing experiments. The
%"wiggle" parameter seems to adequately explain the difference between
%content and spacing experiments.
fitSubject = ...
    (@(subj) fit(SlopeModel('splits', {'subject'}, ...
                            'freeParams', ...
                            {'mu_0', 'beta_0', 'cs', 'beta_summation', ...
                             'saturating_induced', 'wiggle_induced'}, ...
                            'data', extract_subject(subj))));

plotFit = (@(model, varargin) ...
           plotModel(model, 'plotBy', ...
                     {facet_row_spacing{:}, 'col', 'exp_type', varargin{:}}));

%so far we're only not making sense of ML and TL (and I'm not sure
%they aren't random button pushers.) Some of the subjects are not
%fitting with as high of a slope as the data suggests; motion-energy
%effects?
Mpbm = fitSubject('pbm'); plotFit(Mpbm, 'newFigures', true);
Mml = fitSubject('ml'); plotFit(Mml, 'newFigures', true);
Mnj = fitSubject('nj'); plotFit(Mnj, 'newFigures', true);
Mtl = fitSubject('tl'); plotFit(Mtl, 'newFigures', true);
Mjb = fitSubject('jb'); plotFit(Mjb, 'newFigures', true);
Mns = fitSubject('ns'); plotFit(Mns, 'newFigures', true);
Mmc = fitSubject('mc'); plotFit(Mmc, 'newFigures', true);
Mje = fitSubject('je'); plotFit(Mje, 'newFigures', true);
Mcj = fitSubject('cj'); plotFit(Mcj, 'newFigures', true);
Mas = fitSubject('as'); plotFit(Mas, 'newFigures', true);
