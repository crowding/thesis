% here is what we will do. Let's work on a subset of the data, to begin with.

%I'll work with unfolded data, here, since it allows me to separate an
%absolute "bias" (ML ans SK show an abolute bias) from the nonmonotonic
%effects of direction content (similar to those found in Murakami & Shomojo,
%1993, Bex and Dakin 2010).

load('data.mat');
data = doRename(data, false); %no folding
%start with my data, to begin with.
subset = data(strcmp(data.subject, 'pbm') & ...
              strcmp(data.exp_type, 'spacing'), :);

%To start off, we'll just fit a psychometric function: constant slope,
%constant bias. If the subjects perceive the envelope motion
%viridically this should work.
M = SlopeModel();
M.splits = {};
M.freeParams = {'beta_0', 'mu_0'};
M = M.fit(subset);
plotModel(M);
%plot those fits

%Now, it's clear that the slope is generally too shallow, but it might
%be more clear what's going on if I plot residuals binned over dx,
%conditional on spacing.
resid_base = M.residuals({'spacing'}, 'dx', 25);
%                    figure    x     y                row       col
plotResiduals(resid_base, 'const_', 'dx', 'pearson_resid', 'const_', 'const_', ...
              'spacing', 'n_obs');
%             color,     size

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
%spacing and increases when critical spacing decreases below critical.

%So we want seneitivity to be asymptotic when spacing is large, and
%close to zero when spacing is small. However, id would ot make sense
%to have sensitivity ever go negative. So I'll use a sigmoid: The
%seneitivity will be multiplied by 2 - 2/(1+exp(-cs/s)), where "s" is
%the spacing and "cs" is a model parameter. This has the
%attractive feature of only requiring one additional parameter. Here
%"cs" can be interpreted as a "critical spacing" number; it
%corresponds to the spacing at which threshold doubles (which is a
%frequently used empirical measure of critical spacing in crowding
%literature.) A CS of 0 eliminates dependency on spacing, which makes
%the model nicely nestedm, and the model does reasonable thengs even
%if you let CS go negative, which makes fminsearch happier.

M.freeParams = {'beta_0', 'mu_0', 'cs'}
M = M.fit();
plotModel(M);
%Nice, you see how the slope is allowed to change.

resid_slopechange = M.residuals({'spacing'}, 'dx', 25);
plotResiduals(resid_slopechange, 'const_', ...
              'dx', 'pearson_resid', ...
              'const_', 'const_', ...
              'spacing', 'n_obs');

%These residuals now look less patterned, as a function of spacing. So
%the slope dependency looks like a win. But looking at the raw fits,
%there is still something going on with direction content. Let's look at
%the residuals conditioned on content and spacing.
resid_content = M.residuals({'content', 'spacing'}, 'dx', 10000);
plotResiduals(resid_content, 'const_', ...
              'content', 'pearson_resid', ...
              'const_', 'const_', ...
              'spacing', 'n_obs');

%The horizontal axis of this plot is direction content, while colors
%indicate residuals against the model fit thus far.

%Note that at this point the model does not care at all about
%direction content; as far as the model is concerned there is none.
%We expect that as small spacings, the responses will be determined by
%direction content, and we see this: the small spacings (cool colors,
%here) have a positive slope. There are residual "ccw" responses when
%direction content is "ccw" and residual "cw" responses when direction
%content is "cw."

%Here's another way of looking at this, this time putting spacing on
%the x axis using colors for the direction content:
plotResiduals(resid_content, 'const_', ...
              'spacing', 'pearson_resid', ...
              'const_', 'const_', ...
              'content', 'n_obs');

%Pretty cute, huh? The influence of direction content plain just plain
%_reverses_ at some spacing.

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

%If you do this the atraighforward way, it just works out to adding a
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
%around shallow slopes) but looking at residuals brings out the
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

resid_dx = M.residuals({}, 'dx', 50);
plotResiduals(resid_with_switch, 'const_', ...
              'dx', 'pearson_resid', ...
              'const_', 'const_', ...
              'const_', 'n_obs');

%Not sure if there's anything there with that last one.

%I don't seem much of a pattern left over either way

%How about if we see how this works for the content-based
%experiment. Note that the direction-content experiment doesn't have
%data to localize the critical spacing, so I will leave it fixed at
%the value we already know:

M_content = M;
content_subset = data(strcmp(data.subject, 'pbm') & strcmp(data.exp_type, ...
                                                  'content'), :);

%So the dependency on direction content also needs to interact with
%spacing.

%%{ I'm not sure of any of the following
Sketching this out, we see that adding the direction content signal
%cue-combination style requires two additional parameters: the
%variability of the direction content signal, and the proportionality
%constant of the direction content signal.

%Now, we don't have to write the equations in terms of the cue
%combination model, just use any two parameteters that capture the
%same behavior.  combination model. So instead of explicitly using
%variance and constant-of-direction-constant-dependence as model
%parameters, I'll pick two parameters that are capable of capturing
%the same behavior, but will be nicer to do fminsearch or regression
%with.

%Now, of you sketch out the behavior we're postulating, you'll find
%that in completely crowded conditions, response rate will be
%determined by the ratio of the carrier signal to the carrier
%variability, while in completely uncrowded situations, any effect of
%the carrier content will be weighted by the proportion of carrier
%signal variability to envelope signal variability. Working through
%the math, this just makes for a carrier sensitivity term determining
%the response.
%%}


%First, with just one coefficient for the direction content. This is
%what straight up cue combination would have happen.

%Huh, that hardly changed a thing, did it?
%You see, direction content is pulling this model in both directions. 


%So what about the part where the sensitivity to carrier moteion
%actually _reverses_? I This has a similarity to the observed
%phemomena (Mareschal, Morgan and Solomon 2010; Murakami and Shimojo
%1993) where assimilation turns into repulsion. I will presume at first
%this is controlled by the same "critical distance" as the sensitivity
%to global motion, anthough this might have to be adjusted.



%But let's parameterize it so that it's more amenable to model nesting
%and such.

%Sneakily, by allowing the this parameterization allows the "local weight" to be
%negative as well as positive.

%Note I have no reason to distinguish variability from
%while our idea of the variability of the
%envelope signal is not constant


%Now, we already have this model of the sensitivity to global
%motion. If we say that the signal from direction content is roughly
%constant in its uncertainty, and the, then, using cue combination
%ideas, the influence of direction content should be proportional to
%its share of the uncertainty.

%Another thing you can see here is that the spacing effect
%looks nearly as large for 0.05 direction content as it did for 1.0
%direction content. It certainly doesn't look 3 times bigger when the
%content is three times bigger, at least. So what's going on there?
%One problem is that 0.05 and 0.15 contents were tested in separate
%sessions, and there may be some adaptation going on.

%So now let's look at what this does with data where direction content
%is varied within the course of an experiment.
subset = data(strcmp(data.subject, 'pbm') & ...
              strcmp(data.exp_type, 'content'), :);
M = SlopeModel()
M.freeParams = {'mu_0', 'beta_0', 'cs'};
M = M.fit(subset);
plotModel(M);
%Okay, that's really bad, it even got the wrong sign on the slopw on
%account of the adaptive psychophysics. We can see that the also needs
%to move up and down in response to the direction content.

%Approaching thsi from a cue combination standoint, The _bias_, that is, the will be proportional to the 

%That effect is harder to see on the raw graphs, because the slope is
%steeper. But Pearson residuals are normalized against predictive
%uncertainty, and under that scaling both signs of the
%content-dependency look equally large.



M.freeParams{end+1} = 'beta_ks';