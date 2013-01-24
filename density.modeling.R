## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
library(ggplot2)
library(plyr)
library(grid)
library(psyphy)
library(ptools)
source("latexing.R")
source("icons.R")
source("scales.R")
source("library.R")
source("model.R")
setup_theme()

## @knitr do-not-run
if (interactive()) {
  quartz()
} else {
  cairo_pdf("density.pdf", onefile=TRUE)
  #on.exit(dev.off(), add=TRUE)
}

## @knitr density-load

load("../modeling/data.Rdata")
load("../modeling/slopeModel.RData")

#this just illustrates the combinations of number and density.

#these are the columns which define each "experiment"
segment.experiment.vars <-
  c("subject", "displacement", "content", "eccentricity")

#within an experiment these are the vars which separate each "stimulus condition"
segment.config.vars <-
  c("spacing", "target_number_shown", "target_number_all")

## we like to plot with folded data, and with the "segment" data we
## spindle the data too, as the stimuli were presented in different
## hemifields. Averaging foldings and hemifields is useful for
## plotting but not as good for modeling.
## "fold" collapses CW and CCW direction contents.
## "spindle" collapses CW and CCW spindling.
labeler <- mkchain(mutate(label=sprintf("%s d=%.2g C=%.2g",
                            toupper(subject), displacement, content)))

extract_segment <- function(df, fold=FALSE, spindle=FALSE)
  chain(df,
        subset(exp_type=="numdensity" & subject %in% names(models)),
        do.rename(folding=fold),
        mkrates(c(  segment.config.vars, segment.experiment.vars
                  , "eccentricity", if(!spindle) "side")),
        mutate(bias=if(fold) 0 else 1,
               sidedness=if(spindle) 0 else 1,
               side=if(spindle) data$side[[1]] else side,
               extent = spacing * target_number_shown),
        labeler)

##We'll be modeling raw data, but plotting folded/spindled. Here's a
##function that  "re-folds" the predictions
mutilate.predictions <-
  function(pred,
           fold=abs(diff(range(sign(pred$displacement)))) > 1,
           spindle=length(unique(pred$side)) > 1) {
    chain(pred,
          refold(fold),
          ddply(c(as.quoted(segment.config.vars),
                  (if (spindle)
                   c(as.quoted(quote(abs(content))),
                     as.quoted(segment.experiment.vars %-% "content"))
                  else as.quoted(segment.experiment.vars)),
                  if (fold) NULL else as.quoted("side")),
                summarize, fit = mean(fit), se.fit = sqrt(sum(se.fit^2))),
          labeler)
  }

bind[segment, segment.folded, segment.folded.spindled] <-
  Map(extract_segment, list(data),
      fold = c(FALSE, TRUE, TRUE),
      spindle = c(FALSE, FALSE,TRUE))

plot.basic <- (ggplot(segment.folded.spindled)
               + proportion_scale
               + spacing_texture_scale[-1]
               + number_color_scale[-1]
               + facet_wrap(~label)
               )



plot(number.plot <- plot.basic
     + aes(x=target_number_shown, group=factor(spacing), linetype=factor(spacing))
     + geom_line())

plot(spacing.plot <- plot.basic
     + aes(x=spacing)
     + geom_line(aes(group = factor(target_number_shown),
                     color = factor(target_number_shown))))

plot(extent.plot <- plot.basic
     + aes(x = extent,
           group = factor(target_number_shown),
           color = factor(target_number_shown),
           fill = factor(target_number_shown))
     + geom_line(linetype=1)
     + geom_line(aes(group = factor(spacing), linetype = factor(spacing),
                     color="black", fill="black")))

## now how can we reconcile this data with the models?

##let's start with a simple regression on the "unfolded" data
##"sidedness" and "bias" are control coefficients, set to 1 always
##unless I'm doing something funny with the predictions.

##do a regression for each subject:
segment.models <-
  ddply_along(
    segment, "subject",
    function(group, dataset) {
      model = glm(
        cbind(n_cw, n_ccw) ~ spacing + displacement
        + target_number_shown
        + factor(side) - 1,
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)),
        data=dataset)
      mutate(group, model=list(model))
    })
#now that it isn't folding we don't have warnings about failing to converge.

##plot fit lines from this...

segment.models.pred <-
  chain(segment.models,
        adply(1, function(row) {
          bind[model=bind[model], ...=group] <- as.list(row)
          data <- match_df(segment, quickdf(group))
          pred <- predict(model, newdata=data,
                          se.fit=TRUE, type="response")
          cbind(data, pred, model=NA) }),
        mutilate.predictions(fold=TRUE))

## TASK 2. for checking-purposes, plot these (local) model fits
## against the data.

prediction.layers <- function(dataset)
  with_args(
    data=dataset,
    mapping=aes(
      y=fit, ymin = fit - se.fit, ymax = fit + se.fit,
      fill=factor(target_number_shown),
      color=factor(target_number_shown)),
    geom_line(linetype=4), geom_ribbon(alpha=0.3, linetype=0))

(spacing.plot + prediction.layers(segment.models.pred))


#okay, that's a start. You'll note that the PBM data doesn't capture
#the two displacements where PBM always responds with the content. Ah,
#it's because the sensitivity to displacement changes with spacing.
#The fits are also rank deficient, because many fits do not sample
#different values of displacement. (or perhaps I should fit on
#unfolded data?! Ah, I think that might be it too. So what might what
#might make sense for that is to make the

## TASK_3. derive some prediction from the already fit model,
## somehow.
unattr <- function(x) {
  attributes(x) <- NULL; x}

segment_model_subject <- function(model, subject) {
  segment.model = model[[1]]
  model <- models[[subject]]
  matched.data <- match_df(segment, data.frame(subject))
  #skip subjects for whom we dont' have data
  if (is.null(model)) return(data.frame())

  #what we need to to is make the main model separate its two
  #different responses to "spacing." There's two "spacing" responses;
  #the one that parameterizes the slope (which I argue should not
  #change, at least in the subject's better hemifield) and another
  #based on summation within the hemifield; and a third based on
  #"induced motion"

  #so here is the model formula as it stands:
  old_formula <-  (
    cbind(n_cw, n_ccw) ~ displacementTerm(spacing, displacement) +
    content + I(content * abs(content)) + I(1/spacing):content +
    bias - 1)

  if (!identical(unattr(model$formula), unattr(old_formula)))
    stop("model formula has changed?")

  #now we've split "content" into "content_local" and
  #"content_global" so let's update to reflect what we think is
  #going on. Also fill in some other columns.
  new_data <-
    mutate(model$data,
           eccentricity = if(!exists("eccentricity")) 20/3 else eccentricity,
           target_number_shown = round(2*pi*eccentricity/spacing),
           target_number_all = target_number_shown,
           content_global = content * target_number_shown,
           content_local = content / spacing,
           extent = eccentricity * spacing * target_number_all)

  #now the simplistic model of induced motion is that it summates over
  #the whole visual field. However, we do not (in the full-circle data)
  #see more induced motion for larger numbers of elements -- my best
  #model just relates induced motion to the direction content. So
  #perhaps it is a normalization happening.

  #So it's a slight stretch to say that changing number will chance
  #induced motion, but that's what I'll do. I'll relate it to the
  #"extent" (space covered) by the stimulus. This has some resonance
  #with thinking of it as a center-surround type of effect.
  new_formula <-  (
    cbind(n_cw, n_ccw) ~ displacementTerm(spacing, displacement)
    + content:extent + I(content * abs(content)):extent
    + content_local
    + bias - 1)

  #Refit the model (this is still to the full-curcle-data. Despite
  #splitting up the variables we should have the same result (so same
  #residual deviance etc.)
  new_model <- update(model,
                      data=new_data,
                      formula=new_formula)
  if (deviance(model) - deviance(new_model) > 2) stop("models not equivalent...")

  #Let's look at what this recast model predicts...
  mutate(segment.model$data,
         content_local = content/spacing,
         content_global = content*target_number_shown,
         ) -> segment.data

  new.pred <- cbind(
    segment.data,
    predict(new_model, newdata=segment.data, type="link", se.fit=TRUE))

  #return segment-from-circle predictions Might return other things in the future.
  list(segment.pred.from.circle=new.pred)
}
alply(segment.models, 1, splat(segment_model_subject)) -> segment.analysis 

#juggle each subject's analysis into joined data frames
chain(segment.analysis,
      do.call(mapply, c(., FUN=list, SIMPLIFY=FALSE)),
      lapply(do.call, what=rbind)) -> segment.analysis
bind[segment.pred.from.circle=segment.pred.from.circle,
     ...=] <- segment.analysis

# graph these predictions...
(spacing.plot + predictions(segment.models.pred))

test <- do.call(mapply, c(segment.model.data, FUN=rbind))
aaply(do.call, rbind)

segment.model.data <- llply(segment.model.data, splat(rbind.fill))



pred <- do.call(pred, 1)


pred <- do.call(pred, 1)

## TASK 4. how to compare and reconcile? the

#another idea would be to cast raw predictions from the models onto
#this experiment.


