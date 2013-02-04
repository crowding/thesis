## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
library(ggplot2)
library(plyr)
library(grid)
library(gnm)
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
  cairo_pdf("density2.pdf", onefile=TRUE)
  #on.exit(dev.off(), add=TRUE)
}

## @knitr density-load

load("../modeling/data.Rdata")
load("../modeling/slopeModel.RData")

#this just illustrates the combinations of number and density.

#these are the columns which define each "experiment" (facet on the
#graph)
segment.experiment.vars <-
  c("subject", "displacement", "content", "eccentricity")
#within an experiment these are the vars which separate each "stimulus
#condition" (data point on the graph)
segment.config.vars <-
  c("spacing", "target_number_shown", "target_number_all")

#label function for each facet
labeler <- function(data) {
  mutate(data, label=sprintf("%s d=%s C=%s",
                 toupper(subject),
                 format(displacement,digits=2),
                 format(content,digits=2)))
 }


## we like to plot with folded data, and with the "segment" data we,
## uh, "spindle" collapsing stimuli presented in different
## hemifields. Averaging foldings and hemifields is useful for
## plotting but not as good for modeling. "fold" collapses CW and CCW
## direction contents.  "spindle" collapses stimulus locations.
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

bind[segment, segment.folded, segment.folded.spindled] <-
  Map(extract_segment, list(data),
      fold = c(FALSE, TRUE, TRUE),
      spindle = c(FALSE, FALSE,TRUE))

#base plot
plot.basic <- (ggplot(segment.folded.spindled)
               + proportion_scale
               + spacing_texture_scale[-1]
               + number_color_scale[-1]
               + facet_wrap(~label)
               )

#plot with x-axis of target number, lines of constant spacing
plot(number.plot <- plot.basic
     + aes(x=target_number_shown, group=factor(spacing), linetype=factor(spacing))
     + geom_line())

#plot with x-axis of target spacing, lines of constant number
plot(spacing.plot <- plot.basic
     + aes(x=spacing)
     + geom_line(aes(group = factor(target_number_shown),
                     color = factor(target_number_shown))))

#plot with x-axis of "extent"
plot(extent.plot <- plot.basic
     + aes(x = extent,
           group = factor(target_number_shown),
           color = factor(target_number_shown),
           fill = factor(target_number_shown))
     + geom_line(linetype=1)
     + geom_line(aes(group = factor(spacing), linetype = factor(spacing)),
                 color="black", fill="black"))

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
        cbind(n_cw, n_ccw) ~ displacement
        + target_number_shown:content
        + content/spacing + content
        + factor(side) - 1,
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)),
        data=dataset)
      mutate(group, model=list(model))
    })

##We'll be modeling raw data, but plotting folded/spindled. Here's a
##function that  "re-folds" the predictions
mutilate.predictions <-
  function(pred,
           fold=abs(diff(range(sign(pred$content)))) > 1,
           spindle=length(unique(pred$side)) > 1) {
    columns <- c(as.quoted(segment.config.vars),
                 as.quoted(segment.experiment.vars),
                 if (spindle) NULL else as.quoted("side"))
    chain(pred,
          refold(fold),
          ddply(columns, summarize,
                fit = mean(fit), se.fit = sqrt(sum(se.fit^2))),
          labeler)
  }

##So extract model predictions for our descriptive model (and re-fold them)
segment.models.pred <-
  chain(segment.models,
        adply(1, function(row) {
          bind[model=bind[model], ...=group] <- as.list(row)
          data <- match_df(segment, quickdf(group))
          pred <- predict(model, newdata=data,
                          se.fit=TRUE, type="response")
          cbind(data, pred, model=NA) }),
        mutilate.predictions(fold=TRUE, spindle=TRUE)
        )
## TASK 2. for checking-purposes, plot these (local) model fits
## against the data.

#ggplot layers to add predictions
prediction.layers <- function(dataset, connect=c("number","spacing"))  {
  connect <- match.arg(connect)
  eval(template(
         with_arg(
           data=dataset,
           mapping=aes(
             y=fit, ymin = fit - se.fit, ymax = fit + se.fit,
             ...(list(
                   number=alist(
                     color=factor(target_number_shown),
                     fill=factor(target_number_shown)),
                   spacing=alist(
                     linetype=factor(spacing)
                     ))[[connect]])),
           geom_line(...(if (connect=="number") list(linetype=3) else list())),
           geom_ribbon(alpha=0.3, linetype=0))))
}

print(spacing.plot + prediction.layers(segment.models.pred))

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
           extent = spacing * target_number_all)

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
    + content + I(content * abs(content))
    + content_local:extent
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
         extent = extent*2
         ) -> segment.data

  new.pred <- cbind(
    segment.data,
    predict(new_model, newdata=segment.data, type="response", se.fit=TRUE))

  #what if we just added a little offset to reconcile the predictions
  #with the data?
  #this doesn't propabage error bars but:
  offset.pred <-
    chain(segment.data,
          cbind(., predict(new_model, newdata=., type="link", se.fit=FALSE)),
          glm(new_model, data=., cbind(n_cw, n_ccw) ~ fit:content,
              family=binomial(link=logit.2asym(g=0.025, lam=0.025))),
          cbind(segment.data, predict(., type="response", se.fit=TRUE)))

  #return segment-from-circle predictions Might return other things in the future.
  list(segment.circle.model = I(list(new_model)),
       segment.pred.from.circle=new.pred,
       offset.pred.from.circle=offset.pred)
}

#analyze each subject as above, then juggle into topical data frames
alply(segment.models, 1, splat(segment_model_subject)) -> segment.analysis
chain(segment.analysis,
      do.call(mapply, c(., FUN=list, SIMPLIFY=FALSE)),
      lapply(do.call, what=rbind)) -> segment.analysis
bind[segment.pred.from.circle=segment.pred.from.circle,
     offset.pred=offset.pred,
     ...=] <- segment.analysis

# graph these predictions...
print(spacing.plot + prediction.layers(mutilate.predictions(segment.pred.from.circle)))

print(number.plot +
      prediction.layers(mutilate.predictions(segment.pred.from.circle),
                        connect="spacing"))



## TASK 4. how to compare and reconcile? Refit a constrained model?
