##{{{ ---------- SETUP AND LOAD DATA -------------------------------------

## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center = TRUE)
library(ggplot2)
library(plyr)
library(grid)
library(gnm)
library(psyphy)
library(ptools)
#source("latexing.R")
#source("icons.R")
source("scales.R")
source("library.R")
source("model.R")
setup_theme()

## @knitr do-not-run
if (interactive()) {
  figure("plot")
} else {
  cairo_pdf("density2.pdf", onefile=TRUE)
  #on.exit(dev.off(), add=TRUE)
}

## @knitr density-load
load("data.RData")
load("slopeModel.RData")

##}}}
##{{{ ---------- COUNT AND LABEL TRIALS ----------------------------------

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
  mutate(data, label = sprintf("%s d=%s C=%s",
                 toupper(subject),
                 format(displacement, digits = 2),
                 format(content, digits = 2)))
}


## we like to plot with folded data, and with the "segment" data we,
## uh, "spindle" collapsing stimuli presented in different
## hemifields. Averaging foldings and hemifields is useful for
## plotting but not as good for modeling. "fold" collapses CW and CCW
## direction contents.  "spindle" collapses stimulus locations.
extract_segment <- function(df, fold=FALSE, spindle=FALSE)
  chain(df,
        subset(exp_type=="numdensity" & subject %in% names(models)),
        do.rename(folding = FALSE), # we handle the folds more comprehensively.
        refold(fold = fold),
        mkrates(c(  segment.config.vars, segment.experiment.vars
                  , "eccentricity", if(!spindle) "side")),
        mutate(bias = if (fold) 0 else 1,
               sidedness = if (spindle) 0 else 1,
               side = if (spindle) NA else side,
               extent = spacing * target_number_shown),
        labeler)

##Aggregate data into counts of CW and CCW responses, with various
##levels of folding/spindling
bind[segment, segment.folded, segment.folded.spindled] <-
  Map(extract_segment, list(data),
      fold = c(FALSE, TRUE, TRUE),
      spindle = c(FALSE, FALSE,TRUE))

##}}}


##{{{ ---------- BASIC PLOTS ----------------------------------------

#base plot

plot.basic <- (ggplot(segment.folded.spindled)
               + proportion_scale
               + spacing_texture_scale[-1]
               + number_color_alt_scale[-1]
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

##We'll be modeling raw data, but plotting folded/spindled. Here's a
##function that "re-folds" the predictions so that they can be plotted
##on a folded plot.
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

#Build ggplot layers to add predictions to a plot.
prediction_layers <- function(dataset, connect = c("number","spacing"))  {
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

predict_from_model_frame <- function(models, newdata, fold=TRUE, spindle=TRUE) {
  ##take a data frame with a list of models, and the variables to
  ##match by, produce predictions for the folding data.
  chain(models,
        adply(1, function(row) {
          bind[model=bind[model], ...=group] <- as.list(row)
          data <- match_df(segment, quickdf(group),
                           on = names(segment) %^% names(group))
          pred <- predict(model, newdata=data,
                          se.fit=TRUE, type="response")
          cbind(data, pred, model=NA)
        }),
        mutilate.predictions(fold=TRUE, spindle=TRUE)
        )
}
##}}}


##{{{ ---------- DESCRIPTIVE MODEL OF NUMBER/DENSITY DATA ----------

## Let's start by descriptively modeling the segment data. We see from
## the graphs that there is a response to changing spacing, and a
## response to changing number of elements. The data in the plots are
## folded, so obviously there is also a response accirding to the
## displacement and/or direction content; but for most subjects there
## is not the data in this experiment alone to distinguish them. So
## modeling content and displacement will result in an improverished
## model.  So what I'll There also isn't strong data for So what I'll
## do is

prev.descriptive.models <- NULL
descriptive.models <- NULL
library(reshape2)

##keep editing and rerunning this chunk until I'm satisfied.
{
  modelsplit <- "subject"
  testModel <- function() {
    ddply_along(
      subset(segment, abs(content) >= 0), modelsplit,
      function(group, dataset) {
        formula <- (  cbind(n_cw, n_ccw) ~
                      content:target_number_shown
                    + content:I(1/spacing)
                    + factor(side) - 1)
        #only include a displacement or content term if the data support it.
        update.if <- function(formula, update.formula) {
          updated <- update(formula, update.formula)
          m <- model.matrix(updated, dataset)
          if (qr(m)$rank == ncol(m)) {
            #cat(unlist(group), deparse(formula) , "->", deparse(updated), "\n")
            updated
          } else {
            #cat(unlist(group), deparse(formula) , "!>", deparse(updated), "\n")
            formula
          }
        }
        formula <- update.if(formula, . ~ . + displacement)
        formula <- update.if(formula, . ~ . + content)
        model <- glm(formula,
                     family=binomial(link=logit.2asym(g=0.025, lam=0.025)),
                     data=dataset)
        #and we need to check that the models are not
        #underconstrained (rank-deficient)
        mutate(group,
               rank.deficient=ncol(model.matrix(model)) - model$rank,
               model=list(model))
      })
  }
  #
  bind[prev.descriptive.models, descriptive.models] <-
    list(descriptive.models, testModel())
  if(any(descriptive.models$rank.deficient > 0)) {
    cat("Rank deficient:\n")
    chain(descriptive.models, subset(rank.deficient > 0, select=modelsplit), print)
  }
  plot(spacing.plot
       + prediction_layers(predict_from_model_frame(descriptive.models, segment), connect="number")
       + labs(title="Descriptive fits"))
  #
  if (!is.null(prev.descriptive.models) && !is.null(descriptive.models)) {
    chain(ldply(list(old = prev.descriptive.models, new = descriptive.models),
                adply, 1, summarize,
                aic = extractAIC(model[[1]])[[2]] ),
          `$<-`(model, NULL),
          melt,
          subset(variable=="aic"),
          acast(.id ~ subject),
          rbind(., change=aaply(., 2, diff)),
          cbind(total=aaply(.,1,sum), .)
          ) -> scores
    print(scores)
  }
}

#can motion energy explain NJ wobbling?
chain(segment.folded.spindled, subset(subject=="nj" && abs(direction_content==0.2))
      , print
      , add_energies
      , subset(target_number_shown == 3 & target_number_all == 9)
      , ddply(., c('target_number_shown', 'target_number_all'),
              mkchain(summarize(n=sum(n), p =), keep_columns(.)))
      , (colwise(fun(length(unique(x)))))()
      ) -> nj.energy

ggplot(nj.energy, aes())

NULL

#For the next step, I need to incorporate realistic spacing. Since we
#know that at wide spacings, there is no change in displacement
#sensitivity with number of elements (no pooling,) then we expect the
#relationship between displacement sensitivity and spacing to be
#unchanged. So let's capture that relation

##The first thing I'll do is capture some descriptive coefficients
##for each subject.

#okay first of all I'm just going to borrow the "critical spacing" and
#the "displacement" coefficients from the model data and then refit
#adding some other coefficients.

##Here's a data frame of the coefficients of my "full circle" models
circle.models <- data.frame(model=I(models), subject=names(models),
                            stringsAsFactors=FALSE)
circle.coefs <- adply(circle.models, 1, function(row) {
  bind[model=bind[model], ...=group] <-as.list(row)
  data.frame(c(coef(model), group), model=NA)
})

##here's a function to plot those predictions

##Now what I'd like is a "descriptive model" that more or less
##captures what we see in the plots above.
informed.models <-
  adply(circle.models, 1, function(row){
    bind[model=bind[model], ...=group] <- as.list(row)
    #
    #skip if there is not corresponding segment data
    if (empty(match_df(as.data.frame(group), segment, names(group))))
      return(data.frame())
    segdata <- merge(group, segment)
    #
    #Let's fix a model taking the position-discrimination as granted.
    #To do that, we'll predict the old model terms over the new data,
    #then use that as an offset in the new model.
    #
    pred <- predict(model, newdata=segdata, type="terms")
    #
    newfit <- glm((cbind(n_cw, n_ccw)
                   ~ offset(pred[,  "displacementTerm(spacing, displacement)"])
                   + bias + content:factor(side) - 1),
                  data = segdata,
                  family = binomial(link=logit.2asym(g=0.025, lam=0.025))
                  )
    c(group, list(model=newfit))
  })

##here's a function to plot those predictions
plot(spacing.plot
     + prediction_layers(predict_from_model_frame(informed.models, segment))
     + labs(title="Fits using spacing/displacement relation from Exp. 1"))

FALSE || {
  ##here's a function to plot those predictions
  plot(spacing.plot
       + prediction_layers(predict_from_model_frame(informed.models, segment)))


  ##okay first of all this can't be a way we model because that
  ##displacement coefficient 

  ##So extract model predictions for our descriptive model (and re-fold them)
  descriptive.models.pred <- predict_from_model_frame(descriptive.models, segment)

  ###graph these predictions
  plot(spacing.plot + prediction_layers(descriptive.models.pred))

  ##I think was can agree that this model captures the behavior of the
  ##"spacing experiment"_in a couple of coefficients. The question is
  ##how to relate that to the full-circle data.
  descriptive.coefs <- adply(descriptive.models, 1, function(row) {
    bind[model=bind[model], ...=group] <- as.list(row)
    data.frame(c(coef(model), group), model=NA)
  })

  ##we'll ignore the side-to-side factors, leaving four interesting
  ##coefficients -- intercepts for displacement and content, and oh,
  ##maybe I should have a coef for displacement * spacing? Would that
  ##even differ?


  ## just to be an ass, let's plot all the coefficients from one set
  ## against all the coefficients from another set.
  circle.interesting <- names(circle.coefs) %-% c("bias", "model")
  descriptive.interesting <- chain(descriptive.coefs,
                                   names,
                                   . %-% grep("^factor", ., value=TRUE),
                                   . %-% c("model"))

  ## make a long-format data frame of each and join
  library(reshape2)
  circle.coefs <- melt(circle.coefs[c(circle.interesting)], "subject")
  descriptive.coefs <- melt(descriptive.coefs[c(descriptive.interesting)], "subject")
  coef.comparison <- merge(circle.coefs, descriptive.coefs, by="subject",
                           suffixes=c(".circle", ".descriptive"))



  (ggplot(subset(coef.comparison, TRUE ))
   + aes(x=value.circle, y=value.descriptive, label=subject)
   + geom_text()
   + facet_grid(variable.descriptive ~ variable.circle, scales="free")
   + scale_x_continuous(expand=c(0.2,0.2))
   + scale_y_continuous(expand=c(0.2,0.2))
   + coord_trans(ytrans=trans_new("asinh", "asinh", "sinh"),
                 xtrans=trans_new("asinh", "asinh", "sinh"))
   )

  ## A thing I could do is plot the coefficients, of both segment
  ## models and classic models.

  segment.models[[1]]

  ## Now here's where I'll go. I'll take the values for critical distance and

  ##}}}

  ##{{{ ---------- (UNFITTED) PREDICTIONS FROM FULL CIRCLE MODEL -------------

  ## now how can we reconcile this data with the models?

  ##let's start with a simple regression on the "unfolded" data
  ##"sidedness" and "bias" are control coefficients, set to 1 always
  ##unless I'm doing something funny with the predictions.

  ## TASK 2. for checking-purposes, plot these (local) model fits
  ## against the data.

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
    attributes(x) <- NULL; x
  }

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

  print(spacing.plot + prediction_layers(segment.models.pred))

  # graph these predictions...
  print(spacing.plot + prediction_layers(mutilate.predictions(segment.pred.from.circle)))

  print(number.plot +
        prediction_layers(mutilate.predictions(segment.pred.from.circle),
                          connect="spacing"))

  ## TASK 4. how to compare and reconcile? Refit a constrained model?

  ## I think these fits are actually very good! But how do I communicate that?

  ##}}}

}
