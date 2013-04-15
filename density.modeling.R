##{{{ ---------- SETUP AND LOAD DATA -------------------------------------

## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4,
        lyx.graphics.center = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(grid)
  library(gnm)
  library(glm2)
  library(binom)
  library(psyphy)
  library(ptools)
  library(reshape2)
})
#source("latexing.R")
#source("icons.R")
source("scales.R")
source("library.R")

datafile <- "data.RData"
modelfile <- "slopeModel.RData"
plotfile <- "density.modeling.pdf"
detailfile <- "density.modeling.detail.pdf"
savefile <- "density.modeling.RData"

main <- function(datafile="data.RData", modelfile="slopeModel.RData",
                 plotfile="density.modeling.pdf",
                 detailfile = "density.modeling.detail.pdf",
                 savefile = "density.modeling.RData",
                 detailplots=!interactive()) {

  setup_theme()
  if (interactive()) {
    dev.new()
    plot.dev <- dev.cur()
    dev.new()
    detail.dev <- dev.cur()
  } else {
    cairo_pdf(plotfile, onefile=TRUE)
    plot.dev <- dev.cur()
    cairo_pdf(detailfile, onefile=TRUE)
    detail.dev <- dev.cur()
  }
  on.exit(dev.off(plot.dev), add=TRUE)
  on.exit(dev.off(detail.dev), add=TRUE)

  bind[data=data, ...=] <- as.list(load2env(datafile))
  bind[models=models, displacementTerm = displacementTerm, ...=] <- as.list(load2env(modelfile))

  #replicate the data mutation from slopeModel.R
  chain(data
        , do.rename(folding=FALSE)
        , mutate(bias=1)
        , ddply("subject",
                function(x) if ("numdensity" %in% x$exp_type) x else data.frame())
        # mutate the displacement to avoid wagon wheel (this will need done anyway)
        , mutate(data, displacement=wrap(displacement, spacing))
        ) -> data

  #ugh! globals...
  models <<- models
  displacementTerm <<- displacementTerm
  #these are the columns which define each "experiment" (facet on the
  #unfolded graph)
  segment.experiment.vars <-
    c("subject", "displacement", "content", "eccentricity")
  segment.experiment.vars <<- segment.experiment.vars
  #within an experiment these are the vars which separate each "stimulus
  #condition" (data point on the graph)
  segment.config.vars <-
    c("spacing", "target_number_shown", "target_number_all")
  segment.config.vars <<- segment.config.vars

  splits <- c(segment.config.vars, segment.experiment.vars)
  splits <<- splits

  ##Aggregate data into counts of CW and CCW responses, with various
  ##levels of folding/spindling
  bind[segment, segment.folded, segment.folded.spindled,
       segment.folded.spindled.mutilated] <-
    Map(extract_segment, list(data),
        fold = c(FALSE, TRUE, TRUE, TRUE),
        spindle = c(FALSE, FALSE, TRUE, TRUE),
        collapse = c(FALSE, FALSE, FALSE, TRUE))

  segment.trials <-
    extract_segment(data, fold=FALSE, spindle=FALSE,
                                    collapse=FALSE, count=FALSE)

  axes.basic <- list(proportion_scale[-1]
                     , spacing_texture_scale[-1]
                     , number_color_alt_scale[-1]
                     , theme(strip.text=element_text(size=8))
                     )
  #
  plot.basic <- (ggplot(segment.folded.spindled) + aes(y=p)
                 + axes.basic)
  #
  plot.wrap <- list(facet_wrap(~label))
  #
  by.number <- list(aes(x=target_number_shown),
                    geom_line(aes(group=factor(spacing),
                                  linetype=factor(spacing))),
                    geom_point(aes(group=factor(spacing),
                                  linetype=factor(spacing)), size=1))
  #
  by.spacing <- list(  aes(x=spacing), labs(x="Spacing")
                     , geom_line(aes(  group = factor(target_number_shown)
                                     , color = factor(target_number_shown)))
                     , geom_point(aes( group = factor(target_number_shown),
                                       color = factor(target_number_shown)),
                                  size=1))
  #
  by.extent <- list(aes(x = extent,
                        group = factor(target_number_shown),
                        color = factor(target_number_shown),
                        fill = factor(target_number_shown))
                    , geom_line(linetype=1)
                    , geom_line(aes(group = factor(spacing),
                                    linetype = factor(spacing)),
                                color="black", fill="black")
                    , geom_point(size=1))

  #plot with x-axis of target number, lines of constant spacing
  plot.number <- plot.basic + by.number + plot.wrap
  #
  #plot with x-axis of target spacing, lines of constant number
  plot.spacing <- plot.basic + by.spacing + plot.wrap
  plot.spacing <<- plot.spacing
  #
  #plot with x-axis of "extent"
  plot.extent <- plot.basic + by.extent + plot.wrap

  ## Let's start by descriptively modeling the segment data. We see from
  ## the graphs that there is a response to changing spacing, and a
  ## response to changing number of elements. The data in the plots are
  ## folded, so obviously there is also a response accirding to the
  ## displacement and/or direction content; but for most subjects there
  ## is not the data in this experiment alone to distinguish them. So
  ## modeling content and displacement will result in an improverished
  ## model.

  dev.set(plot.dev)
  plot(plot.spacing %+% segment.folded.spindled.mutilated
       + errorbars(segment.folded.spindled.mutilated)
       + labs(title="Number/spacing raw data"))

  #For the next step, I need to incorporate realistic spacing. Since we
  #know that at wide spacings, there is no change in displacement
  #sensitivity with number of elements (no pooling,) then we expect the
  #relationship between displacement sensitivity and spacing to be
  #unchanged. So let's extract the Rrelationships to "critical spacing"
  #and "displacement" from the slope model.

  #okay first of all I'm just going to borrow the "critical spacing" and
  #the "displacement" coefficients from the model data and then refit
  #adding some other coefficients.

  ##As a starting move, we'll use the circle model to "inform" simpler
  ##models.  That is, make predictions using the interesting nonlinear
  ##term that was fit in the circle model, then play around with the
  ##residuals.
  ##Here's a data frame of the coefficients of my "full circle" models
  circle.models <- data.frame(model=I(models), subject=names(models),
                              stringsAsFactors=FALSE)

  basic.informed.models <-
  adply(circle.models, 1, function(row) {
    bind[model=bind[model], ...=group] <- as.list(row)
    #skip if there is not corresponding segment data
    if (empty(match_df(as.data.frame(group), segment, names(group))))
      return(data.frame())
    newdata <- merge(group, segment)
    #Let's fix a model taking the position-discrimination as granted.
    #To do that, we'll predict the old model terms over the new data,
    #then use that as an offset in the new model.
    newfit <- basic_inform_model(model, newdata)
    quickdf(c(group, list(model=I(list(newfit)))))
  })

  ##here's a function to plot those predictions
  dev.set(detail.dev)
  plot(plot.spacing
       + prediction_layers(predict_from_model_frame(basic.informed.models))
       + labs(title="Predictions of spacing-causes-collapse model",
              x=paste("spacing", sep="\n",
                "(ignore the stratification with target number, that is not in this model.",
                "The point is that we got slope of response~spacing right without even trying)")
              ))

  dev.set(plot.dev)
  plot((plot.spacing %+% segment.folded.spindled.mutilated)
       + prediction_layers(predict_from_model_frame(
         basic.informed.models, fold=TRUE, spindle=TRUE, collapse=TRUE))
       + errorbars(segment.folded.spindled.mutilated)
       + labs(title="Predictions of spacing-causes-collapse model",
              x=paste("Spacing", sep="\n")))

  ##That is really cool. this gets the slope with respect to "spacing"
  ##pretty much right, with only offset terms. The slope of every line here
  ##is determined by exp. 1 and only the y-intercept is being adjusted
  ##to exp.2.

  ##Now what does that success actually tell us? It's telling us how
  ##the "content" sensitivity trades off with the spacing
  ##sensitivity. Slope with respect to spacing is an odd metric
  ##though, as it's drawing off the nonlinear term of the model.

  ##Let's think about what that means. We've captures the slope of
  ##lines of constant target number.  In the descriptive model, these
  ##slopes are determined by the term (content:I(1/spacing))

  ## First let's compare this "spacing-causes-collapse" model to a "number-causes-collapse" model.
  ##Rhetorically, we want to say that the data are more consistent
  ##with a collapse with respect to spacing than they are with a
  ##collapse of spacing with respect to number. So let's tweak the
  ##models to respond to number rather than spacing (but otherwise
  ##equal.)

  ## What we will end up finding is that the spacing model based on
  ## Experiment 1 predicts slope a lot better than the number model.

  ##Let's pull the same trick, but pretend that it's elemnent number
  ##that causes sensitivity collapse and not spacing. What should we
  ##see?
  number.informed.models <- adply(circle.models, 1, function(row) {
    bind[model=bind[model], ...=group] <- as.list(row)
    if (empty(match_df(as.data.frame(group), segment, names(group))))
      return(data.frame())
    newdata <- merge(group, segment)
    newfit <- number_inform_model(model, newdata, number.factor=2)
    quickdf(c(group, list(model=I(list(newfit)))))
  })

  dev.set(detail.dev)
  plot(plot.spacing
       + prediction_layers(predict_from_model_frame(number.informed.models))
       + labs(title="Predictions from number-causes-collapse model", x="spacing"))

  dev.set(plot.dev)
  plot(plot.spacing %+% segment.folded.spindled.mutilated
       + prediction_layers(predict_from_model_frame(
         number.informed.models,
         fold=TRUE, spindle=TRUE, collapse=TRUE))
       + labs(title="Predictions of number-causes-collapse model")
       + errorbars(segment.folded.spindled.mutilated))

  dev.set(detail.dev)
  plot(plot.extent %+% segment.folded.spindled.mutilated
       + prediction_layers(predict_from_model_frame(
         basic.informed.models,
         fold=TRUE, spindle=TRUE, collapse=TRUE))
       + labs(title=paste0("Extent plot. note increasing spacing/number\n",
                "effects at larger extents?\n",
                "(predictions are from spacing-causes-collapse model)"))
       + errorbars(segment.folded.spindled.mutilated))

  ## Compare the predictions for NJ and PBM in particular versus the
  ## spacing-causes-collapse predictions.

  # Now let's try to expand the model to describe the data we really see.
  ## This "descriptive models" is really a bit of data smoothing I'm
  ## applying, just to show the. However, it may provide a basis for
  ## comparing the number-model to the spacing-model as well.
  descriptive.models <- make_descriptive_models(segment)
  dev.set(detail.dev)
  informed.models <- ldply(descriptive.models$subject, function(subj) {
    sub <- data.frame(subject=subj, stringsAsFactors=FALSE)
    subject.data <- match_df(segment, sub, on="subject")
    bind[descriptive.model, informed.model, circle.model] <-
      lapply(list(descriptive.models, basic.informed.models, circle.models),
             function(x) match_df(x, sub)[[1,"model"]])
    #
    #here's the a plot of one subject's data without any folding and spindling
    unfolded.prediction.plot <-
      (ggplot(subset(predict_from_model(descriptive.model), content != 0))
       + axes.basic + by.spacing
       + prediction_layers(connect="number") + aes(y=fit)
       + facet_grid(content ~ side ~ displacement, labeller=pretty_strip))
    #
    #if (interactive()) figure("source")
    plot(unfolded.prediction.plot
         + labs(title=sprintf("Descriptive fits for subject %s, unfolded",
                  toupper(sub$subject))))
    #
    #if(interactive()) figure("compare")
    recast.model <- do_recast(circle.model)
    informed.recast.model <- inform_model(recast.model, recast_data(subject.data))
    plot(unfolded.prediction.plot
         %+% subset(predict_from_model(informed.recast.model), content != 0)
         + labs(title=paste0("Displacement model + global content",
                  " sum, subject ", toupper(sub$subject))))
    #
    data.frame(subject=subj, model=I(list(informed.recast.model)), row.names=subj)
  })

  #finally, plot "collapsed" predictions.
  collapsed.predictions <-
    predict_from_model_frame(informed.models,
                             fold=TRUE, spindle=TRUE, collapse=TRUE)
  dev.set(plot.dev)
  plot((plot.spacing %+% segment.folded.spindled.mutilated)
       + prediction_layers(collapsed.predictions)
       + errorbars(segment.folded.spindled.mutilated)
       + labs(title="Model fits (Experiment 1 model + global motion-energy)",
              x=paste("Spacing", sep="\n")))

  save(file=savefile, list=ls())
}

errorbars <- function(segment, x.axis="spacing") {
  ddply(segment, .(label), here(summarize),
        y = 0.5,
        x = max(
          switch(x.axis,
                 number=target_number_shown,
                 spacing=spacing,
                 extent=extent)),
        ymax = 0.5 + binom_se(min(n_obs), 0.5),
        ymin = 0.5 - binom_se(min(n_obs), 0.5)) -> errorbar
  with_arg(data = errorbar
           , inherit.aes = FALSE
           , mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax)
           , show_guide = FALSE, geom_errorbar(width=0.2)
           #, geom_point(size = 4, shape = "+")
           )
}

## we like to plot with folded data, and with the "segment" data we,
## uh, "spindle" collapsing stimuli presented in different
## hemifields. Averaging foldings and hemifields is useful for
## plotting but not as good for modeling. "fold" collapses CW and CCW
## direction contents.  "spindle" collapses stimulus locaions.
extract_segment <- function(df, fold=FALSE, spindle=FALSE, collapse=FALSE,
                            count=TRUE)
  chain(df,
        subset(exp_type=="numdensity" & subject %in% names(models)),
        do.rename(folding = FALSE), # we handle the folds more comprehensively.
        refold(fold = fold),
        if(count) mkrates(., c(  segment.config.vars, segment.experiment.vars
                               , "eccentricity", if(!spindle) "side")) else .,
        mutate(bias = if (fold) 0 else 1,
               sidedness = if (spindle) 0 else 1,
               side = if (spindle) NA else side,
               extent = spacing * target_number_shown),
        if(collapse) collapse(.) else .,
        labeler)

#Build ggplot layers to add predictions to a plot.
prediction_layers <- function(data, connect = c("number","spacing"))  {
  connect <- match.arg(connect)
  eval(template(
         with_arg(
           ...(if (missing(data)) list() else list(data=quote(data))),
           mapping=aes(
             y=fit, ymin = fit - se.fit, ymax = fit + se.fit,
             ...(switch(connect,
                        number=alist(
                          color=factor(target_number_shown),
                          fill=factor(target_number_shown)),
                        spacing=alist(
                          linetype=factor(spacing))))),
           geom_line(...(if (connect=="number") list(linetype=3) else list())),
           geom_ribbon(alpha=0.3, linetype=0))))
}

descriptive_model <- function(dataset) {
  formula <- (  cbind(n_cw, n_ccw)
              ~ content:I(1/spacing)
              - 1)
  #Some of our subjects were tested at multiple
  #displacements/contents, others not. So teh
  #displacement/content coefficient only makes sense to include
  #if the data support it:
  #cat("descriptive model", unique(dataset$subject), "\n")
  update.if <- function(formula, update.formula) {
    updated <- update(formula, update.formula)
    m <- model.matrix(updated, dataset)
    # cat(as.character(update.formula),
    #     kappa(model.matrix(formula, dataset)), '->', kappa(m), "\n")
    if (qr(m)$rank == ncol(m) && kappa(m) < 100) updated else formula
  }
  formula <- update.if(formula, . ~ . + content:target_number_shown)
  formula <- update.if(formula, . ~ . + factor(side))
  formula <- update.if(formula, . ~ . + content:factor(side))
  formula <- update.if(formula, . ~ . + displacement)
  formula <- update.if(formula, . ~ . + content)
  model <- suppress_matching_warnings(
    "truncated",
    glm2(
      formula,
      family=binomial(link=logit #logit.2asym(g=0.025, lam=0.025)
        ),
      data=dataset,
      start=descriptive_starting_values(formula, dataset)
      , maxit=100
#      , trace=TRUE
      ))
  model
}

suppress_matching_warnings <- function(pattern, expr)  {
  withCallingHandlers(
    expr,
    warning=function(w) {
      if (isTRUE(grepl(pattern, conditionMessage(w)))) {
        invokeRestart("muffleWarning")
      } else {
        warning(w)
        invokeRestart("muffleWarning")
      }
    })
}

descriptive_starting_values <- function(formula, dataset) {
  names <- colnames(model.matrix(formula, dataset))
  vapply(names, switch, 0,
         "content:target_number_shown" = 0,
         "content:displacement" = 0,
         "displacement" = 10,
         "content" =
           if (unique(dataset$subject) == "pbm") 0 else -4,
         "content:I(1/spacing)" =
           if (unique(dataset$subject) == "pbm") 10 else -10,
         0)
}

make_descriptive_models <- function(segment) {
  ddply_along(
    subset(segment, abs(content) >= 0), "subject",
    function(group, dataset) {
      model <- descriptive_model(dataset)
      #and we need to check that the models are not
      #underconstrained (rank-deficient)
      mutate(group,
             rank.deficient=ncol(model.matrix(model)) - model$rank,
             model=list(model))
    })
}

basic_inform_model <- function(model, newdata=model$data) {
  offs <- predict(model, newdata=newdata, type="terms")
  whichterm <- grep("displacementTerm", colnames(offs))
  newdata$offs <- rowSums(offs[,whichterm, drop=FALSE])
  #newdata$pred <- pred[,  whichterm] #rowSums(pred) #does this break earlier data?
  newfit <- glm(cbind(n_cw, n_ccw)
                ~ offset(offs)
                + content:factor(side) - 1
                , data = newdata
                , family = binomial(link=logit.2asym(g=0.025, lam=0.025)))
}

number_inform_model <- function(model, newdata=model$data, number.factor=1) {
  num.model <- numberize_model(model)
  #set number.factor to 2 if you think "number of targets in a hemifield"
  #is going to be more important than total number.
  newdata <- recast_data(newdata, number.factor=number.factor)
  offs <- predict(num.model, newdata=newdata, type="terms")
  whichterm <- grep("displacementTerm", colnames(offs))
  newdata$offs <- rowSums(offs[,whichterm, drop=FALSE])
  newfit <- glm(cbind(n_cw, n_ccw)
                ~ offset(offs)
                + content:factor(side) - 1
                , data = newdata
                , family=num.model$family)
}

numberize_model <- function(model, newdata=model$data) {
  #requires a new column in the dataset "number_shown_as_spacing"
  new.data <- recast_data(model$data)
  which.term <- grep("displacementTerm", labels(terms(model)))
  new.terms <- drop.terms(terms(model), which.term, keep.response=TRUE)
  new.formula <- formula(new.terms)
  new.formula <-
    update(new.formula,
           . ~ . + displacementTerm(number_shown_as_spacing, displacement,
                                    start=c(cs=4, beta_dx=14)))
  new.model <- gnm(data=new.data, formula=new.formula,
                      family=model$family)
}

do_recast <- function(model) {
  #One thing we need to to is make the main model separate its two
  #different responses to "spacing." There's two "spacing" responses;
  #the one that parameterizes the slope (which I argue should not
  #change, at least in the subject's better hemifield) and another
  #based on summation within the hemifield; and a third based on
  #"induced motion"

  #now we've split "content" into "content_local" and
  #"content_global" so let's update to reflect what we think is
  #going on. Also fill in some other columns.
  new_data <- recast_data(model$data)

  #now the simplistic model of induced motion is that it summates over
  #the whole visual field. However, we do not (in the full-circle data)
  #see more induced motion for larger numbers of elements -- my best
  #model just relates induced motion to the direction content. So
  #perhaps it is a normalization happening.

  #So it's a slight stretch to say that changing number will change
  #induced motion, but that's what I'll do. I'll relate it to the
  #"extent" (space covered) by the stimulus. This has some resonance
  #with thinking of it as a center-surround type of effect.
  new_formula <- update(model$formula, . ~ .
                        - I(1/spacing:content)
                        + content_global)

  #Refit the model (this is still to the full-circle-data. Despite
  #splitting up the variables we should have the same result (so same
  #residual deviance etc.)
  new_model <- update(model,
                      data=new_data,
                      formula=new_formula, family=model$family)
  if (deviance(model) - deviance(new_model) > 2) {
    cat("models not equivalent...\n")
    print(c(new=extractAIC(new_model), old=extractAIC(model)))
  }
  new_model
}

inform_model <- function(model, newdata=model$data) {
  pred <- predict(model, newdata=newdata, type="terms")
  newdata$pred <- rowSums(pred) #does this break earlier data?
  formula <- (cbind(n_cw, n_ccw)
                   ~ offset(pred)
                   + content:factor(side)                   + content_global
                   )
  suppress_matching_warnings(
    "truncate",
    newfit <- glm2(
      formula
      , data = newdata
      , family = binomial(link=logit #logit.2asym(g=0.025, lam=0.025)
          )
#      , start = inform_model_start(formula, newdata)
#      , trace=TRUE
      ))
}

inform_model_start <- function(formula, dataset) {
  content_factor <- switch(unique(dataset$subject),
                           pbm=10, 0)
  names <- colnames(model.matrix(formula, dataset))
  vapply(names, switch, 0,
         "content:factor(side)left" = content_factor,
         "content:factor(side)right" = content_factor,
         "content:factor(side)bottom" = content_factor,
         "content:factor(side)top" = content_factor,
         0)
}


FALSE && {
  informed.model.descriptions <- adply(basic.informed.models, 1, function(row) {
    bind[model=bind[model], ...=group] <- as.list(row)
    responses <- predict(model, type="response")
    dataset <- model$data
    dataset <- mutate(dataset, p=responses)
    if ("n_cw" %in% names(model$data)) {
      dataset <- mutate(dataset, p=responses, n_obs = n_ccw + n_cw,
                        n_cw = n_obs * p, n_ccw = n_obs * (1 - p))
    }
    description <- descriptive_model(dataset)
    quickdf(c(group, list(model=I(list(description)))))
  })

  informed.coefs <- adply(informed.model.descriptions, 1, function(row) {
    bind[model=bind[model], ...=group] <-as.list(row)
    data.frame(c(coef(model), group), model=NA)
  })

  descriptive.coefs <- adply(descriptive.models, 1, function(row) {
    bind[model=bind[model], ...=group] <-as.list(row)
    data.frame(c(coef(model), group), model=NA)
  })

  dev.set(detail.dev)
  plot(  plot.spacing
       + prediction_layers(predict_from_model_frame(informed.model.descriptions)))

  ## comparo <- merge(
  ##   descriptive.models, informed.model.descriptions,
  ##   by = ((names(informed.model.descriptions) %^% names(descriptive.models))
  ##         %-% "model"),
  ##   suffixes = c(".descriptive", ".informed"))

  ## just to be an ass, let's plot all the coefficients from one set
  ## against all the coefficients from another set.
  informed.interesting <- chain(names(informed.coefs),
                                . %-% c("bias", "model"),
                                . %-% grep("^factor", ., value=TRUE),
                                . %-% c("model"))
  descriptive.interesting <- chain(descriptive.coefs,
                                   names,
                                   . %-% grep("^factor", ., value=TRUE),
                                   . %-% c("model"))

  ## make a long-format data frame of each and join
  library(reshape2)
  circle.coefs <- melt(informed.coefs[c(informed.interesting)], "subject")
  descriptive.coefs <- melt(descriptive.coefs[c(descriptive.interesting)], "subject")
  coef.comparison <- merge(circle.coefs, descriptive.coefs, by="subject",
                           suffixes=c(".circle", ".descriptive"))

  (ggplot(subset(coef.comparison, TRUE ))
   + aes(x=value.circle, y=value.descriptive, label=subject)
   + geom_text()
   + facet_grid(variable.descriptive ~ variable.circle, scales="free")
   + scale_x_continuous(expand=c(0.2,0.2))
   + scale_y_continuous(expand=c(0.2,0.2))
   ## + coord_trans(ytrans=trans_new("asinh", "asinh", "sinh"),
   ##               xtrans=trans_new("asinh", "asinh", "sinh"))
   )

  ##Let's put that on more solid footing.Calculate "what that slope
  ##would be" for these predictions, by fitting out descriptive model to
  ##the predictions.

  ## desc.coef <- t( mapply %<<% descriptive.models %()% list(function(model, ...) {
  ##   c(list(...), obs=model$coefficients[["content:I(1/spacing)"]])
  ## }))

  ## pred.coef <- t( mapply %<<% basic.informed.models %()% list(function(model, ...) {
  ##   #predict the spacing~content slope. How? fit the descriptive model to it.
  ##   ifit <- cbind(model$data, fit=predict(model, type="link"))
  ##   newdata <- mutate(ifit, n_cw = fit*n, n_ccw = (1-fit)*n)
  ##   c(list(...), pred=descriptive_model(newdata)$coefficients[["content:I(1/spacing)"]])
  ## }))

  ##merge(as.data.frame(desc.coef), as.data.frame(pred.coef), on="subject")
}

run_as_command()
