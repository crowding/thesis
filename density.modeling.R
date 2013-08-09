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
  library(vadr)
  library(reshape2)
})
#source("latexing.R")
source("icons.R")
source("scales.R")
source("library.R")
source("slopeModel.R")
source("density_library.R")

datafile <- "data.RData"
modelfile <- "slopeModel.RData"
plotfile <- "density.modeling.pdf"
savefile <- "density.modeling.RData"

main <- function(datafile="data.RData", modelfile="slopeModel.RData",
                 plotfile="density.modeling.pdf",
                 savefile = "density.modeling.RData") {

  setup_theme()
  if (!interactive()) {
    cairo_pdf(plotfile, onefile=TRUE)
    plot.dev <- dev.cur()
  }
  on.exit(dev.off(plot.dev), add=TRUE)

  bind[data=data, ...=] <- as.list(load2env(datafile))
  bind[models=models, ...=] <- as.list(load2env(modelfile))

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

  ## Let's start by descriptively modeling the segment data. We see from
  ## the graphs that there is a response to changing spacing, and a
  ## response to changing number of elements. The data in the plots are
  ## folded, so obviously there is also a response accirding to the
  ## displacement and/or direction content; but for most subjects there
  ## is not the data in this experiment alone to distinguish them. So
  ## modeling content and displacement will result in an improverished
  ## model.

  print(plot.spacing %+% segment.folded.spindled.mutilated
       + errorbars(segment.folded.spindled.mutilated)
       + labs(title="Number/spacing raw data"))

  print(plot.extent %+% segment.folded.spindled.mutilated
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
  ##Here's a data frame of my "full circle" models
  circle.models <- data.frame(model=I(models), subject=names(models),
                              stringsAsFactors=FALSE)


  ##Let's think about what that means. We've captured the slope of
  ##lines of constant target number.  In the descriptive model, these
  ##slopes are determined by the term (content:I(1/spacing))

  ## First let's compare this "spacing-causes-collapse" model to a
  ## "number-causes-collapse" model.  Rhetorically, we want to say that
  ##the data are more consistent with a collapse with respect to
  ##spacing than they are with a collapse of spacing with respect to
  ##number. So let's tweak the models to respond to number rather than
  ##spacing (but otherwise equal.)

  ## What we will end up finding is that the spacing model based on
  ## Experiment 1 predicts slope a lot better than the number model.

  ##Let's pull the same trick, but pretend that it's elemnent number
  ##that causes sensitivity collapse and not spacing. What should we
  ##see? In fact, let's pull all four combinations_(with a fidge for
  ##guessing the field size)

  quad.conditions <- (expand.grid(carrier.local=c(TRUE, FALSE),
                                  envelope.local=c(TRUE, FALSE)))

  bind[models=quad.models, predictions=quad.predictions] <-
    recast_all_models(circle.models, conditions=quad.conditions,
                      carrier.field.guess=2, envelope.field.guess=2,
                      inform=TRUE,
                      inform.data=segment,
                      inform.fmla = . ~ .)

  print(
    condition_prediction_plot(quad.predictions,
                              segment.folded.spindled.mutilated,
                              conditions=quad.conditions))

  tall.colormap.args <- list(
    theme(legend.position="bottom",
          axis.text.x=element_text(size=rel(0.7), angle=90, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=rel(0.7), angle=90, hjust=0.5, vjust=0.5),
          panel.background=element_rect(fill="gray50"),
          strip.text.x=element_text(size=rel(0.7)),
          strip.text.y=element_text(size=rel(0.7))))

    # and plot this with a color scale too, if you care.
  print(
    condition_prediction_colormap_plot(
      quad.predictions, segment.folded.spindled.mutilated,
      conditions=quad.conditions,
      circle.properties = c(size=2, weight=0.2))
    + tall.colormap.args
    + labs(title="Predictions from Experiment 1 for Experiment 2"))

  print(
    condition_prediction_plot(
      quad.predictions, segment.folded.spindled.mutilated,
      match=data.frame(subject=c("nj", "pbm")),
      conditions=quad.conditions,
      orientation="over"))

  short.colormap.args <- list(theme(
    legend.position="bottom",
    panel.background=element_rect(fill="gray50")))

    # and plot this with a color scale too, if you care.
  print(
    condition_prediction_colormap_plot(
      quad.predictions, segment.folded.spindled.mutilated,
      conditions=quad.conditions,
      circle.properties = c(size=3, weight=0.5),
      orientation="over",
      match=data.frame(subject=c("nj", "pbm")))
    + short.colormap.args
    + labs(title="Predictions from Experiment 1 for Experiment 2"))

  ##That is really cool. envelope-local, content-global
  ##gets the slope with respect to "spacing"
  ##pretty much right, with only offset terms. The sense of slope with
  ##respoect to spacing and with respect to number is mostly in the
  ##right direction.

     ##Now what does that success actually tell us? It's telling us how
  ##the "content" sensitivity trades off with the spacing
  ##sensitivity. Slope with respect to spacing is an odd metric
  ##though, as it's drawing off the nonlinear term of the model.

  bind[models=adj.models, predictions=adj.predictions] <-
    recast_all_models(
      circle.models, conditions=quad.conditions,
      carrier.field.guess=2, envelope.field.guess=2,
      inform=TRUE, inform.data=segment,
      inform.fmla = . ~ . + content + content_global)

  print(
    condition_prediction_plot(
      adj.predictions, segment.folded.spindled.mutilated,
      conditions=quad.conditions)
    + labs(title="Predictions adjusted by carrier + other hemifield"))

  print(
    condition_prediction_plot(
      adj.predictions, segment.folded.spindled.mutilated,
      conditions=quad.conditions)
    + labs(title="Predictions adjusted by carrier + other hemifield"))

  # and plot this with a color scale too, if you care.
  print(
    condition_prediction_colormap_plot(
      adj.predictions, segment.folded.spindled.mutilated,
      conditions=quad.conditions,
      circle.properties = c(size=2, weight=0.2))
    + labs(title="Predictions adjusted by carrier + other hemifield")
    + tall.colormap.args)

  # and plot this with a color scale too, if you care.
  print(
    condition_prediction_colormap_plot(
      adj.predictions, segment.folded.spindled.mutilated,
      conditions=quad.conditions,
      circle.properties = c(size=3, weight=0.5),
      orientation="over",
      match=data.frame(subject=c("nj", "pbm")))
    + short.colormap.args
    + labs(title="Predictions adjusted by carrier + other hemifield"))

  # Now let's try to expand the model to describe the data we really see.
  ## This "descriptive models" is really a bit of data smoothing I'm
  ## applying, just to show the. However, it may provide a basis for
  ## comparing the number-model to the spacing-model as well.
  descriptive.models <- make_descriptive_models(segment)

  #tediously, print the desctiptive models against the informed ones...
  informed.models <-
    subset(adj.models, carrier.local==FALSE & envelope.local==TRUE)

  ((Map %<<% modelmerge(informed.models, descriptive.models,
                        c(".informed", ".descriptive")))
   (f = function(model.informed, model.descriptive, subject, ...) {
     #here's the a plot of one subject's data without any folding and spindling
     unfolded.prediction.plot <-
       (ggplot(subset(predict_from_model(model.descriptive), content != 0))
        + axes.basic + by.spacing
        + density_prediction_layers(connect="number") + aes(y=fit)
        + facet_grid(content ~ side ~ displacement, labeller=pretty_strip))
     #
     #if (interactive()) figure("source")
     print(unfolded.prediction.plot
           + labs(title=sprintf("Descriptive fits for observer %s, unfolded",
                    toupper(subject))))
     #
     #if(interactive()) figure("compare")
         #
    print(unfolded.prediction.plot
         %+% subset(predict_from_model(model.informed), content != 0)
         + labs(title=paste0("Displacement model + global content",
                  " sum, observer ", toupper(subject))))
   }))

   save(file=savefile, list=ls())
}

recast_all_models <- function(
  model.df,
  carrier.field.guess = 2,
  envelope.field.guess = 3,
  conditions = (expand.grid(carrier.local=c(TRUE, FALSE),
                            envelope.local=c(TRUE, FALSE))),
  inform = TRUE,
  inform.data = predict.data,
  predict.data = segment,
  inform.fmla = . ~ . + sign(content)) {
  models <-
    chain(
      model.df
      , merge(conditions, by=c(), type="full")
      , (Map %<<% .)(
        function (model, carrier.local, envelope.local, ...) {
#          if (carrier.local)
          group <- list(...)
          group.data <- merge(group, inform.data)
          if (empty(group.data)) return(data.frame())
          model <- flex_recast_model(
            model
            , carrier.local=carrier.local, envelope.local=envelope.local
            , carrier.factor=carrier.field.guess
            , envelope.factor=envelope.field.guess
            , inform = inform
            , inform.data = group.data
            , inform.fmla = inform.fmla)
          quickdf(list(
            ..., model=list(model)
            , carrier.local=carrier.local, envelope.local=envelope.local))
        })
      , rbind %()% .
      , asisify)
  if (!inform) {
    inform.data <- recast_data(inform.data, number.factor=envelope.field.guess,
                               carrier.factor=carrier.field.guess)
    prediction_function <- curr(predict_from_model_frame, newdata=inform.data)
  } else {
    prediction_function <- predict_from_model_frame
  }
  predictions <- chain(
    models,
        ddply(c("carrier.local", "envelope.local"),
          prediction_function,
          fold=TRUE, spindle=TRUE, collapse=TRUE)
    )
  list(models=models, predictions=predictions)
}

condition_prediction_plot <- function(predictions, data, match, conditions,
                                      orientation = c("down", "over")) {
  orientation <- match.arg(orientation)
  facet.fmla <- switch(orientation,
                       down=subject ~ carrier.local + envelope.local,
                       over=carrier.local + envelope.local ~ subject)
  if (!is.missing(match)) {
    data <- merge(data, match)
    predictions <- merge(predictions, match)
  }
  chain(conditions
        , expanded.data=merge(data)
        , plot.spacing %+% .
        , +facet_grid(facet.fmla, labeller=condition_facet_labeller)
        , +density_prediction_layers(predictions)
        , +labs(title="Predictions for Experiment 2 from Experiment 1")
        , +errorbars(expanded.data,
                     facet=c("carrier.local", "envelope.local", "subject")))
}

condition_prediction_colormap_plot <- function(
  predictions, data, match, conditions, orientation = c("down", "over"),
  circle.properties=c()) {
  orientation <- match.arg(orientation)
  facet.fmla <- switch(orientation,
                       down=subject ~ carrier.local + envelope.local,
                       over=carrier.local + envelope.local ~ subject)
  if (!missing("match")) {
    data <- merge(data, match)
    predictions <- merge(predictions, match)
  }
  expanded.data <- merge(conditions, data)
  (ggplot(predictions)
   + aes(y=factor(target_number_shown),
         x = factor(spacing,
           levels=sort(unique(spacing)),
           labels=format(sort(unique(spacing)), digits=2)))
   + scale_fill_gradientn("P(Response CW)", colours=hex(colorful.gradient))
   + scale_y_discrete("Element number")
   + scale_x_discrete("Spacing")
   + geom_tile(aes(fill=fit), color=NA)
   + facet_grid(facet.fmla, labeller=condition_facet_labeller)
   + (geom_circle %<<% circle.properties)(data=data,
                        aes(fill=p), color="green", linetype="11")
   + no_grid
   + theme(aspect.ratio=1)
 )
}

condition_facet_labeller <- function(var, value) {
  switch(
    var,
    carrier.local=
    paste("Carrier", ifelse(value, "local", "global")),
    envelope.local=
    paste("Envelope", ifelse(value, "local", "global")),
    subject=paste("Observer", toupper(value)))
}

modelmerge <- function(
  x, y,
  suffixes=paste0(".", c(deparse(substitute(x)), deparse(substitute(y))))) {
  force(suffixes)
  merge(x, y, by=intersect(names(x), names(y)) %-% c("model"), suffixes=suffixes)
}

descriptive_model <- function(dataset) {
  formula <- (  cbind(n_cw, n_ccw)
              ~ content:I(1/spacing)
              - 1)
  #Some of our subjects were t ested at multiple
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

flex_recast_model <- function(model,
                              carrier.local=FALSE,
                              envelope.local=TRUE,
                              carrier.factor=2, envelope.factor=2,
                              inform=FALSE,
                              inform.fmla=.~.+content,
                              inform.data=model$data) {
  #One thing we need to to is make the main model separate its two
  #different responses to "spacing." There's two "spacing" responses;
  #the one that parameterizes the slope (which I argue should not
  #change, at least in the subject's better hemifield) and another
  #based on summation within the hemifield; and a third based on
  #"induced motion"
  model.data <- recast_data(model$data,
                            envelope.factor=envelope.factor,
                            carrier.factor=carrier.factor)
  #now we've split "content" into "content_local" and
  #"content_global" so let's update to reflect what we think is
  #going on.
  fmla <- model$formula
  i.carrier.term <- grep("content.*/spacing", labels(terms(fmla)))
  new.terms <- drop.terms(terms(fmla), i.carrier.term, keep.response=TRUE)
  fmla <- formula(new.terms)
  if (carrier.local) {
    fmla <- update(fmla, .~.+ content_local)
  } else {
    fmla <- update(fmla, . ~ . + content_global)
  }
  if (!envelope.local) {
    fmla <- as.formula(substituteDirect(fmla, alist(spacing=number_shown_as_spacing)))
  }
  #Refit the model (this is still to the full-circle-data.) Despite
  #splitting up the variables we should have the same result (so same
  #residual deviance etc.)
  new.model <- gnm(data=model.data, formula=fmla, family=model$family)
  #deviance should be equal....
  if (deviance(model) - deviance(new.model) > 2) {
    cat("models not equivalent...\n")
    print(c(new=extractAIC(new_model), old=extractAIC(model)))
  }
  if(inform) {
    new.data <- recast_data(inform.data,
                            envelope.factor=envelope.factor,
                            carrier.factor=carrier.factor)
    new.data$pred <- predict(new.model, newdata=new.data, type="link")
    fmla <- cbind(n_cw, n_ccw) ~ offset(pred)
    fmla <- update(fmla, inform.fmla)
    new.model <- glm(fmla, data=new.data, #family=binomial(link=logit)
                     family=new.model$family
                     )
  }
  new.model
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

  plot(  plot.spacing %+% segment.folded.spindled.mutilated
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
