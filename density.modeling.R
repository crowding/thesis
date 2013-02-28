##{{{ ---------- SETUP AND LOAD DATA -------------------------------------

## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(grid)
  library(gnm)
  library(psyphy)
  library(ptools)
  library(reshape2)
})
#source("latexing.R")
#source("icons.R")
source("scales.R")
source("library.R")
source("model.R")

datafile <- "data.RData"
modelfile <- "slopeModel.RData"
plotfile <- "density.modeling.pdf"

main <- function(datafile="data.RData", modelfile="slopeModel.RData", plotfile="density.modeling.pdf") {

  setup_theme()
  if (interactive()) {
    figure("plot")
  } else {
    cairo_pdf(plotfile, onefile=TRUE)
    on.exit(dev.off(), add=TRUE)
  }

  bind[data=data, ...=] <- as.list(load2env(datafile))
  bind[models=models, ...=] <- as.list(load2env(modelfile))

  #these are the columns which define each "experiment" (facet on the
  #unfolded graph)
  segment.experiment.vars <-
    c("subject", "displacement", "content", "eccentricity")
  #within an experiment these are the vars which separate each "stimulus
  #condition" (data point on the graph)
  segment.config.vars <-
    c("spacing", "target_number_shown", "target_number_all")

  ##Aggregate data into counts of CW and CCW responses, with various
  ##levels of folding/spindling
  bind[segment, segment.folded, segment.folded.spindled,
       segment.folded.spindled.mutilated] <-
    Map(extract_segment, list(data),
        fold = c(FALSE, TRUE, TRUE, TRUE),
        spindle = c(FALSE, FALSE, TRUE, TRUE),
        collapse = c(FALSE, FALSE, FALSE, TRUE))

  axes.basic <- list(proportion_scale[-1]
                     , spacing_texture_scale[-1]
                     , number_color_alt_scale[-1]
                     , theme(strip.text=element_text(size=8))
                     )

  plot.basic <- (ggplot(segment.folded.spindled) + aes(y=p)
                 + axes.basic)

  plot.wrap <- list(facet_wrap(~label))

  by.number <- list(aes(x=target_number_shown, group=factor(spacing),
                        linetype=factor(spacing)),
                    geom_line())

  by.spacing <- list(  aes(x=spacing)
                     , geom_line(aes(  group = factor(target_number_shown)
                                     , color = factor(target_number_shown))))

  by.extent <- list(aes(x = extent,
                        group = factor(target_number_shown),
                        color = factor(target_number_shown),
                        fill = factor(target_number_shown))
                    , geom_line(linetype=1)
                    , geom_line(aes(group = factor(spacing),
                                    linetype = factor(spacing)),
                                color="black", fill="black"))

  #plot with x-axis of target number, lines of constant spacing
  plot.number <- plot.basic + by.number + plot.wrap

  #plot with x-axis of target spacing, lines of constant number
  plot.spacing <- plot.basic + by.spacing + plot.wrap

  #plot with x-axis of "extent"
  plot.extent <- plot.basic + by.extent + plot.wrap


  ## Let's start by descriptively modeling the segment data. We see from
  ## the graphs that there is a response to changing spacing, and a
  ## response to changing number of elements. The data in the plots are
  ## folded, so obviously there is also a response accirding to the
  ## displacement and/or direction content; but for most subjects there
  ## is not the data in this experiment alone to distinguish them. So
  ## modeling content and displacement will result in an improverished
  ## model.  So what I'll There also isn't strong data for So what I'll
  ## do is

  descriptive.models <- make_descriptive_models(segment)

  plot(plot.spacing
       + prediction_layers(predict_from_model_frame(descriptive.models, segment), connect="number")
       + labs(title="Descriptive fits",
              x="spacing\n(details unimportant, just capturing behavior for comparison)"))

  plot(plot.spacing %+% segment.folded.spindled.mutilated)

  #For the next step, I need to incorporate realistic spacing. Since we
  #know that at wide spacings, there is no change in displacement
  #sensitivity with number of elements (no pooling,) then we expect the
  #relationship between displacement sensitivity and spacing to be
  #unchanged. So let's extract the relationships to "critical spacing"
  #and "displacement" from the slope model.

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

  ##As a starting move, we'll use the circle model to "inform" simpler
  ##models.  That is, make predictions using the interesting nonlinear
  ##term that was fit in the circle model, then play around with the
  ##residuals.

  basic.informed.models <-
  adply(circle.models, 1, function(row){
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
  plot(plot.spacing
       + prediction_layers(predict_from_model_frame(basic.informed.models))
       + labs(title="Exp. 1 model of displacement/spacing + offset by 'content'",
              x=paste("spacing", sep="\n",
                "(ignore the stratification with target number, that is not in this model.",
                "The point is that we got slope of response~spacing right without even trying)")
              ))

  ##That is really cool. this gets the slope with respect to "spacing"
  ##pretty much right, with only offset terms. The slope of every line here
  ##is determined by exp. 1 and only the y-intercept is being adjusted
  ##to exp.2.

  ##Now what does that success actually tell us? It's telling us how
  ##the "content" sensitivity trades off with the spacing
  ##sensitivity. Slope with respect to spacing is an odd metric
  ##though, as it's drawing off the nonlinear term of the model.

  ##Let's think about what that means. We've captures the slope of lines of constant target number.
  ##In the descriptive model, these slopes are determined by the term (content:I(1/spacing))

  subj <- "nj"

  informed.models <- ldply(descriptive.models$subject, function(subj) {
    print(subj)
    sub <- data.frame(subject=subj, stringsAsFactors=FALSE)
    subject.data <- match_df(segment, sub)
    bind[descriptive.model, informed.model, circle.model] <-
      lapply(list(descriptive.models, basic.informed.models, circle.models),
             function(x) match_df(x, sub)[[1,"model"]])

    #here's the a plot of one subject's data without any folding and spindling
    unfolded.prediction.plot <-
      (ggplot(subset(predict_from_model(descriptive.model), content != 0))
       + axes.basic + by.spacing
       + prediction_layers(connect="number") + aes(y=fit)
       + facet_grid(content ~ side ~ displacement, labeller=pretty_strip))

    if (interactive()) figure("source")
    plot(unfolded.prediction.plot
         + labs(title=sprintf("Descriptive fits for subject %s, unfolded",
                              toupper(sub$subject))))

    if(interactive()) figure("compare")
    recast.model <- do_recast(circle.model)
    informed.recast.model <- inform_model(recast.model, recast_data(subject.data))
    plot(unfolded.prediction.plot
         %+% subset(predict_from_model(informed.recast.model), content != 0)
         + labs(title=paste0("Displacement model + global content",
                             " sum, subject ", toupper(sub$subject))))

    data.frame(subject=subj, model=I(list(informed.recast.model)))
  })

  #finally, plot "collapsed" predictions.
  collapsed.predictions <-
    predict_from_model_frame(informed.models,
                             fold=TRUE, spindle=TRUE, collapse=TRUE)

  plot(plot.spacing
       + prediction_layers(collapsed.predictions)
#       + labs (title="Exp. 1 model of displacement/spacing, offset by 'content'",
#        x=paste("spacing", sep="\n")
        )

}

#label function for each facet
labeler <- function(data) {
  if ("displacement" %in% names(data)) {
    mutate(data, label = sprintf("%s d=%s C=%s",
                   toupper(subject),
                   format(displacement, digits = 2),
                   format(content, digits = 2)))
  } else {
    mutate(data, label = sprintf("%s", toupper(subject)))
  }
}

## we like to plot with folded data, and with the "segment" data we,
## uh, "spindle" collapsing stimuli presented in different
## hemifields. Averaging foldings and hemifields is useful for
## plotting but not as good for modeling. "fold" collapses CW and CCW
## direction contents.  "spindle" collapses stimulus locations.
extract_segment <- function(df, fold=FALSE, spindle=FALSE, collapse=FALSE)
  chain(df,
        subset(exp_type=="numdensity" & subject %in% names(models)),
        do.rename(folding = FALSE), # we handle the folds more comprehensively.
        refold(fold = fold),
        mkrates(c(  segment.config.vars, segment.experiment.vars
                  , "eccentricity", if(!spindle) "side")),
#        nrow)
        mutate(bias = if (fold) 0 else 1,
               sidedness = if (spindle) 0 else 1,
               side = if (spindle) NA else side,
               extent = spacing * target_number_shown),
        if(collapse) collapse(.) else .,
        labeler)

##We'll be modeling raw data, but plotting folded/spindled. Here's a
##function that "re-folds" the predictions so that they can be plotted
##on a folded plot.
mutilate.predictions <-
  function(pred,
           fold=abs(diff(range(sign(pred$content)))) > 1,
           spindle=length(unique(pred$side)) > 1,
           collapse=FALSE) {
    columns <- c(as.quoted(segment.config.vars),
                 as.quoted(segment.experiment.vars),
                 if (spindle) NULL else as.quoted("side"))
    chain(pred,
          refold(fold),
          ddply_keeping_unique_cols(columns, summarize,
                fit = mean(fit), se.fit = sqrt(sum(se.fit^2))),
          if(collapse) collapse(.) else .,
          labeler)
  }

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

predict_from_model_frame <- function(models, newdata,
                                     fold=TRUE, spindle=TRUE, collapse=FALSE) {
  ##take a data frame with a list of models, and the variables to
  ##match by, produce predictions for the folding data.
  newdata_missing <- missing(newdata)
  chain(descriptive.models,
        adply(1, function(row) {
          bind[model=bind[model], ...=group] <- as.list(row)
          if (newdata_missing) {
            predict_from_model(model)
          } else {
            predict_from_model(model,
                               match_df(newdata, quickdf(group),
                                        on = names(newdata) %^% names(group)))
          }
        }),
        mutilate.predictions(fold=fold, spindle=spindle, collapse=collapse))
}

collapse <- function(data) {
  #collapses different sides and direction contents together (as for
  #most subjects in this experiment these don't matter.)
  args <- dots(
    chain(data, subset(abs(content) > 0 & displacement/sign(content) > -0.45)),
    segment.config.vars %v% segment.experiment.vars %-% c("displacement", "content"),
    summarize)

  if ("n" %in% names(data)) {
    if ("fit" %in% names(data)) {
      args <- args %__% dots(fit = mean(fit*n)/sum(n),
                             se.fit = sqrt(sum((se.fit^2)*n)/sum(n)))
    }
    args <- args %__% dots(n_ccw = sum(n_ccw), n_cw = sum(n_cw),
                           n = sum(n), p = n_cw/n)
  } else {
    if ("fit" %in% names(data))
      args <- args %__% dots(fit = mean(fit), se.fit = sqrt(mean(se.fit^2)))
  }
  ddply_keeping_unique_cols %()% args
}

predict_from_model <- function(model, newdata=model$data) {
  pred <- predict(model, newdata=newdata, se.fit=TRUE, type="response")
  newdata[(names(newdata) %in% names(pred))] <- list()
  cbind(newdata, pred, model=NA)
}

descriptive_model <- function(dataset) {
  formula <- (  cbind(n_cw, n_ccw) ~
              content:target_number_shown
              + content:I(1/spacing)
              + content:factor(side)
              + factor(side) - 1)
  #Some of our subjects were tested at multiple
  #displacements/contents, others not. So teh
  #displacement/content coefficient only makes sense to include
  #if the data support it:
  update.if <- function(formula, update.formula) {
    updated <- update(formula, update.formula)
    m <- model.matrix(updated, dataset)
    if (qr(m)$rank == ncol(m)) updated else formula
  }
  formula <- update.if(formula, . ~ . + displacement)
  formula <- update.if(formula, . ~ . + content)
  glm(formula,
      family=binomial(link=logit.2asym(g=0.025, lam=0.025)),
      data=dataset)
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
  pred <- predict(model, newdata=newdata, type="terms")
  newdata$pred <- pred[,  "displacementTerm(spacing, displacement)"] #rowSums(pred) #does this break earlier data?
  newfit <- glm(cbind(n_cw, n_ccw)
                ~ offset(pred)
                + content:factor(side) - 1
                , data = newdata
                , family = binomial(link=logit.2asym(g=0.025, lam=0.025))
                )
}



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

recast_data <- function(data) {
  mutate(data,
         eccentricity = if(!exists("eccentricity")) 20/3 else eccentricity,
         target_number_shown = (if(!exists("target_number_shown"))
                                round(2*pi*eccentricity/spacing) else
                                target_number_shown),
         target_number_all = target_number_shown,
         content_global = content * target_number_shown,
         content_local = content / spacing,
         extent = spacing * target_number_shown)
}

do_recast <- function(model) {
  #One thing we need to to is make the main model separate its two
  #different responses to "spacing." There's two "spacing" responses;
  #the one that parameterizes the slope (which I argue should not
  #change, at least in the subject's better hemifield) and another
  #based on summation within the hemifield; and a third based on
  #"induced motion"
  old_formula <- (
    cbind(n_cw, n_ccw) ~ displacementTerm(spacing, displacement) +
    content + I(content * abs(content)) + I(1/spacing):content +
    bias - 1)

  if (!identical(unattr(model$formula), unattr(old_formula)))
    stop("model formula has changed?")

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
  new_formula <-  (
    cbind(n_cw, n_ccw) ~ displacementTerm(spacing, displacement)
    + content + I(content * abs(content)) #induced motion ?
    + content_global # summation? in place of I(1/spacing:content)
    + bias - 1)

  # (content/spacing)*(spacing*target_number_shown
  #Refit the model (this is still to the full-circle-data. Despite
  #splitting up the variables we should have the same result (so same
  #residual deviance etc.)
  new_model <- update(model,
                      data=new_data,
                      formula=new_formula)
  if (deviance(model) - deviance(new_model) > 2) {
    cat("models not equivalent...\n")
    print(c(new=extractAIC(new_model), old=extractAIC(model)))
  }
  new_model
}

inform_model <- function(model, newdata=model$data) {
  pred <- predict(model, newdata=newdata, type="terms")
  newdata$pred <- rowSums(pred) #does this break earlier data?
  newfit <- glm(cbind(n_cw, n_ccw)
                ~ offset(pred)
                + content:factor(side)
                + content_global
                , data = newdata
                , family = binomial(link=logit.2asym(g=0.025, lam=0.025))
                )
}


FALSE && {
  #Now that we have that answer, let's fit a combined model...
  combined_model <- function(circle.model, descriptive.model) {
    data <- rbind.fill(recast_data(circle_model$data, recast_data(descriptive_model$data)))

    gnm(formula = (cbind(n_cw, n_ccw) ~ displacementTerm(spacing, displacement)
                   + content + I(content * abs(content))
                   + I(1/spacing):content + bias - 1)
        , family = binomial(link = logit.2asym(g = 0.025, lam = 0.025))
        , data = data)
  }

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
}

