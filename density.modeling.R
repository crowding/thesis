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


##{{{ ---------- PLOT INGREDIENTS ----------------------------------------

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

mutilate <- function(data, fold=TRUE, spindle=TRUE)

prediction_layers(...)

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

predict_from_model_frame <- function(models, newdata, fold=TRUE, spindle=TRUE) {
  ##take a data frame with a list of models, and the variables to
  ##match by, produce predictions for the folding data.
  newdata_missing <- missing(newdata)
  chain(models,
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
        mutilate.predictions(fold=fold, spindle=spindle)
        )
}

predict_from_model <- function(model, newdata=model$data) {
  pred <- predict(model, newdata=newdata, se.fit=TRUE, type="response")
  newdata[(names(newdata) %in% names(pred))] <- list()
  cbind(newdata, pred, model=NA)
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

##keep tweaking the formula until I'm satisfied.
{
  modelsplit <- "subject"
  testModel <- function() {
    ddply_along(
      subset(segment, abs(content) >= 0), modelsplit,
      function(group, dataset) {
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
  #if we did our job with update.if cirrectly above, thsi shouldn't
  #happen, but check it:
  if(any(descriptive.models$rank.deficient > 0)) {
    cat("Rank deficient:\n")
    chain(descriptive.models, subset(rank.deficient > 0, select=modelsplit), print)
  }
  plot(plot.spacing
       + prediction_layers(predict_from_model_frame(descriptive.models, segment), connect="number")
       + labs(title="Descriptive fits"))
  #
  if (!is.null(prev.descriptive.models) && !is.null(descriptive.models)) {
    #while iterating this code chunk and the model formulas, tell me
    #if my changes are improving things.  Did you know, synonyms for
    #"before/after" that are lexicographically ordered are hard to
    #come up w qith.
    chain(ldply(list(ante = prev.descriptive.models, post = descriptive.models),
                adply, 1, summarize, aic = extractAIC(model[[1]])[[2]] ),
          `$<-`(model, NULL),
          melt,
          subset(variable=="aic"),
          acast(.id ~ subject),
          rbind(., change=aaply(., 2, diff)),
          cbind(total=aaply(.,1,sum), .)
          ) -> scores
    print(scores) #lower AIC scores are better
  }
}

#can motion energy explain NJ wobbling?
## chain(segment, subset(subject=="nj" & content==0.4)
##       , add_energies
##       , ddply("displacement", summarize,
##               p=sum(n*p)/sum(n),
##               norm_diff = mean(norm_diff))
##       , ggplot
##       , (. + aes(displacement, p)
##          + geom_point() + geom_line(aes(y=norm_diff/mean(norm_diff)))))
#No, not really. It's just day-to-day variations.

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

inform_model <- function(model, newdata=model$data) {
  pred <- predict(model, newdata=newdata, type="terms")
  newdata$pred <- rowSums(pred) #does this break earlier data?
  newfit <- glm(cbind(n_cw, n_ccw)
                ~ offset(pred)
                + content + factor(side) - 1
                + content_local:extent
                , data = newdata
                , family = binomial(link=logit.2asym(g=0.025, lam=0.025))
                )
}

informed.models <-
  adply(circle.models, 1, function(row){
    bind[model=bind[model], ...=group] <- as.list(row)
    #skip if there is not corresponding segment data
    if (empty(match_df(as.data.frame(group), segment, names(group))))
      return(data.frame())
    newdata <- merge(group, segment)
    #Let's fix a model taking the position-discrimination as granted.
    #To do that, we'll predict the old model terms over the new data,
    #then use that as an offset in the new model.
    newfit <- inform_model(model, newdata)
    quickdf(c(group, list(model=I(list(newfit)))))
  })

##here's a function to plot those predictions
plot(plot.spacing
     + prediction_layers(predict_from_model_frame(informed.models))
     + labs(title="Exp. 1 model of displacement/spacing + offset by 'content'",
            x=paste("spacing", sep="\n",
              "(ignore the stratification with target number, that is not in this model.",
              "The point is that we got slope of response~spacing right without even trying)")
            ))

##what about the spacing predictions?
plot(plot.number
     + prediction_layers(predict_from_model_frame(informed.models), connect="spacing")
     + labs(title="Exp. 1 model of displacement/spacing, + offset by 'content'",
            x=paste("Target number", sep="\n",
              "(The model does not have target number anywhere so you expect this not to work)")))

##That is really cool. this gets the slope with respect to "spacing"
##exactly right, with only offset terms. The slope of every line here
##is determined by exp. 1 and only the y-intercept is being adjusted
##to exp.2.

##Now what does that success actually tell us? It's telling us how the
##"content" sensitivity trades off with the spacing sensitivity. Slope
##with respect to spacing is an odd metric though, as it's the
##nonlinear term of the model.

##Let's think about what that means. We've captures the slope of lines of constant target number.
##In the descriptive model, these slopes are determined by the term (content:I(1/spacing)) 



##  *the displacement sensitivity
##  *the degree to which displacement changes with spacing, 
##   *the degree to which summation within the critical distance influences out
##
##::: sensitivity to displacement is independent of number, which we suspected.
##::: sensitivity to content is dependent on number.

##::: he difference between our descriptive models and the informed
##models is that the descriptive models have a single term every
##descriptive model includes a term content:::target_number_shown
##which should suffice if we can explain it.
##
##The other thing is that the model has already accounted for
##summation -- as we drew predictions based on spacing.

##So the variation with target number must be in the "induced motion
##-- and independent of spacing..."

##Let's start designing a combined model.

##So here's what I'd like to do: define a "subject plot" that shows
##all model fits for the subject and the spacing.

##The graphical technique I will use is to display the distance from the
##combined models to the corrected models as a ribbon.

#POSIT: I'm having a hard time because I should just extract NJ's or
#PBM's data and work on one model at a time, eh?

sub <- data.frame(subject="pbm", stringsAsFactors=FALSE) #"nj"
subject.data <- match_df(segment, sub)
bind[descriptive.model, informed.model, circle.model] <-
  lapply(list(descriptive.models, informed.models, circle.models),
         function(x) match_df(x, sub)[[1,"model"]])

#here's the a plot of one subject's data without any folding and spindling
unfolded.prediction.plot <-
  (ggplot(subset(predict_from_model(descriptive.model), content != 0))
   + axes.basic + by.spacing
   + prediction_layers(connect="number") + aes(y=fit)
   + facet_grid(content ~ side ~ displacement, labeller=pretty_strip))

figure("source")
plot(unfolded.prediction.plot)

#and for comparison here's the "informed" predictions
#unfolded.prediction.plot %+% predict_from_model(informed.model)

#Interesting, the informed model predicts an interaction of
#displacement and spacing. That was hidden by folding? There's not
#evidence for that in the data.  Also note that NJ has an interesting
#bias on one side (which _is_ backed up in the data.)

#Now I'm going to graphically compare the "informed model" to the
#"descriptive model" plot to show (in red) where things are too low
#and (in green) where the prediction is too high.
#If I can fit both PBM and NJ it will be good.
both_predictions <- function(model1, model2, data=model1$data)
  chain(data
        , predict_from_model(model=model1)
        , rename(c(fit="fit1", se.fit="se.fit1"))
        , predict_from_model(model=model2)
        , rename(c(fit="fit2", se.fit="se.fit2"))
        , mutate(number1=target_number_shown, number2=target_number_shown)
        )

#green when upper > lower and red elsewhen
red_green_ribbon <- macro(function(lower, upper, ...) {
  template(
    list(
      with_arg(color=NA, alpha=0.2, ...(list(...))
               , geom_ribbon(  aes(ymin = pmin( .(lower), .(upper) ), ymax=.(upper)
                                   , fill=-target_number_shown
                                   )
                             #                           , fill="green"
                             )
               , geom_ribbon(  aes(ymin = .(upper), ymax = pmax( .(lower), .(upper))
                                   , fill=target_number_shown
                                   )
                             #                           , fill="red"
                             ))
      ))
#  + geom_line(aes(y=))
})

difference_plot <- function(informed.model, descriptive.model) {
  (ggplot(both_predictions(informed.model, descriptive.model))
   + proportion_scale + aes(x=spacing, group=target_number_shown, alpha=0.1)
   + theme(strip.text=element_text(size=8))
   + red_green_ribbon(fit2, fit1)
   + scale_fill_gradientn(  space="Lab", colours=c("cyan", "blue", "black", "red", "yellow")
                          , values=c(0, .3125, 0.5, 0.6875,1)
                          , breaks=c(-8:-3,3:8)
                          )
   + geom_line(aes(y=fit2))
   + facet_grid(content ~ side ~ displacement, labeller=pretty_strip)
   )
}

#difference_plot(informed.model, descriptive.model)

# whatever, we'll improve the color scheme as we come to use it.

# Next step: what we see is that the descriptive model needs to have a
# contribution from "number" that doesn't interact with the spacing.
# The following function recasts the model in these terms.

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
    + content + I(content * abs(content))
    + content_local:extent # the dilemma_these are the same
    + bias - 1)


  (content/spacing)*(spacing*target_number_shown) = content_global
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

figure("compare")
recast.model <- do_recast(circle.model)
informed.recast.model <- inform_model(recast.model, recast_data(subject.data))
unfolded.prediction.plot %+% subset(predict_from_model(informed.recast.model), content != 0)

figure("source")
unfolded.prediction.plot

#currently: there needs to be more contribution of element number with
#less change from spacing.

##With the recast model, recreate the "informed" model.
##okay first of all the "recast" model is doing just the same thing as
##the informed model. but it isn't plotting.

#oh. why the hell doesn't this plot
difference_plot(informed.recast.model, descriptive.model)

FALSE || {

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
  ## that didn't seem to do much.

  ##}}}
  ## I think these fits are actually very good! But how do I communicate that?
}

#fuckit, the answer really is to write a combined model?
combined_model <- function(circle.model, descriptive.model) {
  data <- rbind.fill(recast_data(circle_model$data, recast_data(descriptive_model$data)))

  gnm(formula = (cbind(n_cw, n_ccw) ~ displacementTerm(spacing, displacement)
                 + content + I(content * abs(content))
                 + I(1/spacing):content + bias - 1)
      , family = binomial(link = logit.2asym(g = 0.025, lam = 0.025))
      , data = data)
}

