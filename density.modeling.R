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
setup_theme()

## @knitr do-not-run
if (!interactive()) {
  cairo_pdf("density.pdf", onefile=TRUE)
  #on.exit(dev.off(), add=TRUE)
}

## @knitr density-load

load("../modeling/data.Rdata")
load("../modeling/slopeModel.RData")

#this just illustrates the combinations of number and density.
segment.config.vars <-
  c("spacing", "target_number_shown", "target_number_all")
segment.experiment.vars <-
  c("subject", "displacement", "content", "eccentricity")

## we like to plot with folded data, and with the "segment" data we
## spindle the data too, as the stimuli were presented in different
## hemifields. Averaging foldings and hemifields is useful for
## plotting but not as good for modeling.
## "fold" collapses CW and CCW direction contents.
## "spindle" collapses CW and CCW spindling.
extract_segment <- function(df, fold=FALSE, spindle=FALSE)
  chain(df,
        subset(exp_type=="numdensity" & subject %in% names(models)),
        do.rename(folding=fold),
        mkrates(c(  segment.config.vars, segment.experiment.vars
                  , "eccentricity", if(!spindle) "side")),
        mutate(bias=if(fold) 0 else 1,
               sidedness=if(spindle) 0 else 1,
               side=if(spindle) data$side[[1]] else side,
               extent = spacing * target_number_shown,
               label=sprintf("%s d=%.2g C=%.2g",
                 toupper(subject), displacement, content)))

##We'll be modeling raw data, but plotting folded/spindled. Here's how...
##
mutilate.predictions <- function(pred, fold=FALSE, spindle=FALSE) {

}

bind[segment, segment.folded, segment.spindled] <-
  Map(extract_segment, list(data),
      fold = c(FALSE, TRUE, TRUE),
      spindle = c(FALSE, FALSE,TRUE))

plot.basic <- (ggplot(segment.rates)
               + proportion_scale
               + spacing_texture_scale[-1]
               + number_color_scale[-1]
               + facet_wrap(~label)
               )

plot(extent.plot <- plot.basic
     + aes(x = extent)
     + geom_line(aes(group = factor(spacing), linetype = factor(spacing)))
     + geom_line(aes(group = factor(target_number_shown),
                     color = factor(target_number_shown)),
                 linetype=1))

plot(number.plot <- plot.basic
     + aes(x=target_number_shown)
     + geom_line(aes(group=factor(spacing), color=factor(spacing))))

plot(spacing.plot <- plot.basic
     + aes(x=spacing)
     + geom_line(aes(group = factor(target_number_shown),
                     color = factor(target_number_shown))))

## now how can we reconcile this data with the models?

##let's start with a simple regression on the "unfolded" data (there
##will need to be coefficient to account for side)

segment.models <-
  ddply_along(
    segment.rates.sided, "subject",
    function(group, dataset) {
      model = glm(
        cbind(n_cw, n_ccw) ~ spacing + displacement
        + target_number_shown
        + sidedness:factor(side),
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)),
        data=dataset)
      mutate(group, model=I(list(model)))
    })

##plot fit lines from this...
`%<-%` <- function(x, y) NULL

segment.model.pred <- adply(segment.models, 1,
                            function(row) {
                              bind[model=bind[model], ...=group] <- as.list(row)
                              data <- match_df(segment.rates, quickdf(group))
                              pred <- predict(model, newdata=data,
                                              se.fit=TRUE, type="response")
                              cbind(data, pred, model=NA) })



## TASK 2. for checking-purposes, plot these (local) model fits
## against the data.

(spacing.plot
 + geom_line(data=segment.model.pred,
             aes(y=fit, color=factor(target_number_shown)), linetype=3)
 + geom_ribbon(data=segment.model.pred,
               aes(ymin = fit - se.fit,
                   ymax = fit + se.fit,
                   fill=factor(target_number_shown)),
               alpha=0.3, linetype=0
               ))

#okay, this should fit it But I don't know what's going on with the
#PBm fit. Ah. because sensitivity to displacement changes with
#spacing, is what this simple fit left out.  The fits are also rank
#deficient, because many fits do not sample different values of
#displacement. (or perhaps I should fit on unfolded data?! Ah, I think
#that might be it too.  displacement component. So what might make
#sense for that is to make the

## TASK_3. derive some predictiion from the already fit model,
## somehow.

segment.model.pred <-
  mapply(models, subject=names(models),
         FUN=function(model, subject) {
           spacing_data <- match_df()
         })

## TASK 4. how to compare and reconcile? the

#another idea ould be to cast raw predictions from the models onto
#this experiment.

