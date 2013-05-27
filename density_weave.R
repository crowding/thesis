## @knitr density-setup
library(knitr)
opts_knit$set(
  stop_on_error=2L)
opts_chunk$set(
  cache.extra=file.info(c(
    "data.RData", "slopeModel.RData",
    "numbers.RData", "density.modeling.RData",
    "latexing.R", "icons.R", "scales.R", "slopeModel.R",
    "density.modeling.R", "density.calibration.R", "contours.R",
    "combined.model.R", "combined.model.RData"))$mtime)
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
library(ggplot2)
library(plyr)
library(grid)
library(ptools)

source("latexing.R")
source("icons.R")
source("scales.R")
source("slopeModel.R")
source("density.modeling.R")
source("density.calibration.R")
source("contours.R")
for (name in ls()) {
  assign(name, get(name), globalenv())
} #coz saved fucntion
setup_theme()
density.example.subjects <- c("pbm", "nj")

## @knitr do-not-run
if (!interactive()) {
  cairo_pdf(commandArgs(trailingOnly=TRUE)[1], onefile=TRUE)
}

## @knitr density-load
load("data.RData")
load("slopeModel.RData")
segment <- chain(  data
                 , subset(exp_type=="numdensity" & subject %in% names(models))
                 , do.rename(folding=TRUE))
load("numbers.RData")
load("density.modeling.RData")

#this just illustrates the combinations of number and density.
segment.config.vars <-
  c("spacing", "target_number_shown", "target_number_all")
segment.experiment.vars <-
  c("subject", "displacement", "content", "eccentricity")
segment.splits <- c(segment.config.vars, segment.experiment.vars)
configurations <- unique(segment[segment.config.vars])
personalizations <- unique(segment[segment.experiment.vars])

## Sanity check: for each personalization, check that all
## configurations are represented.
unmatching <-
  ddply(  personalizations
        , segment.experiment.vars
        , mkchain(  match_df(segment, ., names(.))
                  , unique(.[segment.config.vars])
                  , merge(cbind(., .a=1), cbind(configurations, .a=1),
                          by=names(.))
                  , subset(is.na(.a.x) | is.na(.a.y))
                  ))
 if (!empty(unmatching)) stop("unmatching data")

## @knitr density-conditions
#choose four examples to illustrate changes of number and of density.
#lareg number/tight spacing
#med number/narrow spacing
#med number/wide spacing
#small number/wide spacing
 configs <-
  chain(  configurations
        , summarize(  spacing = sort(unique(spacing))[c(2, 2, 5, 5)]
                    , target_number_shown =
                         sort(unique(target_number_shown))[c(6, 3, 3, 1)]
                    , label=as.character(c(4, 3, 2, 1))
                    , color=rep(TRUE, 4)
                  )
        , merge(configurations, all.y=TRUE)
        , mutate(eccentricity = 20/3))

(ggplot(configs)
      + aes(x=spacing,
            y=target_number_shown,
            fill=color)
      + geometric_shape_scale
      + scale_fill_manual(values=c("gray80"), na.value=NA)
      + scale_x_continuous(breaks=unique(configurations$spacing),
                           labels=function(x) format(x, digits=2))
      + scale_y_continuous("Element number")
      + labs(x="Element spacing (at 6.7\u0080 eccentricity)",
             y="Number of elements",
             title="Stimulus set for Experiment 2")
      + geom_text(aes(label=label), fontface="bold", na.rm=TRUE)
      + theme(legend.position="none"))

## @knitr density-measurements
density.example.dataset <- subset(segment.folded.spindled.mutilated,
                         subject %in% density.example.subjects)
(plot.spacing %+% density.example.dataset
 + theme(aspect.ratio=1)
 + errorbars(density.example.dataset))

##plot with spacing...

## @knitr density-predictions
(condition_prediction_plot(
  quad.predictions, segment.folded.spindled.mutilated,
  match=data.frame(subject=density.example.subjects),
  orientation="over",
  conditions=quad.conditions)
 + theme(aspect.ratio=1))

## and calculate deviances
unadjusted.deviances <- chain(quad.models,
                   ddply(c("subject", names(quad.conditions)),
                         summarize,
                         deviance=vapply(model, chain, 0, extractAIC, `[`(2))),
                   acast(carrier.local ~ envelope.local ~ subject, sum,
                         value.var="deviance", margins=),
                   put(names(dimnames(.)),
                       c("carrier.local", "envelope.local", "subject")))

unadjusted.winner <- chain(unadjusted.deviances,
                           apply(names(quad.conditions), sum),
                           melt(value.name="deviance"),
                           arrange(deviance))

each.unadjusted.winner <- chain(unadjusted.deviances,
                                melt(value.name="deviance"),
                                ddply(., "subject", chain,
                                      arrange(deviance), .[1,]))

adj.deviances <- chain(adj.models,
                   ddply(c("subject", names(quad.conditions)),
                         summarize,
                         deviance=vapply(model, chain, 0, extractAIC, `[`(2))),
                   acast(carrier.local ~ envelope.local ~ subject, sum,
                         value.var="deviance", margins=),
                   put(names(dimnames(.)),
                       c("carrier.local", "envelope.local", "subject")))

adj.winner <- chain(adj.deviances,
                           apply(names(quad.conditions), sum),
                           melt(value.name="deviance"),
                           arrange(deviance))

each.adj.winner <- chain(unadjusted.deviances,
                         melt(value.name="deviance"),
                         ddply(., "subject", chain,
                               arrange(deviance), .[1,]))

## @knitr density-combined-model
source("combined.model.R")
load("combined.model.RData")

## @knitr density-combined-model-plot
source("combined.model.R")
(plot_segment_fit(
  subset(selected.fits, subject %in% density.example.subjects),
  subset(prediction.dataset, subject %in% density.example.subjects),
  subset(segment.folded.spindled.mutilated, subject %in% density.example.subjects))
 + theme(aspect.ratio=1))

## @density-extent-profile
source("combined.model.R")
(make_extent_plots(
  subset(selected.fits, subject %in% density.example.subjects),
  combined.data)
 + theme(aspect.ratio=1))
