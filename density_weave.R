## @knitr density-setup
opts_knit$set(stop_on_error=2L)
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
library(ggplot2)
library(plyr)
library(grid)
library(ptools)
source("latexing.R")
source("icons.R")
source("scales.R")
source("slopeModel.R")
source("library.R")
setup_theme()

density.example.subjects <- c("pbm", "nj")

## @knitr do-not-run
if (!interactive()) {
  cairo_pdf(commandArgs(trailingOnly=TRUE)[1], onefile=TRUE)
}

## @knitr density-load
load("data.Rdata")
load("slopeModel.RData")
segment <- chain(  data
                 , subset(exp_type=="numdensity" & subject %in% names(models))
                 , do.rename(folding=TRUE)
                 )
source("density.modeling.R")
load("numbers.RData")
load("density.modeling.RData")
for (name in ls()) {
  assign(name, get(name), globalenv())
} #coz saved fucntion

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

## @knitr density-predictions
(quad_prediction_plot(match=data.frame(subject=density.example.subjects),
                      orientation="over")
 + theme(aspect.ratio=1))

## @knitr do-not-run
print(ggplot(density.prediction.bins)
      + displacement_scale
      + proportion_scale
      + spacing_color_scale
      + aes(group=spacing)
      + geom_point(size=2)
      + ribbonf(density.prediction.curves)
      + no_grid
      + coord_cartesian(xlim=c(-0.75, 0.75))
      + geom_vline(x=-density.prediction.displacement, linetype="11"))

#combine this with the model predictions corresponding

source("density.calibration.R")

segment.plot.sided.gtables <-
  dlply_along(segment.rates.sided, segment.experiment.vars, joinedplot)
segment.plot.gtables <-
  dlply_along(segment.rates, segment.experiment.vars, joinedplot)
#grid.newpage()
#grid.draw(segment.plot.gtables[[5]])

## @knitr do-not-run
if (!interactive()) {
  lapply(segment.plot.gtables, graphics::plot)
}

## @knitr segment-diagnostics
#here is how much data I have (incl. non-incongruent trials)
print(summary(with(segment.trials, interaction(subject,trial.extra.side))))

## @knitr segment-conditions
(ggplot(segment.properties$tested)
 + aes(radians,trial.extra.nVisibleTargets,
       shape=selected,size=selected)
 + geom_point()
 + scale_x_continuous("Element spacing (e)")
 + scale_y_continuous("No. moving elements", breaks=3:8,labels=3:8)
 + scale_shape_manual(breaks=c(F,T),values = c(42,19))
 + scale_size_manual(breaks=c(F,T),values = c(10,2))
 + theme(legend.position = "none")
 ) -> segment.conditions
print(segment.conditions)

## @knitr segment-rates
chain(segment.trials
      , refold(fold=TRUE)
#      , subset(responseTime >= 0.4 & responseTime <= 0.9)
      , ddply_keeping_unique_cols(
        c(segment.splits, "side"), summarize,
        correct=mean(response), n = length(correct))
      ) -> segment.rates

## @knitr segment-colormap
(ggplot(subset(segment.rates))
 + aes(factor(spacing), factor(target_number_shown), fill=correct)
 + geom_point()
 + geom_tile()
 + scale_fill_gradient("Prop. long-range")
 + facet_grid(subject ~ side)
 + theme(aspect.ratio = 1,
         axis.text.x = element_text(angle=45))
 + scale_x_discrete("Spacing (deg.)",
                    labels=function(x) format(as.numeric(x),digits=2))
 + scale_y_discrete("No. moving elements")
 ) -> segment.colormap
print(segment.colormap)


## @knitr segment-by-spacing
(ggplot(subset(segment.rates))
 + aes(spacing, correct, color=factor(target_number_shown))
 + geom_point()
 + geom_line()
 + facet_grid(subject ~ side)
 + theme(aspect.ratio=1)
 + scale_x_continuous("Element spacing (deg)", breaks=c(2,3,4,5), labels=c(2,3,4,5))
 + scale_y_continuous("Prop. long-range")
 + scale_color_hue("No.\\\\moving\\\\elements")
 ) -> segment.by.spacing

## @knitr segment-by-elements
(ggplot(subset(segment.rates))
 + aes(target_number_shown, correct, color=factor(spacing))
 + geom_point()
 + geom_line()
 + facet_grid(subject ~ side)
 + theme(aspect.ratio=1)
 + scale_x_continuous("No. moving elements", breaks=3:8,labels=3:8)
 + scale_y_continuous("Prop. long-range")
 + scale_color_hue("Element\\\\spacing", labels=function(x) format(as.numeric(x), digits=2))
 ) -> segment.by.elements


## @knitr segment-dualplots
vp <- viewport(x=0,y=1, height=0.5, width=1, just=c("left", "top"))
print(segment.by.spacing, vp=vp)
vp <- viewport(x=0, y=0, height=0.5, just=c("left", "bottom"))
print(segment.by.elements, vp=vp)
