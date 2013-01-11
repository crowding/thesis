## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
library(ggplot2)
library(plyr)
library(grid)
library(ptools)
source("latexing.R")
source("icons.R")
source("scales.R")
source("library.R")
setup_theme()

## @knitr do-not-run
cairo_pdf("density.pdf", onefile=TRUE)
#on.exit(dev.off(), add=TRUE)

## @knitr density-load
load("../modeling/data.Rdata")
load("../modeling/slopeModel.RData")
segment <- chain(  data
                 , subset(exp_type=="numdensity" & subject %in% names(models))
                 , do.rename(folding=TRUE)
                 )

#this just illustrates the combinations of number and density.
segment.config.vars <-
  c("spacing", "target_number_shown", "target_number_all")
segment.experiment.vars <-
  c("subject", "displacement", "content", "eccentricity")
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


##knitr index-plot
#choose four examples to illustrate changes of number and of density.
#lareg number/tight spacing
#med number/narrow spacing
#med number/wide spacing
#small number/wide spacing
segment.examples <-
  chain(  configurations
        , summarize(  spacing = sort(unique(spacing))[c(2, 2, 5, 5)]
                    , target_number_shown =
                         sort(unique(target_number_shown))[c(6, 3, 3, 1)]
                    , label=as.character(c(4, 3, 2, 1))
                    , eccentricity=rep(20/3, 4)
                    , color=rep(TRUE, 4)
                  )
        , merge(configurations))
configs <- merge(segment.examples, configurations, all.y=TRUE)

print(ggplot(configs)
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
             title="Stimulus configurations")
      + geom_text(aes(label=label), fontface="bold", na.rm=TRUE)
      + theme(legend.position="none"))

## @knitr spacing-data-plot

segment.rates <-
  mkrates(  segment
          , c(  segment.config.vars, segment.experiment.vars))
segment.rates.sided <-
  mkrates(  segment
          , c(  segment.config.vars, segment.experiment.vars
              , "side","eccentricity"))

#sanity check:
#I think that in each experiment the number of trials is meant to be
#the same for each condition, at least for each side. Close to the same.
#Sometimes a prematurely terminated experiment means one or two are different.
mapply(rates=list(segment.rates, segment.rates.sided),
       extra = list(c(), "side"), function(rates, extra) {
  all(unlist(dlply(  rates
                   , c(segment.experiment.vars, extra)
                   , mkchain(`$`(n), range, diff))) <= 2) || stop("oops")
})

# Compute a standard error bar over a nominal 50% rate.
binom_se <- function(n, p) sqrt(p*(1-p)/n)

##before any more malarkey, let's make the basic graph I've been
##showing people all along. It's complicated enough to need its own file.
source("density.rawplot.R")

##Now illustrate these conjointly with the matching configurations...
joinedplot <- function(...) {
  cbind(ggplot_gtable(ggplot_build(segment.plot(..., number=TRUE))),
        ggplot_gtable(ggplot_build(segment.plot(..., number=FALSE))),
        size="first")
}

##@knitr skip-for-paper
FALSE && {segment.plot.sided.gtables <-
  dlply_along(segment.rates.sided, segment.experiment.vars, joinedplot)
segment.plot.gtables <-
  dlply_along(segment.rates, segment.experiment.vars, joinedplot)
#grid.newpage()
#grid.draw(segment.plot.gtables[[5]])
}
## @knitr segment-calibration-plot

source("density.calibration.R")

combined_plot <- function(row, segment, full, full.predicted, ...) {
  spacing_breaks <- sort(unique(segment$spacing))
  spacing_range <- range(segment$spacing, full.predicted$spacing)
  upper_plot <- calibration.plot(row=row, full=full, ...,
                                 full.predicted=full.predicted,
                                 spacing_breaks=spacing_breaks,
                                 spacing_range=spacing_range)
  lower_left_plot <- segment.plot(row, segment, number=TRUE,
                                  spacing_breaks=spacing_breaks,
                                  spacing_range=spacing_range)
  browser()
  lower_right_plot <- segment.plot(row, segment, number=FALSE)
  lower_left_table <- ggplot_gtable(ggplot_build(lower_left_plot))
  lower_right_table <- ggplot_gtable(ggplot_build(lower_right_plot))
  lower_table <- cbind(lower_left_table, lower_right_table, size="first")
  #now how to mismatched bind these...I need to insert a column to the left
  #and stretch the middle
  #move one of the legend grobs over...
  upper_table <- ggplot_gtable(ggplot_build(upper_plot))
  move_legend <- upper_table[3,5][[1]][[1]]
  total_table <- chain(upper_table
                       , gtable_filter("^[^g]")
                       , gtable_add_grob(move_legend[,c(1,2,5)], 3, 5)
                       , gtable_add_cols(unit(1, "null"), pos=1)
                       , gtable_add_grob(move_legend[,c(1,4,5)], 3, 2)
                       , gtable_add_cols(rep(unit(1, "null"), 5), pos=5)
                       , `$<-`(., layout, mutate(.$layout, r=ifelse(l==5, r+5, r)))
                       , rbind(lower_table, size="last")
                       )
  total_table
}

## @knitr do-not-run

combined_plots <- lapply(calibration.data, splat(combined_plot))

lapply(combined_plots, plot)

dev.off()

##------------------------------------------------------------
FALSE && {
prefixing.assign('segment.', within(list(),{
  load("Segment.Rdata")
  common.manipulations(environment())
  mutate(  trials
         , flanker.span=abs(sapply(trial.extra.flankerAngle,
                                 splat(`-`)))
         ##unfortunately I screwed up the definition of
         ##"top" and "bottom" in experiment code so I must
         ##swap them here
         , trial.extra.side=`levels<-`(
               factor(trial.extra.side)
             , c('top','left','right','bottom'))
         , spacing = 2*pi*trial.extra.r / trial.extra.nTargets
         ) -> trials
  ##some of this data was captured using variable flanker distances,
  ##which I decided against. select only the trials using a minimum
  ##flanker spacing...
  trials <- subset(trials,
     abs(trial.extra.max.extent - trial.extra.min.extent) < .0001)
}))


## @knitr segment-properties
segment.properties <- with(segment.trials,within(list(),{
  ecc <- unique(trial.extra.r)
  min.flanker.angle <- min(flanker.span)
  max.flanker.angle <- max(flanker.span)
  min.distance <- unique(trial.extra.min.distance * trial.extra.r)
  tested <- mutate(unique(data.frame(trial.extra.nTargets,
                                     trial.extra.nVisibleTargets,
                                     trial.extra.r
                                     )),
                   spacing = 2*pi*trial.extra.r / trial.extra.nTargets,
                   radians = 2*pi/trial.extra.nTargets)
  tested <- merge(tested,
                  data.frame(trial.extra.nTargets=(c(12,21,15,15)),
                             trial.extra.nVisibleTargets = c(5,5,3,6),
                             selected=T),
                  all.x=T)
  tested <- mutate(tested, !is.na(selected))
}))


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
      + opts(legend.position = "none")
      ) -> segment.conditions
print(segment.conditions)


## @knitr segment-rates
pipe(segment.trials
     , subset(responseTime >= 0.4 & responseTime <= 0.9)
     , ddply(c("subject","trial.extra.nTargets"
               , "trial.extra.nVisibleTargets","motionCondition"
               , "trial.extra.side")
             , summarize, correct=mean(correct), n = length(correct)
             , spacing = mean(2*pi*trial.extra.r / trial.extra.nTargets))
     ) -> segment.rates


## @knitr segment-colormap
(ggplot(subset(segment.rates, motionCondition == "incongruent"))
 + aes(factor(spacing), factor(trial.extra.nVisibleTargets), fill=correct)
 + geom_point()
 + geom_tile()
 + scale_fill_gradient("Prop. long-range")
 + facet_grid(subject ~ trial.extra.side)
 + opts(aspect.ratio=1,
        axis.text.x = theme_text(angle=45))
 + scale_x_discrete("Spacing (deg.)", 
                    formatter=function(x) format(as.numeric(x),digits=2))
 + scale_y_discrete("No. moving elements")
 ) -> segment.colormap
print(segment.colormap)


## @knitr segment-by-spacing
(ggplot(subset(segment.rates, motionCondition == "incongruent"))
 + aes(spacing, correct, color=factor(trial.extra.nVisibleTargets))
 + geom_point()
 + geom_line()
 + facet_grid(subject ~ trial.extra.side)
 + opts(aspect.ratio=1)
 + scale_x_continuous("Element spacing (deg)", breaks=c(2,3,4,5), labels=c(2,3,4,5))
 + scale_y_continuous("Prop. long-range", formatter=latex.percent)
 + scale_color_hue("No.\\\\moving\\\\elements")
 ) -> segment.by.spacing


## @knitr segment-by-elements
(ggplot(subset(segment.rates, motionCondition == "incongruent"))
 + aes(trial.extra.nVisibleTargets, correct, color=factor(spacing))
 + geom_point()
 + geom_line()
 + facet_grid(subject ~ trial.extra.side)
 + opts(aspect.ratio=1)
 + scale_x_continuous("No. moving elements", breaks=3:8,labels=3:8)
 + scale_y_continuous("Prop. long-range", formatter=latex.percent)
 + scale_color_hue("Element\\\\spacing",formatter=function(x) format(as.numeric(x), digits=2))
 ) -> segment.by.elements


## @knitr segment-dualplots
vp <- viewport(x=0,y=1, height=0.5, width=1, just=c("left", "top"))
print(segment.by.spacing, vp=vp)
vp <- viewport(x=0, y=0, height=0.5, just=c("left", "bottom"))
print(segment.by.elements, vp=vp) 
}
