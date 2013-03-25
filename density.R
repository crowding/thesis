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
source("slopeModel.R")
setup_theme()

## @knitr do-not-run
if (!interactive()) {
  cairo_pdf("density1.pdf", onefile=TRUE)
  #on.exit(dev.off(), add=TRUE)
}

## @knitr density-load
load("data.Rdata")
load("slopeModel.RData")
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

## @knitr density-rates
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
# I think this is bull though
binom_se <- function(n, p) sqrt(p*(1-p)/n)

# what we are going to do is look at predictions under conditions of
# spacing-collapse and number-collapse. In the top of the figure there
# should be a sample data set varying by spacing;

density.prediction.displacement <- 0.1

#using only a subset of the spacing values here, to be less busy.
density.prediction.bins <- local({
  model <- models$nj
  chain(model$data,
        subset(abs(content) == 0.4 & target_number_shown %in% c(6, 9, 12, 16, 20, 24)),
        bin_along_resid(model, ., "response", splits, "displacement", fold=TRUE))
})

density.prediction.curves <-
  makePredictions(models$nj, density.prediction.bins, fold=TRUE)

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

##------------------------------------------------------------
#below this line is old code that I_don't know right now
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

segment.plot.sided.gtables <-
    dlply_along(segment.rates.sided, segment.experiment.vars, joinedplot)
segment.plot.gtables <-
  dlply_along(segment.rates, segment.experiment.vars, joinedplot)
#grid.newpage()
#grid.draw(segment.plot.gtables[[5]])

if (!interactive) {
  lapply(segment.plot.gtables, plot)
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
