 
## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
library(Cairo)
library(ggplot2)
library(plyr)
library(grid)
library(ptools)
theme_set(theme_bw(12, "Apple Symbols"))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())
source("latexing.R")
use_unicode <- TRUE
source("scales.R")
source("library.R")
CairoFonts(regular="Apple Symbols:style=Regular", bold="Geneva:style=Bold",
           italic="Geneva:style=Italic", bolditalic="Geneva:style=BoldItalic")

## @knitr density-load
load("../modeling/data.Rdata")
load("../modeling/slopeModel.RData")
segment <- subset(data, exp_type=="numdensity" & subject %in% names(models))

#this just illustrates the combinations of number and density.
configurations <- unique(segment[c("target_spacing", "target_number_shown")])
personalizations <- unique(segment[c("subject", "folded_displacement", "folded_direction_content")])

## Sanity check: for each personalization, check that all configurations are represented.
unmatching <-
  ddply(personalizations, .(subject, folded_displacement, folded_direction_content)
        , mkchain(  match_df(segment, ., names(.))
                  , unique(.[c("target_spacing", "target_number_shown")])
                  , merge(cbind(., .a=1), cbind(configurations, .a=1), by=names(.))
                  , subset(is.na(.a.x) | is.na(.a.y))
                  ))
 if (!empty(unmatching)) stop("unmatching data")

(ggplot(configurations)
 + aes(x=target_spacing, y=target_number_shown)
 + geom_point()
 + scale_x_continuous(breaks=unique(configurations$target_spacing),
                      labels=function(x)format(x, digits=2))
 + labs(x="Target spacing (degrees)", y="Number of targets", title="Stimulus configurations")
 )

##TODO: pick four of these to demo. What is already demoed?

## @knitr spacing-data-plots

##For each personalization, collect the "segment" data to compare it against.
alply(personalizations, 1,
      mkchain(  match_df(segment, ., on=names(.))
              , do.rename(folding=TRUE)
              , mkrates(splits=c("spacing", "content", "target_number_shown"))) 
      ) -> personalization.segment.data

chain(  personalizations
      , match_df(segment, ., on=names(.))
      , do.rename(folding=TRUE)
      , mkrates(splits=c("subject", "spacing", "content", "target_number_shown"))
      ) -> segment.data

## and let's plot them with texture and color scales... here express a
## texture and color scale....  how about building miniature
## representations of the stimuli out of dots? That might actually be
## better than using Unicode.

black_circled_number_unicode <- 0x278a:0x2793
white_circled_number_unicode <- 0x2780:0x2789


number_symbol_scale <-
  list(aes(shape=factor(target_number_shown, levels=target_number_shown, ordered=TRUE, )))

#add some binomial confidence intervals.
library(binom)
segment.data <- cbind(segment.data,
                      with(segment.data, binom.confint(n_cw, n, method="wilson") ))
## mutate(  segment_data,
##        , pmin=qbinom(0.05, n, p)/n
##        , pmax=qbinom(0.095, n, p)/n) -> segment.data

##before any more malarkey, let's make the basic graph I've been showing people all along.

number_color_scale <-
  c(  list(aes(  color=factor(target_number_shown)
               , fill=factor(target_number_shown)))
       , with_arg(name="Element\nnumber",
             palette=color_pal,
             labels=prettyprint,
             discrete_scale("fill", "manual"),
             discrete_scale("colour", "manual")
             ))
##
Map(alply(personalizations, 1, identity), personalization.segment.data,
    f=function(row, data){
      chain(data,
            with(binom.confint(n_cw, n, method="wilson", conf.level=0.75)),
            cbind(data,.)) -> data
      (  ggplot(data)
       + aes(x=spacing)
       + aes(y=p, ymin=lower, ymax=upper)
       + geom_errorbar()
       + geom_point()
 #      + spacing_texture_scale
       + proportion_scale
       + number_color_scale
       + geom_line(aes(group=spacing))
       ## + facet_wrap( ~ subject)
       + annotate(  "text", x=0, y=0
                  , label=sprintf("%s, dx=%03f", toupper(row$subject),
                      row$folded_direction_content)
                  )
       )}
    #make a gtable plotting it both ways...
    ) -> segment.raw.plots
segment.raw.plots[[2]]

## @knitr illustrated-customizations

## For this figure we describe a "configuration" from an
## experiment. First we have a spacing-series for a subject, showing a
## particular direction content, with a vertical intercept-line
## showing the particular value of direction content the cirrent data was collected under.

#to take "unfolded" data from a model and make it fold

##for each personalization, collect the "spacing" data
alply(personalizations, 1,
      splat(function(subject, freq, ...) {
            #okay one problem is that the models are phrased in terms of
            #the abs displacement whereas I want to plot the folded
            #displacement (as that is more relevant to the fix that, I
            #force the bias term to zero?)
        matchby <- data.frame(subject, ...)
        chain(  data
              , match_df(matchby, on=names(matchby) %-% "folded_displacement")
              , subset(exp_type %in% c("spacing", "content"))
              , do.rename(folding=TRUE)
              , mkrates
              )
      })) -> personalization.spacing.data

personalization.sampling <-
  seq %call% c(range(data$folded_displacement), len=100)


#Here's the question: my prediction lines are supposed to match the
#data, or the test data??? They are not necessarily equal.
Map %call% c(  personalizations
             , spacing.data=personalization.spacing.data
             , segment.data=personalization.segment.data
             , FUN = function(subject, spacing.data, folded_direction_content, ...) {
               model <- models[[subject]]
               predictions.matching.spacing.data <-
                 expand.grid(displacement=personalization.sampling
                             , spacing = unique(rawdata$spacing)
                             , content=folded_direction_content)
                predictions.matching.collected.data <-
                 expand.grid(displacement=personalization.sampling)


             })


(ggplot(rawdata)
         + proportion_scale + displacement_scale
         + balloon + spacing_color_scale + spacing_texture_scale
         + labs(y="Response proportion",
                title=paste(
                  "Subject ", toupper(subject),
                  ", direction content ",
                  format(matchby$folded_direction_content, digits=3))
                )
          + add_predictions(rawdata, models[[subject]])
         )

##------------------------------------------------------------

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

