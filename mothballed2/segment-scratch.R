require(ggplot2)
require(plyr)
source("common.manipulations.R")
source("latexing.R")
source("programming.R")
theme_set(theme_bw())
theme_update(panel.grid.major=theme_blank(),         
             panel.grid.minor=theme_blank())

prefixing.assign('segment.', within(list(),{
  load("Segment.Rdata")
  common.manipulations(environment())
  trials <- mutate(trials,
                   , flanker.span=abs(sapply(trial.extra.flankerAngle,
                                           splat(`-`)))
                   ##unfortunately I screwed up the definition of
                   ##"top" and "bottom" in experiment code so I must
                   ##swap them here
                   , trial.extra.side=`levels<-`(
                       factor(trial.extra.side)
                       ,c('top','left','right','bottom'))
                   )
}))

segment.properties <- with(segment.trials,within(list(),{
  ecc <- unique(trial.extra.r)
  min.flanker.angle <- min(flanker.span)
  max.flanker.angle <- max(flanker.span)
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
  spacing <- 2*pi*trial.extra.r / trial.extra.nTargets
  spacing.breaks <- rev(unique(trial.extra.nTargets))
  spacing.labels <- format(sort(unique(spacing)),digits=1,nsmall=1)
}))


#now let's plot with a density/colormap.
segment.rates <-
  pipe(segment.trials
       , subset(responseTime >= 0.4 & responseTime <= 0.9)
       , ddply(c("subject","trial.extra.nTargets"
                 , "trial.extra.nVisibleTargets","motionCondition"
                 , "trial.extra.side")
               , summarize, correct=mean(correct), n = length(correct)
               , spacing = mean(2*pi*trial.extra.r / trial.extra.nTargets))
       , mutate( stdev = sqrt(correct*(1-correct)/n))
       )

#now make a plot of them
print(ggplot(segment.properties$tested)
      + aes(radians,trial.extra.nVisibleTargets,
            shape=selected,size=selected)
      + geom_point()
      + scale_x_continuous("Element spacing (e)")
      + scale_y_continuous("No. moving elements", breaks=3:8,labels=3:8)
      + scale_shape_manual(breaks=c(F,T),values = c(42,19))
      + scale_size_manual(breaks=c(F,T),values = c(10,2))
      + opts(legend.position = "none")
      )
                 
#I can straight away imagine several different plots...
print(ggplot(subset(segment.rates, motionCondition=='incongruent'))
      + aes(factor(spacing), factor(trial.extra.nVisibleTargets), fill=correct)
      + geom_tile()
      + facet_grid(subject ~ trial.extra.side)
      + opts(aspect.ratio=1)
      + scale_fill_gradient(name="Proportion\nlong-range", low="black", high="pink", space="Lab")
      + scale_y_discrete("No. moving elements")
      + scale_x_discrete(name="Element spacing (deg.)"
                         , breaks = 1:6
                         ##, labels=segment.properties$spacing.labels))
                         )
      )
##the left and bottom sides have the smallest crowding thresholds, but
##maybe it looks like the left and right sides have shallower slopes
##too.?  This would be consistent with the crowding finding that
##the visual field split is between upper and lower visual fields.

print(ggplot(subset(segment.rates, motionCondition=="incongruent"))
      + aes(log(spacing), correct
            , fill=factor(trial.extra.nVisibleTargets)
            , color=factor(trial.extra.nVisibleTargets)
            , ymin=correct-stdev, ymax=correct+stdev
            )
      + geom_line() + geom_ribbon(alpha=0.3, color=NA)
      + facet_grid(trial.extra.side~subject)
      )

print(ggplot(subset(segment.rates, motionCondition=="incongruent"))
      + aes(log(1/trial.extra.nVisibleTargets), correct
            , fill=factor(trial.extra.nTargets)
            , color=factor(trial.extra.nTargets)
            , ymin=correct-stdev, ymax=correct+stdev
            )
      + geom_line() + geom_ribbon(alpha=0.3, color=NA)
      + facet_grid(trial.extra.side~subject)
      + scale_color_discrete("spacing", breaks=
      )

with(segment.trials, unique(
                            data.frame(trial.extra.nTargets,
                                       trial.extra.nVisibleTargets,
                                       cluster.to.unique(flanker.span))))
geom_tile(aes(x = var))

                   ##serial dependencies in the data....
                   
#some of this data was captured using variable flanker distances, whcih I decided against.
ddply(segment.trials, .(source.file), #for each source file...
      function(x) pipe(x,
                       ddply(.,c("trial.extra.nTargets", "trial.extra.nVisibleTargets"),
                             summarize,n=nrow), #for each 
                       summarize(.,n=mean(n))))

cluster.to.unique(segment.trials$flanker.span)

#pull out the minimum and maximum flnker distance.
with(segment.trials,
     ddply(unique(data.frame(trial.extra.min.extent,
                trial.extra.max.extent,
                trial.extra.min.distance,
                       )))

with(segment.properties,
     ggplot(tested)
     + aes() + geom_point()


trials <- pipe(triggers,
               subset(name="ConcentricTrial/run/startMotion", select=c("trials.i", "next")),
               `colnames<-`(c('i','motionBegun')),
               merge(trials, all.y=TRUE))
trials <- pipe(triggers,
               subset(name %in% c("ConcentricTrial/run/cw", "ConcentricTrial/run/ccw"),
                      select=('trials.i', 'knobTime')),
               `colnames<-`(c('i', 'responseTimestamp')),
               merge(trials, all.y=TRUE))
trials <- mutate(trials, responseTime = responseTimestamp - motionBegun - trial.motion.process.t)

pipe(trials,
     subset(select=c("minResponseTime", "maxResponseTime", "subject")),
     unique)

                                                   
                                                   
pipe(trials,
     subset(trial.version...function == "ConcentricTrial"
            & abs(trial.motion.process.radius - 6.67) < .01
            & !responseTime >= 0.40
            & !responseTime <= 0.90
            , select=c("responseTime", "correct", "responseTimestamp", "motionBegun", "responseInWindow"))
)

