### R code from vignette source 'eyemovements.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: setup
###################################################
setCacheDir("cache")
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
require(ggplot2)
require(plyr)
theme_set(theme_bw())
theme_update(panel.grid.major = theme_blank(),
             panel.grid.minor = theme_blank())
source("common.manipulations.R")
source("latexing.R")
source("programming.R")


###################################################
### code chunk number 2: eyemovements-load
###################################################
prefixing.assign('eye.', within(list(), {
	load("Segment.eyemovements.Rdata")
}))
eye.trials <- subset(eye.trials, !is.na(responseTime))


###################################################
### code chunk number 3: eye-by-direction
###################################################

##Average eye position trace over each time step
eye.averages.raw <-
  mkpipe(subset(select="i")
         , unique
         , merge(eye.traces)
         , subset(!is.nan(eye.x))
         , ddply("eye.t", mean.and.sd)
         )

##Average eye position over each timestep,
##shifted so eye position is zero at stimulus onset
eye.averages.locked.onset <-
    mkpipe(  subset(select="i")
           , unique
           , merge(eye.traces)
           , ddply(  "i"
                   , function(d) {
                       ix <- which(d$eye.t==0)
                       numcolwise(function(x) x - x[ix]) (d)
                     }
                   )
           , subset(!is.nan(eye.x))
           , ddply("eye.t", mean.and.sd)
           )

 ##Within outer groups, compute average eye position over time,
 ##then over subgroups, compute average difference from average eye position.
 eye.averages.subgrouped <-
  function(data, subgroups) {
    pipe(  data
         , subset(select=i)
         , unique
         , merge(eye.traces)
         , ddply("eye.t", numcolwise(mean) )
         ) -> averages
    ddply(  data
          , subgroups,
          , mkpipe(  subset(select=i)
                   , merge(eye.traces)
                   , merge(averages, by="eye.t")
                   , function(data) {
                     browser()
                   }
                   , ddply("eye.t", mean.and.sd)
                   )
          )
  }

  pipe(  eye.trials
       , ddply(  c("trial.extra.side", "subject", "trial.extra.globalDirection")
               , eye.averages.raw)
       , (  ggplot(.)
          + aes(  eye.t, mean.eye.y
                , color=factor(trial.extra.globalDirection)
                , ymin = mean.eye.y - sd.eye.y/sqrt(n.eye.y)
                , ymax = mean.eye.y + sd.eye.y/sqrt(n.eye.y)
                , fill = factor(trial.extra.globalDirection)
                )
          + geom_line()
          + scale_x_continuous(lim=c(-.15, .9))
          + coord_cartesian(ylim=c(-0.5,0.5))
          + facet_grid(subject ~ trial.extra.side)
          + geom_ribbon(alpha = 0.2,linetype=0)
          + geom_vline(x = c(0,0.4), alpha=0.2)
          + opts(aspect.ratio=1)
          + scale_y_continuous(breaks=c(-0.5, 0,0.5))
          )
       )

  pipe(  eye.trials
       , subset(motionCondition == "congruent" & trial.extra.globalDirection == 1)
       , ddply(  c("trial.extra.side", "subject", "trial.extra.nVisibleTargets")
               , eye.averages.locked.onset)
       , (  ggplot(.)
          + aes(  eye.t, mean.eye.y
                , ymin = mean.eye.y - sd.eye.y/sqrt(n.eye.y)
                , ymax = mean.eye.y + sd.eye.y/sqrt(n.eye.y)
                , color= trial.extra.nVisibleTargets
                , fill = trial.extra.nVisibleTargets
                , group= trial.extra.nVisibleTargets
                )
          + geom_line()
          + scale_x_continuous(lim=c(-.15, .9))
          + coord_cartesian(ylim=c(-0.5,0.5))
          + facet_grid(subject ~ trial.extra.side)
          + geom_ribbon(alpha = 0.2,linetype=0)
          + geom_vline(x = c(0,0.4), alpha=0.2)
          + opts(aspect.ratio=1)
          + scale_y_continuous(breaks=c(-0.5, 0,0.5))
          )
       )

#what if you stack up ALL traces not just means...
pipe(eye.trials
     , subset(motionCondition == "incongruent" & trial.extra.globalDirection == 1)
     , mutate(target.spacing = factor(target.spacing))
     , subset(select=c("target.spacing", "i", "subject", "trial.extra.side"))
     , merge(eye.traces)
     , ddply(c("i"), function(d) {
       ix <- which(d$eye.t==0)
       cbind(numcolwise(function(x) x - x[ix]) (d)
             ,catcolwise(identity)(d))
     })
     , (  ggplot(.)
        + aes(eye.t, eye.y, group=i, color=factor(target.spacing))
        + geom_line()
        + scale_x_continuous(lim=c(-.15, .9))
        + coord_cartesian(ylim=c(-3,3))
        + facet_grid(subject ~ trial.extra.side)
        + scale_color_brewer(type="seq")
        )
     )


         pipe(  eye.trials
       , ddply(  c("trial.extra.side", "subject", "trial.extra.globalDirection")
               , eye.averages.raw)
       , (  ggplot(.)
          + aes(  eye.t, mean.eye.y
                , color=factor(trial.extra.globalDirection)
                , ymin = mean.eye.y - sd.eye.y/sqrt(n.eye.y)
                , ymax = mean.eye.y + sd.eye.y/sqrt(n.eye.y)
                , fill = factor(trial.extra.globalDirection)
                )
          + geom_line()
          + scale_x_continuous(lim=c(-.15, .9))
          + coord_cartesian(ylim=c(-0.5,0.5))
          + facet_grid(subject ~ trial.extra.side)
          + geom_ribbon(alpha = 0.2,linetype=0)
          + geom_vline(x = c(0,0.4), alpha=0.2)
          + opts(aspect.ratio=1)
          + scale_y_continuous(breaks=c(-0.5, 0,0.5))
          )
       )
