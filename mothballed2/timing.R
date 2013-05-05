require(plyr)
require(ggplot2)

raw.response.rates <-
  ddply(subset(trials, !is.na(responseTime)
               & ((subject == "pbm")
                  | (responseTime > 0.55 & responseTime < 1.05))),
        .(subject, trial.motion.process.radius,motionCondition)
        summarize, p=mean(correct), n=length(correct)
        )


        
               

                                        #the response time plot...
load("Constant.Rdata")
source("common.manipulations.R")

common.manipulations(environment())

response.time.subset <- subset(trials,
                               abs(trial.motion.process.radius - 6.67) < .01
                               & motionCondition=="incongruent")
response.time.spacings <- with(response.time.subset, unique(target.spacing))

require(ggplot2)
(ggplot(response.time.subset)
 + aes(log(target.spacing), responseTime, color=factor(correct))
 + facet_grid(~subject)
 + geom_jitter(width=0.005, height=0)
 + geom_point(size=0.05)
 + theme_bw()
 + opts(panel.grid.major=theme_blank(),
        panel.grid.minor=theme_blank())
 + coord_cartesian(ylim=c(0.25, 1.25))
 + scale_color_discrete(name="Reported direction", breaks=c(FALSE,TRUE),labels=c("With short-range", "With long-range"))
 + scale_y_continuous("Response time (s)", breaks=c(0,0.25,0.5,0.75,1))
 + scale_x_continuous(name="Target spacing (degrees)", breaks=log(c(1,2,4,8)),labels=c(1,2,4,8))
 + opts(aspect.ratio=1)
 )

## try binning the trials according to response time, across all eccentricities
require(plyr)
trials <- mutate(trials,
                 responseTime.b05 = cut(responseTime, breaks=seq(0.3, 1, by=0.05)),
                 responseTime.b10 = cut(responseTime, breaks=seq(0.3, 1.05, by=0.25))
                 )

## now compute rates correct per response time, per subject
rt.rates <- ddply(subset(trials, motionCondition=="incongruent"),
                  .(responseTime.b05, subject),
                  summarize, p = mean(correct), n = length(correct))

#if we include this graph, add error bars (based on arcsin transform?)
(ggplot(subset(rt.rates, n >= 10 & !is.na(responseTime.b05)))
 + aes(responseTime.b05, p, color=subject, group=subject)
 + geom_line())

rt.spacing.rates <- ddply(subset(trials, motionCondition=="incongruent"),
                          .(responseTime.b10, subject, trial.extra.nTargets),
                          summarize, p = mean(correct), n = length(correct))

#I think this is not working because the 

(ggplot(subset(rt.spacing.rates, n >= 10 & !is.na(responseTime.b10)))
 + aes(trial.extra.nTargets, p, color=responseTime.b10, group=responseTime.b10)
 + facet_grid(~subject)
 + opts(aspect.ratio=1)
 + geom_line())

#New tack. How far backward can we push the RT window, until things start to go wrong?

varying.cutoffs <-
  ddply(data.frame(rt.cutoff=seq(0.2, 1, by=0.05)),
        .(rt.cutoff),
        with,
        ddply(subset(trials, responseTime >= rt.cutoff & responseTime < rt.cutoff + 0.30),
              .(subject, motionCondition,
                trial.extra.nTargets),
              summarize,
              p = mean(correct), n = length(correct)))

(ggplot(subset(varying.cutoffs, n > 10))
 + aes(rt.cutoff, p, color=trial.extra.nTargets, group=trial.extra.nTargets)
 + geom_line()
 + facet_grid(motionCondition~subject)
 )

