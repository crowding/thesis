### R code from vignette source 'problem.Rnw'
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
### code chunk number 2: titrate-load
###################################################
prefixing.assign('titrate.', within(list(),{
  load("Titrate.Rdata")
  common.manipulations(environment())
}))


###################################################
### code chunk number 3: titrate-process-basic
###################################################
## the original data format is kinda messy
## extract columns that are in terms of an absolute "left" vs. "right",
## as well as folded over the local contrast
mutate(titrate.trials,
       abs.localDirectionContrast =
         sign(trial.extra.localDirection) * trial.extra.directionContrast,
       abs.displacement =
         sign(trial.extra.globalDirection) * trial.extra.globalVScalar
         * trial.extra.r * trial.motion.process.dt,
       abs.response = -result.response,
       spacing=pipe(target.spacing, ordered, factor(.,
         labels=format(as.numeric(levels(.)),
           digits=3,nsmall=2))),
       ) -> titrate.trials

#too much sign confusion here...
mutate(titrate.trials,
       folded.localDirectionContrast = abs(abs.localDirectionContrast),
       folded.displacement =
          ifelse(abs.localDirectionContrast != 0,
                 abs.displacement * sign(abs.localDirectionContrast),
                 abs.displacement * sign(abs.displacement)),
       folded.response =
         ifelse(abs.localDirectionContrast != 0,
                abs.response * sign(abs.localDirectionContrast),
                abs.response * sign(abs.displacement))
       ) -> titrate.trials

##there were some experiments where I tried contrast as the ordinate, forget those.
titrate.contrast.trials <- subset(titrate.trials, motionCondition=="local")
titrate.trials <- subset(titrate.trials, motionCondition != "local")

#Contrast exploration sessions are those days which tested only at one eccentricity.
pipe(titrate.trials
     , subset(select=c("source.file", "trial.extra.radius"))
     , ddply("source.file", c(x=length(unique(trial.extra.radius) > 1)))
     )


              )
              {
                titrate.crit.distance.trials <- subset(titrate.trials, 
                titrate.trials <- titrate.trials[,!.]
              }
                                 

#note that the "agreement" is with the direction of carrier, not envelope motion.
pipe(  titrate.trials
     , subset(responseTime >= 0.4 & responseTime <= 0.9)
     , ddply(c(  "folded.displacement"
                   , "folded.localDirectionContrast"
                   , "subject"
                   , "target.spacing"
                   , "spacing"
                   , "trial.motion.process.radius")
             , summarize
             , n.agree = sum(sign(folded.response) == 1)
             , n.disagree = sum(sign(folded.response) == -1)
             , n = sum(sign(abs.response) != 0)
             )
     , mutate(p.agree = n.agree / n, p.disagree = n.disagree / n)
     ) -> titrate.folded.rates


###################################################
### code chunk number 4: titrate-what-tested
###################################################
pipe(  titrate.trials
     , titrate.tested <- ddply(.,  c(  "subject", "trial.extra.nTargets"
                                     , "abs.localDirectionContrast")
                               , nrow)
     , (  ggplot(.) + geom_point()
        + aes(trial.extra.nTargets, abs(abs.localDirectionContrast), size=n.sessions)
        + facet_wrap(~ subject)
        + scale_x_continuous("No. of targets")
        + scale_y_continuous("Carrier motion contrast")
        + scale_size("No. of sessions")
        )
     , print
     )


###################################################
### code chunk number 5: titrate-pbm-example
###################################################
##note that I'm fitting separate GLMs when I could have fit GLM
##crossed with contrast and spacing. This is because I sampled my grid
##nonuniformly and this makes it easier to pull out sensitivity and
##slope numbers...

data.with.fits <- function(model) {
  pred <- predict(model, type="response", se.fit=TRUE)
  with(pred, data.frame(model$data, fit, se.fit))
}

pipe(  titrate.folded.rates
     , subset(  subject == "pbm"
              & folded.localDirectionContrast == 0.15
              & abs(trial.motion.process.radius-20/3) < 0.01)
     , ddply("spacing",
             mkpipe(  glm(  formula=cbind(n.agree, n.disagree) ~ folded.displacement
                          , family=binomial(link=logit))
                    , data.with.fits
                    ))
     , mutate(spacing = factor(spacing, levels=unique(spacing)))
) -> test

pipe(test,
     , (  ggplot(.) + aes(folded.displacement, p.agree, color=spacing, fill=spacing)
        + geom_point(aes(size=n)) + scale_area()
        + geom_line(aes(y=fit))
        + geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit), color=0, alpha=0.3)
        + geom_vline(x = 0, alpha=0.5)
        )
     , print
     )

###################################################
### code chunk number 6: titrate-arbitrary-pse
###################################################
print(qplot(0,0))

pipe(  titrate.folded.rates
     , subset(subject== "pbm" & folded.localDirectionContrast == 0.15)
     , dlply("spacing"
             , mkpipe(  glm (  formula=cbind(n.agree, n.disagree) ~ folded.displacement
                             , family = binomial(link=logit))
                      , predict(  data.frame(folded.displacement=seq(-0.2, 0, length.out=5))
                                , type="response", se.fits=T)
                    )
             )
        )
             


###################################################
### code chunk number 7: titrate-slopes
###################################################
ddply(titrate.folded.rates
      ,  c("target.spacing", "subject", "folded.localDirectionContrast")
      , mkpipe(  glm(  formula = cbind(n.agree, n.disagree) ~ folded.displacement
                     , family=binomial(link=logit) )
               , summary
               , `$`(coef)
               , data.frame(slope = .[2,1], bias = .[1,1],
                            slope.sd = .[2,2], bias.sd = .[1,2])
               )
      ) -> titrate.slopes


###################################################
### code chunk number 8: titrate-slopes-fig
###################################################
print(  ggplot(subset(titrate.slopes, folded.localDirectionContrast > 0.01))
      + aes(  target.spacing, slope
            , color=factor(folded.localDirectionContrast)
            , group=factor(folded.localDirectionContrast))
      + geom_point()
      + geom_line()
      + facet_wrap(~ subject)
      + geom_errorbar(aes(  ymin = slope - slope.sd
                          , ymax = slope + slope.sd))
      + scale_color_discrete("carrier contrast")
      + ylim(-10, 50)
      )


###################################################
### code chunk number 9: titrate-bias-fig
###################################################
print(  ggplot(subset(titrate.slopes, folded.localDirectionContrast > 0.01))
 + aes(target.spacing, bias, color=factor(folded.localDirectionContrast), group=factor(folded.localDirectionContrast))
 + geom_point() + geom_line()
 + facet_wrap(~ subject)
 + geom_errorbar(aes(ymin = bias - bias.sd, ymax = bias + bias.sd))
 + scale_color_discrete("carrier contrast")
 + ylim(-10, 10)
 )


