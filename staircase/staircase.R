library(shiny)
library(psyphy)
library(ptools)
library(plyr)
library(ggplot2)
library(binom)

#here is a pSubject and a randomnumber generator
pse <- 0.4
threshold <- 0.2
pSubject <- function(input) { pnorm(input, pse, threshold) }
rSubject <- function(input) { runif(length(input)) < pSubject(input) }

# here is a function to obtain fitted values and confidence intervals on a GLM
predict.conf.glm <- function(model, newdata=model$data, conf.level=.9, ...) {
  linkinv <- fit$family$linkinv
  chain(newdata
        , cbind(., predict(  model, newdata=newdata
                           , se.fit=TRUE, type="link"))
        , mutate(  upper = linkinv(qnorm(0.5+conf.level/2, mean=fit, sd=se.fit))
                 , lower = linkinv(qnorm(0.5-conf.level/2, mean=fit, sd=se.fit))
                 , fit = linkinv(fit)))
}

# We know how to fit a psychometric function to binary response data.
# Let's illustrate that process. Let's pretend we have a yes/no task
# where the subject's behavior is described by a cumulative normal
# function, described by a "pse" and a "threshold" (which are the mean
# and SD of the cumulative normal.)

# Say we do method of constant stimulus on 20 points
n <- 20
trialPoints <- seq(0, 1, 0.1)
MCSdata <- data.frame(stimulus = rep(points, n), response = rSubject(strength))

# Fit a curve against the data
fit <- glm(response ~ stimulus, MCSdata, family=binomial(link=probit))

# Count the data points and give binomial confidence intervals for each
conf.level <- 0.9
measuredPoints <- chain(
    data.frame(stimulus, response)
  , ddply(  "stimulus", summarize
          , n = length(response)
          , y = sum(as.logical(response))
          , response = mean(response))
  , cbind(., with(., binom.confint(y, n, conf.level, methods="exact"))))

# Also generate data to draw a curve showing the model fit
curvePoints <- chain(
    data.frame(stimulus = seq(0, 1, length = 101))
  , mutate(response = pSubject(stimulus))
  , predict.conf.glm(fit, newdata = ., conf.level = conf.level))

#plot data, fitted curve, and confidence intervals, with "real" function in red
(ggplot(measuredPoints)
 + aes(stimulus, response, ymin=lower, ymax=upper)
 + geom_pointrange()
 + with_arg(data=curvePoints
            , geom_line(color="red")
            , geom_line(aes(y=fit), linetype=2)
            , geom_ribbon(alpha=0.3, color=NA, fill="black")))

# This shows you the results of using "asymototic" methods to
# determine the confidence interval.

# A slightly more low-level way to show a confidence interval is a
# parameteric bootstrap.

library(arm) #library associated with Gelman and Hill's excellent book
             #on regression modeling

# What this does is use the covariance matrix fo the fitted model coefficients
# vcov(fit) to generate random coefficients.


# An even more low level way to extract confidence intervals do it is
# to use the fitted model to simulate the subject's yes/no responses,
# and use those si the subject's response, and use those yes/no
# responses
simulateMCS <- function(model,
                        sim.data=model$data,
                        pred.data=model$data,
                        nSims=200) {
  ldply(seq_len(nSims), function(sim) {
    sim.responses <- predict(model, newdata=sim.data, type="response")
    new.data <- mutate(sim.data, response=runif(nrow(sim.data)) < sim.responses)
    newfit <- update(model, data=new.data)
    mutate(pred.data,
           response=predict(newfit, newdata=pred.data, type="response"), sim=sim)
  })
}

theme_set(theme_bw())

#From this we can make a spaghetti plot of simulated fits
simulated <- simulateMCS(fit, MCSdata, curvePoints)
plot(sphagetti <- ggplot(simulated)
     + aes(stimulus, response, group=sim)
     + geom_line(alpha=0.05))

#if we want confidence intervals, we can derive them thusly
simulated.confint <- ddply(  simulated, "stimulus", summarize
                           , lower = quantile(response, 0.5 - conf.level/2)
                           , upper = quantile(response, 0.5 + conf.level/2)
                           , response = 0)

(ggplot(simulated) + aes(x=stimulus, y=response, ymin=lower, ymax=upper)
 + geom_ribbon(data=simulated.confint, color=NA, alpha=0.3, fill="green")
 + geom_line(alpha=0.05, aes(group=sim)))

# Does this confidence interval look similar to the one we got from asymptotic methods?
(ggplot(rbind.fill(mutate(curvePoints, method="asymptotic"),
                   mutate(simulated.confint, method="simulation")))
 + aes(stimulus, response, ymin=lower, ymax=upper, fill=method, color=NA, alpha=0.3)
 + geom_ribbon(alpha=0.5, color=NA))

#They're slightly different but they both seem like valid confidence
#intervals. We could verify the coverage of the asymtotic intervals
#against the stimulated fits, too:
coverage <- chain( simulated
  , merge(curvePoints, by="stimulus", suffixes=c(".sim", ".asym"))
  , ddply("stimulus", summarize
          , coverage = mean(
              response.sim > lower.asym
              & response.sim < upper.asym)))

print(mean(coverage$coverage))

(ggplot(coverage) + aes(stimulus, coverage)
 + coord_cartesian(ylim=c(0,1)) + geom_line())
#That seems like the asymptotic does a good job for method of constant
#stimulus data.

# This last methods of simulation can answer the question, "how does a
# staircase affect accuracy?"

# The question I want to answer is, how does the use of a staircase
# affect accuracy? In a staircase the independent variable is not
# longer independent. How does this affect the confidence we may place
# on model fits?
FALSE && {
  staircase <- function(n, rSubject, staircaseFunction, init) {
    strength = init
    strength = Reduce()
  }

  simSubject <- function(model)
    function(trial) {
      runif(nrow(trial)) < predict(model, newdata=trial, type="response")
    }

  #this function returns a function that returns a new stimulus for every response
  makeStaircase <- function(stimulusSet, index=1, responseFunction,
                            nUp = 2, nDown=1, resetting=FALSE) {
    function (response) {}
    historyLength <- max(nUp, nDown)
    history <- logical(hisotryLength)
    list(
      nextStimulus=function(response) {
        if (historyLength())
        })
    histLength
    history <- numeric(max(nUp, nDown))
    nextStim <- function(response) {
      if(response) {
        if (response[])
        } else {

        }
    }
    nextTrial <- function(response)
    }
}
# another thing I want to do here is write the "fair" way of binning
# staircase data (given model fit, derive from a sum of residuals e.g.)

