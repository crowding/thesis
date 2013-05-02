library(shiny)
library(psyphy)
library(ptools)
library(plyr)
library(ggplot2)
library(binom)

theme_set(theme_bw())

#here is a pSubject and a randomnumber generator
pse <- 0.4
threshold <- 0.2
pSubject <- function(input) { pnorm(input, pse, threshold) }
fSubject <- with %<<% dots(runif(length(response)) < pnorm(stimulus, pse, threshold ))
rSubject <- function(input) { runif(length(input)) < pSubject(input) }

# here is a function to obtain fitted values and confidence intervals on a GLM
predict.conf.glm <- function(model, newdata=model$data, conf.level=.9, ...) {
  linkinv <- model$family$linkinv
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
MCSdata <- mutate(data.frame(stimulus = rep(trialPoints, n)),
                  response = rSubject(stimulus))

# Fit a curve against the data
fit <- glm(response ~ stimulus, MCSdata, family=binomial(link=probit))

# Count the data points and give binomial confidence intervals for each
conf.level <- 0.9
measuredPoints <- chain(
    MCSdata
  , ddply(  "stimulus", summarize
          , n = length(response)
          , y = sum(as.logical(response))
          , response = mean(response))
  , cbind(., with(., binom.confint(y, n, conf.level, methods="exact"))))

# Also generate data to draw a curve showing the model fit
curvePoints <- chain(
  data.frame(stimulus = seq(0, 1, length = 101)),
  mutate(response = pSubject(stimulus)))

#plot data, fitted curve, and confidence intervals, with "real" function in red
(ggplot(measuredPoints)
 + aes(stimulus, response, ymin=lower, ymax=upper)
 + with_arg(data=predict.conf.glm(fit, curvePoints)
            , geom_line(color="red")
            , geom_line(aes(y=fit), linetype=2)
            , geom_ribbon(alpha=0.3, color=NA, fill="black"))
  + geom_pointrange())

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

simulated <- simulateMCS(fit, MCSdata, curvePoints)

#if we want confidence intervals, we can derive them thusly.
simulated.confint <- ddply(  simulated, "stimulus", summarize
                           , lower = quantile(response, 0.5 - conf.level/2)
                           , upper = quantile(response, 0.5 + conf.level/2)
                           , response = 0)

#From this we can make a spaghetti plot of simulated fits
plot(ggplot(simulated)
     + aes(stimulus, response)
     + geom_line(alpha=0.05, aes(group=sim))
     + with_arg(data=simulated.confint, linetype=2, color="red",
                geom_line(aes(y=lower)), geom_line(aes(y=upper)))
     )

# Does this confidence interval look similar to the one we got from
# asymptotic methods?
(ggplot(rbind.fill(mutate(predict.conf.glm(fit, curvePoints), method="asymptotic"),
                   mutate(simulated.confint, method="simulation")))
 + aes(stimulus, response, ymin=lower, ymax=upper, fill=method, color=NA, alpha=0.3)
 + geom_ribbon(alpha=0.5, color=NA))

#They're slightly skewed (which is allowed) but they both seem like
#similar coverage. We could verify the coverage of the
#asymtotic intervals against the stimulated fits, too:
coverage <- chain(
    simulated
  , merge(predict.conf.glm(fit, curvePoints), by="stimulus", suffixes=c(".sim", ".asym"))
  , ddply("stimulus", summarize
          , coverage = mean(
              response.sim >= lower & response.sim < upper)))

print(mean(coverage$coverage))
(ggplot(coverage) + aes(stimulus, coverage)
 + coord_cartesian(ylim=c(0,1)) + geom_line()
 + geom_hline(y = conf.level, linetype=2))

# That seems like the asymptotic does a good job for method of
# constant stimulus data. But what about staircases? The problemith
# staircases (or other adaptive methods) is that the samples are no
# longer independent. We model the staircase function

# We can get a confidence interval for each subject. Here's a trial
# generator. it generates a trial, asks the subject (model) for a
# random response, updates internal state, and returns that trial (as a row)

# First, implmeent a staircase. Here's a function that returns a
# staircase function.
staircase <- function(set, index = 1, nUp, nDown, resetting=FALSE) {
  len <- length(set)
  histLen <- max(nUp, nDown)
  history <- logical(histLen)
  function(up) {
    history[] <<- c(as.logical(up), history[-histLen])
    if (isTRUE(all(history[1:nDown]))) {
      if (resetting) history[] <- NA
      index <<- max(1, index-1)
    } else if (isTRUE(all(!history[1:nUp]))) {
      index <<- min(len, index+1)
      if (resetting) history[] <- NA
    }
    set[index]
  }
}

#Given a staircase, above, a list of trials, and a model to predict
#subject's response probabilities, this function uses a staircase to
#fill out a simulated stimulus strength and response.
generateData <- function(trials, predictor, generator,
                      stimulus="stimulus", response="response") {
  stim <- generator(NA) #initialize
  adply(trials, 1, function(trial) {
    trial[[1, stimulus]] <- stim
    resp <- runif(1) < predictor(trial)
    stim <<- generator(resp)
    trial[[1, response]] <- resp
    trial
  })
}

#Let's use our simulated subject and stimulate some staircase data.
#We'll use a 1-up-1-down staircase and 100 trials

#Here's the parameters of the staircase we'll use.
nTrials <- 200
staircaseArgs <- dots(set=seq(0, 1, 0.05), nUp=1, nDown=1, index=10)

#Generate a single experiment, and fit a psychometric function to it
staircaseTrials <-
  generateData(generator=staircase %()% staircaseArgs,
               data.frame(stimulus=numeric(nTrials), response=logical(nTrials)),
               fSubject)
staircaseFit <- glm(response ~ stimulus, binomial(probit), data=staircaseTrials)

#we can plot the data and the fit...
(ggplot(predict.conf.glm(staircaseFit, curvePoints) #predict.conf.glm(fit, curvePoints)
        )
 + aes(y=fit, x=stimulus)
 + geom_line(linetype=2)
 + geom_line(aes(y=response), color="red")
 + geom_ribbon(alpha=0.3, color=NA, fill="black", aes(ymax=upper, ymin=lower))
 + geom_point(data=ddply(staircaseTrials, "stimulus", summarize,
                n=length(response), response=mean(response)),
              aes(size=n, y=response))
 + coord_cartesian(ylim=c(0,1), xlim=c(0,1)))

## Note the tradeoff in using a staircase: it comes to a better idea of
## the PSE, but has a much worse idea about slope.

#Now let's take this fit and use it to stimulate and fit new
#staircases. Analogous to SimulateMCS above, this function
#re-simulates all the data from an experiment, collecting a number of
#simulated staircases.
#
#Note that each stimulation needs to have the staircase start from
#scratch; therefore the argument doesnt' take a staircase function,
#but a function that constructs staircase functions.
simulateStaircase <- function(model, nSims, constructor,
                              output=c("sim", "pred", "model"),
                              sim.data=model$data, pred.data=model$data, ...) {
  output <- match.arg(output)
  out <- lapply(seq_len(nSims), function(sim) {
    new.data <- generateData(sim.data,
                             function(trial) predict(model, newdata=trial,
                                                     type="response"),
                             generator=constructor())
    if (output != "sim") newfit <- update(model, data=new.data)
    switch(output,
           model=newfit,
           sim=mutate(new.data, sim=sim),
           pred=mutate(pred.data,
             response=predict(newfit, newdata=pred.data, type="response"),
             sim=sim))
  })
  switch(output, model=list, sim=rbind, pred=rbind) %()% out
}

staircase.sim <-
  simulateStaircase(staircaseFit, 200, output="pred", pred.data=curvePoints,
                    constructor=staircase %<<% staircaseArgs)

#And confidence intervals
staircase.confint <- ddply(  staircase.sim, "stimulus", summarize
                           , lower = quantile(response, 0.5 - conf.level/2)
                           , upper = quantile(response, 0.5 + conf.level/2)
                           , response = 0)

#From this we can make a spaghetti plot of simulated fits
plot(ggplot(staircase.sim)
     + aes(stimulus, response)
     + geom_line(alpha=0.05, aes(group=sim))
     + with_arg(data=staircase.confint, linetype=2, color="red",
                geom_line(aes(y=lower)), geom_line(aes(y=upper))))

# Does this confidence interval look similar to the one we got from
# asymptotic methods?
(ggplot(rbind.fill(mutate(predict.conf.glm(staircaseFit, curvePoints), method="asymptotic"),
                   mutate(staircase.confint, method="simulation")))
 + aes(stimulus, response)
 + geom_ribbon(alpha=0.5, aes(ymin=lower, ymax=upper, fill=method, color=method), linetype=2)
 + geom_line(data=staircase.sim, aes(group=sim), color="black", alpha=0.05))

#Oh that's weird, what is the coverage of the asymptotic fit?
coverage <- chain(
    simulated
  , merge(predict.conf.glm(fit, curvePoints), by="stimulus", suffixes=c(".sim", ".asym"))
  , ddply("stimulus", summarize
          , coverage = mean(
              response.sim >= lower & response.sim < upper)))

print(mean(coverage$coverage))
(ggplot(coverage) + aes(stimulus, coverage)
 + coord_cartesian(ylim=c(0,1)) + geom_line()
 + geom_hline(y = conf.level, linetype=2))

## Huh, the coverage is fine. So why is it skewed though?

## If you care about both the threshold and the slope of a psychometric
## function, you can use two interleaved staircases -- say, a
## 1-up-3-down and a 3-down-1-up.

## Here's a function that takes two staircase functions, and returns a
## function that alternates between them.
interleave <- function(...) {
  functions <- list(...)
  which <- 1 #which function is awaiting input
  #but we delay giving input until we get back to that staircase.
  storage <- rep(list(NA), length(functions))
  function(up) {
    storage[[which]] <<- up
    which <<- (which %% length(functions)) + 1
    functions[[which]](storage[[which]])
  }
}

##We can resimulate using two staircases:
interleavedTrials <-
  generateData(generator=(interleave(
    staircase(set=seq(0, 1, 0.05), nUp=4, nDown=1, index=10),
    staircase(set=seq(0, 1, 0.05), nUp=1, nDown=4, index=10))),
               data.frame(stimulus=numeric(nTrials), response=logical(nTrials)),
               fSubject)

interleavedFit <- glm(response~stimulus, binomial(probit), interleavedTrials)

(ggplot(predict.conf.glm(interleavedFit, curvePoints))
 + aes(y=fit, x=stimulus)
 + geom_line(linetype=2)
 + geom_line(aes(y=response), color="red")
 + geom_ribbon(alpha=0.3, color=NA, fill="black", aes(ymax=upper, ymin=lower))
 + geom_point(data=(ddply(interleavedTrials, "stimulus", summarize,
                          n=length(response), response=mean(response))),
              aes(size=n, y=response))
 + coord_cartesian(ylim=c(0,1), xlim=c(0,1)))

#Here the two staircases cluster at different points of the
#psychometric function.

interleave.sim <-
  simulateStaircase(interleavedFit, 200, output="pred", pred.data=curvePoints,
                    constructor=function() interleave(
                      staircase(set=seq(0, 1, 0.05), nUp=4, nDown=1, index=10),
                      staircase(set=seq(0, 1, 0.05), nUp=1, nDown=4, index=10)))

#And confidence intervals
interleaved.confint <- ddply(  interleave.sim, "stimulus", summarize
                             , lower = quantile(response, 0.5 - conf.level/2)
                             , upper = quantile(response, 0.5 + conf.level/2)
                             , response = 0)

#From this we can make a spaghetti plot of simulated fits
plot(ggplot(interleave.sim)
     + aes(stimulus, response)
     + geom_line(alpha=0.05, aes(group=sim))
     + with_arg(data=interleaved.confint, linetype=2, color="red",
                geom_line(aes(y=lower)), geom_line(aes(y=upper))))

# Does this confidence interval look similar to the one we got from
# asymptotic methods?
(ggplot(rbind.fill(mutate(predict.conf.glm(interleavedFit, curvePoints), method="asymptotic"),
                   mutate(interleaved.confint, method="simulation")))
 + aes(stimulus, response, ymin=lower, ymax=upper, fill=method, color=NA, alpha=0.3)
 + geom_ribbon(alpha=0.5, color=NA))
