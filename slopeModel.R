suppressPackageStartupMessages({
library(plyr)
library(ggplot2)
library(ptools)
library(psyphy)
library(gnm)
})
theme_set(theme_bw())
use_unicode=FALSE
source("scales.R")

match_df <- function(...) suppressMessages(plyr::match_df(...))

do.rename <- function(data, folding=TRUE) {
  replacements <- if (folding) {
    c(folded_direction_content="content",
      folded_displacement="displacement",
      folded_response_with_carrier="response",
      target_spacing="spacing"
      )
  } else {
    c(abs_direction_content="content",
      abs_displacement="displacement",
      abs_response_cw="response",
      target_spacing="spacing"
      )
  }
  rename(data, replacements)
}

rates <- mkchain(ddply(c("displacement", "content",
                         "spacing", "subject", "exp_type"),
                       summarize, n = length(response), p = mean(response)),
                 arrange(desc(n)))

seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)

`%++%` <- function(x, y) paste(x, y, sep="")

#Our real model has nonlinear terms in the spacing-dependent
#sensitivities to displacement and to direction content.
#This devfines the response to displacement.
displacementTerm = function(spacing, displacement) {
  #for an intro to "gnm" custom terms see:
  #http://statmath.wu.ac.at/research/friday/resources_WS0708_SS08/gnmTalk.pdf
  #
  #"variables" are in data. "predictors" are model parameters we introduce.
  spacing = substitute(spacing)
  displacement = substitute(displacement)
  list(  predictors = alist(cs=1, beta_dx=1) #not sure why "1"
       , variables = list(spacing, displacement)
       , term = function(predLabels, varLabels) {
         t = parse(text=c(predLabels, varLabels));
         names(t) = c("cs", "beta_dx", "spacing", "displacement")
         deparse(bquote( (2 - 2/(1 + exp( -.(t$cs) / .(t$spacing) ) ) )
                        * .(t$beta_dx) * .(t$displacement) ))
       }
       , start = function(predictors) {
         #spacing, displacement starting points.
         predictors[c(1,2)] <- c(3, 14)
         predictors
       }
       )
}
class(displacementTerm) <- "nonlin"

bmakePredictions <-
  function(model, data=model$data,
           splitting_vars=c("subject", "content", "exp_type", "spacing"),
           ordinate = "displacement"
           ) {
    r <- range(data[ordinate])
    sampling <- merge(unique(data[splitting_vars]),
                      quickdf(structure(list(seq(r[1], r[2], length=100)),
                                        names=ordinate)),
                      by.x=c(), by.y=c())
    pred <- predict(model, newdata=sampling, type="response", se.fit=TRUE)
    cbind(sampling, pred[1:2])
    #one can also do a prediction from the terms...
  }

plotPredictions <- function(...) {
  predictions = makePredictions(...)
  with_arg(data=predictions,
           geom_ribbon(color="transparent", alpha=0.2,
                       aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
           geom_line(aes(y=fit)))
}

#then there are the two components to the "wobble", with a fixed scale
#(for now); they are the
wobble1 <- function(x) 1./(exp(x)+1)
wobble3 <- function(x) exp(x)*(exp(x)-1)/(exp(x)+1)^3
wobbleScale <- 0.3

plot_page <- function(model, subject) {
  rates <- match_df(rates, data.frame(subject=subject), on="subject")
  print((ggplot(rates)
         + displacement_scale
         + proportion_scale
         + content_color_scale
         + facet_spacing_experiment
         + plotPredictions(model, rates)
         + balloon
         + labs(title = "Data and model fits for subject " %++% subject)
         ))
}

main <- function(
          infile = "data.RData"
          ) {

  load(infile, e <- new.env())

  chain(e$data
        , do.rename(folding=FALSE)
        , subset(exp_type %in% c("spacing", "content"))
        , match_df(., subset(count(., "subject"), freq>2000), on="subject")
        ) -> data

  rates = rates(data)

  #replicate the first slope change model...
  model <- gnm(  (  response
                  ~ displacementTerm(spacing, displacement)
                  + content
                  + I(1/spacing):content
                  )
               , family=binomial(link=logit.2asym(g=0.025, lam=0.025))
               , data=subset(data, subject=="pbm")
               )

  plot_page(model, subject="pbm")

  ## %the "induced motion* comes in two types, which I'll fit with a
  ## %logistic plus the third derivative of a logistic (no real
  ## %justification here other than hte combination looks loke
  ## %the data.) The scale parameter of the logit is chosen
  ## %arbitrarily, not fit.
  ## logit = @(x) 1./(exp(x)+1);
  ## logit3 = @(x) exp(x)*(exp(x)-1)./(exp(x)+1).^3;
  ## induced = (p.saturating_induced .* logit(p.induced_scale.*data.content) ...
  ##            + p.wiggle_induced .* logit(p.induced_scale.*data.content));

  ## bias = p.mu_0 ...
  ##        + p.beta_induced.*data.content ...
  ##        + p.beta_summation.*summation.*data.content ...
  ##        + induced;

  ## prelink = bias + p.beta_0.*data.dx.*sens + p.beta_small.*data.dx.*(1-sens);

  ## logit = @(x)0.98./(1+exp(-x)) + 0.01;
  ## prob = logit(prelink);
  #plot the data similarly to how I've been plotting it all along.

  l_ply(unique(data$subject), plot_page)
}
