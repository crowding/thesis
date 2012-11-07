suppressPackageStartupMessages({
library(plyr)
library(ggplot2)
library(ptools)
library(psyphy)
library(gnm)
library(grid)
})
theme_set(theme_bw())
use_unicode=TRUE
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
                       summarize,
                       n = length(response), p = mean(response),
                       n_cw = sum(response), n_ccw = sum(!response)),
                 arrange(desc(n)))

seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)

`%++%` <- function(x, y) paste(x, y, sep="")

#Our model has one term nonlinear in the spacing-dependent
#sensitivity to displacement and to direction content.
#This defines the response to displacement.
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

makePredictions <-
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
wobble1 <- function(dx) {x <- dx / wobbleScale; 1./(exp(x)+1) }
wobble3 <- function(dx) {x <- dx / wobbleScale; exp(x)*(exp(x)-1)/(exp(x)+1)^3 }
wobbleScale <- 0.3

plot_fit <- function(model, subject=model$data$subject[1]) {
  rates <- match_df(model$data,
                    data.frame(subject=subject), on="subject")
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

plot_contours <- function(model, subject) {
  #make contour plots of the model prediction as a function of delta-x
  #and spacing for a fixed (20%) direction content.
}

main <- function(
          infile = "data.RData"
          outfile = "slopeModel.RData"
          ) {

  load(infile, e <- new.env())

  #use only subjects whose calculations exceed 2...
  chain(e$data
        , do.rename(folding=FALSE)
        , subset(exp_type %in% c("spacing", "content"))
        , match_df(., subset(count(., "subject"), freq>2000), on="subject")
        ) -> data

  rates = rates(data)

  #fit models to each subject.
  models <- dlply(rates, "subject", function(chunk) {
    cat("fitting subject " %++% chunk$subject[1] %++% "\n")
    gnm(  (  cbind(n_cw, n_ccw)
           ~ displacementTerm(spacing, displacement)
           + I(wobble1(content))
           + I(wobble3(content))
           + I(1/spacing):content
           )
        , family=binomial(link=logit.2asym(g=0.025, lam=0.025))
        , data=chunk
        )
  })

  save(models, file="slopeModel.RData")

  #plot the models
  cairo_pdf("slopeModel.pdf", onefile=TRUE)
  mapply(models, names(models), FUN=function(model, name) {
    cat("plotting subject " %++% name %++% "\n")
    plot_fit(model)
    grid.newpage()
  })
  dev.off()

  plot_fit(model, subject=subject)
  grid.newpage()
  plot_contours(model)
  dev.off()
  model
  #make some contour plots

  #draw some random simulations of coefficients.
  #plot where the coefficients of the models are (with confidence ellipses?)

  #compare that with the cieff
}
