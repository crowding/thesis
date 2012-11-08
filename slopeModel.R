suppressPackageStartupMessages({
library(plyr)
library(ggplot2)
library(ptools)
library(psyphy)
library(gnm)
library(grid)
library(extrafont)
})
theme_set(theme_bw())
use_unicode=TRUE
source("scales.R")
loadfonts()
quartz.options(family="MS Gothic")
pdf.options(family="MS Gothic")

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

mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type")) {
  chain(data,
        ddply(splits,
              summarize,
              n = length(response), p = mean(response),
              n_cw = sum(response), n_ccw = sum(!response)),
        arrange(desc(n)))
}

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

main <- function(infile = "data.RData", outfile = "slopeModel.RData") {

  load(infile, e <- new.env())

  #use only subjects whose calculations exceed 2...
  chain(e$data
        , do.rename(folding=FALSE)
        , subset(exp_type %in% c("spacing", "content"))
        , match_df(., subset(count(., "subject"), freq>2000), on="subject")
        ) -> data

  rates = mkrates(data)

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
  cairo_pdf("slopeModel.pdf", onefile=TRUE, family="MS Gothic")
  mapply(models, names(models), FUN=function(model, name) {
    cat("plotting subject " %++% name %++% "\n")
    plot_fit(model)
  })
  dev.off()

  #how about some contour plots.
  cairo_pdf("contours.pdf", onefile=TRUE, family="MS Gothic")
  mapply(models, names(models), FUN=function(model, name) {
    plot_contours(model)
  })
  dev.off()

  #make some contour plots

  #draw some random simulations of coefficients.  plot where the
  #coefficients of the models are (with confidence ellipses?)

  #compare that with the coefficients gleaned empirically from
  #number/density data.

  #make the point that there is no pooling, comparing slopes of 2, 4,
  #and 6 elements on screen. Somehow also average the data across
  #subjects, too?
}

plot_contours <- function(model) {
  #make a contour plot with displacement on the x-axis and spacing on
  #the y-axis.

  displacement_sampling <- expand.grid(
                spacing=seq_range(pmin(range(model$data$spacing), 10), length=100),
                displacement = seq_range(range(model$data$displacement), length=100),
                content = 0.1
                )
  displacement_sampling$pred <-
    predict(model, newdata=displacement_sampling, type="response")

  print(ggplot(displacement_sampling)
   + aes(x = displacement, y = spacing, z=pred)
   + geom_vline(x=0, color="gray50", linetype="11", size=0.2)
   + decision_contour
   + displacement_scale_nopadding + y_nopadding
   + no_grid
   + annotate("text", label=toupper(model$data$subject[[1]]),
              x=max(displacement_sampling$displacement),
              y=min(displacement_sampling$spacing),
              hjust=1.2, vjust=-0.5)
   )

  #now to show the (trickier) relation of direction content to spacing
  content_sampling <-
    expand.grid(
      spacing      = seq_range(pmin(range(model$data$spacing), 10), length=100),
      content      = seq_range(range(model$data$content), length=100),
      displacement = 0,
      subject = model$data$subject[[1]]
      )
  content_sampling$pred <-
    predict(model, newdata=content_sampling, type="response")

  print(ggplot(content_sampling)
   + aes(x = content, y = spacing, z = pred)
   + geom_contour(size=0.2, color="gray70", breaks=seq(0,1,0.02))
   + decision_contour
   + y_nopadding
   + scale_x_continuous(name="Direction content",labels=add_arrows, expand=c(0,0))
   + no_grid
   + annotate("text", label=toupper(model$data$subject[[1]]),
              x=max(content_sampling$content), y=min(content_sampling$spacing),
              hjust=1.2, vjust=-0.5)
   )

  #might be interesting to plot deviance over these coordinates. color
  #coded deviance plot?

  #uncertainty inside/outside
  uncertainty_sampling <- expand.grid(
        spacing = c(2, 3, 5, 10, 20),
        displacement = seq_range(range(model$data$displacement), length=100),
        content = 0.1
  )
  pred = predict(model, newdata=uncertainty_sampling, type="response", se.fit=TRUE)
  uncertainty_sampling <- cbind(uncertainty_sampling, pred)

  print(ggplot(uncertainty_sampling)
   + displacement_scale
   + spacing_texture_scale
   + ribbon
   + proportion_scale
   + annotate("text", label=toupper(model$data$subject[[1]]),
              x=max(uncertainty_sampling$displacement),
              y=0,
              hjust=1.2, vjust=-0.5)
   + no_grid
   + geom_vline(x=0, color="gray50", linetype="11", size=0.2)
   )

}

plot_curves <- function(models) {
  #plot interesting curves from each model (one per subject.)

  allData <- ldply(models, `[[`, "data")
  #the effect of direction content at wide spacing (set dx = 0, spacing = 10)

  #plot uncertainty as a function of spacing for all subjects.
  
  #show drop in slope (for one subject)
  
}
