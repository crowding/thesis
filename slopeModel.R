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
#suppressMessages(loadfonts())
#quartz.options(family="MS Gothic")
#pdf.options(family="MS Gothic")

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
  chain(
    rename(data, replacements)
    , mutate(bias=if(folding) 0 else 1))
}

mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type", "bias")) {
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
           splitting_vars=c("subject", "content", "exp_type", "spacing", "bias"),
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

plot_fit <- function(model, subject=model$data$subject[1]) {
  rates <- match_df(model$data,
                    data.frame(subject=subject, stringsAsFactors=FALSE),
                    on="subject")
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

  #use only subjects whose number of trials exceed 2000
  chain(e$data
        , subset(exp_type %in% c("spacing", "content"))
        , do.rename(folding=FALSE)
        , match_df(., subset(count(., "subject"), freq>2000), on="subject")
        ) -> data

  rates = mkrates(data)

  #fit models to each subject.
  models <- dlply(rates, "subject", function(chunk) {
    cat("fitting subject " %++% chunk$subject[1] %++% "\n")
    gnm(  (  cbind(n_cw, n_ccw)
           ~ displacementTerm(spacing, displacement)
           + content
           + I(content*abs(content))
           + I(1/spacing):content
           + bias - 1 #allows us to switch off bias to predict unfolded data
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

  plot_curves(models)
}

plot_contours <- function(model) {
  #make a contour plot with displacement on the x-axis and spacing on
  #the y-axis.

  displacement_sampling <- expand.grid(
                spacing=seq_range(pmin(range(model$data$spacing), 10), length=100),
                displacement = seq_range(range(model$data$displacement), length=100),
                content = 0.1,
                bias=1
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
      subject = model$data$subject[[1]],
      bias=1
      )
  content_sampling$pred <-
    predict(model, newdata=content_sampling, type="response")

  print(ggplot(content_sampling)
   + aes(x = content, y = spacing, z = pred)
   + geom_contour(size=0.2, color="gray70", breaks=seq(0,1,0.02))
   + decision_contour
   + y_nopadding
   + scale_x_continuous(name="Direction content",labels=newline_arrows, expand=c(0,0))
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
        content = 0.1,
        bias=1
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

  content_sampling <-
    expand.grid(
      spacing = c(2, 3, 5, 10, 20),
      content = seq(-1, 1, length=100),
      displacement = 0,
      bias = 1)
  pred = predict(model, newdata=content_sampling, type="response", se.fit=TRUE)
  content_sampling <- cbind(content_sampling, pred)

  print(ggplot(content_sampling)
        + content_x_scale
        + proportion_scale
        + spacing_texture_scale
        + ribbon
        + annotate("text", label=toupper(model$data$subject[[1]]),
                   x=max(model$data$content),
                   y=0,
                   hjust=1.2, vjust=-0.5)
        )

}

plot_curves <- function(models, prefix="../writing/inset_") {
  #plot interesting curves from each model (one per subject.)
  #plot curves in a way that is informative when included in a figure file.

  allData <- ldply(models, `[[`, "data")

  m <- models$pbm

  # first we show the decay of spacing with sensitivity. We mark it
  # with the critical spacing.
  # another way to get at this is to run a prediction.
  sensitivity_data <-
    expand.grid(spacing=seq(0, 10, len=200), content=0, bias=1)

  sensitivity_data <- cbind(sensitivity_data,
                                s=(predict(m,  data.frame(sensitivity_data, displacement=0.5) )
                                   -predict(m, data.frame(sensitivity_data, displacement=-0.5))))
  sensitivity_plot <- (
    ggplot(sensitivity_data)
    + aes(x=spacing, y=s)
    + geom_line()
    + scale_x_continuous(limits=c(0, 10),
                         breaks = c(0, m$coefficients["cs"], 10),
                         labels = c(0, "cs", 10))
    + scale_y_continuous("sensitivity", limits= c(0, m$coefficients["beta_dx"]),
                         breaks = c(0, m$coefficients["beta_dx"]),
                         labels = c("0", "\u03B2\u2080")
                         )
    )
  cairo_pdf(file=paste(prefix, "sensitivity.pdf", sep=""), width=3, height=2, family="MS Gothic")
  print(sensitivity_plot)
  dev.off()

  #another way to get at this is to run a prediction.
  bias_all_data <-
    expand.grid(spacing=seq(0, 10, len=200), content=0.2, displacement=0)
  bias_all_data <- cbind(bias_all_data, p=predict(m, bias_all_data, type="response")
                         )
  bias_all_plot <- (
    ggplot(bias_all_data)
    + proportion_scale
    + aes(x=spacing)
    + geom_line()
    + geom_ribbon(aes(ymin=pmin(0.5, p), ymax=p), color=NA, fill="green", alpha=0.5)
    + geom_ribbon(aes(ymin=p, ymax=pmax(0.5,p)), color=NA, fill="red", alpha=0.5)
    )
  cairo_pdf(file=paste(prefix, "all_bias.pdf", sep=""), width=3, height=2, family="MS Gothic")
  print(bias_all_plot)
  dev.off()

  #another way to look at it is to plot the bias (instead of the scaled plot)
  bias_all_data$localbias <- predict(m, bias_all_data)
  bias_all_plot2 <- (
    ggplot(bias_all_data)
    + aes(y = localbias)
    + scale_y_continuous("Bias", labels=replace_arrows)
    + aes(x=spacing)
    + geom_line()
    + geom_ribbon(aes(ymin=pmin(0, localbias), ymax=localbias), color=NA, fill="green", alpha=0.5)
    + geom_ribbon(aes(ymin=localbias, ymax=pmax(0,localbias)), color=NA, fill="red", alpha=0.5)
    + coord_cartesian(ylim=c(-5,20))
    )
  print(bias_all_plot2)

  #and for my last trick some plot of the distant bias
  wide_content_data <-
    expand.grid( content = seq(-1, 1, len=200), spacing=10, displacement=0)
  wide_content_data <- cbind(wide_content_data,
                             losalbias=predict(m, wide_content_data)
                                      - predict(m, mutate(wide_content_data, content=-content)))
  wide_content_plot <- (
    ggplot(wide_content_data)
    + geom_line()
    + aes(y = localbias)
    + scale_y_continuous("Bias", labels=replace_arrows)
    + aes(x=content)
    + scale_x_continuous(name="Direction content (at 10 degrees spacing)",labels=newline_arrows, expand=c(0,0))
    + geom_ribbon(aes(
                    ymin=ifelse(content>0, pmin(0, localbias), 0),
                    ymax=ifelse(content>0, 0,             pmax(0, localbias))),
                  color=NA, fill="red", alpha=0.5)
    + geom_ribbon(aes(
                    ymin=ifelse(content>0, pmax(0, bias), 0),
                    ymax=ifelse(content>0, 0,             pmin(0, localbias))),
                  color=NA, fill="green", alpha=0.5)
    )
  cairo_pdf(file=paste(prefix, "wide_content.pdf", sep=""), width=3, height=2, family="MS Gothic")
  print(wide_content_plot)
  dev.off()

  #the effect of direction content at wide spacing (set dx = 0, spacing = 10)

  #plot sensitivity as a function of spacing for all subjects.
  #just going to...

}

as.names <- function(names, value=missing.value()) {
  x <- replicate(length(names), value)
  names(x) <- names
  x
}

extract.nonlin.function <- function(nonlin.term) {
  eval(template(function( .a=...(nonlin.term$predictors),
                          .b=...(as.names(nonlin.term$variables)))
    {
       .(parse(text=nonlin.term$term(names(nonlin.term$predictors),
                        nonlin.term$variables ))[[1]])
    }
    ), parent.frame())
}

run_as_command()
