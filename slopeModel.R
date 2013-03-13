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
source("library.R")
#suppressMessages(loadfonts())
#quartz.options(family="MS Gothic")
#pdf.options(family="MS Gothic")

match_df <- function(...) suppressMessages(plyr::match_df(...))
seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)
`%++%` <- paste0
splits <- c("subject", "content", "exp_type", "spacing", "bias")

makePredictions <-
  function(model, data=model$data,
           splits=c("subject", "content", "exp_type", "spacing", "bias"),
           ordinate = "displacement",
           ordinate.values = plot.displacement
           ) {
    r <- range(data[ordinate])
    sampling <- merge(unique(data[splits]),
                      quickdf(structure(list(ordinate.values),
                                        names=ordinate)),
                      by.x=c(), by.y=c())
    #chunk up predictions... gnm has some O(n^2) memory usage bullshit going on
    #with making predictions on binary data, so predict 1000 trials at a time
    #(tested faster than 500 or 2000)
    .nrow <- nrow(sampling); .chunksize <- 1000
    sampling <- mutate(sampling,  .chunk=floor(seq_len(.nrow)/.chunksize))
    pred <- ddply(sampling, ".chunk",
                  function(chunk) {
                    pred <- doPredict(model, newdata=chunk, type="response", se.fit=TRUE)
                    cbind(chunk, pred[1:2])
                  })
  }

plotPredictions <- function(...) {
  predictions = makePredictions(...)
  with_arg(data=predictions,
           geom_ribbon(color="transparent", alpha=0.2,
                       aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
           geom_line(aes(y=fit)))
}

#perhaps make this go using predict_from_model_frame
plot_fit <- function(model, subject=model$data$subject[1],
                     style=c("binned", "bubble")) {
  style <- match.arg(style)
  subdata <- match_df(model$data,
                      data.frame(subject=subject, stringsAsFactors=FALSE),
                      on="subject")
  switch(style, bubble = {
    plotdata <- subdata
  }, binned = {
    plotdata <- bin_along_resid(model, subdata,
                                "response", splits, "displacement")
  })
  print((ggplot(plotdata)
         + displacement_scale
         + proportion_scale
         + content_color_scale
         + facet_spacing_experiment
         + plotPredictions(model)
         + geom_point()
         + (switch(style, bubble=balloon, binned=geom_point()))
         + labs(title = "Data and model fits for subject " %++% subject)
         ))
}

infile <- "data.Rdata"
grid <- "motion_energy.csv"
outfile <- "slopeModel.RData"
plot <- "slopeModel.pdf"

main <- function(infile = "data.RData", grid = "motion_energy.csv",
                 outfile = "slopeModel.RData", plot="slopeModel.pdf") {

  load(infile, envir = e <- new.env())
  motion.energy <- add_energies(read.csv(grid))

  bind[plot.displacement, plot.content, plot.spacing] <- (
    chain(motion.energy, subset(grid==TRUE),
          mutate(spacing=target_number_all * 2*pi/eccentricity),
          .[c("displacement", "content", "spacing")],
          lapply(unique), lapply(sort)))
  plot.displacement <<- plot.displacement
  plot.content <<- plot.content
  plot.spacing <<- plot.spacing

  #use only subjects whose number of trials exceed 2000
  chain(e$data
        , subset(exp_type %in% c("spacing", "content"))
        , do.rename(folding=FALSE)
        , mutate(bias=1)
        , match_df(., subset(count(., "subject"), freq>2000), on="subject")
        , attach_motion_energy(motion.energy)
        ) -> data

  #count trials in each condition. While keeping motion energy
  #information, this kills speed. Might do binning instead.

  #Our model has one term nonlinear in the spacing-dependent
  #sensitivity to displacement and to direction content.
  #This defines the response to displacement.
  displacementTerm <<- (nonlinearTerm(cs, beta_dx)(spacing, displacement)
                       ((2 - 2/(1+exp(-cs/spacing))) * beta_dx * displacement))
  formula <- (  response
              ~ displacementTerm(spacing, displacement,
                                 start=c(cs=4, beta_dx=14))
              + content
              + I(content*abs(content))
              + I(1/spacing):content
              + bias - 1 #"bias" set to 0 to predict folded data (but a
              # better way to predict folded data is to predict both
              # directions then average)
              )
  family <- binomial(link=logit.2asym(g=0.025, lam=0.025))

  #fit models to each subject.
  models <- dlply(data, "subject", function(chunk) {
    cat("fitting subject ", chunk$subject[1], "\n")
    gnm(formula, family=family, data=chunk)
  })
  model.frame <- data.frame(model = I(models), subject=names(models))

  save(model.frame, models, displacementTerm, formula, family,
       plot.displacement, plot.content, plot.spacing, file=outfile)

  #try using motion energy (and normalized motion energy) instead of
  #direction content in the model...
  motion.energy.formula <-
    update(formula, . ~ . + I(contrast_diff/spacing))

  motion.energy.models <- dlply(data, "subject", function(chunk) {
    cat("fitting subject ", chunk$subject[1],  "\n")
    motion_energy_model(gnm(motion.energy.formula, family=family, data=chunk),
                        motion.energy)
  })
  motion.energy.model.frame <-
    data.frame(model=I(motion.energy.models), subject=names(models))
  #compare the models by residual deviance (smaller is better)
  #Positive numbers are better here.
  ddply(merge(model.frame, motion.energy.model.frame,
              by="subject", suffixes=c(".content", ".energy")), "subject",
        summarize, difference = (extractAIC(model.content[[1]])[[2]]
                                 - extractAIC(model.energy[[1]])[[2]]))

  #well, that's weird, it made things worse according to the
  #AIC. You know what, maybe if I plot the fit it will shed some light.
  plot_fit(motion.energy.model.frame[["pbm", "model"]])

  #plot the models
  cairo_pdf(plot, onefile=TRUE, family="MS Gothic")
  (mapply %<<% model.frame)(function(model, subject) {
    cat("plotting subject ", as.character(subject), "\n")
    plot_fit(model)
    #tryCatch(plot_fit(model), error=function(x) warning(x))
  })
  dev.off()

  #how about some contour plots. Everyone hated these.
  cairo_pdf("contours.pdf", onefile=TRUE, family="MS Gothic")
  (mapply %<<% model.frame)(function(model, subject) {
    cat("plotting contours for", as.character(subject), "\n")
    #plot_contours(model)
    tryCatch(plot_contours(model), error=function(x) warning(x))
  })
  dev.off()

  #this is where we might make some 3d plots. Or contour plots along
  #different axes.

  #draw some random simulations of coefficients.  plot where the
  #coefficients of the models are (with confidence ellipses?)

  #compare that with the coefficients gleaned empirically from
  #number/density data.

  #make the point that there is no pooling, comparing slopes of 2, 4,
  #and 6 elements on screen. Somehow also average the data across
  #subjects, too?

  #and this makes interesting plots that show us about the model
  #properties???
  plot_curves(models)
}


plot_contours <- function(model) {
  #make a contour plot with displacement on the x-axis and spacing on
  #the y-axis.

  displacement_sampling <- expand.grid(
                spacing = plot.spacing,
                displacement = plot.displacement,
                content = 0.1,
                bias = 1
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
      spacing      = plot.spacing,
      content      = plot.content,
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
        + scale_x_continuous(name="Direction content",
                             labels=newline_arrows, expand=c(0,0))
        + no_grid
        + annotate("text", label=toupper(model$data$subject[[1]]),
                   x=max(content_sampling$content),
                   y=min(content_sampling$spacing),
                   hjust=1.2, vjust=-0.5))

  #might be interesting to plot deviance over these coordinates. color
  #coded deviance plot?

  #uncertainty inside/outside
  uncertainty_sampling <- expand.grid(
    spacing = match_to(c(2, 3, 5, 10, 20), plot.spacing),
        displacement = plot.displacement,
        content = 0.1,
        bias=1
  )
  pred <- predict(model, newdata=uncertainty_sampling,
                  type="response", se.fit=TRUE)
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
      spacing = match_to(c(2, 3, 5, 10, 20), plot.spacing),
      content = plot.content,
      displacement = 0,
      bias = 1)
  pred = predict(model, newdata=content_sampling, type="response", se.fit=TRUE)
  content_sampling <- cbind(content_sampling, pred)

  print(ggplot(content_sampling)
        + content_scale
        + proportion_scale
        + spacing_texture_scale
        + ribbon
        + annotate("text", label=toupper(model$data$subject[[1]]),
                   x=max(model$data$content),
                   y=0,
                   hjust=1.2, vjust=-0.5)
        )

}

example_plots <- function(model.frame) {

}

illustrative_plots <- function(model.frame) {
  #illustrate sensitivity changes with spacing...

  #we want an x-axis: spacing, a y-axis: sensitivity
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

  sensitivity_data <-
    cbind(sensitivity_data,
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
  cairo_pdf(file=paste(prefix, "sensitivity.pdf", sep=""),
            width=3, height=2, family="MS Gothic")
  print(sensitivity_plot)
  dev.off()

  #another way to get at this is to run a prediction.
  bias_all_data <-
    expand.grid(spacing=seq(0, 10, len=200), content=0.2, displacement=0, bias=1)
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
  cairo_pdf(file=paste(prefix, "all_bias.pdf", sep=""),
            width=3, height=2, family="MS Gothic")
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
  cairo_pdf(file=paste(prefix, "all_bias2.pdf", sep=""), width=3, height=2, family="MS Gothic")
  print(bias_all_plot2)
  dev.off()

  #and for my last trick some plot of the distant bias
  wide_content_data <-
    expand.grid( content = seq(-1, 1, len=200), spacing=10, displacement=0, bias=1)
  wide_content_data <- cbind(wide_content_data,
                             localbias=predict(m, wide_content_data)
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
