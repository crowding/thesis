## Rebooted contour plots, Use shading, with a colormap. Also maybe
## mess with 3d plotting.

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

infile <- "slopeModel.RData"
grid <- "motion_energy.csv"
outfile <- "contours.pdf"

main <- function(infile = "slopeModel.RData", grid = "motion_energy.csv",
                 outfile = "contours.pdf") {
  load(infile)
  motion.energy <- add_energies(read.csv(grid))

  bind[plot.displacement, plot.content, plot.spacing] <- (
    chain(motion.energy, subset(grid==TRUE),
          mutate(spacing=target_number_all * 2*pi/eccentricity),
          .[c("displacement", "content", "spacing")],
          lapply(unique), lapply(sort)))
  plot.displacement <<- plot.displacement
  plot.content <<- plot.content
  plot.spacing <<- plot.spacing

  #how about some contour plots. Everyone hated these.
  if (!interactive()) {
    cairo_pdf(outfile, onefile=TRUE)
    on.exit(dev.off(), add=TRUE)
  }
  (mapply %<<% model.df)(function(model, subject) {
    cat("plotting contours for", as.character(subject), "\n")
    #plot_contours(model)
    tryCatch(plot_contours(model), error=function(x) warning(x))
  })
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
