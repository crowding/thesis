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
source("slopeModel.R")

infile <- "slopeModel.RData"
grid <- "motion_energy.csv"
outfile <- "contours.pdf"

main <- function(infile = "slopeModel.RData", grid = "motion_energy.csv",
                 outfile = "contours.pdf") {

  load(infile)
  motion.energy <- add_energies(read.csv(grid))

  bind[displacement.sampling, content.sampling, spacing.sampling] <- (
    chain(motion.energy, subset(grid==TRUE),
          mutate(spacing = eccentricity * 2 * pi / target_number_all),
          .[c("displacement", "content", "spacing")],
          lapply(unique), lapply(sort)))
  displacement.sampling <<- displacement.sampling
  content.sampling <<- content.sampling
  spacing.sampling <<- spacing.sampling

  #how about some contour plots. Everyone hated these.
  if (!interactive()) {
    cairo_pdf(outfile, onefile=TRUE)
    on.exit(dev.off(), add=TRUE)
  }
}

model <- subset(model.df, subject=="pbm")[[1, "model"]]

plot_contours <- function(model, motion.energy) {
  # we want three contour plots along our three axes --
  # spacing, displacement and direction content --
  # maybe even put it in 3d with the other one.

  # we also need to bin from 3d into 2d showing residuals for each case.
  # so we need to decide where the bins are placed...

  nominal.eccentricity <- take_nearest(20/3, motion.energy$eccentricity)
  # first, displacement versus spacing.

  wide.spacing <- take_nearest(2*pi*nominal.eccentricity/6, spacing.sampling)
  narrow.spacing <-take_nearest(2*pi*nominal.eccentricity/20 , spacing.sampling)
  
  grids <- within(list(), {
    displacement_spacing <- expand.grid(
      spacing = spacing.sampling,
      displacement = displacement.sampling,
      content = 0)
    spacing_content <- expand.grid(
      spacing = spacing.sampling,
      content = content.sampling,
      displacement = 0)
    content_displacement_wide <- expand.grid(
      spacing = wide.spacing,
      content = content.sampling,
      displacement = displacement.sampling)
    content_displacement_narrow <- expand.grid(
      spacing = narrow.spacing,
      content = content.sampling,
      displacement = displacement.sampling)
  })

  #cook in additional fields that the model may need
  #because matlab's idea of 20/3 is different from R's....
  grids <- lapply(
    grids,
    mkchain(
      mutate(
        eccentricity = (
          if (exists("eccentricity")) eccentricity else nominal.eccentricity),
        bias = 1,
        target_number_all = (
          if (exists("target_number_shown")) target_number_all
          else round(2*pi*eccentricity / spacing))),
      recast_data,
      attach_motion_energy(motion.energy),
      mutate(., pred = predict(model, newdata=., type="response"))))

  print(ggplot(grids$displacement_spacing)
        + aes(x = displacement, y = spacing, z = pred)
        + decision_contour
        + geom_vline(x = 0, color = "gray50", linetype = "11", size = 0.2)
        + displacement_scale_nopadding + y_nopadding
        + no_grid
        + annotate("text", label=toupper(model$data$subject[[1]]),
                   x=Inf, y=-Inf, hjust=1.2, vjust=-0.5))

  print(ggplot(grids$spacing_content)
        + aes(x = content, y = spacing, z = pred)
        + decision_contour
        + geom_vline(x = 0, color = "gray50", linetype = "11", size = 0.2)
        + content_scale_nopadding + y_nopadding
        + no_grid
        + annotate("text", label=toupper(model$data$subject[[1]]),
                   x=Inf, y=-Inf, hjust=1.2, vjust=-0.5))

  print(ggplot(grids$content_displacement_narrow)
        + aes(x = content, y = displacement, z = pred)
        + decision_contour
        + geom_vline(x = 0, color = "gray50", linetype = "11", size = 0.2)
        + content_scale_nopadding + y_nopadding
        + no_grid
        + annotate("text", label=toupper(model$data$subject[[1]]),
                   x=Inf, y=-Inf, hjust=1.2, vjust=-0.5))

  print(ggplot(grids$content_displacement_wide)
        + aes(x = content, y = displacement, z = pred)
        + decision_contour
        + geom_vline(x = 0, color = "gray50", linetype = "11", size = 0.2)
        + content_scale_nopadding + y_nopadding
        + no_grid
        + annotate("text", label=toupper(model$data$subject[[1]]),
                   x=Inf, y=-Inf, hjust=1.2, vjust=-0.5))

  #todo:_add "raw" (binned, munged) data, change colors, make 3d plot
}
