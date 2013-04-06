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
  motion.energy <- chain(grid, read.csv, add_energies)

  #how about some contour plots. Everyone hated these.
  if (!interactive()) {
    cairo_pdf(outfile, onefile=TRUE)
    on.exit(dev.off(), add=TRUE)
  }
}

if (exists("model.df")) model <- subset(model.df, subject=="pbm")[[1, "model"]]

plot_contours <- function(model, motion.energy) {
  # we want three contour plots along our three axes --
  # spacing, displacement and direction content --
  # maybe even put it in 3d with the other one.

  # we also need to bin from 3d into 2d showing residuals for each case.
  # so we need to decide where the bins are placed...

  #because 20/3 in the dataset is different form R's idea of 20/3....
  nominal.eccentricity <- take_nearest(20/3, motion.energy$eccentricity)

  # grid coordinates to sample on
  is.motion.energy <- "motion_energy_models" %in% class(model)
  if ("motion_energy_model" %in% class(model)) {
    ##in motion energy models we can only evaluate stimuli whose
    ##motion energies have been precomputed
    bind[displacement.sampling, content.sampling, spacing.sampling] <- (
      chain(motion.energy, subset(grid==TRUE),
            mutate(spacing = eccentricity * 2 * pi / target_number_all),
            .[c("displacement", "content", "spacing")],
            lapply(unique), lapply(sort)))
    ##geom_tile since we are not guaranteed even spacing. Is there a
    ##way to draw this interpolated?
    geom <- geom_tile()
  } else {
    bind[displacement.sampling, content.sampling, spacing.sampling] <- (
      chain(motion.energy, subset(grid==TRUE),
            mutate(spacing = eccentricity * 2 * pi / target_number_all),
            .[c("displacement", "content", "spacing")],
            lapply(range), lapply(seq_range, length=50), lapply(sort)))
    #why doesn't interpolate seem to have an effect?
    geom <- layer(geom="raster", geom_params=list(interpolate=TRUE))
  }

  content.bins <- unique(round_any(content.sampling, 0.2))
  displacement.bins <- unique(round_any(content.sampling, 0.1))
  spacing.bins <- unique(round_any(spacing.sampling, 2))

  wide.spacing <- take_nearest(2*pi*nominal.eccentricity/6, spacing.sampling)
  narrow.spacing <- take_nearest(2*pi*nominal.eccentricity/20 , spacing.sampling)

  bind[grids, bins] <- Map(
    spacing = list(spacing.sampling, spacing.bins),
    displacement = list(displacement.sampling, displacement.bins),
    content = list(content.sampling, content.bins),
    fun(list(
      displacement_spacing = expand.grid(
        spacing = spacing,
        displacement = displacement,
        content = 0),
      spacing_content = expand.grid(
        spacing = spacing,
        content = content,
        displacement = 0),
      content_displacement_wide = expand.grid(
        spacing = wide.spacing,
        content = content,
        displacement = displacement),
      content_displacement_narrow = expand.grid(
        spacing = narrow.spacing,
        content = content,
        displacement = displacement)
      )))

  xvars <- c("displacement", "spacing", "displacement", "displacement")
  yvars <- c("spacing", "content", "content", "content")
  xscales <- list(displacement_scale_nopadding,
                  spacing_scale_x_nopadding,
                  displacement_scale_nopadding,
                  displacement_scale_nopadding)
  yscales <- list(spacing_scale_y_nopadding,
               content_scale_y_nopadding,
               content_scale_y_nopadding,
               content_scale_y_nopadding)

  spacing.threshold <- 4.5
  filters <- list(
    identity,
    identity,
    subset %<<% dots(spacing >= spacing.threshold),
    subset %<<% dots(spacing < spacing.threshold))

  annotations <- with_arg(
    x=Inf, y=-Inf, geom="text", vjust=-0.5, hjust=1.05, size=3,
    color="gray50",
    annotate(label="Carrier = 0"),
    annotate(label="Envelope = 0"),
    annotate(label=sprintf("Spacing = %.2g \n(binning >= %.2g)",
               wide.spacing, spacing.threshold)),
    annotate(label=sprintf("Spacing = %.2g \n(binning < %.2g)",
               narrow.spacing, spacing.threshold)))

  #cook in additional fields that the model may need
  bind[grids, bins] <- lapply(
    list(grids, bins), lapply,
    mkchain(
      mutate(
        eccentricity = (
          if (exists("eccentricity")) eccentricity else nominal.eccentricity),
        bias = 1,
        target_number_all = (
          if (exists("target_number_all")) target_number_all
          else round(2*pi*eccentricity / spacing))),
      recast_data,
      if (is.motion.energy) attach_motion_energy(., motion.energy) else .,
      mutate(., pred = predict(model, newdata=., type="response"))))

  subject <- unique(model$data$subject)
  plot.tables <- Map(
    grid=grids, bin=bins, xscale=xscales, yscale = yscales, fig=2:5,
    xvar=xvars, yvar=yvars, anno=annotations, filt=filters,
    f = function(grid, bin, xscale, yscale, fig, xvar, yvar, anno, filt) {
      ##we also need "actual data" binned along the
      ##missing variable. This function snaps the data to
      ##grid lines, while computing an "average" that
      ##preserves the residual.
      binned_data <- bin_grid_resid(
        model, bin, data=filt(model$data), coords=c(xvar, yvar))
      the.plot <- (
        ggplot(grid)
        + xscale + yscale
        + decision_contour
        + no_grid
        + geom
        + geom_point(data=binned_data, shape=21, color="blue",
                     aes(size=n_obs, fill=bound_prob(p)))
        + anno
        + scale_size_area("N", breaks=c(20, 50, 100, 200, 500))
        + labs(title="foo"
               ## , title=sprintf(
               ##   "Model fit and data for subject %s", toupper(subject))
               ))
      ggplot_gtable(ggplot_build(the.plot))
    })

  #Stuff four plots in one, using the legend from one of them.
  gt <- chain(
    gtable(widths = unit.c(unit(c(1,1), "null"), gtable_width(plot.tables[[1]][,5])),
           heights = unit.c(unit(1, "lines"), unit(c(1, 1), "null"))),
    gtable_add_grob(plot.tables[[2]][-1:-2,-5], 2, 1),
    gtable_add_grob(plot.tables[[3]][-1:-2,-5], 2, 2),
    gtable_add_grob(plot.tables[[1]][-1:-2,-5], 3, 1),
    gtable_add_grob(plot.tables[[4]][-1:-2,-5], 3, 2),
    gtable_add_grob(plot.tables[[1]][,5], 2, 3, 3))

  grid.newpage()
  grid.draw(gt)

  #let's also make a 3d plot to serve as a key.
  library(rgl)

}
