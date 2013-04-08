## Rebooted contour plots, Use shading, with a colormap. Also maybe
## mess with 3d plotting.

suppressPackageStartupMessages({
  library(plyr)
  library(ggplot2)
  library(ptools)
  library(psyphy)
  library(gnm)
  library(grid)
  library(rgl)
  library(fields)
  library(reshape2)
  library(gtable)
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
                 outlist = "contours/contours.list") {

  load(infile)
  motion.energy <- chain(grid, read.csv, add_energies)

  outlist.conn <- file(outlist, "w")
  on.exit(close(outlist.conn), add=TRUE)

  #since our 3d plots can't plot directly into a multipage PDF, we'll
  #have to list separate files
  (Map %<<% model.df)(f = function(model, subject, ...) {
    subject <- as.character(subject)
    cat("plotting subject ", subject, "\n")
    pdf(pdf.file <- replace_extension(outlist, "pdf",
                                     paste0("_", subject, "_2d")))
    on.exit(dev.off(), add=TRUE)
    open3d(windowRect=c(100L, 100L, 1024L, 512L))
    on.exit(rgl.close(), add=TRUE)
    plot_contours(motion.energy=motion.energy, model=model, subject=subject, ...)
    rgl.postscript(rgl.file <- replace_extension(outlist, "pdf",
                                                 paste0("_", subject, "_3d")))
    writeLines(c(pdf.file, rgl.file), outlist.conn)
  })
}

plot_contours <- function(model, subject, motion.energy, outlist, ...) {
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
    stop("extract sampling form motion energy limited to range of stimuli")
    ##Use geom_tile since we are not guaranteed even spacing. Is there a
    ##way to draw this interpolated?
    geom <- geom_tile()
  } else {
    bind[displacement.sampling, content.sampling, spacing.sampling] <- (
      chain(model$data,
            .[c("displacement", "content", "spacing")],
            lapply(range), lapply(seq_range, length=50), lapply(sort)))
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
  xscales <- list(displacement_scale,
                  spacing_scale_x,
                  displacement_scale,
                  displacement_scale)
  yscales <- list(spacing_scale_y,
               content_scale_y,
               content_scale_y,
               content_scale_y)

  spacing.threshold <- 4.5
  filters <- list(
    identity,
    identity,
    here(subset) %<<% dots(spacing >= spacing.threshold),
    here(subset) %<<% dots(spacing < spacing.threshold))

  annotations <- with_arg(
    x=Inf, y=Inf, geom="text", vjust=1, hjust=1.2, size=3,
    color="gray50",
    annotate(label="\nCarrier = 0"),
    annotate(label="\nEnvelope = 0"),
    annotate(label=sprintf("\nSpacing = %.2g \n(binning >= %.2g)",
               wide.spacing, spacing.threshold)),
    annotate(label=sprintf("\nSpacing = %.2g \n(binning < %.2g)",
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
        + no_grid
        + geom
        + decision_contour
        + geom_point(data=binned_data, shape=21, color="blue",
                     aes(size=n_obs, fill=bound_prob(p)))
        + anno
        + scale_size_area("N", breaks=c(20, 50, 100, 200, 500))
        + theme(panel.border=element_blank())
        + labs(title="foo"))
      ggplot_gtable(ggplot_build(the.plot))
    })

  #Stuff four plots in one, using the legend from one of them.
  titleGrob <- textGrob(label=sprintf("Subject %s", toupper(subject)),
                        gp=gpar(fontsize=18))
  gt <- chain(
    gtable(widths = (unit.c(unit(c(1,1), "null"),
                            gtable_width(plot.tables[[1]][,5]))),
           heights = (unit.c(2 * grobHeight(titleGrob),
                             unit(c(1, 1), "null")))),
    gtable_add_grob(titleGrob, 1, 1, 1, 2),
    gtable_add_grob(plot.tables[[2]][-1:-2,-5], 2, 1),
    gtable_add_grob(plot.tables[[3]][-1:-2,-5], 2, 2),
    gtable_add_grob(plot.tables[[1]][-1:-2,-5], 3, 1),
    gtable_add_grob(plot.tables[[4]][-1:-2,-5], 3, 2),
    gtable_add_grob(plot.tables[[1]][,5], 2, 3, 3))
  grid.newpage()
  grid.draw(gt)

  #let's also make a 3d plot to serve as a key.
  plot_3d_grids(model=model, grids=grids, bins=bins)
}

plot_3d_grids <- function(model, grids, ...) {
  #turn data frames in "grids" into matrices for plotting
  bind[x, y, z, value] <- zip(lapply(grids, matrixify), collate=list)
  rgl.clear()
  bg3d(color="gray80")
  par3d(scale=c(10, 6, 1.5))
  view3d(130, 15, 40, 1)
  #For each plane...
  Map(x=x, y=y, z=z, value=value, function(x, y, z, value) {
    #draw the colormapped plane
    colors = hsv(s=0,v=bound_prob(value))
    surface3d(x, y, z, color=colors, lit=FALSE)
    #compute and draw contour lines on the surface
    clines <- contourLines(z=value, levels=seq(0.1,0.9,0.2))
    lapply(clines, splat(function(level, xi, yi) {
      obj <- list(x = seq(0, 1, length=nrow(x)), y=seq(0, 1, length=ncol(y)))
      lineX <- interp.surface(c(obj, list(z=x)), cbind(xi, yi))
      lineY <- interp.surface(c(obj, list(z=y)), cbind(xi, yi))
      lineZ <- interp.surface(c(obj, list(z=z)), cbind(xi, yi))
      lines3d(lineX, lineY, lineZ, color="blue", lit=FALSE, lwd=0.5)
    }))
    #and outline the edges of the plane
    s = dim(x)
    indices <- cbind(c(1, s[1], s[1], 1, 1), c(1, 1, s[2], s[2], 1))
    lines3d(x[indices], y[indices], z[indices], color="gray50", lit=FALSE)
  })
  #add axes
  axis3d("x--", nticks=5, expand=1)
  axis3d("y++", expand=1)
  axis3d("z-+", expand=1)
  mtext3d("Envelope motion", "x--", 1, at=0.4, adj=1.2)
  mtext3d("Carrier strength", 'y++', 1, at=-1.5, adj=0)
  mtext3d("Spacing", 'z-+', 3, at=15)
  ##maybe we want to compute null (PSE) surface...
}

matrixify <- function(grid) {
  lapply(c("displacement", "content", "spacing", "pred"),
         function(var) {
           drop.dims(acast(grid, displacement ~ content ~ spacing,
                           value.var = var))
         })
}

drop.dims <- function(a) {
  `[` %()% (dots(a)
            %__% replicate(length(dim(a)), missing_value())
            %__% list(drop=TRUE))
}

#too slow to use...
rgl.grob <- function(...) {
  the.postscript <- tempfile(fileext=".eps")
  rgl.postscript(the.postscript)
  the.xml <- tempfile(fileext=".xml")
  PostScriptTrace(file=the.postscript, outfilename=the.xml)
  pictureGrob(readPicture(the.xml))
}
