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
source("icons.R")
source("library.R")
source("slopeModel.R")

infile <- "slopeModel.RData"
grid <- "motion_energy.csv"
outlist <- "contours/contours.list"
fold <- TRUE
dev.fun <- dev.new

main <- function(infile = "slopeModel.RData", grid = "motion_energy.csv",
                 outlist = "contours/contours.list", fold=c(FALSE, TRUE),
                 dev.fun= cairo_pdf %<<% dots(width=8, height=6)) {

  out <- match.fun(dev.fun)
  fold <- if(is.logical(fold)) fold[[1]] else match.arg(fold)
  load(infile)
  motion.energy <- chain(grid, read.csv, add_energies)

  outlist.conn <- file(outlist, "w")
  on.exit(close(outlist.conn), add=TRUE)

  #since our 3d plots can't plot directly into a multipage PDF, we'll
  #have to list separate files
  bind[model=bind[model], subject, ...=] <- as.list(model.df[1,])
  open3d(windowRect=c(100L, 100L, 768L, 512L))
  on.exit(rgl.close(), add=TRUE)
  (Map %<<% model.df)(f = function(model, subject, ...) tryCatch({

    subject <- as.character(subject)
    cat("plotting subject ", subject, "\n")
    dev.fun(pdf.file <- replace_extension(outlist, "pdf",
                                      paste0("_", subject, "_2d")))
    on.exit(dev.off(), add=TRUE)
    plot_contours(motion.energy=motion.energy, model=model, subject=subject,
                  fold=fold, ...)
    rgl.postscript(
      fmt="pdf",
      (rgl.file <- replace_extension(outlist, "pdf",
                                      paste0("_", subject, "_3d"))))
    writeLines(c(pdf.file, rgl.file), outlist.conn)
  }, error=warning))
}

plot_contours <- function(model, subject, motion.energy, outlist,
                          fold=FALSE, ...) {
  # we want three contour plots along our three axes --
  # spacing, displacement and direction content --
  # maybe even put it in 3d with the other one.

  # we also need to bin from 3d into 2d showing residuals for each case.
  # so we need to decide where the bins are placed...

  #because 20/3 in the dataset is different form R's idea of 20/3....
  nominal.eccentricity <- take_nearest(20/3, motion.energy$eccentricity)

  roundings <- c(0.2, 0.2, 2)
  # id coordinates to sample on
  is.motion.energy <- "motion_energy_models" %in% class(model)
  if ("motion_energy_model" %in% class(model)) {
    ##in motion energy models we can only evaluate stimuli whose
    ##motion energies have been precomputed
    stop("extract sampling form motion energy limited to range of stimuli")
    ##Use geom_tile since we are not guaranteed even spacing. Is there a
    ##way to draw this interpolated?
    geom <- geom_tile()
  } else {
    bind[displacement.sampling, content.sampling, spacing.sampling] <- chain(
      model$data,
      .[c("displacement", "content", "spacing")],
      lapply(range),
      Map(f=round_any, roundings),
      Map(f=pmax, list(-Inf, if(fold) 0 else -Inf, -Inf)),
      Map(f=function(x,r) x+c(-0.50*r, 0.50*r), roundings),
      lapply(seq_range, length=51), lapply(sort))
   bind[content.bins, displacement.bins, spacing.bins] <- chain(
      model$data,
      refold(fold=fold),
      .[c("content", "displacement", "spacing")],
      Map(f=round_any, roundings),
      lapply(unique))
    geom <- layer(geom="raster", geom_params=list(interpolate=TRUE))
  }

  wide.spacing <- take_nearest(2*pi*nominal.eccentricity/6, model$data$spacing)
  narrow.spacing <- take_nearest(2*pi*nominal.eccentricity/20, model$data$spacing)

  bind[grids, bins] <- Map(
    spacing = list(spacing.sampling, spacing.bins),
    displacement = list(displacement.sampling, displacement.bins),
    content = list(content.sampling, content.bins),
    fun(list(
      displacement_spacing = expand.grid(
        spacing = spacing,
        displacement = displacement,
        content = .Machine$double.xmin), # to avoid folding on displacement
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
    here(subset) %<<% dots(spacing >= spacing.threshold),
    here(subset) %<<% dots(spacing < spacing.threshold))

  annotations <- with_arg(
    x=Inf, y=Inf, geom="text", vjust=1.3, hjust=1.2, size=3,
    color="black",
    annotate(label="Carrier = 0"),
    annotate(label="Envelope = 0"),
    annotate(label=sprintf("Spacing = %.2g \n(using trials >= %.2g)",
               wide.spacing, spacing.threshold)),
    annotate(label=sprintf("Spacing = %.2g \n(using trials < %.2g)",
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
      mutate(., pred = folding_predict(model, newdata=., type="response", fold=fold))))

  plot.tables <- Map(
    grid=grids, bin=bins, xscale=xscales, yscale = yscales, fig=2:5,
    xvar=xvars, yvar=yvars, anno=annotations, filt=filters,
    f = function(grid, bin, xscale, yscale, fig, xvar, yvar, anno, filt) {
      ##we also need "actual data" binned along the
      ##missing variable. This function snaps the data to
      ##grid lines, while computing an "average" that
      ##preserves the residual.
      binned_data <- bin_grid_resid(
        model, bin, data=filt(model$data), coords=c(xvar, yvar), fold=fold)
      the.plot <- (
        ggplot(grid)
        + xscale + yscale
        + no_grid
        + geom
        + decision_contour
        + geom_circle(data=binned_data, linetype="12", color="white", weight=0.35,
                      aes(size=n_obs, fill=bound_prob(p)))
        + anno
        + scale_size_area("N", breaks=c(20, 50, 100, 200, 500))
        + labs(title="foo")
        + guides(size=guide_legend("N",
                   override.aes=list(colour="black")))
        + theme(aspect.ratio=1))
      ggplot_gtable(ggplot_build(the.plot))
    })
  #
  #Stuff four plots in one, using the legend from one of them.
  titleGrob <- textGrob(label=sprintf("Observer %s", toupper(subject)),
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
  grid.draw(gt)

  #let's also make a 3d plot to serve as a key.
  plot_3d_grids(model=model, grids=grids, bins=bins)
}

plot_3d_grids <- function(model, grids, fold=FALSE, ...) {
  #turn data frames in "grids" into matrices for plotting
  bind[x, y, z, value] <- zip(lapply(grids, matrixify), collate=list)
  rgl.clear()
  bg3d(color="gray80")
  par3d(scale=c(6, 6, 1.5))
  view3d(220, 20, 45, 1)
  #For each plane...
  Map(x=x, y=y, z=z, value=value, function(x, y, z, value) {
    #draw the colormapped plane
    colors = gradient_n_pal(decision_colors, decision_values)(value)
    surface3d(-x, y, z, color=colors, lit=FALSE)
    #compute and draw contour lines on the surface
    clines <- contourLines(z=value, levels=seq(0.1,0.9,0.2))
    lapply(clines, splat(function(level, xi, yi) {
      obj <- list(x = seq(0, 1, length=nrow(x)), y=seq(0, 1, length=ncol(y)))
      lineX <- interp.surface(c(obj, list(z=x)), cbind(xi, yi))
      lineY <- interp.surface(c(obj, list(z=y)), cbind(xi, yi))
      lineZ <- interp.surface(c(obj, list(z=z)), cbind(xi, yi))
      lines3d(-lineX, lineY, lineZ, color="#DDDD00", lit=FALSE, lwd=0.5)
    }))
    #and outline the edges of the plane
    s = dim(x)
    indices <- cbind(c(1, s[1], s[1], 1, 1), c(1, 1, s[2], s[2], 1))
    lines3d(-x[indices], y[indices], z[indices], color="gray50", lit=FALSE)
  })
  #add axes
  xat <- chain(par3d("bbox"), .[1:2], pretty(3), .[c(-1, -length(.))])
  axis3d("x--", at=xat, labels=as.character(-xat), expand=1)
  yat <- chain(par3d("bbox"), .[3:4], pretty(3), .[c(-1, -length(.))])
  axis3d("y-+", expand=1, at=yat)
  zat <- chain(par3d("bbox"), .[5:6], pretty(3), .[c(-1, -length(.))])
  axis3d("z++", expand=1, at=zat)

  mtext3d("Envelope motion", "x--", 0, at=xat[[1]], adj=-0.2)
  mtext3d("Carrier strength", 'y-+', 2, at=yat[1] - 0.5, adj=0.5)
  mtext3d("Spacing", 'z++', 3, at=10)
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

run_as_command()
