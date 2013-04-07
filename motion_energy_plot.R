library(ptools)
library(plyr)
library(R.devices)
source("library.R")
source("scales.R")
library(reshape2)
setup_theme()

relabel_energy <- function(data) {
  mutate(data, variable = revalue(variable, c(
                 energy_cw = paste("CW energy", circleright),
                 energy_ccw = paste("CCW energy", circleleft),
                 energy_diff = "Difference",
                 energy_total = "Total",
                 norm_diff = "Normalized")))
}

revalue_named = function(value, x, ...) {
  revalue(value, names(x), x, ...)
}

#Break into segments, connecting each column of an array.
segmentify <- function(x) {
  nrow <- dim(x)[1]
  noth <- length(x) / nrow
  ix <- sort(c(seq_len(nrow), seq(from=2, length=nrow-2)))
  ix <- (matrix(rep(ix, noth), ncol=noth)
         + matrix(rep(nrow*seq(0, length=noth), length(ix)), ncol=noth, byrow=TRUE))
  x[as.vector(ix)]
}

main <- function(motion_energy_file="motion_energy.csv",
                 outfile="energy/motion_energy.list",
                 menergy=read.csv(motion_energy_file)
                 ) {
  fout = file(outfile, "w")
  on.exit(close(fout), add=TRUE)

  chain(menergy,
        subset(grid==1),
        add_energies
        ) -> energy_measurements

  pdffile = replace_extension(outfile, "pdf")
  print(pdffile)
  cairo_pdf(pdffile, onefile=TRUE)
  writeLines(pdffile, fout)
  on.exit(dev.off(), add=TRUE)

  #assert that everything we're plotting isn't from the bicontrast experiment
  with(energy_measurements,
       all(abs(content_cw + content_ccw - 0.5) < 0.01)) || stop("no")

  #plot 5 measures
  energy <- chain(energy_measurements,
                  melt(measure.vars=
                       c("energy_cw", "energy_ccw", "energy_diff",
                         "energy_total", "norm_diff")))

  #this plot is a confusing overlapping mess, but that's kind of the
  #point: direction content doesnt' predict energy.
  #if (interactive()) figure("content_energy")
  plot(ggplot(relabel_energy(subset(energy, target_number_all == 6)),
              aes(y=value, x=content,
                  color=displacement, group=displacement))
       + geom_line(alpha=0.3, size=1.5)
       + facet_grid(variable~., scales="free_y")
       + labs(title="Motion energy as a function of direction content",
              y="Energy")
       + displacement_scale_continuous_waterline
       + content_scale)

  #okay, now this illustrates the problem...
  #if(interactive()) figure("displacement_energy")
  plot(ggplot(relabel_energy(subset(energy, target_number_all == 6)),
              aes(y=value, x=displacement,
                  color=content, group=content))
       + geom_line()
       + facet_grid(variable~., scales="free_y")
       + labs(title="Motion energy as a function of displacement", y="Energy")
       + content_scale_continuous_waterline
       + geom_vline(x=0, alpha=0.3, size=0.5)
       + displacement_scale)

  #The second plot sort of gets at it, but let's also make a surface plot!

  library(rgl)
  axis_vars <- c("displacement", "content")
  axes <- displacement ~ content
  value_vars <- c("energy_cw", "energy_ccw", "energy_diff", "norm_diff")
  grid <- Map(curr(acast, subset(energy_measurements, target_number_all==6), axes),
              value.var=c(axis_vars, value_vars))

  with(grid, {
    Map(file=c("norm_energy", "energy"), data=list(norm_diff, energy_diff),
        title=c("Normalized motion energy", "Motion energy"),
        subtitle=c("(cw - ccw)/(cw + ccw)", "cw - ccw"),
        zlab="Energy",
        f=function(file, data, title, subtitle, zlab) {
          open3d(antialias=4, windowRect=c(0L, 44L, 897L, 772L))
          on.exit(rgl.close(), add=TRUE)
          bg3d(color="gray90")
          persp3d(displacement, content, data,
                  col=mapColors(data, energy.colors),
                  box=FALSE, axes=FALSE, xlab="", ylab="", zlab="")
          segments3d(x=segmentify(displacement), y=par3d("bbox")[3],
                     z=segmentify(data), add=TRUE,
                     col=mapColors(segmentify(data), energy.colors))
          segments3d(x=par3d("bbox")[2], y=segmentify(t(content)),
                     z=segmentify(t(data)), add=TRUE,
                     col=mapColors(segmentify(t(data)), energy.colors))
          axis3d(c("x--", "y+-", "x++"))
          bbox3d(col=c("white", "black"), marklen=30)
          title3d(c(subtitle, title), line=c(2,3))
          mtext3d("Envelope motion", 'x--', 2)
          mtext3d("Carrier strength", 'y+-', 2)
          mtext3d("Spacing", 'z++', 2)
          fdir <- dirname(outfile)
          fname <- file.path(fdir,paste0(file,".html"))
          writeWebGL(fdir, fname)
          writeLines(fname, fout)
        })
  })

  invisible(NULL)
}

run_as_command()
