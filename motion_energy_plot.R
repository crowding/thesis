library(ptools)
library(plyr)
library(R.devices)
source("library.R")
source("scales.R")
library(reshape2)


add_energies <- function(data){
  cw_cols <- grep("(.*)_cw\\.\\d$", names(data), value=TRUE)
  ccw_cols <- grep("(.*)_ccw\\.\\d$", names(data), value=TRUE)
  channel_cols <- grep("\\.\\d*$", names(data), value=TRUE)
  chain(data,
        mutate(
          total_e = rowSums(data[c(cw_cols, ccw_cols)]),
          energy_cw = rowSums(data[cw_cols]) / max(total_e),
          energy_ccw = rowSums(data[ccw_cols]) / max(total_e),
          energy_diff = energy_ccw - energy_cw,
          energy_total = total_e / max(total_e),
          norm_diff = energy_diff / energy_total),
        drop_columns(c(cw_cols, ccw_cols, "total_e")),
#        drop_columns(drop_cols),
        rename(c(abs_displacement="displacement",
                 abs_direction_content="content"))
        )
}

relabel_energy <- function(data) {
  mutate(data, variable = revalue(variable, c(
                 energy_cw = paste("CW Energy", circleright),
                 energy_ccw = paste("CCW energy", circleleft),
                 energy_diff = "CW - CCW",
                 energy_total = "CW + CCW",
                 norm_diff = "(CW-CCW)/(CW+CCW)")))
}

revalue_named = function(value, x, ...) {
  revalue(value, names(x), x, ...)
}

main <- function(motion_energy_file="motion_energy.csv",
                 menergy=read.csv(motion_energy_file),
                 outfile="motion_energy.pdf"
                 ) {

  chain(menergy,
        subset(grid==1),
        add_energies
        ) -> energy_measurements

  if (!interactive()) {
    cairo_pdf(outfile)
    on.exit(dev.off(), add=TRUE)
  }

  #assert that everything we're plotting isn't from the bicontrast experiment
  with(energy_measurements,
       all(abs(content_cw + content_ccw - 0.5) < 0.01)) || stop("no")

  #plot 5 measures
  energy <- chain(energy_measurements,
                  melt(measure.vars=
                       c("energy_cw", "energy_ccw", "energy_diff",
                         "energy_total", "norm_diff")),
                  relabel_energy)

  if (interactive()) figure("content_energy")
  plot(ggplot(subset(energy, target_number_all == 6),
              aes(y=value, x=content,
                  color=displacement, group=displacement))
       + geom_line()
       + facet_grid(variable~., scales="free_y")
       + labs(title="Motion energy as a function of direction content",
              y="Energy")
       + displacement_scale_continuous_waterline
       + content_scale)

  #okay, now this illustrates the problem...
  if(interactive()) figure("displacement_energy")
  plot(ggplot(subset(energy, target_number_all == 6),
              aes(y=value, x=displacement,
                  color=content, group=content))
       + geom_line()
       + facet_grid(variable~target_number_all, scales="free_y")
       + labs(title="Motion energy as a function of displacement", y="Energy")
       + scale_color_gradientn(colours=c("yellow", "red", "black", "blue", "cyan"),
                               values=c(0, 0.4999,0.5,0.5001,1))
       + geom_vline(x=0, alpha=0.3)
       )
  if (interactive)
  
  #let's make a surface plot too...

  #weird as fuck, why is it scaling up and down so much? maybe I need
  #a surface plot. Or maybe the calculation is off?

  library(rgl)

  invisible(NULL)}

