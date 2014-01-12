library(plyr)
library(vadr)
library(gtable)
library(ggplot2)
source("stan_predictor.R", local=TRUE)
source("library.R")

#we might want to plot collapsed residuals
modelfile="models/d_soft_local_c_hemi_e_A.fit.RData"

main <- function(outfile="residuals.pdf",
                 ...) {
  modelfiles <- c(...)
  theme_set(theme_bw())
  cairo_pdf(outfile, width=18, height=12, onefile=TRUE)
  on.exit(dev.off())
  lapply(modelfiles, plot_model)
  NULL
}

plot_model <-  function(modelfile) {

  print(modelfile)
  model <- load_stanfit(modelfile)
  if (all(model$data$full_circle)) return()
  model_name <- paste0(model$model_name, " ", modelfile)
  model <- predictable(model)

  # what are the values that differ in the data...
  title <- textGrob(label=model_name)
  table <- gtable(heights=unit.c(2*grobHeight(title),
                                 unit(c(1,1), "null")),
                  widths=unit(1, "null"))
  table <- gtable_add_grob(table, title, 1, 1)

  qeply({
    resids <- pearson_resids(model,
                             split=c("full_circle", ..(split)),
                             fold=TRUE)
    grob <- ggplotGrob(
      ggplot(resids,
             aes(spacing, pearson_resid, size=n_obs,
                 ..(aes), weight=n_obs))
      + geom_point()
      + scale_size_area()
      + facet_grid(~full_circle, scales="free_x", space="free_x",
                   labeller = (function(var, value)
                               ifelse(value, "Full circle", "Segment")))
      + coord_trans(xtrans="sqrt", limy=c(-15, 15))
      + .(coloring)
      + geom_smooth(fill=NA))
    table <- gtable_add_grob(table, grob, .(facet), 1)
    NULL
  })(
    facet = c(2,3),
    split = list(c("spacing", "target_number_all"), alist("spacing", "subject")),
    aes = list(NULL, alist(color = factor(subject))),
    coloring = list(NULL, scale_color_brewer("Observer", type="qual", palette=3))
    )
  grid.newpage()
  grid.draw(table)

  segdata <- chain(model$data, subset(!full_circle),
                   mutate(extent = ifelse(full_circle, 2*pi * radius,
                                          spacing * target_number_shown)-1))
  segresids <- chain(
    segdata,
    pearson_resids(model=model, fold=TRUE,
                   split=c("spacing", "target_number_shown")))

  #residuals aggregated over observers
  constants <- list(
    geom_point(),
    scale_size_area(),
    theme(aspect.ratio=1),
    geom_hline(y=0, linetype="11", size=0.5, alpha=0.5),
    geom_point(),
    coord_cartesian(ylim=c(-5, 5))
    )

  table <- chain(
    gtable(widths=unit(c(1,1,1), "null"),
           heights=unit.c(grobHeight(title)*2, unit(c(1,1), "null"))),
    gtable_add_grob(title, t=1, l=1, r=3),
    gtable_add_grob(t=2, l=1, ggplotGrob(
      ggplot(segresids, aes(spacing, pearson_resid))
      + number_color_alt_scale
      + constants
      + geom_line())),
    gtable_add_grob(t = 2, l = 2, ggplotGrob(
      ggplot(segresids, aes(target_number_shown, pearson_resid))
      + spacing_texture_scale
      + constants
      + geom_line())),
    gtable_add_grob(t=2, l = 3, ggplotGrob(
      ggplot(segresids, aes(extent, pearson_resid))
      + constants
      + geom_line(aes(linetype=factor(spacing)))
      + geom_line(aes(color=factor(target_number_shown)))
      + spacing_texture_scale[[2]]
      + number_color_alt_scale[2:3])))
  qeply({
    resids <- pearson_resids(model, segdata, split=c("subject", ".(xcoord)"))
    grob <- ggplotGrob(
      ggplot(resids, aes(.(xcoord), pearson_resid, color=factor(subject)))
      + geom_point() + geom_line()
      + scale_color_brewer("Observer", type="qual", palette=3)
      + constants
      + theme(aspect.ratio=1))
    table <- gtable_add_grob(table, grob, 3, .(col))
    NULL
  })(
    xcoord=alist(spacing, target_number_shown, extent),
    col=1:3)
  grid.newpage()
  grid.draw(table)

  #plot the individual link columns aggregated over the scale of the residuals?
}

run_as_command()

if(FALSE) {
  #this actually makes no sense to do, since we're targeting distractions
  # as well as the problems with link vs. binary residual
  # and this makes weirdly weighted data.

  split <- c("full_circle", "spacing", "subject")
  resids <- pearson_resids(model, split=split, fold=TRUE)
  link_cols <- str_match_matching(names(resids), "link_.*")[,1]
  plotdata <- melt(resids, id.vars=split, measure.vars=link_cols)

  (ggplot(plotdata,
          aes(spacing, value,
              color=factor(variable)))
   + geom_line()
   + geom_point()
#   + scale_size_area()
   + facet_grid(subject~full_circle, scales="free", space="free_x",
                labeller = (function(var, value)
                            if(var == "full_circle")
                            ifelse(value, "Full circle", "Segment")
                            else as.character(value)))
   + coord_trans(xtrans="sqrt")
   + scale_color_brewer("Term", type="qual", palette=2)
   )

}
