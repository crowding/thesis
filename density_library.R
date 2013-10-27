#functions for spindling/folding/collapsing and plotting density data.
source("scales.R")
source("icons.R")

#place a single error bar in the right half of each fact.
errorbars <- function(segment, x.axis="spacing", facet="label") {
  ddply(segment, facet, here(summarize),
        y = 0.5,
        x = max(
          switch(x.axis,
                 number=target_number_shown,
                 spacing=spacing,
                 extent=extent)),
        ymax = 0.5 + binom_se(min(n_obs), 0.5),
        ymin = 0.5 - binom_se(min(n_obs), 0.5)) -> errorbar
  with_arg(data = errorbar
           , inherit.aes = FALSE
           , mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax)
           , show_guide = FALSE, geom_errorbar(width=0.2)
           #, geom_point(size = 4, shape = "+")
           )
}

#these are the columns which define each "experiment" (facet on the
#unfolded graph)
segment.experiment.vars <-
   c("subject", "displacement", "content", "eccentricity")
#within an experiment these are the vars which separate each "stimulus
#condition" (data point on the graph)
segment.config.vars <-
  c("spacing", "target_number_shown", "target_number_all")
splits <- c(segment.config.vars, segment.experiment.vars)
segment.splits <- c(segment.config.vars, segment.experiment.vars)

## we like to plot with folded data, and with the "segment" data we,
## uh, "spindle" collapsing stimuli presented in different
## hemifields. Averaging foldings and hemifields is useful for
## plotting but not as good for modeling. "fold" collapses CW and CCW
## direction contents.  "spindle" collapses stimulus locaions.
extract_segment <- function(df, fold=FALSE, spindle=FALSE, collapse=FALSE,
                            count=TRUE, subjects=names(models),
                            splits = c(segment.config.vars,
                                       segment.experiment.vars,
                                       "eccentricity"))
  chain(df,
        subset((
            if (exists("exp_type"))
            exp_type=="numdensity"
            else (target_number_shown < target_number_all)
            ) & subject %in% subjects),
        do.rename(folding = FALSE), # we handle the folds more comprehensively.
        refold(fold = fold),
        if(count) mkrates(., c(splits, if(!spindle) "side")) else .,
        mutate(bias = if (fold) 0 else 1,
               sidedness = if (spindle) 0 else 1,
               side = if (spindle) NA else side,
               extent = spacing * target_number_shown),
        if(collapse) collapse(.) else .,
        labeler)

c.uneval <- function(...) structure(c(NULL, ...), class=class(..1))

geom_nd <- function(mapping, ...)
    list(geom_line(mapping=mapping, ...),
         geom_point(size=3, color="white", mapping=mapping, ...),
         geom_point(size=1, mapping=mapping, ...),
         geom_numdensity(fill=NA, tick_in=0, tick_out=1.5,
                         size=1.5, weight=0.3, color="gray20", linetype=1,
                         mapping = c(aes(number=target_number_shown,
                                         eccentricity=eccentricity, spacing=spacing),
                                     mapping)),
         geom_point(size=1, mapping=mapping, ...))

#a number of ggplot bits
axes.basic <- list(proportion_scale[-1]
                   , spacing_texture_scale[-1]
                   , number_color_alt_scale[-1]
                   , theme(strip.text=element_text(size=8))
                   )
#
plot.basic <- (ggplot(data.frame()) + aes(y=p)
               + axes.basic)
#
plot.wrap <- list(facet_wrap(~label))
#
by.number.nd <- list(aes(x=target_number_shown),
                     geom_nd(aes(group=factor(spacing),
                                 linetype=factor(spacing))))
by.number <- list(aes(x=target_number_shown,
                      group=factor(spacing),
                      linetype=factor(spacing),
                      shape=factor(spacing)),
                   geom_line(), geom_point(color="white", size=2.5, shape=16),
                  geom_point(size=2.5),
                  #geom_text(label="5", size=1.5, fontface="bold"),
                  labs(x="Element number"),
                  spacing_shape_scale)
#
by.spacing.nd <- list(  aes(x=spacing), labs(x="Spacing")
                   , geom_nd(aes(  group = factor(target_number_shown)
                                 , color = factor(target_number_shown))))
by.spacing <- list(aes(x=spacing,
                       group=factor(target_number_shown),
                       label=target_number_shown,
                       color=factor(target_number_shown)),
                   geom_line(), geom_point(color="white", size=3),
                   geom_text(size=2.5),
                   labs(x="Spacing"))
#
by.extent.nd <- list(  aes(x = extent)
                     , geom_line(  aes(  group = factor(spacing)
                                       , linetype = factor(spacing))
                                 , color="black", fill="black")
                     , geom_nd(aes(  group = factor(target_number_shown)
                                   , color = factor(target_number_shown)
                                   , fill = factor(target_number_shown))))
by.extent <- list(aes(x=extent),
                  geom_line(aes(  group = factor(target_number_shown)
                                , color = factor(target_number_shown)
                                , fill = factor(target_number_shown))),
                  spacing_shape_scale,
                  geom_line(aes(  group = factor(spacing)
                                , linetype = factor(spacing)
                                , shape = factor(spacing))),
                  geom_point(color="white", size=2.5, shape=16),
                  geom_point(size=2.5),
                  geom_text(aes(label=factor(target_number_shown)),
                            color="black",
                            size=1.5, fontface="bold"),
                  geom_text(aes(label=factor(target_number_shown),
                                color=factor(target_number_shown)),
                            size=1.5, fontface="bold", alpha=0.5))
#
#plot with x-axis of target number, lines of constant spacing
plot.number.nd <- plot.basic + by.number.nd + plot.wrap
plot.number <- plot.basic + by.number + plot.wrap
#
#plot with x-axis of target spacing, lines of constant number
plot.spacing.nd <- plot.basic + by.spacing.nd + plot.wrap
plot.spacing <- plot.basic + by.spacing + plot.wrap
plot.spacing.nd <<- plot.spacing.nd
plot.spacing <<- plot.spacing
#
#plot with x-axis of "extent"
plot.extent.nd <- plot.basic + by.extent.nd + plot.wrap
plot.extent <- plot.basic + by.extent + plot.wrap

#Build ggplot layers to add predictions to a plot.
density_prediction_layers <- function(data, connect = c("number","spacing"))  {
  connect <- match.arg(connect)
  eval(template(
         with_arg(
           ...(if (missing(data)) list() else list(data=quote(data))),
           mapping=aes(
             y=fit, ymin = fit - se.fit, ymax = fit + se.fit,
             ...(switch(connect,
                        number=alist(
                          color=factor(target_number_shown),
                          fill=factor(target_number_shown)),
                        spacing=alist(
                          linetype=factor(spacing))))),
           geom_line(...(if (connect=="number") list(linetype="11") else list())),
           geom_ribbon(alpha=0.3, linetype=0))))
}
