## This file makes some overly ornate plots...

library(ggplot2)
library(plyr)
library(grid)
library(ptools)
library(psyphy)
source("latexing.R")
source("icons.R")
source("scales.R")
source("library.R")
source("density.rawplot.R")
setup_theme()

#let's have a colorbar for the line fits and a spa
joinedplot <- function(...) {
  cbind(ggplot_gtable(ggplot_build(segment.plot(..., x.axis="spacing"))),
        ggplot_gtable(ggplot_build(segment.plot(..., x.axis="number"))),
        size="first")
}

# now make a plot out of all that.
calibration.plot <-
  function(row, full, full.predicted, segment.predicted, ...,
           spacing_breaks = sort(unique(segment$spacing)),
           spacing_range = range(segment$spacing, full.predicted$spacing)
           ) {
    #we'll have to manually limit the scales for spacing and color, etc.
    # because the range and values differ in each dataset.
    model_spacing_breaks <- unique(full$spacing)
    number_breaks <- sort(round(2*pi*row$eccentricity/spacing_breaks), dec=TRUE)
    spacing_labels <- format(spacing_breaks, digits=2)
    (ggplot(full)
     + proportion_scale + displacement_scale
     + aes(color=spacing, fill=spacing, group=spacing)
     + geom_point(aes(size=n_obs), show_guide=FALSE)
     + scale_size_area()
     + with_arg(name="Spacing",
                trans=log_trans(),
                breaks=spacing_breaks,
                minor_breaks=model_spacing_breaks,
                limits=spacing_range,
                colours = muted(c("blue", "cyan", "yellow", "red"), l=70, c=180),
                labels = spacing_labels,
                scale_color_gradientn(),
                scale_fill_gradientn())
     + continuous_scale("number", "identity", identity_pal(),
                        breaks=number_breaks, labels=spacing_labels)
     + identity_scale(continuous_scale("spacing", "identity", identity,
                                       breaks=spacing_breaks, labels=spacing_labels))
     + geom_vline(x=row$displacement, linetype="11", color="gray50")
     + geom_vline(x=0, size=0.1)
     + geom_hline(y=0.5, size=0.1)
     + labs(title = (sprintf("%s, \u0394x = %s, C = %s",
                             toupper(row$subject),
                             format(row$displacement, digits = 2),
                             format(row$content, digits = 2))))
     + geom_line(data=full.predicted, aes(y=fit), show_guide=FALSE)
     + geom_ribbon(data=full.predicted,
                   aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit),
                   color="transparent", alpha=0.25)
     + geom_numdensity(data=segment.predicted, fill="transparent",
                       aes(y=fit, number=round(2*pi*eccentricity/spacing),
                           eccentricity=eccentricity, spacing=spacing))
     + guides(colour=guide_legend("Spacing"),
              spacing=guide_legend(
                "Spacing", override.aes=list(eccentricity=20/3)),
              number=guide_legend("Spacing"),
              fill=guide_colorbar("Spacing"),
              size="none"
              )
     + theme(legend.box="horizontal"))
  }
 
combined_plot <- function(row, segment, full, full.predicted, ...) {
  spacing_breaks <- sort(unique(segment$spacing))
  spacing_range <- range(segment$spacing, full.predicted$spacing)
  upper_plot <- calibration.plot(row=row, full=full, ...,
                                 full.predicted=full.predicted,
                                 spacing_breaks=spacing_breaks,
                                 spacing_range=spacing_range)
  lower_left_plot <- segment.plot(row, segment, x.axis="number",
                                  spacing_breaks=spacing_breaks,
                                  spacing_range=spacing_range)
  lower_right_plot <- segment.plot(row, segment, x.axis="spacing")
  lower_left_table <- ggplot_gtable(ggplot_build(lower_left_plot))
  lower_right_table <- ggplot_gtable(ggplot_build(lower_right_plot))
  lower_table <- cbind(lower_left_table, lower_right_table, size="first")
  #now how to mismatched bind these...I need to insert a column to the left
  #and stretch the middle
  #move one of the legend grobs over...
  upper_table <- ggplot_gtable(ggplot_build(upper_plot))
  move_legend <- upper_table[3,5][[1]][[1]]
  total_table <- chain(upper_table
                       , gtable_filter("^[^g]")
                       , gtable_add_grob(move_legend[,c(1,2,5)], 3, 5)
                       , gtable_add_cols(unit(1, "null"), pos=1)
                       , gtable_add_grob(move_legend[,c(1,4,5)], 3, 2)
                       , gtable_add_cols(rep(unit(1, "null"), 5), pos=5)
                       , `$<-`(., layout, mutate(.$layout, r=ifelse(l==5, r+5, r)))
                       , rbind(lower_table, size="last")
                       )
  total_table
}



`%unless_empty%` <- function(a,b) if (empty(a)) b else a

datafile <- "data.Rdata"
modelfile <- "slopeModel.RData"
output <- "density.calibration.pdf"

main <- function(datafile="data.Rdata", modelfile="slopeModel.RData",
                 output="density.calibration.pdf") {
  load(datafile)
  load(modelfile)
  cairo_pdf(output, onefile=TRUE)
  on.exit(dev.off, add=TRUE)

  segment <- chain(  data
                   , subset(exp_type=="numdensity" & subject %in% names(models))
                   , do.rename(folding=TRUE)
                   )

  segment.config.vars <-
    c("spacing", "target_number_shown", "target_number_all")
  segment.experiment.vars <-
    c("subject", "displacement", "content", "eccentricity")
  configurations <- unique(segment[segment.config.vars])
  personalizations <- unique(segment[segment.experiment.vars])

  segment.rates <-
    mkrates(  segment
            , c(  segment.config.vars, segment.experiment.vars))
  segment.rates.sided <-
    mkrates(  segment
            , c(  segment.config.vars, segment.experiment.vars
                , "side","eccentricity"))

  get_model <- function(...) {
    index <- data.frame(...)
    match_df(model.df, index,
             on=(names(index) %^% names(model.df)))$model[[1]]
  }

  ## For this figure we describe a "configuration" from an
  ## experiment. First we have a spacing-series for a subject, showing a
  ## particular direction content, with a vertical intercept-line
  ## showing the particular value of direction content the current data
  ## was collected under.

  ##First problem is that we need to to take "unfolded" data from a model
  ##and make it fold, just to show the model.

  ##collect the spacing data for each experiment, then collect model predictions.
  displacement.sampling <- chain(data, do.rename,
                                 subset(select="displacement"),
                                 unique,
                                 mutate(bias=1))

  #collect predictions for all of these
  #!# bind[calibration.row, calibration.data, calibration.predictions] <-

  calibration.dataset <- chain(
    data,
    subset(exp_type %in% c("spacing", "content")),
    do.rename(folding=TRUE),
    mkrates(splits=c(
              "displacement", "content", "spacing",
              "target_number_all", "subject", "eccentricity")),
    mutate(bias=0))

  #construct model predictions along each experiment condition
  #we are comparing "segment" data with earlier collected "calibration" data.

  #each calibration plot needs experiment condition, raw matching data,
  #raw fitted data, calibration
  calibration.data <-
    dlply_along(
      segment.rates,
      segment.experiment.vars,
      function(row, segment.data) {
        raw.calibration.data <-
          match_df(calibration.dataset, row,
                   on=segment.experiment.vars %-% "displacement")
        predicted.data <-
          chain(
            displacement.sampling
            , merge(row[!(names(row) %in% "displacement")])
            #Ah, here's a puzzle, which spacing values do I use?
            , merge(data.frame(
              spacing = unique(
                (raw.calibration.data %unless_empty% segment.rates)
                $spacing),
              bias=1))
            , cbind(., predict(get_model(row), newdata=.,
                               type="response", se.fit=TRUE)))
        #now I also predict the spacing-dependent response for model stimuli.
        predicted.segment.response <-
          chain(row,
                merge(data.frame(  spacing=unique(segment.rates$spacing)
                                 , bias=0)),
                cbind(., predict(get_model(row), newdata=.,
                                 type="response", se.fit=TRUE)))
        list(row=row, segment=segment.data, full=raw.calibration.data,
             full.predicted=predicted.data,
             segment.predicted = predicted.segment.response)
      })

  combined_plots <- lapply(calibration.data, splat(combined_plot))

  if (!interactive()) {
    invisible(lapply(combined_plots, plot))
  } else {
    plot(combined_plots[[5]])
  }
}

run_as_command()
