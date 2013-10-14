## @knitr density-setup
library(knitr)
opts_knit$set(
  stop_on_error=2L)
opts_chunk$set(
  cache.extra=file.info(c(
    "data.RData", "slopeModel.RData",
    "numbers.RData", "density.modeling.RData", "contours.R",
    "latexing.R", "icons.R", "scales.R", "slopeModel.R",
    "density.modeling.R", "density.calibration.R", "contours.R",
    "combined.model.R", "combined.model.RData", "SlopeModel.fit.RData"))$mtime)
options(width = 70, useFancyQuotes = FALSE, digits = 4,
        lyx.graphics.center=TRUE)
options(encoding="")

suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(grid)
  library(vadr)
  source("latexing.R")
  source("icons.R")
  source("scales.R")
  source("slopeModel.R")
  source("density.modeling.R")
  source("density.calibration.R")
  source("contours.R")
})
for (name in ls()) {
  assign(name, get(name), globalenv())
} #coz saved function
setup_theme()
density.example.subjects <- c("pbm", "nj")

## @knitr do-not-run
if (!interactive()) {
  cairo_pdf(commandArgs(trailingOnly=TRUE)[1], onefile=TRUE)
}

## @knitr density-load
suppressMessages({
  circle.model <- load_stanfit("SlopeModel.fit.RData")
})

circle.data <- chain(load2env("data.RData")$data,
                     do.rename(fold=TRUE),
                     subset(subject %in% unique(circle.model$fits$subject)
                            & exp_type %in% c("content", "spacing")),
                     mkrates(circle.model$splits))

segment <- chain(load2env("data.RData")$data
                 , subset(exp_type=="numdensity"
                          & subject %in% unique(circle.model$fits$subject))
                 , do.rename(folding=TRUE))

load("numbers.RData")
#load("slopeModel.RData")
load("density.modeling.RData")

#this just illustrates the combinations of number and density.
segment.config.vars <-
  c("spacing", "target_number_shown", "target_number_all")
segment.experiment.vars <-
  c("subject", "displacement", "content", "eccentricity")
segment.splits <- c(segment.config.vars, segment.experiment.vars)
configurations <- unique(segment[segment.config.vars])
personalizations <- unique(segment[segment.experiment.vars])

## Sanity check: for each personalization, check that all
## configurations are represented.
unmatching <-
  ddply(  personalizations
        , segment.experiment.vars
        , mkchain(  match_df(segment, ., names(.))
                  , unique(.[segment.config.vars])
                  , merge(cbind(., .a=1), cbind(configurations, .a=1),
                          by=names(.))
                  , subset(is.na(.a.x) | is.na(.a.y))
                  ))
 if (!empty(unmatching)) stop("unmatching data")

## @knitr density-conditions
#choose four examples to illustrate changes of number and of density.
#lareg number/tight spacing
#med number/narrow spacing
#med number/wide spacing
#small number/wide spacing
configs <-
    chain( configurations
         , summarize( spacing = sort(unique(spacing))[c(2, 2, 5, 5)]
                    , target_number_shown =
                        sort(unique(target_number_shown))[c(6, 3, 3, 1)]
                    , label=as.character(c(4, 3, 2, 1))
                    , color=rep(TRUE, 4)
                    )
         , merge(configurations, all.y=TRUE)
         , mutate(eccentricity = 20/3))

( ggplot(configs)
 + aes(x=spacing,
       y=target_number_shown,
       fill=color)
 + geom_numdensity(aes(number=target_number_shown,
                       eccentricity=eccentricity, spacing=spacing),
                   size=5, tick_in=0)
 + identity_scale(continuous_scale("spacing", "spacing",
                                   identity, name="Spacing"))
 + scale_fill_manual(values=c("gray80"), na.value=NA)
 + scale_x_continuous(breaks=unique(configurations$spacing),
                      labels=function(x) format(x, digits=2),
                      expand=c(0, 0.4)
                      #        trans="reciprocal"
                      )
 + scale_y_continuous(expand=c(0,0.4))
 + labs(x=expression(paste("Element spacing (at 6.7"*degree, " eccentricity)")),
        y="Element number",
        title="Stimulus set for Experiment 2")
 + geom_text(aes(label=label), fontface="bold", na.rm=TRUE, size=5)
 + theme(legend.position="none", aspect.ratio=1))

## @knitr density-specifications

roundings <- c(0.2, 0.1, 2)

samplings <- chain(
  segment,
  .[c("displacement", "content", "spacing")],
  lapply(range),
  Map(f=round_any, roundings),
  Map(f=pmax, list(-Inf, 0, -Inf)),
  Map(f=function(x,r) x+c(-0.50*r, 0.50*r), roundings),
  lapply(seq_range, length=51), lapply(sort))

specification.subject <- "nj"
example <- chain(
  personalizations,
  subset(subject=="nj" & displacement==-0.1),
  put(rownames(.), NULL)
  )
example.set <- chain(
  example,
  mutate(pred=0.5),
  merge(configurations, by=c()),
  put(.["target_number_shown"], NULL),
  mutate(pred=0.5),
  unique)

elide.dim <- "displacement"

field.data <- chain(
  example,
  drop_columns(names(samplings) %-% elide.dim),
  cbind(expand.grid %()% put(samplings[elide.dim], NULL)),
  put(.[elide.dim], example[elide.dim]),
  mutate(target_number_all = 2 * pi * eccentricity / spacing,
         target_number_shown = target_number_all),
  cbind(., pred=folding_predict(
    ., model=predictable(circle.model), fold=TRUE,
    type="response")))

geom <- layer(geom="raster", geom_params=list(interpolate=TRUE))

(ggplot(field.data, aes(spacing, content, fill=pred))
 + spacing_scale_x_nopadding + content_scale_y_nopadding
 + geom + decision_contour
 + theme(aspect.ratio=1)
 + geom_numdensity(
   data=example.set,
   fill=NA, color="gray30", alpha=0.5,
   tick_in=0, tick_out=1,
   aes(spacing=spacing,
       number=target_number_all,
       eccentricity=eccentricity))
 + annotate(label=sprintf("paste('  ',Delta*x==%g)", example$displacement),
            "text", -Inf, Inf, parse=TRUE, hjust=0, vjust=1.4))

## @knitr density-measurements
density.example.dataset <- subset(segment.folded.spindled.mutilated,
                                  subject %in% density.example.subjects)
(plot.spacing %+% density.example.dataset
 + theme(aspect.ratio=1)
 + errorbars(density.example.dataset))

##plot with spacing...

## @knitr density-predictions
##Here's how I should make this more sensible, plot "actual" as another row.

(condition_prediction_plot(
  quad.predictions, segment.folded.spindled.mutilated,
  match=data.frame(subject=density.example.subjects),
  orientation="over",
  conditions=quad.conditions)
 + theme(aspect.ratio=1))

##

## and calculate deviances
unadjusted.deviances <- chain(quad.models,
                   ddply(c("subject", names(quad.conditions)),
                         summarize,
                         deviance=vapply(model, chain, 0, extractAIC, `[`(2))),
                   acast(carrier.local ~ envelope.local ~ subject, sum,
                         value.var="deviance", margins=),
                   put(names(dimnames(.)),
                       c("carrier.local", "envelope.local", "subject")))

unadjusted.winner <- chain(unadjusted.deviances,
                           apply(names(quad.conditions), sum),
                           melt(value.name="deviance"),
                           arrange(deviance))

each.unadjusted.winner <- chain(unadjusted.deviances,
                                melt(value.name="deviance"),
                                ddply(., "subject", chain,
                                      arrange(deviance), .[1,]))

adj.deviances <- chain(adj.models,
                   ddply(c("subject", names(quad.conditions)),
                         summarize,
                         deviance=vapply(model, chain, 0, extractAIC, `[`(2))),
                   acast(carrier.local ~ envelope.local ~ subject, sum,
                         value.var="deviance", margins=),
                   put(names(dimnames(.)),
                       c("carrier.local", "envelope.local", "subject")))

adj.winner <- chain(adj.deviances,
                           apply(names(quad.conditions), sum),
                           melt(value.name="deviance"),
                           arrange(deviance))

each.adj.winner <- chain(unadjusted.deviances,
                         melt(value.name="deviance"),
                         ddply(., "subject", chain,
                               arrange(deviance), .[1,]))

## @knitr density-combined-model
source("combined.model.R")
load("combined.model.RData")

## @knitr density-combined-model-plot
source("combined.model.R")
(plot_segment_fit(
  subset(selected.fits, subject %in% density.example.subjects),
  subset(prediction.dataset, subject %in% density.example.subjects),
  subset(segment.folded.spindled.mutilated, subject %in% density.example.subjects))
 + theme(aspect.ratio=1))

## @density-extent-profile
source("combined.model.R")
(make_extent_plots(
  subset(selected.fits, subject %in% density.example.subjects),
  combined.data)
 + theme(aspect.ratio=1))
