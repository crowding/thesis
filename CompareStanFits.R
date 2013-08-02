library(scales)
library(grid)
library(rstan)
library(ptools)
library(plyr)
library(reshape2)
library(ggplot2)
source("library.R")
theme_set(theme_bw())
theme_update(panel.grid.major=element_blank(),
             panel.grid.minor=element_blank())

modelfiles <- dir(pattern=".fit.RData")

plotStanFits <- function(modelfiles=dir(pattern=".fit.RData")) {
  fitdata <- llply(modelfiles, get_samples)
  bind[samples=samples, optimized=optimized] <- collect_model_data(fitdata)
  violinPlot(samples, optimized)
  crossPlots(samples, optimized)
}

get_samples <- function(modelfile) {
  #load fit file, extract samples from each fit,
  #remember model name, return named list
  x <- load2env(modelfile)
  mutate(x$fits,
         fit=lapply(fit, as.data.frame),
         model_name=x$model@model_name)
}

unfactor <- colwise(function(x) if (is.factor(x)) as.character(x) else x)

#unpack a fitlist into two long-format data frames
collect_model_data <- function(model.list) {
  #fitlist is a list with one item per model type
  list(samples = ldply(model.list, collect_fit_data, "fit"),
       optimized = ldply(model.list, collect_fit_data, "optimized"))
}

collect_fit_data <- function(model, extract,
                             ignore = c("fit", "optimized") %-% extract) {
  #"model" is a list of lists:
  #"fit": data.frame posterior samples, one per fit
  #"optimized": max-likelihood fit, one per fit
  #and other, arbitrary identifying variables, one per fit
  model[ignore] <- NULL
  identifiers <- names(model) %-% ignore %-% extract
  chain(
    model,
    (Map %<<% .)(list), #zip
    lapply(., mkchain(
      c(., .[[extract]]), put(.[[extract]], NULL),
      data.frame %()% .,
      mutate(., .n=1:nrow(.)), #remember the sample id
      melt(id.vars=c(".n", identifiers)), #long format
      unfactor
      )),
    rbind.fill %()% .
    )
}

violinPlot <- function(samples, optimized) {
  #shows the posterior distributions over each parameter for each subject
  bind[samples, optimized] = shift_likelihoods(samples, optimized, "subject")
  print(
    ggplot(samples)
    + aes(subject, value, fill=model_name, color=model_name)
    + facet_wrap("variable", scales="free")
    + geom_violin( size=0.1, alpha=0.5
                  , position="identity")
    + geom_point(data=optimized, shape=4)
    )
}

#given a long format data frame, match each variable against the other
cross_variables <- mkchain[
    .
  , matchnames=names(samples) %-% c("value")
](
    samples=.
  , unique(.$variable)
  , expand.grid(xvar=., yvar=., stringsAsFactors=FALSE)
  , subset(xvar < yvar)
  , merge(samples, by.x="xvar", by.y="variable")
  , merge( samples
         , by.x = revalue(matchnames, c("variable"="yvar"))
         , by.y = matchnames)
)

#shift likelihoods to bring the maximum likelihood to zero...
shift_likelihoods <- function(samples, optimized, group="subject") {
  optimized.shift <- ddply(optimized, group, function(x) {
    alter(x$value[x$variable=="lp__"], . - max(.))
  })
  samples.shift <- ddply(samples, group, function(x) {
    lpshift__ <- chain(optimized, subset(variable=="lp__"),
                       match_df(., x,
                                on=names(.) %-% c(".n", "value")),
                       max(.$value))
    alter(x$value[x$variable=="lp__"], . - lpshift__)
  })
  list(samples=samples.shift, optimized=optimized.shift)
}

rescale_likelihoods <- function(samples, optimized) {
  #normalize the log probabilities by subtracting the maximum-likelihood value.
  samples <- alter(samples[which(samples$variable == "lp__"), ],
                   mutate(., value = (
                     value - merge( optimized, .
                                  , names(samples) %-% c("value", ".n")
                                  , all.y=TRUE
                                  )$value.x
               )))
  optimized <- alter(optimized[which(optimized$variable == "lp__"), ],
                     mutate(value=0))
  list(samples, optimized)
}

crossPlot <- function(samples, optimized, filter, subsample=500) {
  #shows posterior distributions as density plots over two variables;
  #input is long format data frames
  if (!missing(filter)) {
    samples <- match_df(samples, filter, on=names(filter))
    optimized <- match_df(optimized, filter, on=names(filter))
  }

  matchnames <- names(samples) %-% c("value")

  #subset samples if too many
  if (subsample < Inf) {
    samples <- ddply(samples, matchnames %-% ".n", function(x) {
      if (max(x$.n) > subsample) {
        subset(x, .n %in% sample(max(x$.n), subsample))
      } else x
    })
  }

  bind[samples, optimized] <- rescale_likelihoods(samples, optimized)

  #match samples with every pair of parameters
  cross.samples <- cross_variables(samples)
  cross.optimized <- cross_variables(optimized)

  interaction.aes.expr <- template(
    interaction( ...( lapply( matchnames %-% c(".n", "variable")
                            , as.name) ) )
  )

  print(
    ggplot(cross.samples)
   + aes(x=value.x, y=value.y)
   + geom_hline(y=0, size=0.2, color="gray50")
   + geom_vline(x=0, size=0.2, color="gray50")
   + eval(template(
        aes(color=.(interaction.aes.expr), fill=.(interaction.aes.expr))
     ))
   + scale_color_brewer("Group", type="qual", palette=3)
   + scale_fill_brewer("Group", type="qual", palette=3)
   + facet_grid(yvar ~ xvar, scales="free")
   + stat_density2d(bins=2, geom="polygon", alpha=0.5, size=0.2)
   + geom_point(shape=3, data=cross.optimized, size=2)
   + theme_grey()
   + theme(
       strip.text.x=element_text(size=rel(0.5)),
       strip.text.y=element_text(size=rel(0.5)),
     strip.background=element_blank(),
     axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0,
                              size=rel(0.6), face=2),
     axis.text.y=element_text(angle=0, hjust=1, size=rel(0.5), face=2),
     panel.grid.major=element_line(size=0.2),
     panel.grid.minor=element_blank(),
     panel.margin=unit(1, "mm"),
           aspect.ratio=1)
   + scale_y_continuous(breaks = curr(pretty, n=3))
   + scale_x_continuous(breaks = curr(pretty, n=3))
   + labs(x="", y="", title="FWHM outlines of posterior parameter distributions")
   )
}

crossPlots <- function(samples, optimized) {
  l_ply(unique(samples$model_name), function(x) {
    crossPlot(samples, optimized,
              filter=data.frame(model_name=x, stringsAsFactors=FALSE))
  })
}

main <- function(outfile, ...) {
  cairo_pdf(width=10.5, height=8, outfile, onefile=TRUE)
  on.exit(dev.off, add=TRUE)
  plotStanFits(c(...))
}

try_command <- mkchain(strsplit(" +"), unlist, .[-1:-2], main %()% .)

run_as_command()
