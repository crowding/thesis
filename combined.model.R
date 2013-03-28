## model <- circle.models[["nj", "model"]]
library(ptools)
library(plyr)
library(ggplot2)
library(gnm)
source("density.modeling.R")
source("library.R")
setup_theme()

options(width=130)
infile <- "density.modeling.RData"
outfile <- "combined.model.RData"
plotfile <- "combined.model.pdf"

formula.additions <- list(
  #look, the fact that
  global1 = . ~ . + content_global + I(content_global*full_circle) ,
  #fits acceptably, while
  global2 = . ~ . + content_global + I(content_global*full_circle) ,
  #doesn't and
  global3 = . ~ . + content_global:full_circle
  #also doesn't, is disturbing. Makes me not want to trust the GNM package.
  #I could try fitting these models in Stan...
  )

#we toy with weights to tease the fitting...
segment.weight <- 1
#circle.weight <- 2

fit_combined_model <- function(dataset, relweight=segment.weight) {
  #what really isn't making sense here is why content_local seems to be useless
  #(rank deficient, leading to NA values in the coefficients)

  # Also why does content:factor(side) fit poorly?
  # Why does content_local fit poorly?

  #What's interesting, though is that I get a passable fit using only
  #content_global and ignoring content_local entirely.

  #Also it seems prediction standard errors are NOT_calculated correctly by the
  #GNM package, so forget them. Do a simulation later.

  #dataset <- subset(dataset, full_circle==TRUE)
  fmla.base <- (cbind(n_cw, n_ccw) ~
           displacementTerm(spacing, displacement, start=c(cs=6, beta_dx=10))
           + content + I(content * abs(content))
           + bias:side - 1
           )
  fmla <- update(fmla.base, formula.additions[[1]])
  gnm(  formula=fmla, data=dataset, family=binom.fam
      , weights=with(dataset, ifelse(full_circle, 1, relweight))
      )
}

coef_frame <- function(model.frame) {
  ddply(model.frame, names(model.frame) %-% "model", function(x)
        quickdf(as.list(coef(x$model[[1]]))))
}

#fit a "combined" model to all subject data
fit_combined_models <- function(combined.data) {
  combined.models <- ddply_along(
    combined.data,
    "subject", function(vars, dataset) {
      if (!any(dataset$exp_type == "numdensity")) return(data.frame())
      print(vars)
      #pull in the old stuff and recase the data
      cm <- fit_combined_model(dataset)
      quickdf(c(list(model=I(list(cm)))))
    })
}

main <- function(infile="density.modeling.RData",
                           outfile="combined.model.RData",
                           plotfile="combined.model.pdf") {
  load(infile)

  if (interactive()) {
    while(length(dev.list()) < 3) dev.new()
    plot1 <- dev.list()[2]
    plot2 <- dev.list()[2]
  } else {
    cairo_pdf(plotfile, onefile=TRUE)
    plot.dev <- dev.cur()
    detail.dev <- dev.cur()
    on.exit(dev.off(plot.dev), add=TRUE)
  }

  splits <- c(segment.config.vars, segment.experiment.vars, "exp_type")
  splits <<- splits
  if (!exists("combined.data")) {
    combined.data <- chain(
      data,
      subset(exp_type %in% c("content", "spacing", "numdensity")),
      mutate(side=factor(ifelse(exp_type %in% "numdensity", side, "all"))),
      recast_data,
      mkrates(splits %v% c("exp_type", "side")))
    combined.data <<- combined.data
  }

  combined.models <- fit_combined_models(combined.data)
  # check that coefs are not NA
  coefs <- coef_frame(combined.models)
  print(coefs)

  prediction.dataset <- chain(
    c(segment.config.vars, segment.experiment.vars, "bias", "exp_type", "full_circle", "side"),
    subset(combined.data, exp_type=="numdensity", select=.),
    unique,
    recast_data)

  combined.predictions <-
    predict_from_model_frame(combined.models,
                             newdata=prediction.dataset)

  dev.set(2)
  plot((plot.spacing %+% segment.folded.spindled)
        + prediction_layers(combined.predictions))
  dev.set(3)
  plot((plot.spacing %+% segment.folded.spindled.mutilated
        + prediction_layers(predict_from_model_frame(
          combined.models,
          newdata=prediction.dataset,
          fold=TRUE, spindle=TRUE, collapse=TRUE))))

  save(file=outfile, list=ls())
}

run_as_command()

FALSE && {
  print(do_recast(circle.models$model[[1]])$formula)
  print(informed.models$model[[1]]$formula)
  print(combined.models$model[[1]]$formula)
} #original....

# let's make sure I'm computing "local" and "global" correctly.
FALSE && {
  chain(combined.data,
        subset(select=c("content_local", "content_global",
                 "target_number_all", "target_number_shown", "full_circle")),
        unique, ggplot,
        .+aes(content_local, content_global,
              color=factor(target_number_shown),
              shape=full_circle) +
        facet_wrap(~target_number_all) + geom_point())}

#let's make a prediction on all of it...

#there needs to be a "fudge factor" I think... for full curcle versus
#simple circle.  because global sum only goes over a hemifield.

main()
