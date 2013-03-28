## model <- circle.models[["nj", "model"]]
library(ptools)
library(plyr)
library(ggplot2)
library(gnm)
source("density.modeling.R")
source("library.R")

infile <- "density.modeling.RData"
outfile <- "combined.model.RData"
plotfile <- "combined.model.pdf"

main <- function(infile="density.modeling.RData",
                           outfile="combined.model.RData",
                           plotfile="combined.model.pdf") {
  load(infile)
  if (interactive()) {
    plot1 <- {dev.new(); dev.cur()}
    plot2 <- {dev.new(); dev.cur()}
  } else {
    cairo_pdf(plotfile, onefile=TRUE)
    plot.dev <- dev.cur()
    detail.dev <- dev.cur()
    on.exit(dev.off(plot.dev), add=TRUE)
  }
  splits <<- c(segment.config.vars, segment.experiment.vars, "exp_type")
  if (!exists("combined.data")) {
    combined.data <- chain(
      data,
      subset(exp_type %in% c("content", "spacing", "numdensity")),
      mutate(side=factor(ifelse(exp_type %in% "numdensity", side, "all"))),
      recast_data,
      mkrates(splits %v% c("exp_type", "side")))
  }
  fitmodel <- function(dataset) {
    #dataset <- subset(dataset, full_circle==TRUE)
    fmla <- (cbind(n_cw, n_ccw) ~
             displacementTerm(spacing, displacement, start=c(cs=4, beta_dx=14))
             + content_local + content_global
             + content + I(content * abs(content)) - 1 + bias#:factor(side)
             )
    gnm(formula=fmla, data=dataset, family=binom.)
  }

  #fit a "combined" model to all subject data
  combined.models <- ddply_along(
    combined.data,
    "subject", function(vars, dataset) {
      if (!any(dataset$exp_type == "numdensity")) return(data.frame())
      print(vars)
      #pull in the old stuff and recase the data
      cm <- fitmodel(dataset)
      quickdf(c(list(model=I(list(cm)))))
    })
  # check that coefs are not NA
  coefs <- ddply(combined.models, names(combined.models) %-% "model", function(x)
                 quickdf(as.list(coef(x$model[[1]]))))
  print(coefs)

  prediction.dataset <- chain(
    c(segment.config.vars, segment.experiment.vars, "bias", "exp_type", "side"),
    subset(combined.data, exp_type="numdensity", select=.),
    unique,
    recast_data)

  combined.predictions <-
    predict_from_model_frame(combined.models,
                             newdata=recast_data(prediction.dataset))

  dev.set(2)
  plot((plot.spacing %+% segment.folded.spindled)
       + prediction_layers(combined.predictions))
  dev.set(3)
  plot((plot.spacing %+% segment.folded.spindled.mutilated
        + prediction_layers(predict_from_model_frame(
          combined.models,
          newdata=recast_data(prediction.dataset),
          fold=TRUE, spindle=TRUE, collapse=TRUE))))
  #These are still getting the "global" in the wrong direction too.
  save(file=outfile)
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

#there needs to be a "fudge factor" I think... for full curcle versus simple circle.
#because global sum only goes over a hemifield.

