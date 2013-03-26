## model <- circle.models[["nj", "model"]]

if (!exists("combined.data")) {
  combined.data <- chain(
    data,
    subset(exp_type %in% c("content", "spacing", "numdensity")),
    mutate(side=factor(ifelse(exp_type %in% "numdensity", side, "all"))),
    recast_data,
    mkrates(splits %v% c("exp_type", "side"))
    )
}

## unified_models <- function(model.frame, density.model.frame) {
##   merge(model.frame, spacing.model.frame,
##         by=names(model.frame) %-% names(density.model.frame) %-% "model",
##         suffixes=c(".circle", ".density")) -> total.model.frame
##   Map %<<<-% unify_model %()% total.model.frame
## }

## model.circle <- circle.models[["nj", "model"]]
## model.density <- descriptive.models[["nj", "model"]]

## unify_model <- function(model.circle, model.density, ...) {
##   group <- list(...)
##   mkrates(data, splits)
## }

#we also need some way of plotting the whole thing as a function of
#spacing...

combined_model <- function(dataset) {
  fmla <- (cbind(n_cw, n_ccw) ~
           displacementTerm(spacing, displacement, start=c(cs=4, beta_dx=14))
           + content_global + content_local + bias:factor(side)
           + content + I(content * abs(content)) - 1)
  gnm(formula=fmla, data=dataset, family=binom.fam)
}

#fit a "combined" model to all subject data
combined.models <- ddply_along(combined.data,
                               "subject", function(vars, dataset) {
  if (!any(dataset$exp_type == "numdensity")) return(data.frame())
  print(vars)
  #pull in the old stuff and recase the data
  cm <- combined_model(dataset)
  quickdf(c(list(model=I(list(cm)))))
})

prediction.dataset <- chain(
  c(segment.config.vars, segment.experiment.vars, "bias", "exp_type", "side"),
  segment.trials[.],
  unique,
  recast_data)

combined.predictions <-
  predict_from_model_frame(combined.models,
                           newdata=recast_data(prediction.dataset))

plot((plot.spacing %+% segment.folded.spindled)
     + prediction_layers(combined.predictions))

#let's make a prediction on all of it...
