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

combined_model <- function(dataset) {
  fmla <- (cbind(n_cw, n_ccw) ~
           displacementTerm(spacing, displacement, start=c(cs=4, beta_dx=14))
           + content_global + content_local:full_circle + bias:factor(side)
           + I(content * abs(content)) - 1)
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
#
prediction.dataset <- chain(
  c(segment.config.vars, segment.experiment.vars, "bias", "exp_type", "side"),
  segment.trials[.],
  unique,
  recast_data)
#
combined.predictions <-
  predict_from_model_frame(combined.models,
                           newdata=recast_data(prediction.dataset))
#
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

# let's make sure I'm plotting the stuff reasonably.
chain(combined.data,
      subset(select=c("content_local", "content_global",
               "target_number_all", "target_number_shown", "full_circle")),
      unique, ggplot,
      .+aes(content_local, content_global,
            color=factor(target_number_shown),
            shape=full_circle) +
      facet_wrap(~target_number_all) + geom_point())

#let's make a prediction on all of it...

#there needs to be a "fudge factor" I think... for full curcle versus simple circle.
#because global sum only goes over a hemifield.
