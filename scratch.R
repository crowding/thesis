#scratchpad for testing/debugging code I'm refining elsewhere
source("modeling_resids.R")

print(modelfile)
model <- load_stanfit(modelfile)
if (all(model$data$full_circle)) return()
model_name <- paste0(model$model_name, " ", modelfile)
model <- predictable(model)

resids <- pearson_resids(
  model,
  split=c("full_circle", "spacing", "target_number_all", "subject",
          sgn="sign(content)"),
  fold=TRUE)

unfolded_resids <- pearson_resids(
  model,
  split=c("full_circle", "spacing", "target_number_all", "subject",
          sgn="sign(content)"),
  fold=FALSE)

refold <- function(data, fold=TRUE) {
  if(fold == "over") {
    fold.trial <- TRUE
  } else {
    fold.trial <- with(data, fold & ((content < 0)
                                     | (content == 0 & displacement < 0)))
  }
  fold_trials(data, fold.trial)
}

overfolded_resids <- pearson_resids(
  model,
  split=c("full_circle", "spacing", "target_number_all", "subject",
          sgn="sign(content)"),
  fold="over")

(ggplot(resids,
        aes(x=spacing,
            weight=n_obs,
            color = factor(subject)))
 + geom_point(aes(y=total_obs/n_obs, size=n_obs))
 + geom_line(aes(y=total_pred/n_obs))
 + scale_size_area()
 + facet_grid(.(sgn,full_circle), scales="free", space="free_x",
              labeller = interply(".(..1) = .(..2)"))
 + coord_trans(xtrans="sqrt")
 + scale_color_brewer("Observer", type="qual", palette=3))

#How about if we fold data twice, now, what is different?
data1 <- chain(model, predict(type="terms"))
data2 <- chain(model, predict(type="terms"))
data1 <- numcolwise(identity)(data1)
data2 <- numcolwise(identity)(data2)
sort(vapply(data2 - data1, sd, 0, USE.NAMES=TRUE))

data2 %<~% chain(refold(fold="over"), refold(fold="over"))
