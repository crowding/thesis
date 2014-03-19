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
  library(reshape2)
  library(xtable)
  source("latexing.R")
  source("icons.R")
  source("scales.R")
  source("contours.R")
  source("stan_predictor.R")
  source("density.modeling.R")
  source("density.calibration.R")
})
for (name in ls()) {
  assign(name, get(name), globalenv())
} #coz saved function
setup_theme()
density.example.subjects <- c("pbm", "nj", "ns")

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

i.care.about <- unique(segment$subject)

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

## @knitr density-table

#haskish: encode range as complex numbers.
complex_range <- mkchain(range, complex(real=.[1], imag=.[2]))

formatter <- curr(format, digits=2)

chain(segment.folded,
      subset(abs(content) > 0 & displacement/sign(content) < 0.45 &
             ((displacement < 0.4) | (subject != "tl"))),
      ddply("subject",
            colwise(complex_range, .cols=c("displacement", "content"))),
      colwise(function(x) {
        if (is.complex(x)) {
          ifelse(Re(x)==Im(x),
                 formatter(Re(x)),
                 paste(formatter(Re(x)), formatter(Im(x)), sep=" - "))
        } else toupper(x)
      })(),
      rename(c(subject="Observer",
               displacement="Envelope displacement",
               content="Carrier strength")),
      as.matrix) -> segment.table

print(xtable(segment.table), floating=FALSE)

## @knitr density-measurements
density.example.dataset <- subset(segment.folded.spindled.mutilated,
                                  subject %in% density.example.subjects)
(plot.spacing %+% density.example.dataset
 + theme(aspect.ratio=1)
 + errorbars(density.example.dataset))

density.data.plot <-
    chain(plot.spacing %+% segment.folded.spindled.mutilated,
          + theme(aspect.ratio=1),
          + errorbars(segment.folded.spindled.mutilated),
          ggplotGrob)

##plot with spacing...

## @knitr density-extent-interaction

extent.models <- chain(
  segment,
  subset(abs(content) > 0 & displacement/sign(content) < 0.45),
  ddply("subject", summarize, model=I(list(glm(
    cbind(n_cw, n_ccw) ~ content:(spacing*target_number_shown)+content,
    family=binomial(link=logit))))))

extent.fmlas <- list(
  interaction=(
    cbind(n_cw, n_ccw) ~ content:(spacing*target_number_shown)+content),
  interaction.no.spacing=(
    cbind(n_cw, n_ccw) ~ content:spacing:target_number_shown
    + content:spacing+content),
  interaction.no.number=(
    cbind(n_cw, n_ccw) ~ content:spacing:target_number_shown
    + content:target_number_shown+content),
  full=(
    cbind(n_cw, n_ccw)
    ~ content:(factor(interaction(spacing, target_number_shown)))+content),
  null=(
    cbind(n_cw, n_ccw)
    ~ content:spacing + content:target_number_shown + content),
  full.number=(
    cbind(n_cw, n_ccw)
    ~ content:spacing + content:factor(target_number_shown) + content),
  full.spacing=(
    cbind(n_cw, n_ccw)
    ~ content:factor(spacing) + content:target_number_shown + content),
  null.number=(
    cbind(n_cw, n_ccw) ~ content:spacing + content),
  null.spacing=(
    cbind(n_cw, n_ccw) ~ content:factor(target_number_shown) + content))

extent.models <- chain(
  segment,
  subset(abs(content) > 0 & displacement/sign(content) < 0.45),
  ddply("subject", function(x) {
    quickdf(lapply(
      extent.fmlas,
      function(f)(I(list(glm(f, data=x, family=binomial(link=logit)))))))
  }))

if (FALSE) {
  extent.preds <- chain(
    extent.models,
    rename(c("interaction"="model")),
    .[c("subject", "model")],
    predict_from_model_frame(collapse=TRUE))
  (plot.number %+% segment.folded.spindled.mutilated
   + prediction_layer(extent.preds))
}

lp.frame <- function(frame) {
  colwise(function(col) {
    if(is.list(col) && ("lm" %in% class(col[[1]]))) {
      vapply(col, logLik, 0)
    } else(col)
  })(frame)
}

extent.lp.frame <- lp.frame(extent.models)
extent.lp <- numcolwise(sum)(extent.lp.frame)
pseudo_r2 <- function(frame, null, model, full) {
  (frame[[model]]-frame[[null]]) / (frame[[full]]-frame[[null]])
}

r2.interaction <- pseudo_r2(extent.lp, "null", "interaction", "full")
r2.spacing <- pseudo_r2(extent.lp, "null.spacing", "null", "full.spacing")
r2.number <- pseudo_r2(extent.lp, "null.number", "null", "full.number")

extent.coefs <- chain(
  extent.models,
  .[c("subject", "interaction")],
  rename(c("interaction"="model")),
  adply(., 1, mkchain(.$model[[1]], summary, .$coefficients,
                      as.data.frame, mutate(., coef=rownames(.), model=NA))),
  rename(c(`Pr(>|z|)`="p", "Estimate"="est", `Std. Error`="se", "z value"="z")),
  subset(coef=="content:spacing:target_number_shown"))

extract_pvals <- mkchain(
  x=.,
  drop_columns(extent.models, names(extent.fmlas) %-% .),
  rename(structure("model", names=x)),
  adply(., 1, mkchain(.$model[[1]], summary, .$coefficients,
                      as.data.frame, mutate(., coef=rownames(.), model=NA))),
  rename(c(`Pr(>|z|)`="p", "Estimate"="est", `Std. Error`="se", "z value"="z")),
  drop_columns(c("est", "se", "z", "model")),
  dcast(subject ~ coef, value.var="p"))

extract_pvals("null") -> null.pvals
extract_pvals("interaction.no.spacing")

## @knitr density-pred-modelfiles
##Here's how I should make this more sensible, plot "actual" as another row.

quad.modelfiles <- chain(
  data.frame(carrier.local=c(TRUE, TRUE, FALSE, FALSE),
             envelope.local=c(TRUE,FALSE, TRUE, FALSE)),
  mutate(carrier.name=ifelse(carrier.local, "local", "global"),
         envelope.name=ifelse(envelope.local, "local", "global")), #or make it include hemi??
  mutate(modelfile=paste0("models/d_", envelope.name, "_c_",
                          carrier.name, "_e_sA.fit.RData")))

## @knitr density-pred-load
quad.models <- chain(
  quad.modelfiles,
  mutate(model=I(lapply(modelfile, load_stanfit))))

## @knitr density-predict
unfolded.quad.preds <- chain(
  quad.models,
  adply(1, splat(function(model, ...) {
    pred <- predict(predictable(model[[1]]), segment,
                    type="response", se.fit=TRUE)
    cbind(segment, data.frame(pred), model=NA)
  }), .parallel=TRUE)
  )

quad.preds <- chain(
  unfolded.quad.preds,
  ddply(c("carrier.local", "envelope.local"),
        mutilate.predictions,
        fold=TRUE, spindle=TRUE, collapse=TRUE))

## @knitr density-predict-plot
quad.match <- merge(data.frame(subject=density.example.subjects),
                    rbind.fill(quad.conditions, data.frame(carrier.local=NA)))

(condition_prediction_plot(
  quad.preds,
  data = mutate(segment.folded.spindled.mutilated,
                carrier.local=NA, envelope.local =NA),
  match=quad.match,
  orientation="over",
  conditions=quad.conditions)
 + theme(aspect.ratio=1,
         title = element_text(size=rel(1)))
 + labs(title="Predictions for local and global models"))

all.match <- merge(data.frame(subject=unique(segment$subject)),
                    rbind.fill(quad.conditions, data.frame(carrier.local=NA)))

chain(
  condition_prediction_plot(
    quad.preds,
    data = mutate(segment.folded.spindled.mutilated,
                  carrier.local=NA, envelope.local = NA),
    match=all.match,
    orientation="down",
    conditions=quad.conditions),
  + theme_bw(10),
  + theme(aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()),
  + labs(title="Predictions for local and global models"),
  ggplotGrob) -> all.quad.prediction.plot

if(FALSE) {
  grid.newpage(); grid.draw(all.quad.prediction.plot)
}

## @knitr presentation-density-plots

print(
  condition_prediction_plot(
    quad.preds,
    data = mutate(segment.folded.spindled.mutilated,
                  carrier.local=NA, envelope.local=NA),
    match=quad.match,
    orientation="over",
    conditions=quad.conditions,
    letters=FALSE, presentation=TRUE)
  + theme(aspect.ratio=1,
          title = element_text(size=rel(1)))
  + labs(title="Predictions for local and global models"))

(condition_prediction_plot(
  subset(quad.preds, carrier.local==FALSE & envelope.local == TRUE),
  data=chain(segment.folded.spindled.mutilated,
             mutate(carrier.local=NA, envelope.local=NA),
             subset(TRUE)),
  letters=FALSE, presentation=TRUE,
  match=quad.match,
  orientation="over",
  conditions=quad.conditions)
 + theme(aspect.ratio=1))


a_ply(data.frame(carrier.local=c(FALSE, TRUE, TRUE, FALSE),
                 envelope.local=c(FALSE, TRUE, FALSE, FALSE)),
      1,
      function(match) {
        print(condition_prediction_plot(
          match_df(quad.preds, match),
          data=chain(segment.folded.spindled.mutilated,
                     mutate(carrier.local=NA, envelope.local=NA),
                     subset(FALSE)),
          letters=FALSE, presentation=TRUE,
          match=quad.match,
          orientation="over",
          conditions=quad.conditions)
              + theme(aspect.ratio=1))
      })

## @knitr density-effects

effect.coefs <- mkchain[., breaks=c()](
  mdply(function(model, ...) {
    cbind(as.data.frame(t(coefficients(model[[1]]))), model=NA)
  }),
  melt(id.vars=c("subject", breaks),
       measure.vars=c("content:spacing", "content:target_number_shown")))

model.effects <- mkchain[
  .,
  breaks=c(),
  formula=cbind(n_cw, n_ccw) ~ content:(spacing+target_number_shown)+content + side](
    subset(abs(content) > 0 & displacement/sign(content) < 0.45),
    ddply(c("subject", breaks),
          function(data) {
            model=
            glm(data=data, formula=formula, family=binomial(link="logit"))
            data.frame(model=I(list(model)))}),
  effect.coefs(breaks))

prediction.effects <- model.effects(
  unfolded.quad.preds,
  formula = fit ~ content:(spacing+target_number_shown)+content + side,
  c("carrier.local", "envelope.local"))

observed.effects <- model.effects(segment)

(ggplot(prediction.effects)
 + aes(subject, value)
 + geom_point(position=position_dodge(width=0.9),
              aes(color=interaction(carrier.local, envelope.local)),
              stat="identity")
 + facet_grid(variable~., scales="free_y")
 + geom_point(data=observed.effects))

## @knitr density-j-test

jtest <- function(a, b, data = a$data) {
  #NB. to deal with lapse rate, backconvert rates into binomial link.
  a <- predictable(a)
  b <- predictable(b)
  a.preds <- predict(a, data, type="response")
  b.preds <- predict(b, data, type="response")
  fam <- binomial(link="logit")
  a.link <- fam$linkfun(a.preds)
  b.link <- fam$linkfun(b.preds)
  a.test <- glm(cbind(data$n_cw, data$n_ccw) ~ offset(a.link) + b.link - 1,
                family=binomial(link="logit"))
  b.test <- glm(cbind(data$n_cw, data$n_ccw) ~ offset(b.link) + a.link - 1,
                family=binomial(link="logit"))
  chain(cbind(vadr::alter(summary(a.test)$coefficients, dimnames[[2]], paste0(".a")),
              vadr::alter(summary(b.test)$coefficients, dimnames[[2]], paste0(".b"))),
        .[1,],
        rename(c("Pr(>|z|).a"="p.a", "Pr(>|z|).b"="p.b")))
}

chain(quad.models,
      merge(.,., by=c(), suffixes=c(".a", ".b")),
      subset(modelfile.a < modelfile.b),
      {print(nrow(.)); .},
      mdply(function(model.a, model.b, ...) {
        cat(".")
        x <- jtest(model.a[[1]], model.b[[1]])
      }),
      drop_columns(c("model.a", "model.b"))) -> full.jtests

chain(quad.models,
      merge(.,., by=c(), suffixes=c(".a", ".b")),
      subset(modelfile.a < modelfile.b),
      {print(nrow(.)); .},
      mdply(function(model.a, model.b, ...) {
        cat(".")
        ddply(model.a[[1]]$data, "subject", jtest, a=model.a[[1]], b=model.b[[1]])
      }),
      drop_columns(c("model.a", "model.b")),
      subset(subject %in% i.care.about)) -> subject.jtests

## @knitr density-j-test-interpret

joinby <-  c("carrier.local.a", "envelope.local.a",
             "carrier.local.b", "envelope.local.b")
pr <- c(joinby, c("p.a", "p.b", "Estimate.a", "Estimate.b",
                  "modelfile.a", "modelfile.b"))

mkchain[., data=full.jtests](
  put(names, joinby),
  as.list,
  quickdf,
  join(data),
  .[pr],
  mutate(favorA = p.b / p.a)) -> get_test

get_test(c(TRUE, FALSE, FALSE, TRUE)) #favors A
get_test(c(TRUE, FALSE, TRUE, TRUE)) #favors B?
get_test(c(FALSE, FALSE, TRUE, FALSE)) #favors B

print(get_test(c(FALSE, FALSE, FALSE, TRUE), subject.jtests)) #favors A (carrier local)
print(get_test(c(TRUE, FALSE, TRUE, TRUE), subject.jtests)) #favors B (envelope local)
print(get_test(c(TRUE, FALSE, FALSE, TRUE), subject.jtests)) #favors B (envelope local)

## @knitr density-predict-likelihoods
quad.subject.lls <- chain(
  quad.models,
  adply(., 1, mkchain(.$model[[1]]$fits,
                      mutate(lp__=vapply(optimized, `[[`, 0, "lp__"),
                             model=NA, fit=NA, optimized=NA))),
  subset(subject %in% i.care.about),
  drop_columns(c("model", "fit", "optimized", "modelfile")))

envelope.better <- chain(
  quad.subject.lls,
  subset(carrier.local==FALSE),
  ddply(., c("subject"),
        mkchain(arrange(desc(lp__)), .[1,])))

carrier.better <- chain(
  quad.subject.lls,
  subset(envelope.local==TRUE),
  ddply(., c("subject"),
        mkchain(arrange(desc(lp__)), .[1,])))

quad.max.ll <- chain(
  quad.models,
      adply(1, function(row) { # for each model,
        adply(row$model[[1]]$fits, 1, function(fitrow) {
          #each subject fit w/in each mdoel
          samples <- as.data.frame(fitrow$fit[[1]])
          data.frame(lp__ = max(samples$lp__)) #get max log prob
        })
      }),
      drop_columns(c("model", "optimized"))
      )

quad.sum.lls <- chain(
  quad.max.ll,
  subset(subject %in% segment$subject),
  ddply(c("carrier.local", "envelope.local"),
        numcolwise(sum)),
  arrange(desc(lp__)))

best.subject.lls <- chain(
  quad.subject.lls,
  subset(subject %in% i.care.about),
  ddply("subject", mkchain(
    arrange(desc(lp__)),
    `[`(1,))))

## @knitr density-combined-model-plot
preds <- subset(quad.preds, (carrier.local==FALSE) & (envelope.local==TRUE))

(plot.spacing %+% segment.folded.spindled.mutilated
 + density_prediction_layers(preds, connect="number")
 )

## @density-save
save(file="density-save.RData",
     all.quad.prediction.plot,
     density.data.plot)
