## @knitr do-not-run
library(knitr)

## ----------------------------------------------------------------------
## @knitr results-libraries
library(vadr)
library(plyr)
library(lmtest)
library(reshape2)
opts_chunk$set(cache.extra=file.info(
    c("slopeModel.R", "latexing.R",
      "slopeModel.RData", "motion_energy.csv"))$mtime)
source("slopeModel.R")
source("latexing.R")
theme_set(theme_bw(10))
#in this file I define pretty scales.
use_unicode = TRUE
pdf.options(encoding='ISOLatin2.enc')
#options(width = 60)
#options(encoding="native.enc")=

## @knitr do-not-run
if (!interactive()) {
  cairo_pdf(commandArgs(trailingOnly=TRUE)[1], onefile=TRUE)
}

## ----------------------------------------------------------------------

## @knitr results-loadData
load("slopeModel.RData")
motion.energy <- read.csv("motion_energy.csv")
bind[plot.displacement, plot.content, plot.spacing] <- (
  chain(motion.energy, subset(grid==TRUE),
        mutate(spacing = target_number_all * 2*pi/eccentricity),
        .[c("abs_displacement", "abs_direction_content", "spacing")],
        lapply(unique), lapply(sort)))

## ----------------------------------------------------------------------

#results-functions
#first make models with the linear induced motion
buildModel <- function(modelList, update.arg) {
  update.arg <- substitute(update.arg)
  lapply(modelList, function(model) {
    fmla <- eval(substitute(update(model$formula, arg), list(arg=update.arg)))
    gnm(formula=fmla, family=family(model), data=model$data, verbose=F)
  })
}

mutateModelData <- function(models, ...)
  llply(models, function(model) `$<-`(model,"data", mutate(model$data, ...)))

cbind_predictions <- function(dataset, model, ...)
  cbind(dataset, predict(model, newdata=dataset, ...))

## ----------------------------------------------------------------------
## @knitr results-numbers

#browser()
results.eccentricity <- chain(
  model.df$model, lapply(function(x) x$data$eccentricity),
    do.call(what=c), unique)

bind[narrow.spacings, wide.spacings] <-
    chain(model.df, .$model,
          lapply(`[[`, "data"),
          lapply(`[[`, "spacing"),
          c %()% .,
          x=unique,
          take_nearest(c(2*pi*20/3/6, 2*pi*20/3/20)),
          factor,
          tapply(x,.,identity)
          )

## ----------------------------------------------------------------------
## @knitr results-spacing-collapse

spacing.collapse.plotdata <- ziprbind(Map(
  model = list(models$pbm, models$cj, models$jb),
  match_content = c(0.15, 0.4, 0.2),
  match_number = list(c(4, 20)
    ),
  f = function(model, match_content, match_number) list(
    bubbles = bubbles <- chain(
      model$data,
      subset(abs(content) == match_content),
      subset(target_number_all %in% match_number),
      refold(fold=TRUE),
      mkrates),
    predictions = (makePredictions(model, bubbles,
                                   splits=splits %-% "exp_type", fold=TRUE)))))

print(ggplot(shuffle(spacing.collapse.plotdata$bubbles))
      + displacement_scale
      + proportion_scale
      + coord_cartesian(xlim=c(-0.6, 0.6))
      + spacing_color_scale
      + aes(group=spacing)
      + ribbonf(spacing.collapse.plotdata$predictions)
      + balloon
      # + geom_point(size=3)
      # + binom_pointrange()
      + facet_grid( ~ subject,
                   labeller=function(var, value) sprintf("Observer %s", toupper(value)))
      + no_grid
      + theme(aspect.ratio=1)
      + labs(title=paste(sep="\n","Sensitivity to envelope motion collapses at smaller spacings"))
      + theme(plot.title=element_text(size=rel(1.0)),
              strip.background=element_rect(colour=NA, fill="gray90"))
      + with_arg(inherit.aes=FALSE, show_guide=FALSE,
                 data=ddply(spacing.collapse.plotdata$bubbles,
                    c("subject", "content"),
                   summarize, content = content[1], n_obs = sum(n_obs)),
                 geom_text(aes(
                   label = sprintf("C = %0.2g", content)),
                             x=Inf, y=-Inf, size=3, vjust=-0.5, hjust=1.05),
                 geom_text(aes(
                   label = sprintf("N = %d", n_obs),
                   x = -Inf, y = Inf, size = 6, vjust = 1.5, hjust = -0.2)))
        )

## ----------------------------------------------------------------------
## @knitr results-summation-increases

summation.increases.plotdata <- ziprbind(Map(
  model = models[c("nj","ns","pbm")],
  targnum = list(c(20, 5)),
  strength = list(c(0.1, 0.4)),
  f = function(model, targnum, strength) {
    bubbles <- chain(
      model$data
      , refold(fold=TRUE)
      , subset((exp_type == "content")
               & (target_number_shown %in% targnum)
               & (content %in% strength))
      , mkrates)
    predictions <- makePredictions(model, bubbles, fold=TRUE)
    list(bubbles = bubbles, predictions = predictions)
  }))

print(ggplot(summation.increases.plotdata$bubbles)
      + geom_hline(y=0.5, color="gray")
      + geom_vline(x=0, color="gray")
      + displacement_scale
      + proportion_scale
      + balloon
      + content_color_scale
      + labs(title="Carrier strength has more effect at reduced spacing",
             color="Carrier\nstrength")
      + aes(group=content)
      #+ geom_point(size=2)# + binom_pointrange()
      + ribbonf(summation.increases.plotdata$predictions)
      + no_grid
      + theme(plot.title=element_text(size=rel(1.0)),
              strip.background=element_rect(colour=NA, fill="gray90"),
              aspect.ratio=1)
      + facet_spacing_subject
      + coord_cartesian(xlim=c(-0.35, 0.35))
      + label_count(summation.increases.plotdata$bubbles,
                    c("spacing", "subject"),
                    size=3, y=-Inf, x=-Inf)
      )

## ----------------------------------------------------------------------
## @knitr results-spacing-sensitivity

## Illustrate the claim that the models "slope"
## declines. First, fit models that freely vary slope at each tested
## spacing.
sensitivity.term.label <- ""
sensitivity.models <- adply(model.df, 1, function(row) {
  m <- row[[1,"model"]]
  l <- labels(terms(m))
  toDrop <- grep("displacementTerm", l)
  sensitivity.term.label <<- l[toDrop]
  length(toDrop) == 1 || stop("no")
  dropped <- drop.terms(terms(m), toDrop, keep.response=TRUE)
  f <- formula(dropped)
  environment(f) <- environment(m$formula)
  #??? it should not be generating the "displacement" term (or should
  #leave another out)...
  free <- update(f, . ~ . + (displacement-1):factor(spacing) - displacement)
  flat <- update(f, . ~ . + displacement)
  mutate(row,
         free.model=I(list(glm(formula=free, data=m$data, family=m$family))),
         flat.model=I(list(glm(formula=flat, data=m$data, family=m$family))))
})

idvars <- names(model.df) %-% "model"

model_stats <- function(model.list) {
  id <- model.list[idvars]
  rownames(id) <- NULL
  chain(model.list,
        put(.[idvars], NULL),
        lapply(`[[`, 1),
        sapply(function(m)
               c(df=extractAIC(m)[[1]],
                 aic=extractAIC(m)[[2]],
                 deviance=deviance(m),
                 ll=unclass(logLik(m)),
                 nobs=unclass(attr(logLik(m), "nobs")),
                 ntrials=sum(m$data$n_obs))
               ),
        put(names(dimnames(.)), c("stat", "model.type")),
        melt,
        cbind(id))
}

#not sure what this is good for...
all_model_stats <- function(model_list) {
  x <- alply(model_list, 1, model_stats)
  rbind.fill %()% x
}

sensitivity.stats <- all_model_stats(sensitivity.models)

check <- function(x, ...) if(chain(x, ...)) x else stop("check failed")
browse <- function(x) {browser(); x}

#also run stats on the aggregation of models.
aggregate.sensitivity.stats <- chain(
    sensitivity.stats,
    dcast(list(c(idvars, "model.type"), "stat")),
    ddply(., "model.type", chain,
          subset(., select=names(.) %-% c(idvars)),
          numcolwise(sum)()),
    melt(id.vars="model.type", variable.name="stat"))

#run a log likelihood ratio test on each pair of models...
model_pair_tests <- mkchain[
    ., id = idvars[idvars %in% names(.)],
    tests=c("df", "ll", "ll.x", "aic", "chisq", "lpval",
            "pval", "ntrials", "nobs")](
    dcast(., list(c(id, "model.type"), "stat")),
    merge(.,.,by=id),
    subset(unclass(df.x) < unclass(df.y)),
    check(with(all(nobs.y==nobs.x))),
    mutate(df = df.y - df.x,
           ll = ll.y - ll.x,
           ll.x = ll.x,
           aic = aic.y - aic.x,
           chisq = 2 * abs(ll),
           lpval = pchisq(chisq, round(df), lower.tail=FALSE, log.p=TRUE),
           pval = pchisq(chisq, round(df), lower.tail=FALSE),
           ntrials = ntrials.x,
           nobs = nobs.x),
    .[c(id, "model.type.x", "model.type.y", tests)],
    melt(c(id, "model.type.x", "model.type.y"), variable.name="test"),
    acast(list(id, "test", "model.type.x", "model.type.y")),
    put(names(dimnames(.)), c("id", "test", "model.type.x", "model.type.y"))
    )

sensitivity.tests <- model_pair_tests(sensitivity.stats)
aggregate.sensitivity.tests <- model_pair_tests(aggregate.sensitivity.stats)

#here are some pseudo r-2 measures based on these tables
hosmer_lemeshow_r2 <- function(tests,
                               models=c("flat.model", "model", "free.model")) {
  (  tests[, "ll", models[1], models[2], drop=FALSE]
   / tests[, "ll", models[1], models[3], drop=FALSE])
}

#the axes are group, test stat, model 1, model 2
cox_snell_r2 <- function(tests) {
  1-exp( -2/tests[, "nobs", "flat.model", "model"]
        *  tests[ , "ll", "flat.model", "model"])
}

nagelkerke_r2 <- function(tests) {
  cox_snell_r2(tests) / (1-exp( 2/tests[, "nobs", "flat.model", "model"]
                               * tests[ , "ll.x", "flat.model", "model"]))
}

sensitivity.plot.r2 <-
    chain(sensitivity.tests,
          hosmer_lemeshow_r2,
          melt,
          rename(c(id="subject", value="r2")),
          mutate(., label=qqply("  R"[L]^2 == .(sprintf("%.2f  ", x)))(x=r2)),
          mutate(label=vapply(label, deparse, "")),
          mutate(., observer=chain(subject, toupper, paste("Observer", .))),
          unfactor)

#%plot spacing versus sensitivity?%
sensitivity.plot.data <- rbind.fill %()% (
  Map %<<% sensitivity.models
  %()% list(
    f = function(model, free.model, ...)
    rbind.fill(
      chain(
        free.model$data,
        count(splits %-% c("displacement", "content", "exp_type"), "n_obs"),
        cbind(displacement=1, content=0, bias=0, type="points"),
        cbind_predictions(free.model, type="terms", se.fit=TRUE),
        rename(c("fit.displacement:factor(spacing)"="fit",
                 "se.fit.displacement:factor(spacing)"="se.fit",
                 "freq" = "n_obs"))),
      chain(
        model$data,
        count(splits %-% c("displacement", "content", "exp_type",
                           "spacing", "target_number_all",
                           "target_number_shown"), "n_obs"),
        cbind(displacement=1, content=0, type="curve",
              spacing=seq(2, 20, length=100)),
        cbind_predictions(model, type="terms", se.fit=TRUE),
        rename(structure(names=paste0(c("fit.", "se.fit."), sensitivity.term.label),
                         c("fit", "se.fit")))))))
sensitivity.plot.data$observer <- chain(
    sensitivity.plot.data$subject, toupper, paste("Observer", .))

sensitivity.example.subjects <- c("jb", "pbm", "nj")

sensitivity_plot <- function(data=sensitivity.plot.data, stats=sensitivity.plot.r2){
  (ggplot(subset(data, type=="points"))
   + geom_hline(y=0, color="gray50")
   + aes(x=spacing, y=fit)
   + geom_point(aes(size=n_obs), color="gray50")
   + geom_segment(aes(y=fit-se.fit, yend = fit+se.fit, xend=spacing))
   + scale_size_area(breaks=c(50, 100, 200, 500, 1000))
   + with_arg(data=subset(data, type=="curve"),
              geom_line(),
              geom_ribbon(alpha=0.1, aes(ymin=fit-se.fit, ymax=fit+se.fit)))
   + facet_wrap(~ observer)
   + no_grid
   + spacing_scale_x
   + theme(aspect.ratio=1)
   + geom_text(data=match_df(stats, data, on="subject"),
               aes(x=-Inf, y=Inf, label=label), parse=TRUE,
               vjust=1.5, hjust=0, size=3)
   + labs(title="Sensitivity to envelope motion declines below 3 deg. spacing",
          y=expression(paste("Sensitivity ", beta[Delta*x])),
          size="N"))
}
sensitivity_plot(subset(sensitivity.plot.data,
                        subject %in% sensitivity.example.subjects))

save(file="sensitivity-plot.RData",
     sensitivity.plot.data, sensitivity_plot, sensitivity.tests)


## @knitr do-not-run

hosmer_lemeshow_r2(sensitivity.tests)
hosmer_lemeshow_r2(aggregate.sensitivity.tests)

cox_snell_r2(sensitivity.tests)
cox_snell_r2(aggregate.sensitivity.tests)

nagelkerke_r2(sensitivity.tests)
nagelkerke_r2(aggregate.sensitivity.tests)

sensitivity.tests[, "pval", "flat.model", "model.model"] < 0.05

## @knitr results-sensitivity-content-interactions


## @knitr results-spacing-summation
# ----------------------------------------------------------------------
# show the crossover from summation to induced motion as spacing changes

# The big question is why all of these end up NA when i simply use GNM.
# How does that happen?
whichTerm <- ""
makeSummationModel <- function(m) {
  l <- labels(terms(m))
  toDrop <- grep("1/spacing", l)
  f <- formula(drop.terms(terms(m), toDrop, keep.response=TRUE))
  f <- update(f, . ~ . + (factor(spacing) - 1):content - content)
  m2 <- update(m, family=m$family, formula = f, data=m$data)
  #the problem is for some subsets of data we were completely
  #dominated and can't sensibly fit that subset. Eliminate that subset of data
  #and try again.
  oops <- predict(m2, type="terms")[,"content:factor(spacing)"]
  update(m2, data=m2$data[oops!=0 & abs(oops) < 100,])
}
summation.free.models <- lapply(models, makeSummationModel)

#The next step is to extract coefs and
#curves. Individual coefs and curves from "terms" predictions of each
#model, like so:
summation.plot.data <- rbind.fill %()% Map(
  model=models, free.model=summation.free.models,
  f=function(model, free.model) {
    rbind.fill(
    chain(
      free.model$data,
      count(splits %-% c("displacement", "content", "exp_type"), "n_obs"),
      cbind(displacement=0, content=1, type="points"),
      cbind_predictions(free.model, type="terms"),
      mutate(n_obs = freq,
             fit=`content:factor(spacing)`,
             fit.content=0,
             `fit.content:I(1/spacing)`=0)),
    chain(
      model$data,
      count(splits %-% c("displacement", "content", "exp_type",
                         "spacing", "target_number_all", "target_number_shown")),
      cbind(displacement=0, spacing = seq(2, 20, length=100),
            content=1, type="curve"),
      cbind_predictions(model, type="terms", se.fit=TRUE), 
      mutate(fit=`fit.content:I(1/spacing)` + `fit.content`,
              se.fit=`se.fit.content:I(1/spacing)`))
      )
  })
summation.plot.data$observer <- chain(
  summation.plot.data$subject, toupper, paste("Observer", .))

summation_plot <-
  function(data=summation.plot.data) (
             ggplot(subset(data, type=="points"))
             + spacing_scale_x
             + aes(y=fit, ymin=fit-se.fit, ymax = fit+se.fit)
             + geom_point(aes(size=n_obs), alpha=0.5)
             + scale_size_area()
             + no_grid
             + with_arg(data=subset(data, type=="curve"),
                        geom_line(), geom_ribbon(alpha=0.1))
             + facet_wrap(~observer, scales="free_y")
             + labs(title="Sensitivity to carrier direction is inversely related to spacing"
                    , y=expression(paste("Carrier sensitivity " , beta[S]*M[S](S)))
                    , size="N")
             + geom_hline(y=0, alpha=0.5)
             + theme(aspect.ratio=1))

save(file="summation-plot.RData", summation_plot, summation.plot.data)

summation_plot(subset(summation.plot.data, subject %in% sensitivity.example.subjects))

## @knitr results-no-pooling

results.pooling.sample.size.threshold <- 0

results.pooling.candidates <-
  chain(
  model.df$model,
  lapply(mkchain(
    `$`(data),
    count(c("spacing", "target_number_all", "subject"), "n_obs"))),
    rbind %()% .,
    subset(freq >= results.pooling.sample.size.threshold),
    ddply("subject", subset, all(c(2, 4) %in% target_number_all)),
    subset(target_number_all %in% c(2,4)),
    arrange(subject, desc(spacing)))

results.pooling.data <- chain(
  model.df$model,
  ldply(mkchain(
    `$`(data),
    merge(drop_columns(results.pooling.candidates,"freq")))),
  mutate(target_number_all = C(
           factor(target_number_all,
                  levels=sort(unique(target_number_all))),
           contr.treatment)))

library(glm2)
results.pooling.change <-
  ddply(results.pooling.data, c("subject"), function(set) {
    if (length(unique(abs(set$content))) > 1) {
      fmla <- (cbind(n_cw, n_ccw)
               ~ displacement*target_number_all
               + sign(content):factor(abs(content)))
    } else {
      fmla <- (cbind(n_cw, n_ccw)
               ~ displacement*target_number_all + content)
    }
    chain(set,
          glm2(
            formula = fmla,
            family = binomial(link=logit.2asym(g = 0.025, lam = 0.025))),
          summary, coef, as.data.frame, mutate(., coef=rownames(.)),
          subset(coef=="displacement:target_number_all2"),
          mutate(change=ifelse(Estimate>0, "decreased", "increased"),
                 signif=`Pr(>|z|)` < 0.05))
  })

## @knitr results-induced-modeling
null.models <-
  buildModel(models,
             . ~ . - content - I(content*abs(content)))

#with a linear induced motion response
linear.models <- buildModel(null.models, .~.+content)

#then make models with a 3rd-order component to induced motion what
#might make sens is a madel that responds to leftward and rightward
#components.

second.models <- buildModel(linear.models, .~. + I(abs(content) * content))

#the other idea is "adding up the components differently" with the
#idea that less dominant component matters more than the more dominant
#component. how about square roots?

#and we want some models just fit all direction contents freely, to plot as points.
free.asym.models <- buildModel(mutateModelData(null.models,
                                               fContent = factor(abs(content))),
                               . ~ . + fContent:sign(content)
                               )

## @knitr results-induced-stats
#record some statistics on the likelihood ratio tests of null versus linear and linear versus second-order
suppressMessages(library(lmtest))
pvals <- function(models) sapply(models, `[`, 2, "Pr(>Chisq)")
linear.model.test  <- Map(lrtest, null.models, linear.models)
second.model.test <- Map(lrtest, linear.models, second.models)
complete.model.test <- Map(lrtest, null.models, second.models)
weakest.second.test <- which.max(pvals(second.model.test))
p.level <- 0.05

## @knitr results-induced-additionally
linear.slopes <- sapply(linear.models, function(x) x$coefficients[["content"]])
second.slopes <- sapply(second.models, function(x) x$coefficients[["content"]])
curve.coef.name <- "I(abs(content) * content)"
second.curves <- sapply(second.models, function(x) x$coefficients[["content"]])
carissas.wierd <- names(second.slopes[second.slopes > 0 | second.curves > 0])
#^ that is not the name of the subject. It's a seattle rock band.

## @knitr results-induced-model-plot
# here's a plot showing what the model differences are.
# this probably won't be included in the paper, though.

iunique <- function(x, ...) {
  #indices of unique elements
  i <- seq_along(x)
  setdiff(i, i[duplicated(x)])
}

#curious about how it looks? make a fake plot showing the predicted biases for stimuli with zero displacement
#over direction content at ten degree spacing
spacingData <- function(models, type, newspacing=10, categorical=TRUE) {
  ldply(models, function(model) {
    subject = unique(model$data$subject)
    if (!categorical) {
      newContent <- seq(0,1,length=100)
      fContent = newContent
    } else {
      iu <- iunique(abs(model$data$content))
      newContent = abs(model$data$content[iu])
      fContent = model$data$fContent[iu]
    }
    chain(data.frame(subject=subject, model=type, spacing=newspacing,
            content=newContent, displacement=0, fContent = fContent, bias=0),
          cbind(., folding_predict(model, ., type="response", se.fit=TRUE, fold=TRUE)),
          mutate(observer=paste("Observer", toupper(subject))))
  })
}

nonlinearity_plot <- function(..., spacing=10) {
   curve.fits <-
   mapply(models = list(linear.models, second.models),
           type = c("linear", "2nd order"),
           spacingData, MoreArgs=list(newspacing=spacing, categorical=FALSE),
           SIMPLIFY=FALSE
           )
  curve.fits <- do.call(rbind.fill, curve.fits)
  ##
  cat.fits <- spacingData(free.asym.models, "full", newspacing=spacing,
                          categorical=TRUE)
   ##
  #here is a depiction of what all the alternative induced motion
  #models.
  ##
   (ggplot(subset(curve.fits, ...))
    + content_scale
    + proportion_scale
    + facet_wrap(~observer)
    + aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)
    + with_arg(mapping=aes(color=model, fill=model),
               geom_line(),
               geom_ribbon(color=NA, alpha=0.3))
    + geom_point(data=subset(cat.fits, ...), alpha=0.5)
    + geom_errorbar(data=subset(cat.fits, ...))
#   + coord_cartesian(ylim=c(-10,10))
    + theme(aspect.ratio=1)
    + no_grid
    + labs(
      x="Carrier strength (CW)",
      title=paste( sep="\n",
        "Repulsion as a function of carrier strength",
        "(stationary envelope, 10 deg. spacing)"))
   )
 }

save(file="nonlinearity-plot.RData", spacingData, linear.models,
     second.models, free.asym.models, nonlinearity_plot, iunique)

nonlinearity_plot(subject %in% sensitivity.example.subjects)
##I think the 2nd order plot is totally fine for my purposes.

## Now make a graph to justify the claim that sensitivity declines with spacing.
## do it by altering the formula to remove the displacement term....

## @knitr results-summation-models
##The other thing I_want to do now is look at the behavior at short ranges.

## @knitr do-not-run
dev.off()
