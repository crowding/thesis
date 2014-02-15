## @knitr do-not-run
library(knitr)

## ----------------------------------------------------------------------
## @knitr results-libraries
library(vadr)
library(plyr)
library(lmtest)
library(reshape2)
library(gtable)
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
load("data.RData")
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
    x <- gnm(formula=fmla, family=family(model), data=model$data, verbose=T)
    if (is.null(x)) {
      warning("couldn't fit model for .(..1)" %#% unique(model$data$subject))
      model
    } else x
  })
}

mutateModelData <- function(models, ...)
  llply(models, function(model) `$<-`(model,"data", mutate(model$data, ...)))

cbind_predictions <- function(dataset, model, fold=TRUE, ...)
  cbind(dataset, folding_predict(model, newdata=dataset, fold=fold, ...))

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

# draw a figure illustrating carrier motion sensitivity

lookfor <- list(
  expand.grid(folded_direction_content=c(0.1, 0.2), target_number_shown=c(5, 20))
)
chain(data,
      do.rename(folding=FALSE),
      subset(exp_type %in% c("content", "spacing")),
      subset(subject %in% names(models)),
      ddply("subject", function(d) {
        for (l in lookfor) {
          matches <- match_df(l, d)
          if(nrow(matches) == nrow(l)) return(match_df(d, l))
          NULL
        }
      })) -> summation.data

chain(summation.data,
      ddply(., c("subject", "target_number_shown"),
            mkchain(
              glm(formula = response
                  ~ displacement
                  + content
                  + sign(content),
                  family=binomial(link=logit)),
              list, I, data.frame(model=.)))) -> summation.models

chain(summation.models,
      mdply(., function(model, ...) {
        chain(model[[1]],
              c(coefficients(.), sd=sqrt(diag(vcov(.)))),
              as.array, t, as.data.frame,
              mutate(model=NA))
      })) -> summation.sensitivities

sensitivities.plot <- (
  ggplot(
    summation.sensitivities,
    aes(x = target_number_shown,
        y = content,
        ymin = content - sd.content,
        ymax = content + sd.content,
        group = subject,
        label = toupper(subject)))
  + scale_y_continuous("Carrier sensitivity")
  + scale_x_continuous("Number of elements", breaks=c(5, 20), limits=c(2,23))
  + geom_hline(yintercept=0, color="gray")
  + geom_point()
  + geom_line(color="black")
  + guides(color="none")
  + geom_text(data=subset(summation.sensitivities, target_number_shown==5),
              color="black", hjust=1.3, size=2.5)
  + theme(aspect.ratio=2)
  )

summation.intercepts <- chain(
  summation.models,
  mdply(., function(model, ...) {
    chain(expand.grid(displacement=0,
                      content=c(0.1,0.2)),
          cbind(., folding_predict(model[[1]], fold=TRUE,
                                   newdata=., se.fit=TRUE,
                                   type='response')),
          mutate(model=NA))
  }),
  dcast(subject+target_number_shown+displacement~abs(content), value.var="fit"))

summation.curves <- chain(
  summation.models,
  mdply(., function(model, ...) {
    chain(seq(-0.5, 0.5, len=100),
          expand.grid(displacement=., content=c(0.1, 0.2), model=NA),
          cbind(., folding_predict(model[[1]], fold=TRUE,
                                   newdata=., se.fit=TRUE,
                                   type='response')))
  }))

summation.bubbles <- chain(summation.data, refold(fold=TRUE), mkrates)
summation.examples <- data.frame(subject=c("nj", "pbm"))

summation.plot <- (
  ggplot(data=match_df(summation.curves, summation.examples))
  + aes(x=displacement, y=fit, yend=NA_real_,
        color=factor(content), group=factor(content))
  + (geom_point(
    data=match_df(summation.bubbles, summation.examples),
    aes(size=n_obs, y=p)))
  + scale_size_area("Trials", max_size=3)
  + scale_color_discrete("Carrier strength")
  + scale_y_continuous("P(Response CW)", breaks=c(0,0.5,1))
  + scale_x_continuous("Envelope displacement", breaks=c(-0.25, 0, 0.25),
                       limits=c(-0.4, 0.4))
  + geom_line()
  + (geom_segment(
    data=match_df(summation.intercepts, summation.examples),
    color="black",
    aes(group=NA, xend=displacement, y=`0.1`, yend=`0.2`),
    arrow=arrow(length=unit(0.05, "inches"), type="closed")))
  + facet_grid(subject ~ target_number_shown,
               labeller=function(var, value) {
                 switch(var,
                        target_number_shown=paste(value, "elements"),
                        subject = paste("Observer", toupper(value)))
               })
  + theme(aspect.ratio=1,
          legend.position="top"))

tab <- gtable(widths=unit(c(1.4, 1), "null"), heights=unit(c(1), "null"))
tab %<~% gtable_add_grob(ggplotGrob(summation.plot), 1, 1)
tab %<~% gtable_add_grob(ggplotGrob(sensitivities.plot), 1, 2)
grid.draw(tab)

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
                 nobs=length(m$residuals), #unclass(attr(logLik(m), "nobs")),
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
aggregate_stats <- mkchain(
    dcast(list(c(idvars, "model.type"), "stat")),
    ddply(., "model.type", chain,
          subset(., select=names(.) %-% c(idvars)),
          numcolwise(sum)()),
    melt(id.vars="model.type", variable.name="stat")
    )

aggregate.sensitivity.stats <- aggregate_stats(sensitivity.stats)

#run a log likelihood ratio test on each pair of models...
model_pair_tests <- mkchain[
    ., id = idvars[idvars %in% names(.)],
    tests=c("df", "ll", "ll.x", "aic", "chisq", "lpval",
            "pval", "ntrials", "nobs")
    ](
        dcast(., list(c(id, "model.type"), "stat")),
        merge(.,.,by=id),
        subset(df.x <= df.y & (nobs.y == nobs.x)),
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
        put(names(dimnames(.)), c("id", "test", "model.type.x", "model.type.y")))

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

label_r2 <- mkchain(
    hosmer_lemeshow_r2,
    melt,
    rename(c(id="subject", value="r2")),
    mutate(., label=qqply("  R"[L]^2 == paste(.(sprintf("%.2f", x)), "  "))(x=r2)),
    mutate(label=vapply(label, deparse, "")),
    mutate(., observer=chain(subject, toupper, paste("Observer", .))),
    unfactor
    )

sensitivity.plot.r2 <- label_r2(sensitivity.tests)

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
    sensitivity.plot.data$subject, toupper, paste("Observer", . ))

sensitivity.example.subjects <- c("jb", "pbm", "nj")

error.segment <- geom_segment(aes(y=fit-se.fit, yend = fit+se.fit, xend=spacing))

r2_label <- function(data, labels, x=-Inf, y=Inf, vjust=1.5, hjust=0) {
  geom_text(data=match_df(labels, data, on="subject"),
            aes(label=label), parse=TRUE,
            vjust=vjust, hjust=hjust, size=3, x=x, y=y)
}


xaxis = geom_hline(y=0, size=0.2, linetype="11")

sensitivity_plot <- function(data=sensitivity.plot.data, stats=sensitivity.plot.r2){
  (ggplot(subset(data, type=="points"))
   + aes(x=spacing, y=fit)
   + geom_point(aes(size=n_obs), color="gray50")
   + error.segment
   + scale_size_area(breaks=c(50, 100, 200, 500, 1000))
   + with_arg(data=subset(data, type=="curve"),
              geom_line(),
              geom_ribbon(alpha=0.1, aes(ymin=fit-se.fit, ymax=fit+se.fit)))
   + facet_wrap(~ observer)
   + no_grid
   + spacing_scale_x
   + theme(aspect.ratio=1)
   + r2_label(data, stats)
   + xaxis
   + labs(title="Sensitivity to envelope motion declines at small spacings",
          y=expression(paste("Displacement sensitivity ", beta[Delta])),
          size="N"))
}
sensitivity_plot(subset(sensitivity.plot.data,
                        subject %in% sensitivity.example.subjects))

save(file="sensitivity-plot.RData",
     sensitivity.plot.data, sensitivity_plot, sensitivity.tests,
     sensitivity.plot.r2, error.segment, r2_label, xaxis)


## @knitr do-not-run

hosmer_lemeshow_r2(sensitivity.tests)
hosmer_lemeshow_r2(aggregate.sensitivity.tests)

sensitivity.tests[, "pval", "flat.model", "model"] < 0.05

## @knitr results-sensitivity-content-interactions


## @knitr results-spacing-summation
# ----------------------------------------------------------------------
# show the crossover from summation to induced motion as spacing changes

makeSummationModel <- function(row) {
  bind[model=bind[model], ...=group] <- row
  #eliminate the nonlinear term
  toSwap <- grep("content \\*", labels(terms(model)))
  newdata <- mutate(model$data, content.nonlin = content)
  fmla <- chain(terms(model),
                drop.terms(toSwap, keep.response=TRUE),
                formula,
                update(., . ~ . + I(content.nonlin * abs(content.nonlin))))
  model <- update(model, data=newdata, formula=fmla)
  toDrop <- grep("1/spacing", labels(terms(fmla)))
  flat.fmla <- formula(drop.terms(terms(fmla), toDrop, keep.response=TRUE))
  free.fmla <- update(flat.fmla, . ~ . + (factor(spacing) - 1):content - content)
  bind[flat.model, free.model] <- lapply(
      list(flat.fmla, free.fmla),
      function(fmla) update(model, family=model$family, formula=fmla, data=model$data))
  #the problem is for some subsets of data we were completely
  #dominated by carrier motion and we can't sensibly fit that
  #subset. Eliminate that subset of data and try again.
  oops <- predict(free.model, type="terms")[,"content:factor(spacing)"]
  censored.data = model$data[oops!=0 & abs(oops) < 100,]
  free.model2 <- update(free.model, data=censored.data)
  flat.model2 <- update(flat.model, data=censored.data)
  model2 <- update(model, data=censored.data)
  quickdf(c(group,
            lapply(list(uncensored.flat.model=flat.model,
                        uncensored.model=model,
                        uncensored.free.model=free.model,
                        flat.model=flat.model2,
                        model=model2,
                        free.model=free.model2),
                   mkchain(list, I))))
}
summation.models <- adply(model.df, 1, makeSummationModel)

summation.stats <- all_model_stats(summation.models)
aggregate.summation.stats <- aggregate_stats(summation.stats)

#tests between unconsored flat model and uncensored model are valid, because
#the problem that forces censoring is in the free model.
summation.tests <- model_pair_tests(summation.stats)
aggregate.summation.tests <- model_pair_tests(aggregate.summation.stats)

summation.plot.data <- rbind.fill %()% Map(
  model=summation.models$uncensored.model, free.model=summation.models$free.model,
  f=function(model, free.model) rbind.fill(
    chain(
      #the points
      free.model$data,
      count(splits %-% c("displacement", "content", "exp_type"), "n_obs"),
      cbind(displacement=0, content=1, content.nonlin = 0, type="points"),
      cbind_predictions(free.model, type="link", se.fit=TRUE),
      mutate(n_obs = freq)),
    chain(
      #the lines
      model$data,
      count(splits %-% c("displacement", "content", "exp_type",
                         "spacing", "target_number_all", "target_number_shown")),
      cbind(displacement=0, spacing = seq(1.5, 24, length=100),
            content=1, content.nonlin=0, type="curve"),
      cbind_predictions(model, type="link", se.fit=TRUE))
    )
  )
summation.plot.data$observer <- chain(
  summation.plot.data$subject, toupper, paste("Observer", .))

summation.plot.r2 <- label_r2(summation.tests)

summation_plot <- function(data=summation.plot.data, labels=summation.plot.r2) (
    ggplot(subset(data, type=="points"))
    + spacing_scale_x
    + aes(y=fit)
    + geom_point(aes(size=n_obs), alpha=0.5)
    + error.segment
    + scale_size_area()
    + no_grid
    + with_arg(data=subset(data, type=="curve"),
               geom_ribbon(alpha=0.1,
                           aes(ymin=fit-se.fit, ymax=fit+se.fit)),
               geom_line())
    + facet_wrap(~observer, scales="free_y")
    + r2_label(data, labels, x=Inf, hjust=1)
    + labs(title="Sensitivity to carrier motion is inversely related to spacing"
           , y=expression(paste("Carrier sensitivity " , beta[S] + beta[I][a]))
           , size="N")
    + xaxis
    + theme(aspect.ratio=1))
summation_plot(subset(summation.plot.data,
                      subject %in% sensitivity.example.subjects))

save(file="summation-plot.RData", summation_plot, xaxis,
     summation.plot.data, summation.plot.r2, error.segment, r2_label)

## @knitr do-not-run

summation_plot()

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
repulsion.models = mutate(
    model.df,
    model=lapply(model, function(m)update(m, data=subset(m$data, content!=0))),
    flat.model = buildModel(model, . ~ . - content - I(content*abs(content))),
    linear.model = buildModel(flat.model, .~.+content),
    free.model = buildModel(flat.model,
                            . ~ . + sign(content):factor(abs(content))))

argh <- function(model) {
  newdata <- chain(it=model, predict(type="terms"),
                   .[,1], data.frame(offs=.), cbind(it$data))
  which.term <- chain(model, terms, labels, grep(pattern="displacement"))
  newformula <- chain(model, terms,
                      drop.terms(which.term, keep.response=TRUE), formula,
                      update(., . ~ . + offset(offs)))
  glm(data=newdata, formula=newformula, family=model$family)
}
repulsion.models <- mutate(repulsion.models, free.model=lapply(free.model, argh))

#not used?
narrowed.repulsion.models <- adply(repulsion.models, 1, function(row) {
  ix <- row[idvars]
  row[idvars] <- NULL
  if(any(is.na(coef(row$free.model[[1]])))) {
    model <- row$free.model[[1]]
    chain(model, model.matrix, .[, is.na(coef(model)), drop=FALSE],
          alply(2, `==`, 0), Reduce(f=`&`), .&model$data$content != 0) -> kept
    row <- lapply(row, function(m) {
      model <- m[[1]]
      list(update(model, formula=model$formula, data=model$data[kept,]))
    })
  }
  quickdf(as.list(c(ix, row)))
})



## @knitr results-induced-stats
#record some statistics on the likelihood ratio tests of null versus
#linear and linear versus second-order
#repulsion.stats <- all_model_stats(narrowed.repulsion.models)
repulsion.stats <- all_model_stats(repulsion.models)
aggregate.repulsion.stats <- aggregate_stats(repulsion.stats)
repulsion.tests <- model_pair_tests(repulsion.stats)
aggregate.repulsion.tests <- model_pair_tests(aggregate.repulsion.stats)

linear.model.test <- repulsion.tests[, "pval", "flat.model", "linear.model"]
second.model.test <- repulsion.tests[, "pval", "linear.model", "model"]
complete.model.test <- repulsion.tests[, "pval", "flat.model", "model"]
weakest.second.test <- which.max(second.model.test)
p.level <- 0.05

## @knitr results-induced-additionally
getCoefs <- function(models, whatcoef)
  sapply(models, function(x)
         if (whatcoef %in% names(x$coefficients)) {
           x$coefficients[[whatcoef]]
         } else NA)

linear.slopes <- getCoefs(repulsion.models$linear.model, "content")
second.slopes <- getCoefs(repulsion.models$model, "content")
curve.coef.name <- "I(content * abs(content))"
second.curves <- getCoefs(repulsion.models$model, curve.coef.name)
carissas.wierd <- names(second.slopes[second.slopes > 0 | second.curves < 0])
#^ that is not the name of the subject. It's a seattle rock band.

## @knitr do-not-use


## @knitr results-induced-model-plot

repulsion.plot.data <- rbind.fill %()% (
  Map %<<% repulsion.models
  %()% list(
    f = function(model, linear.model, free.model, ...)
    rbind.fill(
      chain(
          free.model$data,
          refold,
          count(splits %-% c("displacement", "exp_type", "spacing",
                             "target_number_all", "target_number_shown"), "n_obs"),
          cbind(offs=0, displacement=0, spacing=1e6, bias=0,
                type="points", model.type="full"),
          #cbind(., folding_predict(free.model, ., type="response", se.fit=TRUE, fold=TRUE)),
          cbind_predictions(free.model, type="link", se.fit=TRUE),
          rename(c("freq" = "n_obs"))),
      chain(
        linear.model$data,
        count(idvars, "n_obs"),
        merge(data.frame(offs=0, displacement=0, content=seq(0, 1, length=100),
                         type="curve", spacing=1e6, bias=0, model_type="linear"),
              by=c()),
        cbind_predictions(linear.model, type="link", se.fit=TRUE)),
      chain(
        model$data,
        count(idvars, "n_obs"),
        merge(data.frame(offs=0, displacement=0, content=seq(0, 1, length=100),
                         type="curve", spacing=1e6, bias=0,
                         model_type="second order"),
              by=c()),
        cbind_predictions(model, type="link", se.fit=TRUE))
        )))
repulsion.plot.data$observer <- chain(
    repulsion.plot.data$subject, toupper, paste("Observer", .))
repulsion.plot.r2 <- label_r2(repulsion.tests)

repulsion_plot <- function(data=repulsion.plot.data, labels=repulsion.plot.r2) {
   (ggplot(subset(data, type == "points"))
    + aes(x=content)
    + facet_wrap(~observer)
    + aes(y=fit)
    + with_arg(mapping=aes(color=model_type, fill=model_type,
                           ymin=fit-se.fit, ymax=fit+se.fit),
               data=subset(data, type=="curve"),
               geom_line(),
               geom_ribbon(color=NA, alpha=0.3))
    + geom_pointrange(aes(ymin=fit-se.fit, ymax=fit+se.fit),
                      alpha=0.5)
#   + coord_cartesian(ylim=c(-10,10))
    + theme(aspect.ratio=1)
    + no_grid
    + r2_label(data, labels, x=Inf, hjust=1)
    + scale_x_continuous(breaks=c(0,0.5,1))
    + xaxis
    + labs(
        x="Carrier strength",
        color="Model type", fill="Model type",
        y=expression(paste("Carrier repulsion ", beta[I])),
      title=paste( sep="\n",
        "Repulsion as a function of carrier strength")))
 }
repulsion_plot(data=subset(
    repulsion.plot.data, subject %in% sensitivity.example.subjects))

save(file="repulsion-plot.RData", repulsion.plot.data, repulsion.plot.r2,
     repulsion_plot, xaxis, r2_label, error.segment, xaxis)


## @knitr do-not-run
latex.format(hosmer_lemeshow_r2(aggregate.repulsion.tests), digits=2)

## @knitr results-summation-models
##The other thing I_want to do now is look at the behavior at short ranges.

## @knitr do-not-run
dev.off()
