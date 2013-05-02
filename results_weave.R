## ----------------------------------------------------------------------
## @knitr results-libraries
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
        mutate(spacing=target_number_all * 2*pi/eccentricity),
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

results.eccentricity <- chain(
  model.df$model, lapply(function(x) x$data$eccentricity),
  do.call(what=c), unique)

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

## The alternative to ths point cloud is to use binning, as commented out

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
sensitivity.free.models <- lapply(models, function(m) {
  l <- labels(terms(m))
  toDrop <- grep("displacementTerm", l)
  sensitivity.term.label <<- l[toDrop]
  length(toDrop) == 1 || stop("no")
  dropped <- drop.terms(terms(m), toDrop, keep.response=TRUE)
  f <- formula(dropped)
  f <- update(f, . ~ . + displacement:factor(spacing))
  environment(f) <- environment(m$formula)
  glm(formula=f, data=m$data, family=m$family)
})

sensitivity.plot.data <- rbind.fill %()% Map(
  model=models, free.model=sensitivity.free.models,
  f=function(model, free.model) rbind.fill(
    chain(
      free.model$data,
      count(splits %-% c("displacement", "content", "exp_type")),
      cbind(displacement=1, content=0, type="points"),
      cbind_predictions(free.model, type="terms", se.fit=TRUE),
      rename(c("fit.displacement:factor(spacing)"="fit",
               "se.fit.displacement:factor(spacing)"="se.fit"))),
    chain(
      model$data,
      count(splits %-% c("displacement", "content", "exp_type",
                         "spacing", "target_number_all", "target_number_shown")),
      cbind(displacement=1, content=0, type="curve",
            spacing=seq(2, 20, length=100)),
      cbind_predictions(model, type="terms", se.fit=TRUE),
      rename(structure(names=paste0(c("fit.", "se.fit."), sensitivity.term.label),
                       c("fit", "se.fit"))))))

print(ggplot(subset(sensitivity.plot.data, type=="points"))
      + aes(x=spacing, y=fit, ymin=fit-se.fit, ymax = fit+se.fit)
      + geom_pointrange()
      + with_arg(data=subset(sensitivity.plot.data, type=="curve"),
                 geom_line(), geom_ribbon(alpha=0.1))
      + facet_wrap(~subject, scales="free_y")
      + coord_cartesian(xlim=c(0, 10))
      + labs(title="Sensitivity to envelope motion declines below 3 deg. spacing"))

## @knitr results-spacing-summation
# ----------------------------------------------------------------------
# show the crossover from summation to induced motion as spacing changes

# The big question is why all of these end up NA when i simply use GNM.
# How does that happen?
whichTerm <- ""
summation.free.models <- lapply(models, function(m) {
  ofs <- predict(m, type="terms")
  dat <- mutate(m$data, ofs=ofs[, whichTerm <<- grep("displacementTerm", colnames(ofs))])
  f <- drop.terms(terms(m), 1:length(terms(m)), keep.response=TRUE)
  f <- update(f, . ~ . + offset(ofs) + content:factor(spacing) +
              I(content * abs(content)) - 1 + bias)
  glm(family=m$family, data=dat,
      formula = f)
})

#The next step is to extract coefs and
#curves. Individual coefs and curves from "terms" predictions of each
#model, like so:
summation.plot.data <- rbind.fill %()% Map(
  model=models, free.model=summation.free.models,
  f=function(model, free.model) {
    rbind.fill(
    chain(
      free.model$data,
      count(splits %-% c("displacement", "content", "exp_type")),
      cbind(displacement=0, content=1, type="points"),
      mutate(., ofs = predict(newdata=., model, type="terms")[,whichTerm]),
      cbind_predictions(free.model, type="terms", se.fit=TRUE),
      rename(c("fit.content:factor(spacing)"="fit",
               "se.fit.content:factor(spacing)"="se.fit"))),
    chain(
      model$data,
      count(splits %-% c("displacement", "content", "exp_type",
                         "spacing", "target_number_all", "target_number_shown")),
      cbind(displacement=0, spacing = seq(2, 20, length=100),
            content=1, type="curve"),
      cbind_predictions(model, type="terms", se.fit=TRUE),
      mutate(fit=`fit.content:I(1/spacing)` + `fit.content`,
             se.fit=`se.fit.content:I(1/spacing)`)))
  })

print(ggplot(subset(summation.plot.data, type=="points"))
      + aes(x=spacing, y=fit, ymin=fit-se.fit, ymax = fit+se.fit)
      + geom_pointrange()
      + with_arg(data=subset(summation.plot.data, type=="curve"),
                 geom_line(), geom_ribbon(alpha=0.1))
      + facet_wrap(~subject, scales="free_y")
      + labs(title="Sensitivity to carrier direction is inversely related to spacing"
             , y="Summation strength"))

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
                                               fContent = factor(content)),
                               . ~ . + fContent:content
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
      newContent <- seq(-1,1,length=100)
      fContent = newContent
    } else {
      iu <- iunique(model$data$content)
      newContent = model$data$content[iu]
      fContent = model$data$fContent[iu]
    }
    data = data.frame(subject=subject, model=type, spacing=newspacing,
                      content=newContent, displacement=0, fContent = fContent, bias=1)
    cbind(data, predict(model, data, se.fit=TRUE))
  })
}

plotBiasAtSpacing <- function(spacing) {
   curve.fits <-
   mapply(models = list(linear.models, second.models),
           type = c("linear", "2nd"),
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
  (ggplot(curve.fits)
   + content_scale
   + facet_wrap(~subject)
   + aes(y=fit, color=model, fill=model, ymin=fit-se.fit, ymax=fit+se.fit)
   + geom_line()
   + geom_ribbon(color=NA, alpha=0.3)
   + geom_point(data=cat.fits)
   + geom_errorbar(data=cat.fits)
#   + coord_cartesian(ylim=c(-10,10))
   + labs(y="bias",
          title=paste( sep="\n",
            "Induced motion as a function of direction content",
            "(stationary envelope, 10 deg. spacing)"))
   )
 }

print(plotBiasAtSpacing(spacing=10))
##I think the 2nd order plot is totally fine for my purposes.

## Now make a graph to justify the claim that sensitivity declines with spacing.
## do it by altering the formula to remove the displacement term....

## @knitr results-summation-models
##The other thing I_want to do now is look at the behavior at short ranges.
FALSE && {
  #start with the 2nd order model for induced motion, deconstruct its short-range behavior
  null.summation.models <- buildModel(second.models,
                                      . ~ . - I(1/spacing):content)

  ## @knitr results-induced-element-number

  ##the next thing is to justify the 1/x dependence on spacing?
  null.summation.models <- buildModel(models, .~. - content:I(1/spacing))

  free.summation.models <- buildModel(null.summation.models, .~. + content:I(factor(spacing)))

  ##
  #library(gridExtra)
  #do.call(grid.arrange, c(plots, ncol=2))
  #then compare them all.

  #another possible series of functions, motion energy balance of the
  #scene? divided by leftward squared plus rightward squared? Could be a
  #thing... and probably not detailed enough to worry about.
  curve( 2*x / ((1-x)^2+(x+1)^2), -1, 1)
  curve( x - 2*x / ((1-x)^2+(x+1)^2), -1, 1, add=TRUE)
  curve( abs(x) * x / 2, -1, 1, add=TRUE, col="red") 
  curve( x - abs(x) * x / 2, -1, 1, add=TRUE, col="red")
  curve( x^3/2, add=TRUE, col="blue")
  ##it would be really nifty if we could collapse the induced motion
  ##into one thing summation (opponent) versus divisive normalization
  ##(induced)

  #The thing is, this is misleading unless we can also let the other
  #function of dorection content -- the short-range summation --vary
  #freely with direction content.

}

## @knitr do-not-run
dev.off()
