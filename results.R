## ----------------------------------------------------------------------
## @knitr results-libraries
source("slopeModel.R", chdir=TRUE)
source("latexing.R")
theme_set(theme_bw(10))
#in this file I define pretty scales.
use_unicode = TRUE
pdf.options(encoding='ISOLatin2.enc')
#options(width = 60)
#options(encoding="native.enc")=

## @knitr do-not-run
cairo_pdf("results_plots.pdf", onefile=TRUE)

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
## @knitr results-functions
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

ziprbind <- function(l, collector=rbind.fill) Map %<<% dots(f=collector) %()% l

bound_prob <- function(x) pmax(pmin(x, 1), 0)

## ----------------------------------------------------------------------
## @knitr results-spacing-collapse

spacing.collapse.plotdata <- ziprbind(Map(
  model = list(models$pbm, models$cj),
  match_content = c(0.15, 0.4),
  f = function(model, match_content) list(
    bins = bins <- chain(
      model$data,
      subset(content == match_content), #but we really want to combine data and fold?
      #subset(abs(content) == match_content),
      # how to generate folded data from this binning function?
      bin_along_resid(model, ., "response",
                      splits %-% "exp_type", "displacement", bins=4)),
    #also refold these predictions?
    predictions = makePredictions(model, bins, splits=splits %-% "exp_type"))))

#this only includes half the available data (no folding) so the
#scatter can be expected to be improved

print(ggplot(shuffle(spacing.collapse.plotdata$bins))
      + displacement_scale
      + proportion_scale
      + spacing_color_scale
      + aes(group=spacing)
      + ribbonf(spacing.collapse.plotdata$predictions)
      + geom_point(size=3) # + binom_pointrange()
      + facet_grid(subject  ~ ., labeller=function(var, value) toupper(value))
      + no_grid
      + labs(title=paste(sep="\n","Sensitivity collapses at smaller spacings",
                   "*FIXME* half of data missing, add folded data"))
      + with_arg(inherit.aes=FALSE, show_guide=FALSE,
                 data=ddply(spacing.collapse.plotdata$bins,
                    c("subject", "content"),
                   summarize, content = content[1], n_obs = sum(n_obs)),
                 geom_text(aes(
                   label = sprintf("direction content = %0.2g", content)),
                             x=Inf, y=-Inf, size=3, vjust=-0.5, hjust=1.05),
                 geom_text(aes(
                   label = sprintf("N = %d", n_obs),
                   x = -Inf, y = Inf, size=3, vjust = 1.5, hjust = -0.2))))

## The error bars are ugly with the trials spread over so many spacings...
## I wonder if there's a way to bin over displacement and spacing at
## all the direction contents, projecting it onto one direction
## content? That is sort of massive cheating though.

## ----------------------------------------------------------------------
## @knitr results-summation-increases
summation.increases.plotdata <- ziprbind(Map(
  model = models[c("cj","pbm")],
  targnum = list(c(20, 5)),
  f = function(model, targnum) {
    bins <- chain(
      model$data
      , subset(exp_type == "content" & target_number_shown %in% targnum)
      , bin_along_resid(model, ., "response", splits, "displacement"))
    predictions <- makePredictions(model, bins)
    list(bins = bins, predictions = predictions)
  }))

print(ggplot(summation.increases.plotdata$bins)
      + displacement_scale
      + proportion_scale
      + content_color_scale
      + aes(group=content)
      + geom_point(size=2)# + binom_pointrange()
      + ribbonf(summation.increases.plotdata$predictions)
      + no_grid
      + facet_spacing_subject
      + coord_cartesian(xlim=c(-0.75, 0.75))
      + label_count(summation.increases.plotdata$bins,
                    c("spacing", "subject")))

## ----------------------------------------------------------------------
## @knitr results-induced-crossover
spacing.crossover.plotdata <- ziprbind(Map(
  model=list(models$jb, models$nj, models$pbm),
  match_content=c(1.0, 0.4, 0.15),
  f=function(model, match_content) {
    bins <- chain(
      model$data
      , subset((exp_type=="spacing") & (abs(content) == match_content))
      , bin_along_resid(model, ., "response", splits, "displacement"))
    predictions <- makePredictions(model, bins)
    list(bins=bins, predictions=predictions)
  }))

#the suppressWarnings is ther because of stupid encoding things and
#the way that knitr chokes if it doesn't everything about the
#encoding. Gah, someone give me a reproducible example. for knitr +
#latex + unicode that actually works.
#plot example data from subject JB
suppressWarnings(

  print(
    ggplot(spacing.crossover.plotdata$bins)
    + displacement_scale
    + blank_proportion_scale
    + content_color_scale
    + with_arg(data=spacing.crossover.plotdata$predictions,
               geom_ribbon(color="transparent", alpha=0.2,
                           aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
               geom_line(aes(y=fit)))
    + geom_point(size=2)
    + facet_spacing_subject
    + no_grid
    + labs(title=paste0("Carrier motion assimilates at narrow spacing, ",
             "repels at wide spacing\n(rows indicate spacing)"),
           y="Proportion CW responses" )
    + coord_cartesian(xlim=c(-1.1,1.1))
    )

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
  dropped <- drop.terms(terms(m), toDrop, keep.response=FALSE)
  f <- reformulate(labels(dropped), response="response", intercept=FALSE)
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
  glm(family=m$family, data=dat,
      formula = (response ~ offset(ofs) + content:factor(spacing) +
                 I(content * abs(content)) - 1 + bias))
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
      rename(c("fit.content:I(1/spacing)"="fit",
               "se.fit.content:I(1/spacing)"="se.fit"))))})

print(ggplot(subset(summation.plot.data, type=="points"))
      + aes(x=spacing, y=fit, ymin=fit-se.fit, ymax = fit+se.fit)
      + geom_pointrange()
      + with_arg(data=subset(summation.plot.data, type=="curve"),
                 geom_line(), geom_ribbon(alpha=0.1))
      + facet_wrap(~subject, scales="free_y"))

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

## @knite results-sensitivity-modeling



## @knitr results-example-effect

## Next, maybe we can make some surface plots etc.


## @knitr results-example-effect




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
