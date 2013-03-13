## @knitr results-libraries
setwd("../modeling")
suppressWarnings(source("slopeModel.R", chdir=TRUE))
setwd("../writing")
source("latexing.R")
theme_set(theme_bw(10))
#in this file I define pretty scales.
use_unicode = TRUE
pdf.options(encoding='ISOLatin2.enc')
#options(width = 60)
#options(encoding="native.enc")=

## @knitr results-loadData
load("../modeling/slopeModel.RData")
motion.energy <- read.csv("../modeling/motion_energy.csv")
bind[plot.displacement, plot.content, plot.spacing] <- (
  chain(motion.energy, subset(grid==TRUE),
        mutate(spacing=target_number_all * 2*pi/eccentricity),
        .[c("abs_displacement", "abs_direction_content", "spacing")],
        lapply(unique), lapply(sort)))
#


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


## @knitr results-spacing-example
 plotdata <- mapply(
  list(models$jb, models$nj, models$pbm),
  c(1.0, 0.4, 0.15),
  FUN=function(model, cherrypick) {
    bins <- chain(model$data,
                  subset((exp_type=="spacing") & (abs(content) %in% cherrypick)),
                  bin_along_resid(model, ., "response", splits, "displacement"))
    predictions <- makePredictions(model, bins)
    list(bins=bins, predictions=predictions)
  }, SIMPLIFY=FALSE)

bins <- rbind.fill %()% lapply(plotdata, `[[`, "bins")
predictions <- rbind.fill %()% lapply(plotdata, `[[`, "predictions")

#the suppressWarnings is ther because of stupid encoding things and
#the way that knitr chokes if it doesn't everything about the
#encoding. Gah, someone give me a reproducible example. for knitr +
#latex + unicode that actually works.
#plot example data from subject JB
suppressWarnings(print(
  (ggplot(bins)
   #+ displacement_scale
   + aes(x=displacement*10)
   + scale_x_continuous("Envelope motion (degree/sec +CW)")
#   + proportion_scale
   + aes(y=p)
   + scale_y_continuous(breaks=c(0,0.5,1))
   + content_color_scale
   + with_arg(data=predictions,
              geom_ribbon(color="transparent", alpha=0.2,
                          aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
              geom_line(aes(y=fit)))
   + geom_point(size=2)
   + facet_spacing_subject
   #  + labs(title = "Data and model fits for subject " %++% subject)
   + no_grid
   + labs(title=paste0("Carrier motion supports at narrow spacing, ",
            "repels at wide spacing\n(rows indicate spacing)"),
          y="Proportion CW responses" )
   )))

## @knitr results-induced-modeling

null.models <-
  buildModel(models,
             . ~ displacementTerm(spacing, displacement) + I(1/spacing):content)

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
free.asym.models <- buildModel(mutateModelData(models,
                                          fContent = factor(content)),
                               . ~ displacementTerm(spacing, displacement)
                                   + I(1/spacing):content + fContent:content
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
                      content=newContent, displacement=0, fContent = fContent)
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
   + content_x_scale
   + facet_wrap(~subject)
   + aes(y=fit, color=model, fill=model, ymin=fit-se.fit, ymax=fit+se.fit)
   + geom_line()
   + geom_ribbon(color=NA, alpha=0.3)
   + geom_point(data=cat.fits)
   + geom_errorbar(data=cat.fits)
#   + coord_cartesian(ylim=c(-10,10))
   + labs(y="bias",
          title=paste(
            "Estimated bias with stationary envelope\nat ",
            format(spacing, digits=3),
            "\u0080 spacing", sep=""))
   )
 }

plotBiasAtSpacing(spacing=10)
##I think the 2nd order plot is totally fine for my purposes.

## @knite results-sensitivity-modeling

## @knitr results-crossover
#show the crossover from summation to induced motion as spacing changes


#note I_also have to rip out the nonlinear term in order to get many of these to fit?
crossover.free.models <-
  buildModel(models,
             . ~ . - I(1/spacing):content - content - I(content * abs(content)) - bias
             + factor(spacing):content)


free.model <- crossover.free.models[[1]]
#what the hell, man, why isn't that working?  Might have to fake it with GLM instead of GNM.
#The next step is to extract coefs and curves. Like so:

crossover.plot.data <- rbind.fill %()% Map(
  model=models, free.model=crossover.free.models,
  f=function(model, free.model) rbind.fill(
    chain(
      free.model$data,
      count(splits %-% c("displacement", "content", "exp_type")),
      cbind(displacement=0, content=1, type="points"),
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
               "se.fit.content:I(1/spacing)"="se.fit")))))

(ggplot(subset(crossover.plot.data, type=="points"))
 + aes(x=spacing, y=fit, ymin=fit-se.fit, ymax = fit+se.fit)
 + geom_pointrange()
 + with_arg(data=subset(crossover.plot.data, type=="curve"),
            geom_line(), geom_ribbon(alpha=0.1))
 + facet_wrap(~subject, scales="free_y"))

## @knitr results-steepness

## do the same thing with the claim that the models "slope" declines.

slope.models <-
  buildModel(models, )

## @knitr results-example-effect




## Now make a graph to justify the claim that sensitivity declines with spacing.
## do it by altering the formula to remove the displacement term....


## @knitr results-summation-models
##The other thing I_want to do now is look at the behavior at short ranges.

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
