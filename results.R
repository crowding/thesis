## @knitr results-libraries
suppressWarnings(source("../modeling/slopeModel.R", chdir=TRUE))
source("latexing.R")
theme_set(theme_bw(10, "MS Gothic"))
#in this file I define pretty scales.
use_unicode = TRUE
pdf.options(encoding='ISOLatin2.enc')
#options(width = 60)
#options(encoding="native.enc")=

## @knitr results-loadData
load("../modeling/slopeModel.RData")
#

## @knitr results-jb-example
model <- models$jb
rates <- subset(model$data, exp_type=="spacing")
#
#the suppressWarnings is ther because of stupid encoding things and
#the way that knitr chokes if it doesn't everything about the
#encoding. Gah, someone give me a reproducible example. for knitr +
#latex + unicode that actually works.
#plot example data from subject JB
suppressWarnings(print(ggplot(rates)
                       + displacement_scale
                       + proportion_scale
                       + content_color_scale
                       + plotPredictions(model, rates)
                       + geom_point(aes(size=n))
                       + scale_size_area(max_size=4)
                       + facet_spacing_rows
                       #  + labs(title = "Data and model fits for subject " %++% subject)
                       + no_grid
                       ))

## @knitr results-induced-modeling

#first make models with the linear induced motion
buildModel <- function(modelList, update.arg) {
  update.arg <- substitute(update.arg)
  lapply(modelList, function(model) {
    fmla <- eval(substitute(update(model$formula, arg), list(arg=update.arg)))
    gnm(formula=fmla, family=family(model), data=model$data, verbose=F)
  })
}

null.models <-
  buildModel(models,
             . ~ displacementTerm(spacing, displacement) + I(1/spacing):content)

#with a linear induced motion response
linear.models <- buildModel(null.models, .~.+content)

#then make models with a 3rd-order component to induced motion what
#might make sens is a madel that responds to leftward and rightward
#components.

#third.models <- buildModel(linear.models, .~. + I(content^3))
second.models <- buildModel(linear.models, .~. + I(abs(content) * content))
#wobble.models <- buildModel(null.models, .~. + I(wobble1(content)) + I(wobble3(content)))

#the other idea is "adding up the components differently" with the
#idea that less dominant component matters more than the more dominant
#component. how about square roots?

#and we want some models just fit all direction contents freely, to plot as points.
mutateModelData <- function(models, ...)
  llply(models, function(model) `$<-`(model,"data", mutate(model$data, ...)))

free.asym.models <- buildModel(mutateModelData(models,
                                          fContent = factor(content)),
                               . ~ displacementTerm(spacing, displacement)
                                   + I(1/spacing):content + fContent:content
                               )

#record some statistics on the likelihood ratio tests of null versus linear and linear versus second-order

## @knitr results-induced-stats
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

## @knitr results-summation-models
##The other thing I_want to do now is look at the behavior at short ranges.

#start with the 2nd order model for induced motion, deconstruct its short-range behavior
null.summation.models <- buildModel(second.models,
                                    . ~ . - I(1/spacing):content)

linear.summation.models <- buildModel(null.summation.models, . ~ . - I(1/spacing):content)

##the next thing is to justify the 1/x dependence on spacing?
                                      

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
