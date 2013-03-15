suppressPackageStartupMessages({
library(plyr)
library(ggplot2)
library(ptools)
library(psyphy)
library(gnm)
library(grid)
})
theme_set(theme_bw())
use_unicode=TRUE
source("scales.R")

lack_of_pooling <- function() {
  #select from data where number of targets 6 or less. We can use an
  #ANOVA to assert lack of dependence of slope on these large element
  #numbers, and plot the slopes as a bar graph.
  
  chain(  data
        , subset(  exp_type %in% c("content", "spacing")
                 & data$target_number_all <= 6 )
        , mkrates(c("target_number_all", "spacing",
                    "subject", "content", "displacement"))
        #but we only want data where these is more than one density from a
        #subject.
        , ddply("subject", mutate, n_slopes = length(unique(target_number_all)))
        , subset(n_slopes > 1)
        ) -> smallnumber

  smallnumber$content <- factor(smallnumber$content)
  smallnumber$target_number_all <- factor(smallnumber$target_number_all)

  (  ggplot(smallnumber)
   + aes(shape = factor(target_number_all))
   + content_color_scale
   + proportion_scale
   + facet_wrap(~subject)
   + displacement_scale
   + balloon
   )

  #we want to ask, are the slopes any different among these values?
  #we'll use a basic GLM model, letting slope be constant and
  #intercept vary according to the direction content, then use an
  #ANOVA or something to assess dependence of slope on element number.

  model <- glm(  data=smallnumber
               , cbind(n_cw, n_ccw) ~ (displacement + target_number_all + content) : subject
               , family=binomial(link=logit.2asym(g=0.025, lam=0.025))
               )

  model <- glm(  data=smallnumber
               , cbind(n_cw, n_ccw) ~ (displacement:target_number_all + content) : subject
               , family=binomial(link=logit.2asym(g=0.025, lam=0.025))
               )

  #uh, extract the slopes from that. Those are the "displacement" coefficients.
  

  #the basic test is to look at the coefficients of the full regressiopn model
}
