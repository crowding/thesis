#Our model has one term nonlinear in the spacing-dependent
#sensitivity to displacement and to direction content.
#This defines the response to displacement.
enclass <- function(class, x) `class<-`(x, class)

displacementTerm = enclass("nonlin", function(spacing, displacement) {
  #for an intro to "gnm" custom terms see:
  #http://statmath.wu.ac.at/research/friday/resources_WS0708_SS08/gnmTalk.pdf
  #
  #"variables" are in data. "predictors" are model parameters we introduce.
  spacing = substitute(spacing)
  displacement = substitute(displacement)
  list(  predictors = alist(cs=1, beta_dx=1) #not sure why "1"
       , variables = list(spacing, displacement)
       , term = function(predLabels, varLabels) {
         t = parse(text=c(predLabels, varLabels));
         names(t) = c("cs", "beta_dx", "spacing", "displacement")
         deparse(bquote( (2 - 2/(1 + exp( -.(t$cs) / .(t$spacing) ) ) )
                        * .(t$beta_dx) * .(t$displacement) ))
       }
       , start = function(predictors) {
         #spacing, displacement starting points.
         predictors[c(1,2)] <- c(3, 14)
         predictors
       }
       )
})

displacementFixed = function(fixCS, fixBetaDX)
  enclass("nonlin", function(spacing, displacement) {
    spacing = substitute(spacing)
    displacement = substitute(displacement)
    list(  predictors = list() #not sure why "1"
         , variables = list(spacing, displacement)
         , term = function(predLabels, varLabels) {
           t = parse(text=c(predLabels, varLabels));
           names(t) = c("spacing", "displacement")
           deparse(bquote( (2 - 2/(1 + exp( -.(fixCS) / .(t$spacing) ) ) )
                          * .(fixBetaDX) * .(t$displacement) ))
         }
         , start = identity
         )
  })


