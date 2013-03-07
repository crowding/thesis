#Our model has one term nonlinear in the spacing-dependent
#sensitivity to displacement and to direction content.
#This defines the response to displacement.
enclass <- function(class, x) `class<-`(x, class)

## the parts that we need to *specify" for a nonlinear term are: the
## "predictors" might either be constant or expressions, and
## interpreted similarly
## How to use is by currying:
## nonlinearTerm(predictors)(expr_variables)(expr)(data_terms)

## function to build a nonlinear term
nonlinearTerm <- function(..., start=NULL) {
  predictors <- quote_args(...)
  predictors[is.missing(predictors)] <- list(1)
  function(...) {
    variables <- quote_args(...)
    #
    function(expr) {
      expr <- substitute(expr)
      eval(template(structure(
        class="nonlin",
        function(.=...(variables)) {
          ...( lapply(names(variables),
                      function(x) template( .(as.name(x))
                                           <- substitute(.(as.name(x))))))
          list( predictors = alist(...(predictors))
               , variables = alist(...(names(variables)))
               , term = function(predLabels, varLabels) {
                 t = as.list(parse(text=c(predLabels, varLabels)))
                 names(t) <- .(c(names(predictors), names(variables)))
                 deparse(substitute(.(expr), t))
               },
               #we'll bother with "start" for now, but would be base
               ...(if (is.null(start)) NULL else stop("start unsupported"))
               )
        })))
    }
  }
}

displacementTerm <- (nonlinearTerm(cs, beta_dx)(spacing, displacement)
                     ((2 - 2/(1+exp(-cs/spacing))) * beta_dx * displacement))


