#scratchpad for testing/debugging code I'm refining elsewhere

library(vadr)
source_attach <- function(filename) {
  e <- new.env(parent=globalenv())
  source(filename, e)
  attach(e, name=filename)
}
source_attach("library.R")

A <- load2env("diagnose/fit_circle_manual.fit.RData")
B <- load2env("diagnose/fit_circle_subst.fit.RData")

ls(A)

A$stan_predict
B$stan_predict
