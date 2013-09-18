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


a <- function(N) {
  x <- 0
  for (i in 1:N) {
    x <- x + rnorm(1, mean=i)
  }
  x
}
0
#this function explodes when enableJIT(3) is on
#this function explodes when enableJIT(3) is on

options(enableJIT=3)

b <- qe(function(N) .(`{`)(
    .(`<-`)(x,0),
    .(`for`)(
        i,
        .(`:`)(1, N),
        .(`{`)(
            .(`<-`)(
                x,
                .(`+`)(
                    x,
                    .(`rnorm`)(1, mean=i))))),
    x))

microbenchmark(
    a(1000), b(1000))

microbenchmark(
    a(1000), b(1000))


fun <- macro(function(...) {
  args <- list(...)
  anames = all.vars(substitute(x(...)))
  qq(function(.=..(
      qqply(`.(name)`=parent.env(environment())$`.(name)`
            )(name=anames))) {
    ..(args)
  })
})



if(FALSE) {

  grid <- diag(5)
  data <- melt(expand.grid(x=seq(.5, 1.5, len=5), y=seq(0.5, 1.5, len=5)),
               c("x", "y"))
  ref <- list(seq(0, 2, len=5), seq(0, 2, len=5))
  test <- interp.nd(data, grid, ref)
  acast(cbind(data, test), x ~ y)

data <- array(floor(runif(30)*11)/10, c(10,3))
grid <- array(floor(runif(27)*11)/10, c(3,3,3))
ref <- rep(list(seq(0, 1, length=3)), 3)

interp.nd(data, grid, ref)

x <- lennon

data(lennon)
data <- melt(expand.grid(x=seq(0.5, 1.5, len=512), y=seq(0.5, 1.5, len=512)),
             c("x", "y"))
grid <- lennon
ref <- list(seq(0, 2, len=256), seq(0, 2, len=256))
test <- interp.nd(data, grid, ref)
interpolated <- acast(cbind(data, var=test), x ~ y)

range(interpolated)
range(lennon)

}

names(grids$spacing_content)

acast(grids$spacing_content, spacing ~ content, value_var=)

ggplot(grids$spacing_content,
       aes(spacing, content, fill=)) + geom_tile()

ls()

model$model$stan_predict

##  [1] "spacing"                 "content"
##  [3] "displacement"            "eccentricity"
##  [5] "bias"                    "target_number_all"
##  [7] "target_number_shown"     "carrier.factor"
##  [9] "envelope.factor"         "check"
## [11] "content_cw"              "content_ccw"
## [13] "side"                    "content_local"
## [15] "content_global"          "full_circle"
## [17] "extent"                  "number_shown_as_spacing"
## [19] "pred"                    "fullcircle"
## [21] "norm_diff"               "energy_diff"

chain(grids$spacing_content,
      model$model$interpolator(),
      ggplot(),
      aes(spacing, content, fill=content_global)
      +geom_tile())

interpolate.by <-

    chain(interpolating,
          . %in% names(data),
          interpolating[.], #column names
          vapply(count_unmatched_values, 0)
          )

(
   .[.>0],
   sort, rev, names)
