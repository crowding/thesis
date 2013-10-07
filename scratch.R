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

microbenchmark(
    a=mapply(list, a=1:100, b=letters[1:100], MoreArgs=list(c="ha")),
    b=mply(list, c="ha")(a=1:100, b=letters[1:100]))


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

microbenchmark(
    a=do.call(list, as.list(letters)),
    b=qe(list(..(letters))),
    c=list %()% letters)

env.list <- macro(function(...) {
  names <- list(...)
  qq(function(`.(names)`=..(missing_value(length(names)))) environment())
})

env.list.2 <- macro(function(...) {
  names <- list(...)
  qe(function(`.(names)`=..(missing_value(length(names)))) environment())
})

fab <- function(a, b, c, d, e) environment()
fab2 <- function() function(a, b, c, d, e) environment()

#b and f are misleadingly slow when enableJIT(3) -- microbenchmark artifact
microbenchmark(
    a = env.list("a", "b", "c", "d", "e")(1, 2, 3, 4, 5),
    a1 = env.list.2("a", "b", "c", "d", "e")(1, 2, 3, 4, 5),
    b = (function(a,b,c,d,e) environment())(1,2,3,4,5),
    c = fab(1,2,3,4,5),
    c1 = fab2()(1,2,3,4,5),
    d = list2env(list(a=1, b=2, c=3, d=4, b=5)),
    e = list2env(do.call(list, list(a=1,b=2,c=3,d=4,e=5))),
    f = do.call(function(a,b,c,d,e) environment(), list(a=1,b=2,c=3,d=4,e=5)),
    g = do.call(fab, list(a=1,b=2,c=3,d=4,e=5)))


#diagnose something wrong with my predictions, by taking predictions from the
#computed values of Stan.

my.subset <- function(...) {
  result <- base::`[`(...)

  # 
  result
}


chain(d$data, (e$filter_data)(), (e$format_data)(menergy))
