library(grid)
library(ggplot2)
library(plyr)
library(vadr)
source("library.R")
source("stan_predictor.R")

x <- load_stanfit("SlopeModel_segment_adj2_recording.fit.RData")

slotapply <- function(thing, fn, ...) {
  slots <- slotNames(thing)
  names(slots) <- slots
  lapply(slots, function(x) fn(slot(thing, x)), ...)
}

envapply <- function(thing, fn, ...) {
  slots <- ls(thing)
  names(slots) <- slots
  lapply(slots, function(x) fn(get(x, envir=thing)), ...)
}

print(object.size(x$fits), units="Mb")

chain(x, envapply(object.size), unlist, sum, ./2^20)

chain(x$fits$fit,
      lapply(slot, name="sim"),
      lapply(function(x) x$samples),
      lapply(lapply, lapply, length),
      unlist, sum,
      . * 8 / 2^20)

#How does 700mb vectors turn into 3 GB actual.

#get the predictions extracted from stan
stan.predictions <- mply(function(optimized, fit, ...) {
  match <- data.frame(...)
  data <- match_df(x$data, match)
  predictions <- list()
  vadr::qe({
    ..(rev(
      qqply(
        .({
          temp <- parse(text=name)[[1]]
          if (is.call(temp)) {
            temp[[2]] <- call("$", as.name("predictions"), as.name(temp[[2]]))
          }
          temp
        }) <- optimized[[.(index)]])
      (name=names(optimized), index=seq_along(optimized))))
  })
  predictions <- rename(predictions,
                        c(trial_id_stan = "trial_id", fit="response"))
  cbind(data, predictions)
}) %()% x$fits
stan.predictions <- rbind.fill(stan.predictions)

#having extracted the posterior predictions,
#replace the "optimized" with shorter lists
x$fits <- chain(
  x$fits,
  splat(mply(function(optimized, fit, ...) {
    optimized <- quickdf(optimized)
    optimized <- optimized[!grepl("\\[", names(optimized))]
    quickdf(list(optimized=list(optimized), fit=list(fit), ...))
  }))(),
  rbind.fill)

#extract coefficients as well...
coefs <- chain(x$fits,
               invoke(mply(
                 function(fit, ..., optimized) quickdf(c(optimized, list(...))))),
               rbind.fill)

preds <- chain(x,
               Rpreds = predict,
               fields = c("trial_id", "response",
                          grep("link", names(.), value=TRUE)),
               rename(Rpreds[fields],
                      structure(names=., paste0(., ".r"))),
               merge(stan.predictions, ., by.x="trial_id", by.y="trial_id.r"),
               merge(coefs, by="subject"))

chain(preds,
      .[c("response", "response.r", "trial_id", "full_circle")],
      melt(id.vars=c("trial_id", "full_circle"),
           measure.vars=c("response.r", "response")),
      ggplot(aes(value, color=variable)),
      + geom_density(),
      + facet_grid( ~ full_circle),
      + theme(aspect.ratio=2))

ggplot(preds, aes(response, response.r)) +
    geom_density2d() + theme(aspect.ratio=1)

theme_set(theme_gray(9))
chain(preds, ggplot,
      + aes(response, response.r, color=factor(subject)),
      + scale_color_hue(h=c(240,0)),
      + geom_point(size=0.6),
      + guides(colour=guide_legend(override.aes=list(size=3))),
      + theme(aspect.ratio=1)
      )

#what does the difference between R and Stan predictions correlate with?
chain(preds,
      diffs=mutate(diff = link.r - link),
      numcolwise(function(x) cor(x, diffs$diff))(),
      abs, sort, t)

chain(preds, #subset(subject=="pbm"),
      sample.data.frame(1000),
      mutate(diff = link.r - link),
      colwise(function(x) (
        if (is.character(x)) {
          as.numeric(as.factor(x))
        } else if (is.factor(x)) {
          as.numeric(x)
        } else x))(),
      numcolwise(identity)()) -> plotpreds

(ggplot(plotpreds)
 + aes(response, response.r, color=factor(subject))
 + geom_point(size=1)
 + theme_grey(9)
 + theme(aspect.ratio=1))

#so why are there a few very negative links?
#looks like there is a bias thing that's wacko...

#%rem% x[c("trial_id", "diff", "subject")]

chain(melt(plotpreds, id.vars=c("trial_id")),
      merge(plotpreds,
            by="trial_id")) -> diffs

(ggplot(diffs, aes(value, diff))
 + geom_point(size=0.5, alpha=0.5)
 + facet_wrap(~variable, scales="free")
 + aes(color=factor(subject))
 + scale_color_hue(h=c(240,0))
 + scale_y_continuous(breaks=function(x) pretty(x, 3))
 + scale_x_continuous(breaks=function(x) pretty(x, 2))
 + guides(colour=guide_legend(override.aes = list(size=5)))
 + theme_grey(5)
 + theme(strip.background = element_blank(),
         strip.text.x=element_text(colour="grey50"),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         axis.ticks=element_line(size=0.5),
         aspect.ratio=1,
         axis.ticks.margin=unit(0.5, "mm"),
         axis.ticks.length=unit(0.5, "mm")))

#this is interesting -- 100% content actually
#had lower contrast on the forward freq.
chain(preds
      , ggplot
      , + aes(content_cw, content_ccw)
      , + geom_point()
      , + theme(aspect.ratio=1))

chain(preds
      , ggplot
      , + aes(content, link_summation.r, color=target_number_shown)
      , + facet_wrap(~subject, scales="free")
      , + geom_point()
      , + theme(aspect.ratio=1))

chain(preds,
      .[c("trial_id", "content", "subject",
          "link_repulsion", "link_repulsion.r")],
      ggplot,
      + aes(link_repulsion.r, link_repulsion, color=factor(content)),
      + facet_wrap(~subject),
      + geom_point(size=0.5),
      + theme(aspect.ratio=1))

chain(preds,
      .[c("link_nonlinearity", "link_nonlinearity.r",
          "full_circle", "content", "subject")],
      unique,
      melt(measure.vars=c("link_nonlinearity.r", "link_nonlinearity")),
      ggplot(aes(content, value, color=variable)),
      + geom_line(),
      + geom_point(alpha=0.5),
      + facet_wrap( ~ subject),
      + theme(aspect.ratio=1))

#101 unique values of link_repulsion for 9000 data points?
length(unique(preds$link_repulsion.r))

#474 unique values of diff for 9000 data points???
length(unique(with(preds, link.r - link)))

length(unique(plotpreds$bias))

#weird what happens with "grid..."

#So it definitely has something to do with nonlinearity...

#grump! content_cw and content_ccw are mapped weird?
ggplot(preds, aes(content_cw, content_ccw) + geom_point())

#also looks like "content" has some subgroup...
