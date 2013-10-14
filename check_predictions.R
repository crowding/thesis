library(grid)
library(ggplot2)
library(plyr)
library(vadr)
source("library.R")
source("stan_predictor.R")

fitfile <- "models/d_soft_local_c_repulsive_endpoints.fit.RData"

main <- function(outfile="checks/models/d_soft_local_c_repulsive_endpoints.pdf",
                 fitfile="models/d_soft_local_c_repulsive_endpoints.fit.RData",
                 plot=FALSE) {
  x <- load_stanfit(fitfile)
  if (!interactive()) {
    cairo_pdf(outfile, width=9, height=9, onefile=TRUE)
    on.exit(dev.off())
  }
  if ("trial_id_stan[1]" %in% names(x$fits$optimized[[1]])) {
    preds <- get_preds(x)
    plot_simple(preds)
    if (!check_preds(preds)) {
      if (plot) {
        plot_vars(preds)
      }
      message("Mismatched information from Stan and R!")
    }
  } else {
    message("No link information found!")
  }
}

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

#How does 700mb vectors turn into 3 GB actual.

check_preds <- function(preds) {
  with(preds, cor(link, link.r)) > 0.99
}

pause <- function() {
  if (interactive()) {scan()}
}

`%^%` <- intersect

get_preds <- function(x) {
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
            } else {
              temp <- call("$", as.name("predictions"), temp)
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

  Rpreds <- chain(x,
                 Rpreds = predict,
                 fields = c("trial_id", "response",
                            grep("link", names(.), value=TRUE)) %^% names(.),
                 rename(Rpreds[fields],
                        structure(names=., paste0(., ".r"))),
                 merge(stan.predictions, ., by.x="trial_id", by.y="trial_id.r"),
                 merge(coefs, by="subject"))

  Rpreds
}

plot_simple <- function(preds) {

  if (FALSE) {
    chain(preds,
          .[c("response", "response.r", "trial_id", "full_circle")],
          melt(id.vars=c("trial_id", "full_circle"),
               measure.vars=c("response.r", "response")),
          ggplot(aes(value, color=variable)),
          + geom_density(),
          + facet_grid( ~ full_circle),
          + theme(aspect.ratio=2))
    pause()

    print(ggplot(preds, aes(response, response.r)) +
          geom_density2d() + theme(aspect.ratio=1))
    pause()
  }

  theme_set(theme_gray(9))
  chain(preds, ggplot,
        + aes(response, response.r, color=factor(subject)),
        + scale_color_hue(h=c(240,0)),
        + geom_point(size=0.6),
        + guides(colour=guide_legend(override.aes=list(size=3))),
        + theme(aspect.ratio=1),
        print
        )
  pause()

  #what does the difference between R and Stan predictions correlate with?
  chain(preds,
        diffs=mutate(diff = link.r - link),
        numcolwise(function(x) cor(x, diffs$diff))(),
        abs, sort, t)
}

plot_vars <- function(preds) {

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

  chain(melt(plotpreds, id.vars=c("trial_id")),
        merge(plotpreds,
              by="trial_id")) -> diffs

  par(bg="black")
  print(
    ggplot(diffs, aes(value, diff))
    + geom_point(size=0.5, alpha=0.5)
    + facet_wrap(~variable, scales="free")
    + aes(color=factor(subject))
    + scale_color_hue(h=c(270, 630), direction=1, l=65)
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
            axis.ticks.length=unit(0.5, "mm"),
            plot.background=element_rect(fill="black"),
            panel.background=element_rect(fill="gray20"),
            legend.background=element_rect(fill="black"),
            legend.key=element_rect(fill="gray20", colour="black"),
            legend.title=element_text(color="gray50"),
            legend.text=element_text(color="gray50")))
  pause()

  #this is interesting -- 100% content actually
  #had lower contrast on the forward freq.

  if(FALSE) {
    chain(preds,
          .[c("link_nonlinearity", "link_nonlinearity.r",
              "full_circle", "content", "subject")],
          unique,
          melt(measure.vars=c("link_nonlinearity.r", "link_nonlinearity")),
          ggplot(aes(content, value, color=variable)),
          + geom_line(),
          + geom_point(alpha=0.5),
          + facet_wrap( ~ subject),
          + theme(aspect.ratio=1),
          print)
  }

}

run_as_command()
