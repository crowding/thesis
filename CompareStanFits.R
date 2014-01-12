suppressPackageStartupMessages({
  library(digest)
  library(R.cache)
  library(scales)
  library(grid)
  library(rstan)
  library(vadr)
  library(plyr)
  library(reshape2)
  library(ggplot2)
  library(gtable)
  source("library.R")
  source("stan_predictor.R")
})
theme_set(theme_bw())
theme_update(panel.grid.major=element_blank(),
             panel.grid.minor=element_blank())

modelfiles <- dir(pattern=".fit.RData")
modelfiles <- c("models/d_soft_global_c_hemi_e_none.fit.RData",
                "models/d_soft_global_c_hemi_e_A.fit.RData")
outfile <- "test.pdf"

main <- function(outfile, ...) {
  filehashes <- lapply(list(...), digest, file=TRUE)
  if (file.exists(outfile)) {
    outhash <- digest(outfile, file=TRUE)
    x <- loadCache(key=filehashes)
    if (identical(x, outhash)) {
      print("match found")
      touchFile(outfile)
      return(NULL)
    }
  }

  cairo_pdf(width=10.5, height=8, outfile, onefile=TRUE)
  on.exit(dev.off, add=TRUE)
  plotStanFits(c(...))
#  dev.off()

  outhash <- digest(outfile, file=TRUE)
  saveCache(outhash, key=filehashes)
}

`chain<-` <- function(obj, ..., value) {
  f <- mkchain(...)
  environment(f) <- parent.frame()
  f(value)
  #see, with proper environment discovery, this should work...
  #mkchain(...)(value)
}

plotStanFits <- function(modelfiles=dir(pattern=".fit.RData")) {
  resids <- NULL #hmm, this would be a problem with the 'chain' idea.
  bind[fitdata=fitdata,
       chain(resids, rbind.fill)
       ] <- zip(llply(modelfiles, get_samples), list)
  bind[samples=samples, optimized=optimized] <- collect_model_data(fitdata)
  residualPlot(resids)
  violinPlot(samples, optimized)
  crossPlots(samples, optimized)
}

get_samples <- function(modelfile) {
  #load fit file, extract samples from each fit,
  #remember model name, return named list
  x <- load_stanfit(modelfile)
  samples <- mutate(x$fits,
                    fit=lapply(fit, as.data.frame),
                    model_name=x$model@model_name)
  p <- predictable(x)
  resids <- get_resids(p)
  resids %<~% mutate(model_name=x$model@model_name)
  list(fitdata=samples, resids=resids)
}

get_resids <- function(model) {
  split <- data.frame(
    split=I(list("spacing", "spacing", "target_number_shown",
                 "content", "content")),
    fc=c(TRUE, FALSE, FALSE, TRUE, FALSE),
    group=c(1, 2, 3, 4, 5),
    l=c(3, 1, 1, 2, 3),
    r=c(3, 2, 1, 2, 3),
    t=c(1, 1, 2, 2, 2),
    b=c(1, 1, 2, 2, 2),
    rng=c(10, 10, 10, 20, 10),
    legend=c("none", "right", "none", "none", "none"),
    stringsAsFactors=FALSE)
  split %<~% subset(fc %in% unique(model$data$full_circle))
  mdply(split, function(split, fc, ...) {
    chain(model,
          pearson_resids(split=unlist(split),
                         data=subset(model$data, full_circle==fc)),
          cbind(., split=I(rep_len(split, nrow(.))),
                full_circle=fc, ..., stringsAsFactors=FALSE))
  })
}

unfactor <- colwise(function(x) if (is.factor(x)) as.character(x) else x)

`%<<~%` <- macro(function(x, y) {
  if ("." %in% all.names(y)) {
    vadr::qq(.(x) <<- (function(.) .(y))(.(x)))
  }
  else if (is.call(y)) {
    vadr::qq(.(x) <<- .(y[[1]])(.(x), ..(as.list(y)[-1])))
  }
  else {
    vadr::qq(.(x) <<- .(y)(.(x)))
  }
})

residualPlot <- function(residuals) {
  #Each model gets a line for its residuals (smoothed?)and plotting a
  #number of different ways gets a line.
  groups <- unique(residuals$group)
  tab <- gtable(widths=unit(rep(1, 3), "null"),
                heights=unit(rep(1, 2), "null"))
  each_plot <- function(chunk) {
    bind[t=t,b=b,l=l,r=r,
         split=bind[split],legend=legend,rng=rng,
         ...=] <- chunk[1,]
    plot <- (
      ggplot(chunk)
      + (if (split=="spacing")
         (coord_trans(xtrans="sqrt", limy=c(-rng, rng)))
      else (coord_cartesian(ylim=c(-rng, rng))))
      + annotate("segment", y=0, yend=0, x=-Inf, xend=1000, color="gray50")
      + aes_string(x=split[[1]], y="pearson_resid", size="n_obs")
      + aes(color=model_name)
      + geom_point()
      + geom_smooth(fill=NA)
      + scale_size_area(max_size=3)
      + scale_color_brewer("Model", type="qual", palette="Paired")
      + theme(legend.position=legend)
      )
    tab %<<~% gtable_add_grob(ggplotGrob(plot), t=t, l=l, b=b, r=r)
  }
  d_ply(residuals, .(group), each_plot)
  grid.draw(tab)
}

#unpack a fitlist into two long-format data frames
collect_model_data <- function(model.list) {
  #fitlist is a list with one item per model type
  samp <- ldply(model.list, collect_fit_data, "fit")
  optim <- ldply(model.list, collect_fit_data, "optimized")
  list(samples = samp, optimized = optim)
}

collect_fit_data <- function(model, extract,
                             ignore = c("fit", "optimized") %-% extract) {
  #"model" is a list of lists:
  #"fit": data.frame posterior samples, one per fit
  #"optimized": max-likelihood fit, one per fit
  #and other, arbitrary identifying variables, one per fit
  model[ignore] <- NULL
  identifiers <- names(model) %-% ignore %-% extract
  chain(
    model,
    (Map %<<% .)(list),
    lapply(., mkchain(
      vadr::alter(., it[[extract]], .[!grepl('\\[', names(.))]),
      c(., .[[extract]]), vadr:: put(., it[[extract]], NULL),
      data.frame %()% .,
      mutate(., .n=1:nrow(.)), #remember the sample id
                          melt(id.vars=c(".n", identifiers)), #long format
                          unfactor
                          )),
        rbind.fill %()% .)
}

violinPlot <- function(samples, optimized) {
  #shows the posterior distributions over each parameter for each subject
  cat("violin plot\n")
  samples2 <- condition_warn(samples)
  if(nrow(samples2) == 0) error("seriously messed")
  samples <- samples2
  optimized <- condition_warn(optimized)
  bind[samples, optimized] <- shift_likelihoods(samples, optimized, "subject")
  if (nrow(optimized) > 1000) {
    browser()
  }
  print(
    ggplot(samples)
    + aes(subject, value, fill=model_name, color=model_name)
    + facet_wrap("variable", scales="free")
    + geom_violin( size=0.1, alpha=0.5
                  , position="identity")
    + geom_point(data=optimized, shape=4)
      )
}


condition_warn <- function(samples) {
  cond <- with(samples, (!is.finite(value)) | (abs(value) > 1E10))
  if (any(cond)) {
    message("Infinite or nan or large values!")
    message(deparse(substitute(samples)), ":", "\n")
    print(unique(samples[cond,c("variable", "subject", "model_name")]))
  }
  if(sum(!cond) == 0) browser()
  samples[!cond,]
}

#given a long format data frame, match each variable against the other
cross_variables <- mkchain[
    .
  , matchnames=names(samples) %-% c("value")
](
    samples=.
  , unique(.$variable)
  , expand.grid(xvar=., yvar=., stringsAsFactors=FALSE)
  , subset(xvar < yvar)
  , merge(samples, by.x="xvar", by.y="variable")
  , merge( samples
         , by.x = revalue(matchnames, c("variable"="yvar"))
         , by.y = matchnames)
)

#shift likelihoods to bring the maximum likelihood to zero. Delete lower tails...
shift_likelihoods <- function(samples, optimized, group="subject") {
  optimized.shift <- ddply(optimized, group, function(x) {
    vadr::alter(x, it$value[it$variable=="lp__"],
          . - max(.),
          . / abs(quantile(., 0.95)),
          asinh) #asinh transform yo
  })
  samples.shift <- ddply(samples, group, function(x) {
    lpshift__ <- chain(optimized, subset(variable=="lp__"),
                       match_df(., x,
                                on=names(.) %-% c(".n", "value")),
                       max(.$value))
    vadr::alter(x, it$value[it$variable=="lp__"], . - lpshift__)
  })
  list(samples=samples.shift, optimized=optimized.shift)
}

rescale_likelihoods <- function(samples, optimized) {
  #normalize the log probabilities by subtracting the maximum-likelihood value.
  samples %<~% vadr::alter(., it[which(it$variable == "lp__"), ],
                           mutate(., value = (
                             value - merge( optimized, .
                                           , names(samples) %-% c("value", ".n")
                                           , all.y=TRUE
                                           )$value.x
                             )))
  optimized %<~% vadr::alter(it[which(it$variable == "lp__"), ],
                             mutate(value=0))
  list(samples, optimized)
}

crossPlot <- function(samples, optimized, filter, subsample=500) {
  #shows posterior distributions as density plots over two variables;
  #input is long format data frames
  if (!missing(filter)) {
    samples <- match_df(samples, filter, on=names(filter))
    optimized <- match_df(optimized, filter, on=names(filter))
  }

  samples <- condition_warn(samples)
  optimized <- condition_warn(optimized)

  matchnames <- names(samples) %-% c("value")

  #subset samples if too many
  if (subsample < Inf) {
    samples <- ddply(samples, matchnames %-% ".n", function(x) {
      if (max(x$.n) > subsample) {
        subset(x, .n %in% sample(max(x$.n), subsample))
      } else x
    })
  }

  bind[samples, optimized] <- rescale_likelihoods(samples, optimized)

  #match samples with every pair of parameters
  cross.samples <- cross_variables(samples)
  cross.optimized <- cross_variables(optimized)

  interaction.aes.expr <- template(
    interaction( ...( lapply( matchnames %-% c(".n", "variable")
                            , as.name) ) )
  )

  print(
    ggplot(cross.samples)
   + aes(x=value.x, y=value.y)
   + geom_hline(y=0, size=0.2, color="gray50")
   + geom_vline(x=0, size=0.2, color="gray50")
   + eval(template(
        aes(color=.(interaction.aes.expr), fill=.(interaction.aes.expr))
     ))
   + scale_color_brewer("Group", type="qual", palette="Paired")
   + scale_fill_brewer("Group", type="qual", palette="Paired")
   + facet_grid(yvar ~ xvar, scales="free")
   + stat_density2d(bins=2, geom="polygon", alpha=0.5, size=0.2)
   + geom_point(shape=3, data=cross.optimized, size=2)
   + theme_grey()
   + theme(
       strip.text.x=element_text(size=rel(0.5)),
       strip.text.y=element_text(size=rel(0.5)),
     strip.background=element_blank(),
     axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0,
                              size=rel(0.6), face=2),
     axis.text.y=element_text(angle=0, hjust=1, size=rel(0.5), face=2),
     panel.grid.major=element_line(size=0.2),
     panel.grid.minor=element_blank(),
     panel.margin=unit(1, "mm"),
           aspect.ratio=1)
   + scale_y_continuous(breaks = curr(pretty, n=3))
   + scale_x_continuous(breaks = curr(pretty, n=3))
   + labs(x="", y="", title="FWHM outlines of posterior parameter distributions")
   )
}

crossPlots <- function(samples, optimized) {
  l_ply(unique(samples$model_name), function(x) {
    cat("cross plot ", x, "\n")
    crossPlot(samples, optimized,
              filter=data.frame(model_name=x, stringsAsFactors=FALSE))
  })
}

try_command <- mkchain(strsplit(" +"), unlist, .[-1:-2], main %()% .)

run_as_command()
