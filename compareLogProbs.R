library(data.table)
library(ggplot2)
library(scales)
library(gtable)
library(vadr)

source("library.R")

test <- function() {
  files <-
    c("models/d_soft_local_c_global_e_none.fit.RData",
      "models/d_soft_global_c_global_e_none.fit.RData",
      "models/d_soft_hemi_c_global_e_none.fit.RData",
      "models/d_soft_local_c_hemi_e_none.fit.RData",
      "models/d_soft_global_c_hemi_e_none.fit.RData",
      "models/d_soft_hemi_c_hemi_e_none.fit.RData",
      "models/d_soft_windowed_c_hemi_e_none.fit.RData",
      "models/d_soft_local_c_local_e_none.fit.RData",
      "models/d_soft_global_c_local_e_none.fit.RData",
      "models/d_soft_hemi_c_local_e_none.fit.RData",
      "models/d_soft_local_c_windowed_e_none.fit.RData",
      "models/d_soft_hemi_c_windowed_e_none.fit.RData",
      "models/d_soft_local_c_global_e_R.fit.RData",
      "models/d_soft_global_c_global_e_R.fit.RData",
      "models/d_soft_hemi_c_global_e_R.fit.RData",
      "models/d_soft_local_c_hemi_e_R.fit.RData",
      "models/d_soft_global_c_hemi_e_R.fit.RData",
      "models/d_soft_hemi_c_hemi_e_R.fit.RData",
      "models/d_soft_windowed_c_hemi_e_R.fit.RData",
      "models/d_soft_local_c_local_e_R.fit.RData",
      "models/d_soft_global_c_local_e_R.fit.RData",
      "models/d_soft_hemi_c_local_e_R.fit.RData",
      "models/d_soft_local_c_windowed_e_R.fit.RData",
      "models/d_soft_hemi_c_windowed_e_R.fit.RData",
      "models/d_soft_local_c_hemi_e_RE.fit.RData",
      "models/d_soft_local_c_global_e_A.fit.RData",
      "models/d_soft_hemi_c_global_e_A.fit.RData",
      "models/d_soft_local_c_hemi_e_A.fit.RData",
      "models/d_soft_global_c_hemi_e_A.fit.RData",
      "models/d_soft_hemi_c_hemi_e_A.fit.RData",
      "models/d_soft_windowed_c_hemi_e_A.fit.RData",
      "models/d_soft_local_c_local_e_A.fit.RData",
      "models/d_soft_global_c_local_e_A.fit.RData",
      "models/d_soft_hemi_c_local_e_A.fit.RData",
      "models/d_soft_local_c_windowed_e_A.fit.RData",
      "models/d_soft_hemi_c_windowed_e_A.fit.RData",
      "models/d_soft_local_c_hemi_e_AE.fit.RData",
      "models/d_soft_local_c_global_e_RA.fit.RData",
      "models/d_soft_global_c_global_e_RA.fit.RData",
      "models/d_soft_hemi_c_global_e_RA.fit.RData",
      "models/d_soft_global_c_hemi_e_RA.fit.RData",
      "models/d_soft_hemi_c_hemi_e_RA.fit.RData",
      "models/d_soft_local_c_local_e_RA.fit.RData",
      "models/d_soft_global_c_local_e_RA.fit.RData",
      "models/d_soft_hemi_c_local_e_RA.fit.RData",
      "models/d_soft_local_c_windowed_e_RA.fit.RData",
      "models/d_soft_hemi_c_windowed_e_RA.fit.RData")

  plotfile <- "models/compare_lp.pdf"

  main %<<% dots(plotfile=plotfile, pdfout=TRUE) %()% files
}

main <- function(plotfile, ..., pdfout=!interactive()) {
  files <- list(...)
  if (pdfout) {
    cat("Outputting to ", plotfile, "\n")
    cairo_pdf(plotfile, width=18, height=24, onefile=TRUE)
  }
  if (!exists("logprobs")) {
    logprobs <- getprobs(files)
  }
  collapsed_plots(logprobs)
  violin_plots(logprobs)
  if (pdfout) {
    dev.off()
  }
}


test <- function() {
  files <- file.path("models", dir("models", "fit.RData"))
  main %()% c("likelihoods.pdf", files)
}

getprobs <- function(files) {
  ldply(files, getprob)
}

split_array <- function(a, dim) alply(a, dim, identity)

getprob <- function(file) {
  cat(file, "\n")
  e <- load2env(file)
  bind[, displacement, carrier, endpoint] <-
      str_match(file, "d_(.*)_c_([^.]*)_e_([^.]*)")
  chain(
    e$fits,
    invoke(mply(function(..., fit, optimized, .progress="text") {
      cbind(model_name=e$model_name,
            displacement=displacement,
            carrier=carrier,
            endpoint=endpoint,
            rbind.fill(cbind(as.data.frame(fit)["lp__"],
                             class="sample", ...),
                       cbind(as.data.frame(optimized)["lp__"],
                             class="optimized", ...)))
    })),
    invoke(rbind.fill))
}

invoke <- function(x, f) f %()% x

asinh_trans <- function() {
  trans_new("asinh", "asinh", "sinh")
}

`%<~%` <- macro(function(dest, value) {
  qq( .(dest) <- chain(.(dest), .(value)) )
})

label_varvalue <- function(variable, value) {
  paste(variable, "=", value)
}

violin_plots <- function(samples, fix=NULL) {
  #shows the posterior distributions over each parameter for each subject
  plots <- chain(
    cor = c("displacement", "carrier", "subject", "endpoint"),
    expand.grid(gridvar = ., colorvar=., xvar=., fixvar=.,
                stringsAsFactors=FALSE),
    .[apply(., 1, function(x) length(setdiff(cor, x)) == 0), ])

  if (!is.null(fix)) {
    plots <- subset(plots, fixvar %in% names(fix))
  }

  dlply_along(plots, .(fixvar), function(vargroup, varchunk) {
    varchunk$row <- seq_len(nrow(varchunk))
    #generate a page for each value of vargroup

    ddply_along(samples, vargroup$fixvar[[1]], function(fixgroup, samps) {
      fixgroup <- unfactor(fixgroup)
      if (!is.null(fix)
          && !((names(fixgroup)[[1]] %in% names(fix))
               && (fixgroup[[1]] %in% fix[[names(fixgroup)[[1]]]]))) {
        return()
      }
      browser()
      #a row for each in vargroup

      label <- str_match(deparse(fixgroup, control=c()), '\\((.*)\\)')[,2]
      message(label)
      samps <- normalize(samps)
      optimized <- subset(samps, class=="optimized")
      samps <- subset(samps, class=="sample")
      rng <- range(samples$lp__)
      orng <- range(optimized$lp__)

      labelgrob <- textGrob(label=label, gp=gpar(fontsize=18))
      table <- gtable(widths = unit(c(1), "null"),
                      heights = unit.c(unit(1.5, "grobheight", labelgrob),
                                       unit(rep(1, nrow(varchunk)), "null")))
      do.call(qeply({
        the.plot <- (
          ggplot(samps)
          + aes(`.(xvar)`, lp__, fill=`.(colorvar)`, color=`.(colorvar)`)
          + scale_fill_discrete(guide = guide_legend(ncol=2))
          + facet_grid(.~`.(gridvar)`, labeller=label_varvalue)
          + geom_violin(data=samps, size=0.1, alpha=0.5, position="identity")
          + geom_point(data=optimized, shape=4)
          + geom_line(data=optimized, size=0.5, linetype="11", alpha=0.5,
                      aes(group=`.(colorvar)`))
          + coord_trans(ytrans=asinh_trans(), limy = orng * 1.1 + c(-1, 0))
          + scale_y_continuous(breaks=trans_breaks("asinh", "sinh", 5))
          + theme_bw(9)
          + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle=-30, hjust=0, size=rel(0.6))))
        table <- gtable_add_grob(table, ggplotGrob(the.plot), .(row) + 1, 1)
        NULL
      }), varchunk)
      table <- gtable_add_grob(table, labelgrob, 1, 1)
      grid.newpage()
      grid.draw(table)
    })
  })
}

collapsed_plots <- function(samples) {
  #then, again summing over observers

  chain(
    samples,
    data.table,
    .[, lp__ := lp__ - max(lp__), by=subject],
    .[, sample_id := seq_along(lp__),
      by=list(displacement, carrier, endpoint, class, subject)],
    .[, list(lp__ = sum(lp__), over = length(lp__)),
      by=list(displacement, carrier, endpoint, class, sample_id)],
    .[, lp__ := lp__ - max(lp__)],
    .[, list(lp__=max(lp__)),
      by=list(displacement, carrier, endpoint, class)]
    ) -> samps

  plots <- chain(
    cor = c("displacement", "carrier", "endpoint"),
    expand.grid(gridvar = ., colorvar=., xvar=.,
                stringsAsFactors=FALSE),
    .[apply(., 1, function(x) length(setdiff(cor, x)) == 0), ],
    mutate(., row = seq_len(nrow(.))))

  optimized <- subset(samps, class=="optimized")
  samps <- subset(samps, class=="sample")
  rng <- range(samples$lp__)
  orng <- range(optimized$lp__)
  table <- gtable(widths = unit(c(1), "null"),
                  heights = unit(rep(1, nrow(plots)), "null"))

  do.call(qeply({
    the.plot <- (
      ggplot(samps)
      + aes(`.(xvar)`, lp__, fill=`.(colorvar)`, color=`.(colorvar)`)
      + scale_fill_discrete(guide = guide_legend(ncol=2))
      + facet_grid(.~`.(gridvar)`, labeller=label_varvalue)
      + geom_violin(data=samps, size=0.1, alpha=0.5, position="identity")
      + geom_point(data=optimized, shape=4)
      + geom_line(data=optimized, size=0.5, linetype="11", alpha=0.5,
                  aes(group=`.(colorvar)`))
      + coord_trans(ytrans=asinh_trans(), limy = orng * 2 + c(-1, 0))
      + scale_y_continuous(breaks=trans_breaks("asinh", "sinh", 5))
      + theme_bw(9)
      + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle=-30, hjust=0, size=rel(0.6))))
    table <- gtable_add_grob(table, ggplotGrob(the.plot), .(row), 1)
    NULL
  }), plots)

  grid.newpage()
  grid.draw(table)
}

normalize <- function(samples, group=c()) {
  samples <- ddply(samples, c("subject", group), mutate, lp__ = lp__ - max(lp__))
  samples <- mutate(samples, lp__ = lp__ / abs(quantile(lp__, 0.95)))
}

run_as_command()
