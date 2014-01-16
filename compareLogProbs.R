library(data.table)
library(ggplot2)
library(scales)
library(gtable)
library(vadr)

source("library.R")

test <- mkchain[
  .="NineModels.list",
  plotfile="models/compare_lp.pdf"](
    readLines,
    grep(pattern="models/", value=TRUE),
    gsub(".stan.", ".fit.", .),
    main %<<<% dots(plotfile=plotfile, pdfout=TRUE) %()% .
    )

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

getprobs <- function(files) {
  ldply(files, getprob, .parallel=TRUE)
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
      label <- str_match(deparse(fixgroup, control=c()), '\\((.*)\\)')[,2]
      message(label)
      if (!is.null(fix)
          && !((names(fixgroup)[[1]] %in% names(fix))
               && (fixgroup[[1]] %in% fix[[names(fixgroup)[[1]]]]))) {
        return()
      }

      quants <- make_boxes(samps)
      rng <- range(c(quants$q.0, quants$q.100)) + c(-5, 5)
      #a row for each in vargroup

      labelgrob <- textGrob(label=label, gp=gpar(fontsize=18))
      table <- gtable(widths = unit(c(1), "null"),
                      heights = unit.c(unit(1.5, "grobheight", labelgrob),
                                       unit(rep(1, nrow(varchunk)), "null")))
      do.call(qeply({
        the.plot <- (
          ggplot(quants)
          + aes(`.(xvar)`, lp__, fill=`.(colorvar)`, color=`.(colorvar)`,
                middle=q.50, ymin=q.0, ymax=q.100,
                lower=q.25, upper=q.75)
          + scale_fill_discrete(guide = guide_legend(ncol=2))
          + facet_grid(.~`.(gridvar)`)
#          + (stat="identity", size=0.1, alpha=0.2)
          + geom_line(size=0.25, aes(group=`.(colorvar)`))
          + geom_point(size=0.5)
          + coord_trans(ytrans=asinh_trans(), limy = rng)
          + scale_y_continuous(breaks=trans_breaks("asinh", "sinh", 5))
          + theme_bw(9)
          + labs(title=".(gridvar)")
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

#hate "#" in names...
rename_quantile <- sub %<<<% dots(".*?([0-9]+).*", "q.\\1")

make_boxes <- function(samples, collapse=character(0)) {

  samps <- data.table(samples)

  keys <- c("displacement", "carrier", "endpoint", "subject")
  keysi <- c(keys, "sample.id")
  keys_s <- setdiff(keys, collapse)

  optim <- chain(samps,
                 .[class=="optimized"],
                 .[, list(lp__=sum(lp__)), by=keys_s],
                 .[, lp__ := lp__ - max(lp__), by="subject"],
                 .[, lp__ := lp__-max(lp__)-50],
                 )

  samps <- chain(samps,
                 .[class=="sample"],
                 data.table,
                 .[, sample.id := seq_along(lp__), by=keys],
                 .[, list(lp__=sum(lp__)), by=c(keys_s, "sample.id")],
                 merge(optim, by=keys_s),
                 .[, lp__ := lp__.x - max(lp__.x) + lp__.y,
                   by=keys_s],
                 )

  quants <- chain(samps,
                  .[, c(lp__=max(lp__), as.list(quantile(lp__))), by=keys_s],
                  setnames(., names(.), rename_quantile(names(.)))
                  )
}

collapsed_plots <- function(samples) {
  #then, again summing over observers

  quants <- make_boxes(samples, collapse="subjects")

  plots <- chain(
    cor = c("displacement", "carrier", "endpoint"),
    expand.grid(gridvar = ., colorvar=., xvar=.,
                stringsAsFactors=FALSE),
    .[apply(., 1, function(x) length(setdiff(cor, x)) == 0), ],
    mutate(., row = seq_len(nrow(.))))

  rng <- range(c(quants$q.0, quants$q.100)) + c(-5, 5)

  table <- gtable(widths = unit(c(1), "null"),
                  heights = unit(rep(1, nrow(plots)), "null"))
  do.call(qeply({
    the.plot <- (
      ggplot(quants)
      + aes(`.(xvar)`, lp__, fill=`.(colorvar)`, color=`.(colorvar)`,
            middle=q.50, ymin=q.0, ymax=q.100,
            lower=q.25, upper=q.75)
      + scale_fill_discrete(guide = guide_legend(ncol=2))
      + facet_grid(.~`.(gridvar)`)
      ## + geom_boxplot(stat="identity", size=0.1, alpha=0.2,
      ##                position="dodge")
      + geom_line(size=0.25,
                  aes(group=`.(colorvar)`))
      + geom_point(size=0.5)
      + coord_trans(ytrans=asinh_trans(),
                    limy = rng
                    )
      + scale_y_continuous(breaks=trans_breaks("asinh", "sinh", 5))
      + theme_bw(9)
      + labs(title=".(gridvar)")
      + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle=-30, hjust=0, size=rel(0.6))))
    table <- gtable_add_grob(table, ggplotGrob(the.plot), .(row), 1)
    NULL
  }), plots)
  #grid.newpage()
  grid.draw(table)

}

run_as_command()
