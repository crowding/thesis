library(ggplot2)
library(scales)
library(gtable)
library(vadr)

source("library.R")

main <- function(plotfile, ...) {
  files <- list(...)
  logprobs <- getprobs(files)
  if (!interactive()) {
    cairo_pdf(plotfile, width=18, height=18)
    on.exit(dev.off())
  }
  violin_plots(logprobs)
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
  bind[, displacement, carrier] <- str_match(file, "d_(.*)_c_([^.]*)")
  chain(
    e$fits,
    invoke(mply(function(..., fit, optimized, .progress="text") {
      cbind(model_name=e$model_name,
            displacement=displacement, carrier=carrier,
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

violin_plots <- function(samples) {
  #shows the posterior distributions over each parameter for each subject
  samples <- normalize(samples)
  optimized <- subset(samples, class=="optimized")
  samples <- subset(samples, class=="sample")

  table <- gtable(widths = unit(c(1), "null"),
               heights = unit(c(1,1,1,1), "null"))
  qeply({
    the.plot <- (
        ggplot(samples)
        + aes(.(xvar), lp__, fill=.(colorvar), color=.(colorvar))
        + facet_grid(.~.(gridvar))
        + geom_violin( size=0.1, alpha=0.5, position="identity")
        + with_arg(data=optimized, geom_point(shape=4),
                   geom_line(size=0.5, linetype="11", alpha=0.5,
                             aes(group=model_name)))
        + coord_trans(ytrans=asinh_trans())
        + scale_y_continuous(breaks=trans_breaks("asinh", "sinh", 10))
        + theme_bw()
        + theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle=-30, hjust=0, size=rel(0.6))))
    table <- gtable_add_grob(table, ggplotGrob(the.plot), .(row), 1)
    NULL
  })(row=c(1,2,3,4),
     gridvar = alist(displacement, subject, carrier, subject),
     colorvar = alist(carrier, carrier, displacement, displacement),
     xvar = alist(subject, displacement, subject, carrier))
  grid.newpage()
  grid.draw(table)
}

normalize <- function(samples, group=c()) {
  samples <- ddply(samples, c("subject", group), mutate, lp__ = lp__ - max(lp__))
  samples <- mutate(samples, lp__ = lp__ / abs(quantile(lp__, 0.95)))
}

run_as_command()
