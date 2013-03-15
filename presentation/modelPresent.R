
## @knitr preliminaries
suppressPackageStartupMessages({
library(plyr)
library(ggplot2)
library(ptools)
library(psyphy)
library(scales)
library(grid)
})
theme_set(theme_bw())
#in this file I define pretty scales.
use_unicode = FALSE
source("../scales.R")

match_df <- function(...) suppressMessages(plyr::match_df(...))

## @knitr loadData
load("../data.RData")

do.rename <- function(data, folding=TRUE) {
  replacements <- if (folding) {
    c(folded_direction_content="content",
      folded_displacement="displacement",
      folded_response_with_carrier="response",
      target_spacing="spacing"
      )
  } else {
    c(abs_direction_content="content",
      abs_displacement="displacement",
      abs_response_cw="response",
      target_spacing="spacing"
      )
  }
  rename(data, replacements)
}

## @knitr staircase-example
staircase_example <-
  suppressMessages(match_df(data,
           data.frame(subject="pbm", exp_type="spacing",
                      folded_direction_content=0.15, target_number_all=12)))
staircase_example <- do.rename(staircase_example, folding=TRUE)

(ggplot(staircase_example, aes(trial.order, displacement,
                               shape=response,
                               color=response,
                               fill=response))
 + geom_point()
 + labs(title="PBM, $C=0.15$, $\\mathrm{spacing}=3.5^\\circ$",
        x="Trial no.", y="$\\Delta x$")
 + with_arg(name="Response", breaks=c(FALSE, TRUE), labels=c("CW", "CCW"),
            scale_shape_manual(values=c(24, 25)), scale_color_discrete(),
            scale_fill_discrete())
 + theme(legend.position="bottom")
 )


## @knitr example-content
example_content <- match_df(data, data.frame(exp_type="content", subject="pbm"))
example_content <- do.rename(example_content, folding=FALSE)

rates <- mkchain(ddply(c("displacement", "content", "spacing"),
                       summarize, n = length(response), p = mean(response)),
                 arrange(desc(n)))

example_content_fits <-
    glm(data=example_content,
        response ~ displacement*spacing+content+content:spacing,
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)))

seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)

## @knitr plotting-preliminaries


## @knitr example-content-plot
(ggplot(rates(example_content))
 + displacement_scale
 + proportion_scale
 + balloon
 + content_color_scale + facet_spacing_rows
 + labs(title="Effect of direction content at narrow and wide spacing")
 + add_predictions(example_content, example_content_fits)
 )

## @knitr example-spacing
chain(data,
      match_df(data.frame(exp_type="spacing")),
      do.rename(folding=FALSE),
      count(c("subject", "content")),
      arrange(desc(freq)),
      head(10)
      ) -> examples_spacing

example_spacing <-
  chain(data,
        do.rename(folding=TRUE),
        match_df(data.frame(subject="pbm", exp_type="spacing", content=0.15)))

example_spacing_fits <-
    glm(data=example_spacing,
        response ~ displacement + I(content/spacing),
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)))

#spacing_shape_scale fuck it. there
#was a nice sequence but I can't
#make unicode work. It's LyX's
#problem, or the yihui's LyX script.
shape_sequence <- c(16, 17, 15, 18, 9, 8, 13, 1)

#shape_sequence <- -rev(c(-16, -17, 0x2726L, 0x2605L, 0x2736L, 0x2737L, 0x2739L, 0x273AL))
shape_chooser <- function(n) {
  shape_sequence[seq(floor(length(shape_sequence)-n)/2+1, length=n)]
}
spacing_shape_scale <-
  list(
    aes(shape=factor(spacing)),
    discrete_scale("shape", "manual", name="Spacing", palette=shape_chooser, labels=prettyprint)
    )
editGtable <- function(gt, idx, ...) {
  if (is.character(idx))
    idx <- which(gt$layout$name == idx)
  gt$grobs[[idx]] <- editGrob(gt$grobs[[idx]], ...)
  gt
}

## @knitr example-spacing-plot
gg <- ggplotGrob(ggplot(rates(example_spacing))
 + displacement_scale
 + proportion_scale
 + spacing_texture_scale
 + spacing_shape_scale
 + add_predictions(example_spacing, example_spacing_fits)
 + balloon
# + facet_spacing_rows
 )
#gg <- editGtable(gg, "panel", "geom_point.points", grep=TRUE, global=TRUE,
#                 gp=gpar(fontfamily="MS Gothic"))
grid.draw(gg)

