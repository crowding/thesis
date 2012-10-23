
## @knitr preliminaries
suppressPackageStartupMessages({
library(plyr)
library(ggplot2)
library(ptools)
library(psyphy)
library(scales)
})
theme_set(theme_bw())
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
add_arrows <- function(x) {
  first <- which(!is.na(x))[1]
  last <- length(x) + 1 - which(!is.na(rev(x)))[1]
  if (substr(x[[first]], 1, 1) == "-") {
    x[first] <- paste(x[first], "$\\curvearrowleft$")
    x[last] <- paste(x[last], "$\\curvearrowright$")
  } else {
    x[first] <- paste(x[first], "$\\curvearrowright$")
    x[last] <- paste(x[last], "$\\circlearrowright$")
  }
  x
}

replace_arrows <- function(x) {
  first <- which(!is.na(x))[1]
  last <- length(x) + 1 - which(!is.na(rev(x)))[1]
  x[first] <- "$\\curvearrowleft$"
  x[last] <- "$\\curvearrowright$"
  x
}

proportion_scale <-
  list(aes(y=p),
       scale_y_continuous(
         "Response", breaks=c(0, 0.5, 1),
         labels=replace_arrows
         ),
       theme(axis.text.y=element_text(angle=90)))

facet_spacing_rows <-
  facet_grid(spacing ~ .,
               labeller=function(v, value) {
                 paste("$S = ", format(value, digits=3), "^\\circ$", sep="")
               }
             )

balloon <- list(geom_point(aes(size=n))
                , scale_size_area()
                )


displacement_scale <-
  list( aes(x=displacement),
        scale_x_continuous("$\\Delta x$", labels=add_arrows))

#have a custom stats function that shows the predictions with error bars.
add_predictions <- function(data, model) {
  chain(data,
        subset(select=c("content", "spacing")),
        unique,
        merge(data.frame(displacement=seq_range(
                           range(data$displacement),
                           length=100)),
              all.x=TRUE, all.y=TRUE),
        cbind(.,
              predict(model, newdata=.,
                      type="response", se.fit=TRUE)[1:2] )
        ) -> predictions
  with_arg(data=predictions,
           geom_ribbon(color="transparent", alpha=0.2,
                       aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
           geom_line(aes(y=fit)))
}

comp <- function(a, b) function(...) b(a(...))
prettyprint <- function(x) format(as.numeric(x), digits=3)

#I want a discrete color scape derived from a 3-point gradient, so I wrote:
discretize <- function(pal) {
  function(n) {
    pal(seq(0,1,length=n))
  }
}

color_pal <- #function(n) rainbow(n, start=0, end=0.75)
  discretize(gradient_n_pal(
               colours=muted(c("blue", "cyan", "yellow", "red"), l=70, c=180)))

spacing_color_scale <-
  c(
    list(aes(color=factor(spacing),
             fill=factor(spacing))),
    with_arg(name="Spacing",
             palette=color_pal,
             labels=prettyprint,
             discrete_scale("fill", "manual"),
             discrete_scale("colour", "manual")
             )
    )

spacing_texture_scale <-
  list(aes(linetype=factor(spacing)),
       discrete_scale(name="Spacing", "linetype", "linetype_m", function(x) {
         paste(y <- c(1:9, LETTERS)[seq_len(x)], y, sep="")
       }, labels=prettyprint)
       )

content_color_scale <-
  c(
    list(aes(color=factor(content),
             fill=factor(content))),
    with_arg(name="$C$",
             palette=color_pal,
             labels=add_arrows,
             discrete_scale("fill", "manual"),
             discrete_scale("colour", "manual")
             )
    )

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
        do.rename(folding=FALSE),
        match_df(data.frame(subject="pbm", exp_type="spacing")))

example_spacing_fits <-
    glm(data=example_spacing,
        response ~ displacement + I(content/spacing),
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)))

# @knitr example-spacing-plot
 (ggplot(rates(example_spacing))
 + displacement_scale
 + proportion_scale
 + content_color_scale
 + spacing_texture_scale
 + add_predictions(example_spacing, example_spacing_fits)
 + balloon
 + facet_spacing_rows
 )


##Now let's plot all the data from one subject. A nice way to do thsi
##is with spacing in rows.

#we can also do .... plots where we interchange residuals for spacing?
#for spacing plots...

#Ah! Plot residuals on top of each psychometric function!!! That'd be
#a nice summary of the data.

#there should be a way to convert the residual "yes" values into an
#errorbar value.

# I think one focus should be to get
# a model fit on every subject. Then
# to look at number/density.

