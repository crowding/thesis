
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


## @knitr staircase-example
staircase_example <-
  suppressMessages(match_df(data,
           data.frame(subject="pbm", exp_type="spacing",
                      folded_direction_content=0.15, target_number_all=12)))

(ggplot(staircase_example, aes(trial.order, folded_displacement,
                               shape=folded_response_with_carrier,
                               color=folded_response_with_carrier,
                               fill=folded_response_with_carrier))
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

example_content_rates <- ddply(example_content,
                .(abs_displacement, abs_direction_content, target_spacing),
                summarize, n = length(abs_response_cw), p = mean(abs_response_cw))

example_content_fits <-
    glm(data=example_content,
        abs_response_cw ~ abs_displacement*target_spacing+abs_direction_content,
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)))

seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)

## @knitr plotting-preliminaries
add_arrows <- function(x) {
  x[1] <- paste(x[1], "$\\curvearrowleft$")
  x[length(x)] <- paste(x[length(x)], "$\\curvearrowright$")
  x 
}

replace_arrows <- function(x) {
  x[1] <- "$\\curvearrowleft$"
  x[length(x)] <- "$\\curvearrowright$"
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
  facet_grid(target_spacing ~ .,
               labeller=function(v, value) {
                 paste("$S = ", format(value, digits=3), "^\\circ$", sep="")
               }
             )

balloon <- list(geom_point(aes(size=n))
                , scale_size_area()
                )


displacement_scale <-
  list( aes(x=abs_displacement),
        scale_x_continuous("$\\Delta x$", labels=add_arrows))

#have a custom stats function that shows the predictions with error bars.
add_predictions <- function(data, model) {
  chain(data,
        subset(select=c("abs_direction_content", "target_spacing")),
        unique,
        merge(data.frame(abs_displacement=seq_range(
                           range(example_content$abs_displacement),
                           length=100)),
              all.x=TRUE, all.y=TRUE),
        cbind(.,
              predict(model, newdata=.,
                      type="response", se.fit=TRUE)[1:2] )
        ) -> predictions
  with_arg(data=predictions,
           geom_ribbon(color="transparent", alpha=0.2,
                       aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit,
                           fill=factor(abs_direction_content),
                           )),
           geom_line(aes(y=fit)))
}


#I want a discrete color scape derived from a 3-point gradient, so I wrote:
discretize <- function(pal) {
  function(n) {
    pal(seq(0,1,length=n))
  }
}

color_pal <- #function(n) rainbow(n, start=0, end=0.75)
  discretize(gradient_n_pal(
               colours=muted(c("blue", "cyan", "yellow", "red"), l=70, c=180)))

content_color_scale <-
  c(
    list(aes(color=factor(abs_direction_content),
             fill=factor(abs_direction_content))),
    with_arg(name="$C$",
             palette=color_pal,
             labels=add_arrows,
             discrete_scale("fill", "manual"),
             discrete_scale("colour", "manual")
             )
    )

## @knitr example-content-plot
(ggplot(example_content_rates)
 + displacement_scale
 + proportion_scale
 + balloon
 + content_color_scale + facet_spacing_rows
 + labs(title="Effect of direction content at narrow and wide spacing")
 + add_predictions(example_content_rates, example_content_fits)
 )

## @knitr example-spacing
chain(data,
      match_df(data.frame(exp_type="content", exp_type="spacing")),
      count(c("subject", "folded_direction_content")),
      arrange(desc(freq)),
      head(10)
      ) -> examples_spacing

example_spacing <- match_df(data, examples[6,])

#ugly to use here, but
folding <- TRUE

example_spacing_rates <- ddply(example_spacing,
                .(folded_displacement, folded_direction_content, target_spacing),
                summarize, n = length(folded_response_with_carrier), p = mean(folded_response_with_carrier))

example_spacing_fits <-
    glm(data=example_spacing,
        folded_response_with_carrier ~ folded_displacement*target_spacing+folded_direction_content,
        family=binomial(link=logit.2asym(g=0.025, lam=0.025)))

## @knitr example-spacing-plot
