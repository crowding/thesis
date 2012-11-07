library(scales)

if (use_unicode) {
  curveleft = "\u21BA"
  curveright = "\u21BB"
  circleleft = "\u21BA" #"0x27F2L"
  circleright = "\u21BB" #"0x27F3L"
} else {
  curveleft = "$\\curvearrowleft$"
  curveright = "$\\curvearrowright$"
  circleleft = "$\\curvearrowleft$"
  circleright = "$\\curvearrowright$"
}

add_arrows <- function(x) {
  first <- which(!is.na(x))[1]
  last <- length(x) + 1 - which(!is.na(rev(x)))[1]
  if (substr(x[[first]], 1, 1) == "-") {
    x[first] <- paste(x[first], curveleft)
    x[last] <- paste(x[last], curveright)
  } else {
    x[first] <- paste(x[first], curveright)
    x[last] <- paste(x[last], circleleft)
  }
  x
}

replace_arrows <- function(x) {
  first <- which(!is.na(x))[1]
  last <- length(x) + 1 - which(!is.na(rev(x)))[1]
  x[first] <- curveleft
  x[last] <- curveright
  x
}

proportion_scale <-
  list(aes(y=p),
       scale_y_continuous(
         "Response", breaks=c(0, 0.5, 1),
         labels=replace_arrows
         ),
       theme(axis.text.y=element_text(angle=90)))

facet_spacing_experiment <-
  facet_grid(spacing ~ exp_type,
               labeller=function(v, value) {
                 switch(  v
                        , spacing = paste(
                            "$"[!use_unicode],
                            "S = ",
                            format(value, digits=3),
                            if(use_unicode) "\u0080" else "^\\circ",
                            "$"[!use_unicode], sep=""
                            )
                        , exp_type = paste(as.character(value), "experiment"))
               }
             )

facet_spacing_rows <-
  facet_grid(spacing ~ .,
               labeller=function(v, value) {
                 browser()
                 paste("$S = ", format(value, digits=3),
                       if(use_unicode) "\u0080" else "^\\circ$", sep="")
               }
             )

balloon <- list(geom_point(aes(size=n))
                , scale_size_area()
                )

displacement_scale <-
  list( aes(x=displacement),
       scale_x_continuous(if (use_unicode) "\u0394x" else "$\\Delta x$"
                          , labels=add_arrows))

comp <- function(a, b) function(...) b(a(...))

prettyprint <- function(x) format(as.numeric(x), digits=3)

#I want a discrete color scape derived from a 4-point gradient, so I wrote:
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
    with_arg(name=if (use_unicode) "C" else "$C$",
             palette=color_pal,
             labels=add_arrows,
             discrete_scale("fill", "manual"),
             discrete_scale("colour", "manual")
             )
    )

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
