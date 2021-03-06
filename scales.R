library(vadr)
library(scales)
library(ggplot2)
source("colors.R")

if ( !exists("use_unicode") || use_unicode) {
  curveleft = "\u21BA"
  curveright = "\u21BB"
  circleleft = "\u27F2"
  circleright = "\u27F3"
  use_unicode <- TRUE
} else {
  curveleft = "$\\curvearrowleft$"
  curveright = "$\\curvearrowright$"
  circleleft = "$\\curvearrowleft$"
  circleright = "$\\curvearrowright$"
}

combine_arrows <- function(x, combine) {
  first <- which(!is.na(x))[1]
  last <- length(x) + 1 - which(!is.na(rev(x)))[1]
  if (substr(x[[first]], 1, 1) %in% c("-", "0")) {
    x[first] <- combine(circleleft, x[first])
    x[last] <- combine(circleright, x[last])
  } else {
    x[first] <- combine(curveright, x[first])
    x[last] <- combine(circleright, x[last])
  }
  x
}

append_arrows <- function(x) {
  combine_arrows(x, combine=function(x, y) paste(y, x, sep=""))
}

prepend_arrows <- function(x) {
  combine_arrows(x, combine=function(x, y) paste(x, y, sep=" "))
}

newline_arrows <- function(x) {
  combine_arrows(x, combine = function(x, y) paste(y, x, sep="\n"))
}

replace_arrows <- function(x) {
  combine_arrows(x, combine=function(x, y) x)
}

setup_theme <- function(base_size=theme_get()$text$size) {
#  if (interactive()) base_size <- 16
  theme_set(theme_bw(base_size))
  theme_update(
    ## axis.title.x = element_text_with_symbols(
    ##   size=list(c(12, 24)), family=list(c("Myriad Pro", "Apple Symbols")),
    ##   lineheight=list(c(0.9, 1.5)), stack="h")
    ## , axis.title.y = element_text_with_symbols(
    ##     size=list(c(12, 24)), family=list(c("Myriad Pro", "Apple Symbols")),
    ##     lineheight=list(c(0.9, 1.5)), angle=0, stack="v")
    #  axis.text.x = element_text(size=base_size*0.8, family="Apple Symbols", vjust=1)
    #, axis.text.y = element_text(size=base_size*0.8, family="Apple Symbols", hjust=1)
    # ,
    panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    ## , axis.text.x = element_text_with_symbols(
    ##     size=list(c(10, 20)), family=list(c("Myriad Pro", "Apple Symbols")),
    ##     angle=0, stack="v", vjust=1)
    ## , axis.text.y = element_text_with_symbols(
    ##     size=list(c(10, 20)), family=list(c("Myriad Pro", "Apple Symbols")),
    ##     angle=0, stack="h", hjusont=1)
    )
}

blank_proportion_scale <- list(
  aes(y=p),
  scale_y_continuous(
    "P(Response CW)", breaks = c(0, 0.5, 1), labels = c("", "", "")))

proportion_scale <-
  list(aes(y=p),
       scale_y_continuous(
         "P(Response CW)", breaks=c(0, 0.5, 1)
         #, labels=replace_arrows
         ),
       coord_cartesian(ylim=c(-0.1,1.1)),
       theme(axis.text.y=element_text(angle=0)))

facet_labeller <- function(v, value)
  switch(v,
         subject= sprintf("Observer %s", toupper(value)),
         spacing = paste(
           "$"[!use_unicode],
           sprintf("Spacing = %.02f", value),
           #if(use_unicode) "\u0080" else "^\\circ",
           "$"[!use_unicode], sep=""),
         exp_type = paste(as.character(value), "experiment"))

facet_spacing_experiment <-
    list(facet_grid(spacing ~ exp_type,
                    labeller=facet_labeller
                    ),
         theme(axis.text.y = element_text(size=rel(0.7)),
               strip.text.y=element_text(size=rel(0.7))))

facet_spacing_rows <-
    list(facet_grid(spacing ~ ., labeller=facet_labeller),
         theme(axis.text.y = element_text(size=rel(0.7)),
               strip.text.y=element_text(size=rel(0.7))))

facet_spacing_subject <-
    list(facet_grid(spacing ~ subject, labeller=facet_labeller),
         theme(axis.text.y = element_text(size=rel(0.6))))


prediction_layer <- function(predictions)
      with_arg(data=predictions,
           geom_ribbon(color="transparent", alpha=0.2,
                       aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
           geom_line(aes(y=fit)))

balloon <- list(geom_point(aes(size=n_obs))
                , scale_size_area("N")
                )

 if (exists("geom_numdensity")) {
  geometric_shape_scale <-
    list(
      aes(  number=factor(target_number_shown)
          , spacing=spacing
          , eccentricity=eccentricity)
      , geom_numdensity(size=7)
      , identity_scale(continuous_scale("spacing", "spacing",
                                        identity, name="Spacing")))
}

decision_contour <-
  list(
    aes(z = pred, fill = pred),
#    geom_contour(breaks=seq(0.1, 0.9, 0.2), size=0.35, color="white", alpha=0.5),
    geom_contour(breaks=seq(0.1, 0.9, 0.2), size=0.35, linetype="11", color="#FFFF00"),
    scale_fill_gradientn("Responses CW", colours=decision_colors,
                         values=decision_values,
                         breaks = seq(0.1, 0.9, 0.2),
                         rescale=function(x,...) x,
                         oob=squish, space="rgb")
    )

no_padding <- with_arg(expand=c(0,0), scale_x_continuous(), scale_y_continuous())
#no_padding <- list(scale_x_continuous(expand=c(0,0)), scale_y_continuous(expand=c(0,0)))
no_grid <- theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

displacement_scale <-
  list( aes(x=displacement),...=
       scale_x_continuous(
         paste0("Envelope motion") # (\u0394x)
#                          , labels=newline_arrows
#         , labels=function(x) sprintf("%s\u00B0", x)
         ))
displacement_scale_nopadding <- displacement_scale
displacement_scale_nopadding[[2]]$expand <- c(0,0)
y_nopadding <- scale_y_continuous(expand = c(0,0))

comp <- function(a, b) function(...) b(a(...))

prettyprint <- function(x) format(as.numeric(x), digits=3)

pretty_strip <- function(var, value) {
  switch(class(value),
         numeric=format(as.numeric(value), digits=2),
         factor=as.character(value),
         value)
}

#I want a discrete color scale derived from a 4-point gradient, so I wrote:
discretize <- function(pal) {
  function(n) {
    pal(seq(0,1,length=n))
  }
}

sat <- function(x, min, max)
  ifelse(x > min, ifelse(x < max, x, max), min)

clip <- function(x, min, max, replace=NA)
  ifelse(x > min, ifelse(x < max, x, replace), replace)

mapColors <- function(data, map,
                      min=base::min(data), max=base::max(data), limit=sat, ...,
                      sampling=1024) {
  data <- limit(data, min, max, ...)
  colors <- map(sampling)
  index <- round((data - min) / (max-min) * (sampling - 1)) + 1
  colors[index]
}

color_pal <- #function(n) rainbow(n, start=0, end=0.75)
  discretize(gradient_n_pal(
               colours=muted(c("blue", "cyan", "yellow", "red"), l=70, c=180)))

number_color_scale <-
  c(  list(aes(  color=factor(target_number_shown)))
       , with_arg(name="Element\nnumber",
                  palette=discretize(gradient_n_pal(
                    muted(c("cyan", "magenta", "yellow"),
                          l=70, c=180))),
                  labels= prettyprint,
                  discrete_scale("fill", "manual"),
                  discrete_scale("colour", "manual")
                  ))

number_color_alt_scale <-
    c(  list(aes(  color=factor(target_number_shown)))
       , with_arg(name="Element\nnumber",
                  palette=discretize(gradient_n_pal(
                    muted(c("blue", "cyan", "yellow", "red"), l=70, c=180))),
                  labels= prettyprint,
                  discrete_scale("fill", "manual"),
                  discrete_scale("colour", "manual")
                  ))

spacing_color_scale <-
  c(
    list(aes(color=factor(spacing),
             fill=factor(spacing))),
    with_arg(name="Spacing",
             palette=color_pal,
             labels=function(x)sprintf("%0.2f", as.numeric(x)), #\u00B0
             discrete_scale("fill", "manual"),
             discrete_scale("colour", "manual")
             )
    )

spacing_texture_scale <-
  list(aes(linetype=factor(spacing)),
       discrete_scale(name="Spacing", "linetype", "linetype_m", function(x) {
         paste(y <- c(1:9, LETTERS)[seq_len(x)], "1", sep="")
       }, labels=prettyprint)
       )

spacing_shape_scale <-
    list(
      aes(shape=factor(spacing)),
       scale_shape_manual("Spacing", labels=prettyprint, values=c(11, 5, 0, 2, 6, 1))
      )

spacing_scale_x <-
  list( aes(x=spacing),
       scale_x_continuous(name="Spacing"
#                          , labels=newline_arrows
                          ))

spacing_scale_x_nopadding <-
  list( aes(x=spacing),
       scale_x_continuous(name="Spacing", expand=c(0,0)
#                          , labels=newline_arrows
                          ))

spacing_scale_y <-
  list( aes(y=spacing),
       scale_y_continuous(name="Spacing"
#                          , labels=newline_arrows
                          ))

spacing_scale_y_nopadding <-
  list( aes(y=spacing),
       scale_y_continuous(name="Spacing", expand=c(0,0)
#                          , labels=newline_arrows
                          ))

content_color_scale <-
  c(
    list(aes(color=factor(content),
             fill=factor(content))),
    with_arg(name="Carrier\nStrength",
             palette=color_pal,
#             labels=append_arrows,
             discrete_scale("fill", "manual"),
             discrete_scale("colour", "manual")
             )
    )

#a palette with a waterline in the middle
energy_palette <- gradient_n_pal(
                    colours=c("#00DDDD", "blue", "black", "red", "#DDDD00"),
                    values=c(0, 0.4999,0.5,0.5001,1))

energy.colors <- function(n) energy_palette(seq(0, 1, length=n))

content_scale_continuous_waterline <-
  scale_color_gradientn("Carrier\nStrength",
                        colours=c("#00DDDD", "blue", "black", "red", "#DDDD00"),
                        values=c(0, 0.4999,0.5,0.5001,1), limits=c(-1,1),
                        labels=append_arrows
                        )

displacement_scale_continuous_waterline <-
  scale_color_gradientn("Displacement",
                        colours=c("#00DDDD", "blue", "black", "red", "#DDDD00"),
                        values=c(0, 0.4999,0.5,0.5001,1),
                        labels=append_arrows
                        )

label_count <- function(data, group, countvar=n_obs,
                        x = -Inf, y = Inf, size=3, vjust = 1.5, hjust = -0.2)
  geom_text(data=eval(bquote(ddply(
              data, group, summarize,
              .count=sum(.(substitute(countvar)))))),
    inherit.aes=FALSE, show_guide=FALSE,
    aes(label = sprintf("N = %d", .count)),
    x = x, y = y, size= size, vjust = vjust, hjust = hjust)

content_scale <-
  list(aes(x=content),
       scale_x_continuous("Carrier\nStrength",
                          #labels=newline_arrows,
                          expand=c(0,0))
       )
content_scale_nopadding <- content_scale
content_scale_nopadding[[2]]$expand <- c(0,0)

content_scale_y <-
  list( aes(y=content), scale_y_continuous(name="Carrier strength"
                          #, labels=newline_arrows
                          ))

content_scale_y_nopadding <-
  list( aes(y=content),
       scale_y_continuous(name="Carrier strength", expand=c(0,0)
                          #, labels=newline_arrows
                          ))

ribbon <- list(
            geom_ribbon(color="transparent", alpha=0.2,
                        aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
            geom_line(aes(y=fit)))

ribbonf <- function(data)
  with_arg(data=data,
           geom_ribbon(color="transparent", alpha=0.2,
                       aes(y=fit, ymin=fit-se.fit, ymax=fit+se.fit)),
           geom_line(aes(y=fit)))

binom_pointrange <- function(...)
  geom_pointrange(aes(ymin = binom_se_lower(n_obs, bound_prob(p)),
                      ymax = binom_se_upper(n_obs, bound_prob(p))),
                  ...)


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
