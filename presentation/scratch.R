load("../data.RData")
library(ggplot2)
library(plyr)
library(ptools)
library(psyphy)
theme_set(theme_bw())


add_arrows <- function(x) {
  x[1] <- paste(x[1], "$\\curvearrowleft$")
  x[length(x)] <- paste(x[length(x)], "$\\curvearrowright$")
  x
}

proportion_scale <-
  list(aes(y=p),
       scale_y_continuous(
         "Response", breaks=c(0, 0.5, 1),
         labels=add_arrows
         ),
       theme(axis.text.y=element_text(angle=90)))

#3 point gradient for this scale(red-gray-blue)
content_color_scale <-
  c(
    list(aes(color=factor(abs_direction_content), fill=factor(abs_direction_content))),
    with_arg(name="Direction\\\nContent",
             labels=add_arrows,
             scale_fill_discrete(),
             scale_color_discrete()
             )
    )

facet_spacing_rows <-
  facet_grid(target_spacing ~ .,
               labeller=function(v, value) {
                 paste("$S = ", format(value, digits=3), "^\\circ$", sep="")
               }
             )

balloon <- list(geom_point(), aes(size=n), scale_size_area())

displacement_scale <-
  list( aes(x=abs_displacement),
        scale_x_continuous("$\\Delta x$", labels=add_arrows))

(ggplot(example_content_rates)
 + displacement_scale
 + proportion_scale
 + balloon
 + content_color_scale + facet_spacing_rows
 + labs(title="Effect of direction content at narrow and wide spacing")
 )

# + with_arg(data=predictions,
#            geom_ribbon(color=NULL, alpha=0.3,
#                        aes(ymin=fit+se.fit, ymax=fit+se.fit,
#                            fill=factor(abs_direction_content))),
#            geom_line(aes(y=fit)))
