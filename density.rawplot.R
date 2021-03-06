library(gtable)
segment.plot <- function(row, segment, x.axis="spacing",
                         spacing_breaks=unique(segment$spacing),
                         spacing_range=range(segment$spacing), ...) {
  #use template to provide a version for spacing and a version for number
  eval(template({
    errorbar <-
      ddply(segment, .(if("side" %in% names(segment)) "side" else c()), summarize,
            y = 0.5,
            x = max(.(
              switch(x.axis,
                     number=quote(target_number_shown),
                     spacing=quote(spacing),
                     extent=quote(extent)))),
            ymax = 0.5 + binom_se(min(n), 0.5),
            ymin = 0.5 - binom_se(min(n), 0.5))
    (ggplot(segment)
     + proportion_scale
     + geom_point(size = 6, color = "white")
     + geom_numdensity(
         aes(size = eccentricity,
             eccentricity = eccentricity,
             number = factor(target_number_shown),
             spacing = spacing,
             center = {
               if (exists("side"))
                 c(top = pi/2, bottom = 3 * pi/2, left = pi, right = 0)[side]
               else pi/2}
             ))
     + identity_scale(
         continuous_scale(
           "spacing", "spacing", identity,
           breaks=sort(unique(segment$spacing)),
           labels=function(x) format(x, digits=2)))
     + discrete_scale("number", "identity", identity_pal())
     + scale_size("Eccentricity", guide="none")
     #   + geom_text(size = 2.5, aes(label = target_number_shown), fontface = "bold")
     + annotate("text",
                size = 3, xmin = -Inf, xmax = Inf, x = Inf, y = -Inf,
                vjust = -0.5, hjust = 1.1,
                label = (sprintf("%s, \u0394x = %s, C = %s",
                                 toupper(row$subject),
                                 format(row$displacement, digits = 2),
                                 format(row$content, digits = 2))))
     + with_arg(data = errorbar
                , inherit.aes = FALSE
                , mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax)
                , show_guide = FALSE, geom_linerange(),
                geom_point(size = 4, shape = "+"))
     + list(
         ...(if ("side" %in% names(segment)) alist(facet_wrap(~side, nrow=2)) else list()),
         ...({
           ## The rest changes conditional on number or density scale
           switch(
             x.axis,
             number = alist(
               aes(x=target_number_shown, color=spacing),
               scale_x_continuous("Element number"),
               geom_line(
                 aes(group=spacing), show_guide=FALSE),
               continuous_scale(
                 "colour", "colour",
                 trans=log_trans(),
                 gradient_n_pal(
                   colours = muted(c("blue", "cyan", "yellow", "red"),
                                   l=70, c=180)),
                 breaks=spacing_breaks,
                 limits=spacing_range,
                 labels=function(x) format(x, digits=2)),
               guides(colour=guide_legend("Spacing (deg.)"),
                      spacing=guide_legend("Spacing (deg.)",
                        override.aes=list(eccentricity=20/3)),
                      number="none"),
               theme(legend.position="left")),
             spacing = alist(
               aes(x=spacing),
               scale_x_continuous("Spacing (deg.)", breaks=unique(segment$spacing),
                                  labels=function(x) format(x, digits=2)),
               number_color_scale,
               geom_line(
                 aes(group=target_number_shown),
                 show_guide=FALSE),
               guides(colour=guide_legend("Element\nnumber"),
                      number=guide_legend("Element\nnumber"),
                      spacing="none")
               ## theme(axis.text.y=element_blank(),
               ##       axis.ticks.y=element_blank(),
               ##       axis.title.y=element_blank())
               ),
             extent = alist(
               x=spacing * target_number_shown,
               scale_x_continuous("Spacing (deg.)", breaks=unique(segment$spacing),
                                  labels=function(x) format(x, digits=2)),
               number_color_scale,
               geom_line(
                 aes(group=target_number_shown),
                 show_guide=FALSE),
               guides(colour=guide_legend("Element\nnumber"),
                      number=guide_legend("Element\nnumber"),
                      spacing="none")
               )
             )
         }))
     )
  }))
}
