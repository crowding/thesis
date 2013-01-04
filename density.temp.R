bind[arow, adata] <- dlply_along(segment.rates.sided, segment.experiment.vars, list)[[2]]

segment.plot <- function(row, data, number=FALSE) {
  eval(template({
    errorbar <-
      ddply(data, "side", summarize,
            y = 0.5,
            x = max(.(
              if (number) quote(target_number_shown) else quote(spacing))),
            ymax = 0.5 + binom_se(min(n), 0.5),
            ymin = 0.5 - binom_se(min(n), 0.5))
    (ggplot(data)
     + proportion_scale
     + geom_point(size = 6, color = "white")
     + geom_numdensity(
         aes(size = eccentricity,
             number = factor(target_number_shown),
             spacing = spacing/eccentricity * 20/3,
             center = {
               if (exists("side"))
                 c(top = pi/2, bottom = 3 * pi/2, left = pi, right = 0)[side]
               else pi/2}
             ))
     + continuous_scale(
         "spacing", "spacing",
         function(x) {print(x); x},
         breaks=unique(data$spacing),
         #     labels=function(x) format(x, digits=2),
         rescaler = function(x, from) {
           ##THIS_IS_GETTING_AN_NA_ARGUMENT in scale_map.continuous
           #if (any(is.na(x))) stop("WHAT")
           rescale(x, from = from, to = from/20 * 3)
         })
     + discrete_scale("number", "identity", identity_pal())
     + scale_size("Eccentricity", guide="none")
     + list(...({
       #choose the rest conditional on number or density scale
       if (number)
         alist(
           aes(x=target_number_shown, color=spacing),
           scale_x_continuous("Element number"),
           geom_line(
             ##color="black",
             aes(group=spacing), show_guide=FALSE),
           continuous_scale("colour", "colour",
                            gradient_n_pal(muted(c("cyan", "magenta", "yellow"),
                                                 l=70, c=180)),
                            breaks=unique(data$spacing),
                         #   labels=function(x) format(x, digits=2)
                            ),
           guides(colour=guide_legend("Spacing (deg.)"),
                  spacing=guide_legend("Spacing (deg.)"),
                  number="none"
                  )
           )
       else
         alist(
           aes(x=spacing),
           scale_x_continuous("Spacing (deg.)", breaks=unique(data$spacing),
                              labels=function(x) format(x, digits=2)),
           number_color_scale,
           geom_line(
             aes(group=target_number_shown),
             show_guide=FALSE),
           guides(colour=guide_legend("Element\nnumber"),
                  number=guide_legend("Element\nnumber"),
                  spacing="none"))
     }))
     #   + geom_text(size = 2.5, aes(label = target_number_shown), fontface = "bold")
     + facet_wrap(~side)
     + annotate("text",
                size = 3, xmin = -Inf, xmax = Inf, x = Inf, y = -Inf,
                vjust = -0.5, hjust = 1.1,
                label = (sprintf("%s, \u0394x = %s, C = %s",
                                 toupper(row$subject),
                                 format(row$displacement, digits = 2),
                                 format(row$content, digits = 2))))
     + with_arg(data = errorbar,
                inherit.aes = FALSE
                , mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax)
                , show_guide = FALSE, geom_linerange(),
                geom_point(size = 4, shape = "+"))
     + theme(legend.position=.(if (number) "left" else "right"))
     )}))
}
dev.set(2); grid.newpage(); dev.flush()
print(segment.plot(arow, adata, number=FALSE))
dev.set(3); grid.newpage(); dev.flush()
print(segment.plot(arow, adata, number=TRUE))
