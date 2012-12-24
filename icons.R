## In this file we define some custom graphical objects for ggplot2.
library(plyr)
library(grid)
library(ggplot2)
library(scales)
library(proto)

#we also use these private symbols from ggplot2
Geom <- ggplot2:::Geom
ggname <- ggplot2:::ggname
.pt <- ggplot2:::.pt
identity_scale <- ggplot2:::identity_scale

## The "numdensity" geom draws a circle with some evenly spaced
##
## tickmarks sticking out of it.
## It accepts the following aesthetics:
## "number" -- the number of ticks.
## "size" -- the diameter of the circle, in mm.
## "center" the angle (in radians measured CCW from rightward) that the ticks stick out at.
## "tick_in" the length by which the ticks stick "in" to the circle (in mm)
## "tick_out" the length by which hte ticks stick out.
## "spacing" the angle (in radians, assuming you use an identity scale) between tickmarks.
## "in most cases youant to covary spacing and number.
## "weight" -- the thickness of the tickmarks.

## "color" and "fill" have their usual meanings.

geom_numdensity <- function(mapping=NULL, data=NULL,
                            stat="identity", position="identity", ...) {
  GeomNumdensity$new(mapping=mapping, data=data,
                     stat=stat, position=position, ...)
}

GeomNumdensity <- proto(Geom, {
  objname <- "numdensity"
  default_stat <- function(.) StatIdentity
  default_aes <- function(.) aes(  colour = "black", size=4, weight=0.5, linetype=1
                                 , tick_in=1, tick_out = 2, alpha = NA, shape=16, fill=NA
                                 , center = pi/2)
  required_aes <- c("x", "y", "spacing", "number")

  circles_and_ticks <- function(munched) {
    # legend defaults; use of "with" below will fall back on these for legends.
    x <- 0.5; y <- 0.5; spacing <- 2*pi/9; number <- 3;
    tickmarked <- adply(  munched, 1, function(row) {
      angles <- with(row, seq(  center - spacing/2 * (number - 1)
                              , by=spacing, length=number))
      with(row, data.frame(
                  begin.x = size/2 * cos(angles) * (size - tick_in) / size
                  , begin.y = size/2 * sin(angles) * (size - tick_in) / size
                  , end.x = size/2 * cos(angles) * (size + tick_out) / size
                  , end.y = size/2 * sin(angles) * (size + tick_out) / size
                  ))
    })
    gList(
      with(munched,
           circleGrob(
             x = unit(x, "native"), y = unit(y, "native")
             , r = unit(size/2, "mm")
             , gp = gpar(
                 col = colour
                 , fill = NA ###alpha(munched$fill, munched$alpha)
                 , lwd = weight * .pt
                 , lty = linetype
                 )))
      , with(tickmarked,
             segmentsGrob(
               x0 = unit(x, "native") + unit(begin.x, "mm")
               , y0 = unit(y, "native") + unit(begin.y, "mm")
               , x1 = unit(x, "native") + unit(end.x, "mm")
               , y1 = unit(y, "native") + unit(end.y, "mm")
               , gp = gpar(
                   col = colour
                   , lwd = weight * .pt
                   , lty = 1
                   , lineend = "square"
                   ))))
  }

  draw <- function(., data, scales, coordinates, ...) {
    munched <- coord_transform(coordinates, data, scales)

    ## draw the major glyphs out of minor glyphs.  the minor glyphs are
    ## IRL scaled, not data scaled, so there are separate 'offset.x'
    ## and 'offset.y (which have units of mm.)
    ## Otherwise, what defaults do we use though?

    (ggname)(.$my_name(), gTree(.$my_name(),
                                children = circles_and_ticks(munched)))
   }

  #We use geom_point's guide for now, but we will need a custom guide...
  guide_geom <- function(.) "numdensity"

  draw_legend <- function(., data, ...) {
    #draw the legend symbol...
    #here we have the legen
    data <- aesdefaults(data, .$default_aes(), list(...))

    #Each time, we are drawing one "thing"
    gTree(  gp = gpar(
              col = data$colour
              , lwd = data$weight * .pt
              , fill = alpha(data$fill, data$alpha)
              , lty = 1
              , lineend = "square")
          , children=circles_and_ticks(data)
          )
  }})

custom_geom_demo <- function() {
  #here's an example that uses these.
  library(ptools)
  theme_set(theme_bw(12, "Apple Symbols"))
  theme_update(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
  use_unicode <<- TRUE
  source("scales.R")
  library(binom)

  load("../modeling/data.Rdata")
  load("../modeling/slopeModel.RData")
  segment <- subset(data, exp_type=="numdensity" & subject %in% names(models))

  #this just illustrates the combinations of number and density.
  configurations <- unique(segment[c("target_spacing", "target_number_shown")])
  personalizations <- unique(segment[c("subject", "folded_displacement",
                                       "folded_direction_content")])

  ##For each personalization, collect the "segment" data to compare it against.
  alply(personalizations, 1,
        mkchain(  match_df(segment, ., on=names(.))
                , do.rename(folding=TRUE)
                , mkrates(splits=c(
                            "spacing", "content", "target_number_shown")))
        ) -> personalization.segment.data

  number_color_scale <-
    c(  list(aes(  color=factor(target_number_shown)) )
      , with_arg(name="Element\nnumber",
                 palette=color_pal,
                 labels=prettyprint,
                 discrete_scale("colour", "manual")
                 ))
  ##
  Map(alply(personalizations, 1, identity), personalization.segment.data,
      f=function(row, data){
        chain(data,
              with(binom.confint(n_cw, n, method="wilson", conf.level=0.75)),
              cbind(data,.)) -> data
        (  ggplot(data)
         + aes(x=spacing)
         + aes(y=p, ymin=lower, ymax=upper)
         + aes(spacing=spacing/20*3, center=pi/2, number=target_number_shown, )
         + identity_scale(continuous_scale(c("spacing", "number"),
                                           "identity", identity_pal()))
         #       + geom_errorbar()
         + geom_numdensity()
         #      + spacing_texture_scale
         + proportion_scale
         + number_color_scale
         + geom_line(aes(group=target_number_shown))
         ## + facet_wrap( ~ subject)
         + annotate(  "text", x=mean(range(data$spacing)), y=0
                    , label=sprintf("%s, dx=%03f", toupper(row$subject),
                        row$folded_direction_content)
                    , vjust=0, lineheight=6
                    )
         )}
      #make a gtable plotting it both ways...
      ) -> segment.raw.plots

  print(segment.raw.plots[[5]])
}
