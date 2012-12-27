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
    x <- 0.5; y <- 0.5; spacing <- pi/9; number <- 3;
    tickmarked <- adply(  munched, 1, function(row) {
      angles <- with(row, seq(  center - spacing/2 * (as.numeric(number) - 1)
                              , by=spacing, length=as.numeric(number)))
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
                 , fill = alpha(munched$fill, munched$alpha)
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

    ## The glyphs are IRL scaled, not data scaled, so there are
    ## separate 'offset.x' and 'offset.y (which have units of mm.)
    ## Otherwise, what defaults do we use though?

    ggname(.$my_name(), gTree(.$my_name(),
                              children = circles_and_ticks(munched)))
   }

  guide_geom <- function(.) "numdensity"

  draw_legend <- function(., data, ...) {
    #draw the legend symbol...
    data <- aesdefaults(data, .$default_aes(), list(...))
    gTree(  gp = gpar(
              col = data$colour
              , lwd = data$weight * .pt
              , fill = alpha(data$fill, data$alpha)
              , lty = 1
              , lineend = "square")
          , children=circles_and_ticks(data)
          )
  }})

## here is a theme object for mixing different fonts using a fallback
## unicode symbols into the display of
## a text label.

#

#' This is an extension of element_text which can specify multiple
#' fonts to use, for example a symbol font to fall back on if a glyph
#' is not available in the main font.
#'
#' @param family font family. This may be a vector, in which case it
#' specifies a list of font faces in priority order.
#' @param face font face ("plain", "italic", "bold", "bold.italic")
#' @param colour text colour
#' @param size text size (in pts)
#' @param hjust horizontal justification (in [0, 1])
#' @param vjust vertical justification (in [0, 1])
#' @param angle angle (in [0, 360])
#' @param lineheight line height
#' @param color an alias for \code{colour}
#' @param fallback A function (glyph, family, face) that returns FALSE
#' if the specified font should not be used.
#' @param symbolsize
#' @export
element_text_with_symbols <- function(family = NULL, face = NULL, colour = NULL,
  size = NULL, hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
  color = NULL, fontcheck=NULL) {

  if (!is.null(color))  colour <- color
  structure(
    list(family = family, face = face, colour = colour, size = size,
      hjust = hjust, vjust = vjust, angle = angle, lineheight = lineheight, fontcheck=renderable_char),
    class = c("element_text_with_symbols", "element_text", "element")
  )
}

renderable <- function(label, family, face, size) {
  lapply(strsplit(label, ""), vapply, renderable_char, FALSE, family, face, size)
}

renderable_char <- function(char, family, face, size) {
  is.character(char) && (nchar(char) == 1) || stop("Must be called on individual chars.")
  length(charToRaw(char)) == 1
}

#return the label, split into
#family, face and size are
pbutfirst <- function(x) if (length(x) > 1) x[-1] else x

renderable_split <- function(label, family=NULL, face=NULL, size=NULL,
                             fontcheck = renderable_char) {
  renderable_p = vapply(strsplit(label, "")[[1]], renderable_char, FALSE,
    family[[1]], face[[1]], size[[1]])

  #each "TRUE" gets the first family element. Each "FALSE" get s
  renderable_splits = mutate(
    as.data.frame(unclass(rle(renderable_p)))
    , ends = cumsum(lengths)
    , starts = ends-lengths+1
    )

  mdply(renderable_splits,
        function(lengths, values, ends, starts) {
          if (values || all(vapply(list(family, face, size), length, 0) == 1)) {
            data.frame(labels=substr(label, starts, ends),
                 family=family[[1]], face=face[[1]], size=size[[1]])
          } else {
            renderable_split(substr(label, starts, ends),
                             family=pbutfirst(family),
                             face=pbutfirst(face),
                             size=pbutfirst(size))
          }
        })
}

family <- c("Myriad", "Apple Symbols")
face <- "plain"
size <- c(12, 24)

#' @S3method element_grob element_text_with_symbols
element_grob.element_text_with_symbols <-
  function(element, label = "", x = NULL, y = NULL,
           family = NULL, face = NULL, colour = NULL, size = NULL,
           hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
           default.units = "npc", ...)
 {

  vj <- vjust %||% element$vjust
  hj <- hjust %||% element$hjust

  angle <- angle %||% element$angle
  if (is.null(angle)) {
    stop("Text element requires non-NULL value for 'angle'.")
  }
  angle <- angle %% 360
  
  if (angle == 90) {
    xp <- vj
    yp <- hj
  } else if (angle == 180) {
    xp <- 1 - hj
    yp <- vj
  } else if (angle == 270) {
    xp <- vj
    yp <- 1 - hj
  }else {
    xp <- hj
    yp <- vj
  }

  x <- x %||% xp
  y <- y %||% yp

  #split the string up into substrings based on the criterion
    if (length(label) != 1) {
    stop("Multiple labels are specified, this is not supported.")
  }

  #divide the string up into renderable and non-renderable characters.
  
  if(renderable) {}

  #then make a tree of text grobs that attach to each others' ends

  # The gp settings can override element_gp
  gp <- gpar(fontsize = size, col = colour,
    fontfamily = family, fontface = face,
    lineheight = lineheight)
  element_gp <- gpar(fontsize = element$size, col = element$colour,
    fontfamily = element$family, fontface = element$face,
    lineheight = element$lineheight)

  textGrob(
    label, x, y, hjust = hj, vjust = vj,
    default.units = default.units,
    gp = modifyList(element_gp, gp),
    rot = angle, ...
  )
}

#

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
      , with_arg(  name="Element\nnumber"
                 , palette=color_pal
                 , labels=prettyprint
                 , discrete_scale("colour", "manual")
                 ))
  ##
  Map(alply(personalizations, 1, identity), personalization.segment.data,
      f=function(row, data){
        chain(data,
              with(binom.confint(n_cw, n, method="wilson", conf.level=0.75)),
              cbind(data,.)) -> data
        rad <- 20/3
        (  ggplot(data)
         + aes(x=spacing)
         + aes(y=p, ymin=lower, ymax=upper)
         + aes(spacing=spacing, number=factor(target_number_shown))
         + continuous_scale("spacing", "spacing", identity,
                            , rescaler=function(x, from) {
                              rescale(x , from = from , to = from/rad)
                            }
                            , name="Spacing")
         + identity_scale(discrete_scale("number", "identity",
                                         identity_pal(), name="Element\nnumber"))
         #       + geom_errorbar()
         + geom_numdensity()
         #      + spacing_texture_scale
         + proportion_scale
         + number_color_scale
         + geom_line(aes(group=target_number_shown), show_guide=FALSE)
         ## + facet_wrap( ~ subject)
         + annotate(  "text", x=mean(range(data$spacing)), y=0
                    , label=sprintf("%s, dx=%03f", toupper(row$subject)
                        , row$folded_direction_content)
                    , vjust=-0.5, size=6
                    )
         )}
      #make a gtable plotting it both ways...
      ) -> segment.raw.plots

  print(segment.raw.plots[[5]])
}
