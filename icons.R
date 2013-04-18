## In this file we define some custom graphical objects for ggplot2.
options(encoding = "UTF-8")
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
`%||%` <- ggplot2:::`%||%`

###Until ggplot2 merges pull request #748, monkey patch it.
if (packageVersion("ggplot2") < package_version('0.9.4')) {
  e <- environment(ggplot2:::combine_elements)
  unlockBinding("combine_elements", e)
  evalq(envir = e, {
    combine_elements <- function(e1, e2) {

      # If e2 is NULL, nothing to inherit
      if (is.null(e2))  return(e1)

      # If e1 is NULL, or if e2 is element_blank, inherit everything from e2
      if (is.null(e1) || inherits(e2, "element_blank"))  return(e2)

      # If e1 has any NULL properties, inherit them from e2
      n <- vapply(e1[names(e2)], is.null, logical(1))
      e1[names(n[n])] <- e2[names(n[n])]

      # Calculate relative sizes
      if (is.rel(e1$size)) {
        e1$size <- e2$size * unclass(e1$size)
      }

      e1
    }
  })
}

## A circle geom (that allows better control over line weight and fill
## than geom_point...)
geom_circle <- function(mapping=NULL, data=NULL,
                        stat="identity", position="identity", ...) {
  GeomCircle$new(mapping = mapping, data = data,
                 stat = stat, position = position, ...)
}

GeomCircle <- proto(Geom, {
  objname <- "circle"
  default_stat <- function(.) StatIdentity
  default_aes <- function(.)
    aes(colour = "black", fill=NA, size = 4, weight = 0.5, linetype = 1,
        alpha = NA, shape = 16, fill = NA)
  required_aes <- c("x", "y")

  draw <- function(., data, scales, coordinates, ...) {
    munched <- coord_transform(coordinates, data, scales)
    ggname(.$my_name(), circles_grob(munched))
  }

  circles_grob <- function(munched) {
    x <- 0.5; y <- 0.5
    with(munched
         , circleGrob(
             x = unit(x, "native"), y = unit(y, "native")
           , r = unit(size/2, "mm")
           , gp = gpar(
               col = alpha(munched$colour, munched$alpha)
               , fill = alpha(munched$fill, munched$alpha)
               , lwd = weight * .pt
               , lty = linetype
               , lineend = "butt"
               )))
  }

  guide_geom <- function(.) "circle"

  draw_legend <- function(., data, ...) {
    #draw the legend symbol...
    data <- aesdefaults(data, .$default_aes(), list(...))
    circles_grob(data)
  }
})

demo_circle <- function() {

  circles <- expand.grid(
    a=1:5,
    b=1:5,
    c=1:5,
    d=1:5,
    e=1:5
    )
  (ggplot(circles)
   + aes(x=a, y=b, color=c, fill=d, size=a, linetype=factor(b))
   + facet_grid(c ~ d)
   + geom_circle())

}


## 

## The "numdensity" geom draws a circle with some evenly spaced
## tickmarks sticking out of it.
##
## It accepts the following aesthetics:
## "number" -- the number of ticks.
## "size" -- the diameter of the circle, in mm.
## "center" the angle (in radians measured CCW from rightward) that the ticks
## cluster around.
## "tick_in" the length by which the ticks stick "in" to the circle (in mm)
## "tick_out" the length by which the ticks stick out.
## "spacing" the angle (in radians, assuming you use an identity scale)
## between tickmarks.
## "eccentricity" a scaling factor for "size"
## "weight" -- the thickness of the tickmarks.
## "color" and "fill" have their usual meanings.

geom_numdensity <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity", ...) {
  GeomNumdensity$new(mapping = mapping, data = data,
                     stat = stat, position = position, ...)
}

GeomNumdensity <- proto(Geom, {
  objname <- "numdensity"
  default_stat <- function(.) StatIdentity
  default_aes <- function(.)
    aes(colour = "black", size = 4, weight = 0.5, linetype = 1, tick_in = 1,
        tick_out = 2, alpha = NA, shape = 16, fill = NA, center = pi/2,
        eccentricity = 1)
  required_aes <- c("x", "y", "spacing", "number")

  circles_and_ticks <- function(munched) {
    # legend defaults; use of "with" below will fall back on these for legends.
    x <- 0.5; y <- 0.5; spacing <- pi/9; number <- 3;
    tickmarked <- adply(  munched, 1, function(row) {
      angles <- with(row, seq(
        center - spacing/eccentricity/2 * (as.numeric(number) - 1),
        by=spacing/eccentricity, length=as.numeric(number)))
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
  }
})

## The Geom_envelope draws a Gabor function, oriented horizontally.
## Aesthetics:
## 'sigma' how may sigmas the icon includes in either direction (default t0)
## 'width', 'height' size of the Gabor icon, in mm
## 'linetype', 'alpha', 'colour' as with geom_segment;
## 'samples' the number of points along which to sample,
## 'freq' the frequency of the osciallation ( in sigma units )
## 'phase' the phase of the oscillation
## 'center' shifts the envelope

geom_gabor <- function(mapping=NULL, data=NULL,
                        stat="identity", position="identity", ...) {
  GeomGabor$new(mapping = mapping, data = data,
                state = stat, position = position, ...)
}

GeomGabor <- proto(Geom, {
  objname <- "gabor"
  default_stat <- function(.) StatIdentity
  default_aes <- function(.)
    aes(width = 5, height = 5, size = 1,
        linetype = 1, sigma = 2, freq = 0, phase = 0,
        center = 0,
        colour = "black", samples=30, alpha=NA)
  required_aes <- c("x", "y")

  gabor_lines <- function(munched) {
    # legend defaults; use of "with" below will fall back on these for legends.
    x <- 0.5; y <- 0.5

    .id <- 0
    gabors <- adply(
      munched, 1, function(x) {
        coo = with(x, seq(-sigma, sigma, length = samples))
        merge(
          mutate(x, .id = .id <<- .id + 1),
          with(x, data.frame(
            xx = coo/sigma,
            yy = (exp(-(coo-center)^2)
                  * Re(exp(complex(
                    imaginary = (phase + freq * (coo - center)))))))),
          by=c())
      })
    with(gabors, polylineGrob(
      x = unit(xx * width/2, "mm") + unit(x, "native"),
      y = unit(yy * height/2, "mm") + unit(y, "native"),
      id = .id,
      gp = (gpar(col = alpha(colour, alpha), lwd = size * .pt, lty=linetype))
      ))
  }

  draw <- function(., data, scales, coordinates, ...) {
    munched <- coord_transform(coordinates, data, scales)
    ggname(.$my_name(), gabor_lines(munched))
  }

  guide_geom <- function(.) ("gabor")

  draw_legend <- function(., data, ...) {
    data <- aesdefaults(data, .$default_aes(), list(...))
    gabor_lines(data)
  }
})

# a quiver geom. Quiver plots can also be done with geom_segment, but
# these are different in that they are scales to physical units rather
# than graph units.
geom_quiver <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", arrow = NULL, lineend = "butt",
                        length = 5, ...) {
  GeomQuiver$new(mapping = mapping, data = data, stat = stat,
                 position = position, arrow = arrow, lineend = lineend,
                 length = length, ...)
}

GeomQuiver <- proto(Geom, {
  objname <- "quiver"
  default_state <- function(.) StatIdentity
  default_aes <- function(.)
    aes(length=5, size=0.5, linetype=1, colour="black", alpha=NA,
        xstart=0, ystart=0)
  required_aes <- c("x", "y", "xlen", "ylen")
  guide_geom <- function(.) "quiver"

  draw <- function(., data, scales, coordinates, arrow = NULL,
                   lineend = "butt", linetype=1, ...) {
    munched <- coord_transform(coordinates, data, scales)
    quiver_segments(.$my_name(), munched, linetype=linetype,
                    lineend=lineend, arrow=arrow, ...)
  }

  quiver_segments <- function(name, munched, linetype=1,
                              lineend="butt", arrow=NULL, ...) {
    #these defaults apply when drawing guides
    x <- 0.5; y <- 0.5; xlen <- 1; ylen <- 0; length <- 3; lineend = "butt"
    segments <- with(munched, segmentsGrob(
      x0 = unit(x, "native") + unit(xstart * length, "mm"),
      y0 = unit(y, "native") + unit(ystart * length, "mm"),
      x1 = unit(x, "native") + unit((xstart + xlen) * length, "mm"),
      y1 = unit(y, "native") + unit((ystart + ylen) * length, "mm"),
      gp = gpar(col = alpha(colour, alpha), lwd = size * .pt,
        lty = linetype, lineend = lineend),
      arrow = arrow
      ))
    ggname(name, segments)
  }

  draw_legend <- function(., data, ...) {
    data <- aesdefaults(data, .$default_aes(), list(...))
    quiver_segments(.$my_name(), data)
  }
})

quiverDemo <- function() {
  data <- chain(
    expand.grid(x = seq(-5, 5, 0.25), y = seq(-5, 5, 0.25)),
    mutate(dx = sin(y),
           dy = cos(x))
    )
  (ggplot(data) + aes(x=x, y=y, xlen=dx, ylen=dy, xstart=-dx/2, ystart=-dy/2) +
   geom_quiver(size=0.25, arrow=arrow(length=unit(1, "mm"))))
}

iconsDemo <- function() {
  #build some ways to depict "envelope" and "carrier" motion
  chain(data.frame(envelope.shift = c(-1, 0, 1),
                   envelope = c("left", "null", "right")),
        merge(data.frame(phase.shift = c(-1, 0, 1),
                         carrier = c("left", "null", "right"))),
        merge(data.frame(freq = c(0, 0, 3*pi/2),
                         linetype=c("11", "11", "1"),
                         phase = c(0, pi, 0))),
        subset(ifelse(carrier != "null", freq != 0, TRUE)),
        subset(ifelse(envelope != "null", freq == 0, TRUE))
        ) -> carriers

  chain(carriers,
        mutate(xstart = envelope.shift - phase.shift/2,
               xlen = phase.shift + envelope.shift),
        subset(phase==0 & xlen != 0)) -> arrows

  plot(ggplot(carriers)
       + aes(x = carrier, y = envelope, freq = freq, linetype = linetype,
             group=freq, phase = phase)
       + geom_gabor(sigma = 2, size = 0.5,
                    width = 15, height = 15, samples = 50, show_guide=FALSE)
       + geom_quiver(data=arrows, aes(xlen = xlen, ylen=0, xstart=xstart),
                     size=2, colour="white", lineend="square", linetype=1,
                     arrow=arrow(length=unit(2, "mm")), show_guide=FALSE)
       + geom_quiver(data=arrows, aes(xlen = xlen, ylen=0, xstart=xstart),
                     size=1, colour="red", lineend="square", linetype=1,
                     arrow=arrow(length=unit(2, "mm")),
                     show_guide=FALSE)
       + theme_bw()
       + theme(panel.grid.major = element_blank()
               , panel.grid.minor = element_blank(), aspect.ratio=1)
       )
}



## here is a theme object for mixing different fonts using a fallback
## unicode symbols into the display of
## a text label.

#

#' This is an extension of element_text which can specify multiple
#' fonts to use, for example a symbol font to fall back on if a glyph=
#' is not available in the main font. For example, using a Unicode
#' symbol from one font but the rest of text from another font.
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
#' @param stack -- "h" or "v" to stack chunks horizontally or vertically.
#' @export
element_text_with_symbols <-
  function(
    family = NULL, face = NULL, colour = NULL,
    size = NULL, hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
    stack=NULL, fontcheck=renderable_char, color=NULL) {

  if (!is.null(color))  colour <- color
  structure(
    list(family = family, face = face, colour = colour, size = size,
      hjust = hjust, vjust = vjust, angle = angle, lineheight = lineheight,
          stack = stack, fontcheck = fontcheck),
    class = c("element_text_with_symbols", "element_text", "element")
  )
}

#Decide if a character is renderable in the given face. Or it
#should. I just punt and say that things outside ASCII are not
#renderable.
renderable_char <- function(char, family, face, size) {
  is.character(char) && (nchar(char) == 1) || stop("Must be called on individual chars.")
  length(charToRaw(char)) == 1
}

#given s string, and vectors of fonts to try, split the string up into
#chunks returning a data frame.
renderable_split <- function(label, family=NULL, face=NULL, size=NULL,
                             lineheight=NULL, fontcheck = renderable_char) {
  pbutfirst <- function(x) {
    if (length(x) > 1) x[-1] else x
  }
  pfirst <- function(x) if (is.null(x)) list(x) else x[[1]]

  renderable_p = vapply(strsplit(label, "")[[1]], fontcheck, FALSE,
    family[[1]], face[[1]], size[[1]])

  #each "TRUE" gets the first font spec.
  #Each "FALSE" we try over again with the second, if it exists.
  renderable_splits = mutate(
    as.data.frame(unclass(rle(renderable_p)))
    , ends = cumsum(lengths)
    , starts = ends-lengths+1
    )

  mdply(renderable_splits,
        function(lengths, values, ends, starts) {
          lengths <- vapply(list(family, face, size, lineheight), length, 0)
          if (values
              || all(lengths <= 1)) {
            #no more fonts to try
            quickdf(list( labels=substr(label, starts, ends)
                         , family=pfirst(family)
                         , face=pfirst(face)
                         , size=pfirst(size)
                         , lineheight=pfirst(lineheight)))
          } else {
            #descend through the font list
            renderable_split(  substr(label, starts, ends)
                             , family=pbutfirst(family)
                             , face=pbutfirst(face)
                             , size=pbutfirst(size)
                             , lineheight=pbutfirst(lineheight)
                             )
          }
        })
}

#' @S3method element_grob element_text_with_symbols
#' Draws a text label with potentially mixed fonts.
element_grob.element_text_with_symbols <- function(
           element, label = "", x = NULL, y = NULL,
           family = NULL, face = NULL, colour = NULL, size = NULL,
           hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
           fontcheck = renderable_char, stack=NULL,
           default.units = "npc", ...) {
  firsts <- function(x) if (is.list(x=x)) sapply(x, `[`, 1) else x

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
  } else {
    xp <- hj
    yp <- vj
  }

  x <- x %||% xp
  y <- y %||% yp

  family <- family %||% element$family
  face <- face %||% element$face
  size <- size %||% element$size
  fontcheck <- fontcheck %||% element$fontcheck
  lineheight <- lineheight %||% element$lineheight

  if (length(label) > 1) {
    ## hack to support multiple labels.
    args <- list(label = label, x = x, y = y,
                 family = family, face = face, colour = colour, size = size,
                 hjust = hjust, vjust = vjust, angle = angle,
                 lineheight = lineheight, fontcheck = fontcheck,
                 stack = stack, default.units = default.units, ...)
    listy <- vapply(args, is.recursive, TRUE)
    nully <- vapply(args, is.null, TRUE)
    args[listy | nully] <- lapply(args[listy | nully], list)

    children <- do.call("mapply",
                        c(element=list(list(element)), args,
                          FUN=element_grob.element_text_with_symbols,
                          SIMPLIFY=FALSE))
    return(gTree(children=do.call("gList", children)))
  } else {
    if (is.list(family)) family <- family[[1]]
    if (is.list(face)) face <- face[[1]]
    if (is.list(size)) size <- size[[1]]
    if (is.list(lineheight)) lineheight <- lineheight[[1]]
  }

  gp <- gpar(fontsize = size, col = colour,
             fontfamily = family, fontface = face,
             lineheight = lineheight)
  element_gp <- gpar(fontsize = firsts(element$size), col = element$colour,
                     fontfamily = firsts(element$family),
                     fontface = firsts(element$face),
                     lineheight = firsts(element$lineheight))

  if (length(label) == 0) {
    textGrob(
      label, x, y, hjust = hj, vjust = vj,
      default.units = default.units,
      gp = modifyList(element_gp, gp),
      rot = angle, ...
      )
  } else {
    stringsplits <- renderable_split(label, family, face,
                                     size, lineheight, fontcheck)
    stack <- stack %||% element$stack
    offset.x <- 0
    offset.y <- 0

    if (dim(stringsplits)[1] > 1) {
      #plot with multiple fonts in the same frame.  what we should do
      #is pack a frame, and something something alignment...  the
      #angle, horisontal aliignment and vertical alignment should
      #become peroperties of how the frame is drawn...
      children <-
        mlply(stringsplits[c("labels", "family", "face", "size", "lineheight")],
              function(labels, family, face, size, lineheight, ...) {
                gp <- gpar(fontsize = size,
                           fontfamily = family, fontface = face,
                           lineheight = lineheight)
                textGrob(
                  labels, x=0, y=0,
                  hjust=ifelse(stack=='h', 0, hj),
                  vjust=ifelse(stack=='v', vj, 0),
                  default.units = default.units,
                  gp = modifyList(element_gp, gp), ...
                  )
              })

      #need to work in hj and vj in there somehow
      fr = frameGrob(vp=viewport(angle=angle, just=c(hj, vj)))
      for (ch in children) {
        fr = switch(stack
          , h=packGrob(fr, ch, side="right")
          , v=packGrob(fr, ch, side="bottom")
          , stop("'stack' must be 'h' or 'v'"))
      }
      fr
    } else {
      textGrob(
        stringsplits$labels, x, y, hjust = hj, vjust = vj,
        default.units = default.units,
        gp = modifyList(element_gp, gp),
        rot = angle, ...
        )
    }
  }
}


unicode_demo <- function() {
  curveleft <- "\u21BA"
  curveright <- "\u21BB"
  circleleft <- "\u27F2"
  circleright <- "\u27F3"

  combine_arrows <- function(x, combine) {
    first <- which(!is.na(x))[1]
    last <- length(x) + 1 - which(!is.na(rev(x)))[1]
    if (substr(x[[first]], 1, 1) == "-") {
      x[first] <- combine(circleleft, x[first])
      x[last] <- combine(circleright, x[last])
    } else {
      x[first] <- combine(curveright, x[first])
      x[last] <- combine(circleright, x[last])
    }
    x
  }

  append_arrows <- function(x) {
    combine_arrows(x, combine = function(x, y) paste(y, x, sep = ""))
  }

  prepend_arrows <- function(x) {
    combine_arrows(x, combine = function(x, y) paste(x, y, sep = " "))
  }

  newline_arrows <- function(x) {
    combine_arrows(x, combine = function(x, y) paste(y, x, sep = "\n"))
  }

  theme_set(theme_bw(12, "Myriad Pro"))
  theme_update(
    axis.title.x = element_text_with_symbols(
      size = list(c(12, 24)), family = list(c("Myriad Pro", "Apple Symbols")),
      lineheight = list(c(0.9, 1.5)), stack = "h")
    , axis.title.y = element_text_with_symbols(
        size = list(c(12, 24)), family = list(c("Myriad Pro", "Apple Symbols")),
        lineheight = list(c(0.9, 1.5)), angle = 0, stack = "v")
    , axis.text.x = element_text(size = 10, family = "Apple Symbols", vjust = 1)
    , axis.text.y = element_text(size = 10, family = "Apple Symbols", hjust = 1)
    ## , axis.text.x = element_text_with_symbols(
    ##     size = list(c(10, 20)), family = list(c("Myriad Pro", "Apple Symbols")),
    ##     angle = 0, stack = "v", vjust = 1)
    ## , axis.text.y = element_text_with_symbols(
    ##     size = list(c(10, 20)), family = list(c("Myriad Pro", "Apple Symbols")),
    ##     angle = 0, stack = "h", hjust = 1)
      )
  (qplot(x = -10:10, y = (-10:10)^2, geom = "line",
         xlab = "foo\u21BBbaz", ylab = "abla\u21BAedio")
   ## + scale_x_continuous(labels = append_arrows)
   ## + scale_y_continuous(labels = prepend_arrows)
   + scale_x_continuous(labels = newline_arrows)
   + scale_y_continuous(labels = prepend_arrows)
   )
}

#
custom_geom_demo <- function() {
  #here's an example that uses these.
  library(ptools)
  theme_set(theme_bw(12, "Myriad Pro"))
  theme_update(  panel.grid.major = element_blank()
               , panel.grid.minor = element_blank()
               )
  use_unicode <<- TRUE
  source("scales.R")
  library(binom)

  load("../modeling/data.Rdata")
  load("../modeling/slopeModel.RData")
  segment <- subset(data, exp_type == "numdensity" & subject %in% names(models))

  #this just illustrates the combinations of number and density.
  configurations <- unique(segment[c("target_spacing", "target_number_shown")])
  personalizations <- unique(segment[c("subject", "folded_displacement",
                                       "folded_direction_content")])

  ##For each personalization, collect the "segment" data to compare it against.
  alply(personalizations, 1,
        mkchain(  match_df(segment, ., on = names(.))
                , do.rename(folding = TRUE)
                , mkrates(splits = c(
                            "spacing", "content", "target_number_shown")))
        ) -> personalization.segment.data

  number_color_scale <-
    c(  list(aes(  color = factor(target_number_shown)) )
      , with_arg(  name = "Element\nnumber"
                 , palette = color_pal
                 , labels = prettyprint
                 , discrete_scale("colour", "manual")
                 ))
  ##
  Map(alply(personalizations, 1, identity), personalization.segment.data,
      f = function(row, data){
        chain(data,
              with(binom.confint(n_cw, n, method = "wilson", conf.level = 0.75)),
              cbind(data,.)) -> data
        rad <- 20/3
        (  ggplot(data)
         + aes(x = spacing)
         + aes(y = p, ymin = lower, ymax = upper)
         + aes(spacing = spacing, number = factor(target_number_shown))
         + continuous_scale("spacing", "spacing", identity,
                            , rescaler = function(x, from) {
                              rescale(x , from = from , to = from/rad)
                            }
                            , name = "Spacing")
         + identity_scale(discrete_scale("number", "identity",
                                         identity_pal(), name = "Element\nnumber"))
         #       + geom_errorbar()
         + geom_numdensity()
         #      + spacing_texture_scale
         + proportion_scale
         + number_color_scale
         + geom_line(aes(group = target_number_shown), show_guide = FALSE)
         ## + facet_wrap( ~ subject)
         + annotate(  "text", x = mean(range(data$spacing)), y = 0
                    , label = sprintf("%s, dx=%03f", toupper(row$subject)
                        , row$folded_direction_content)
                    , vjust = -0.5, size = 6
                    )
         )}
      #make a gtable plotting it both ways...
      ) -> segment.raw.plots

  print(segment.raw.plots[[5]])
}
