  
## @knitr density-setup
options(width = 70, useFancyQuotes = FALSE, digits = 4, lyx.graphics.center=TRUE)
library(Cairo)
library(ggplot2)
library(plyr)
library(grid)
library(ptools)
theme_set(theme_bw(12, "Apple Symbols"))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())
source("latexing.R")
use_unicode <- TRUE
source("scales.R")
source("library.R")
CairoFonts(regular="Apple Symbols:style=Regular", bold="Geneva:style=Bold",
           italic="Geneva:style=Italic", bolditalic="Geneva:style=BoldItalic")

geom_numdensity <- function(mapping=NULL, data=NULL,
                            stat="identity", position="identity", ...) {
  GeomNumdensity$new(mapping=mapping, data=data,
                     stat=stat, position=position, ...)
}

library(proto)
GeomNumdensity <- proto(ggplot2:::Geom, {
  objname <- "numdensity"
  default_stat <- function(.) StatIdentity
  default_aes <- function(.) aes(  colour = "black", size=5, weight=0.5, linetype=1
                                 , tick_in=1, tick_out = 3, alpha = NA, shape=16, fill=NA)
  #We use geom_point's guide for now, but we will need a custom guide...
  guide_geom <- function(.) "point"
  required_aes <- c("x", "y", "spacing", "number", "center")


  draw <- function(., data, scales, coordinates, ...) {
    munched <- coord_transform(coordinates, data, scales)

    ## draw the major glyphs out of minor glyphs.  the minor glyphs are
    ## IRL scaled, not data scaled, so there are separate 'offset.x'
    ## and 'offset.y (which have units of mm.)
    expanded <-
      adply(  munched, 1, function(row) {
        angles <- with(row, seq(  center - spacing/2 * (number - 1)
                                , by=spacing, length=number))
        data.frame(  angle = angles
                   , offset.x = row$size/2 * cos(angles)
                   , offset.y = row$size/2 * sin(angles))
      })
    tickmarked <- mutate(  expanded
                         , end.x = offset.x * (size+tick_out/2)/size
                         , end.y = offset.y * (size+tick_out/2)/size,
                         , offset.x = offset.x * (size-tick_in/2)/size
                         , offset.y = offset.y * (size-tick_in/2)/size)
    (ggplot2:::ggname)(.$my_name(), gTree(.$my_name(), children = gList(
        circleGrob(
          x = unit(munched$x, "native"), y = unit(munched$y, "native")
          , r = unit(munched$size/2, "mm")
          , gp = gpar(
              col = expanded$colour
              , fill = NA ###alpha(munched$fill, munched$alpha)
              , lwd = expanded$weight * (ggplot2:::.pt)
              , lty = expanded$linetype
              ))
      , segmentsGrob(
           x0 = unit(tickmarked$x, "native") + unit(tickmarked$offset.x, "mm")
           , y0 = unit(tickmarked$y, "native") + unit(tickmarked$offset.y, "mm")
           , x1 = unit(tickmarked$x, "native") + unit(tickmarked$end.x, "mm")
           , y1 = unit(tickmarked$y, "native") + unit(tickmarked$end.y, "mm")
           , gp = gpar(
               col = expanded$colour
               , lwd = tickmarked$weight * ggplot2:::.pt
               , lty = 1
               , lineend = "square"
               )))))

   }
})

##

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
       + ggplot2:::identity_scale(continuous_scale(c("spacing", "number"),
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

 
