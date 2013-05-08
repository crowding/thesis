library(ggplot2)

#(beginnings of some functions for constructing subtle and wonderful colormaps.

#mess with color scales.....
library(colorspace)

#conversion of cone-space to CIE space
lms2xyz = function(lms) {
  if(!(is.matrix(lms))) {
    lms <- matrix(lms, nrow=1)
  }
  lms <- t(lms)
  xyz <- solve(
    matrix(c(
      .7328, .4296, -0.1624,
      -.7036, 1.6975, .0061,
      .0030, -.0136, .9834),
           nrow=3, byrow=TRUE),
    lms)
  XYZ(t(xyz))
}
col2RGB <- function(col) RGB(t(col2rgb(col)) / 255)
setY0 <- mkchain(as("XYZ"), coords, `[<-`( , 2, 0), XYZ)

#grayRamp <- as(col2RGB(c("gray10", "gray75")), "XYZ")

#some isoluminant directions in RGB space
redtint <- chain(c(1,-1,0), lms2xyz, setY0)
bluetint <- chain(c(0,0,1), lms2xyz, setY0)
greentint <- chain(c(-1,1,0), lms2xyz, setY0)
yellowtint <- chain(c(0,0,-1), lms2xyz, setY0)

pushRGB <- function(colors, direction, pointwise=FALSE) {
  #adjust some colors to a direction, pushing as far as display will allow.
  #do this in RGB_space because
  #RGB is a linear space (unlike sRGB) where [0,1] bounds the gamut
  #if given "pointwise"_will push each color separately; otherwise
  #pushes all colors by the given direction.
  coo <- coords(as(colors, "RGB"))
  direction <- coords(as(direction, "RGB"))
  lower.headroom1 <- coo
  lower.headroom2 <- lower.headroom1
  for (i in 1:dim(lower.headroom1)[2])
    lower.headroom2[,i] = coo[,i] + direction[,i]
  upper.headroom1 <- (1-coo)
  upper.headroom2 <- upper.headroom1
  for (i in 1:dim(upper.headroom1)[2])
    upper.headroom2[,i] = coo[,i] - direction[,i]
  howmuch <- function(a, b) a/(a-b)
  p_max_pos <- function(a,b) {
    #the maximum of the positive numbers...
    pmax(pmax(a, 0, na.rm=TRUE), pmax(b, 0, na.rm=TRUE), na.rm=TRUE)
  }
  scale <- p_max_pos(howmuch(lower.headroom1, lower.headroom2),
                     howmuch(upper.headroom1, upper.headroom2))
  scale <- Reduce(pmin,
                  lapply(1:dim(scale)[[2]], function(x)scale[,x,drop=FALSE]))
  if(!pointwise) {
    scale <- min(scale[is.finite(scale) & scale >= 0], na.rm=TRUE)
  }
  out <- coo
  for (i in 1:dim(coo)[2]) out[,i] <- coo[,i] + direction[,i] * scale
  RGB(pmin(pmax(out, 0, na.rm=TRUE), 1, na.rm=TRUE))
}

low <- 0.10; med <- 0.75; high <- 0.79
decision_colors <- c(hex(pushRGB(RGB(c(low, med), c(low, med), c(low, med)),
                                 bluetint)),
                     hex(RGB(high, high, high)),
                     hex(pushRGB(RGB(c(med, low), c(med, low), c(med, low)),
                                 redtint)))
## print(as(col2RGB(decision_colors), "XYZ"))
decision_values <- c(0, 0.5*((med-low)/(high-low)),
                     0.5,
                     1 - 0.5*((med-low)/(high-low)), 1)
## (ggplot(melt(volcano), aes(x=Var1, y=Var2, fill=value)) + geom_raster()
##  +     scale_fill_gradientn(colours=decision_colors,
##                             values=decision_values,
##                             breaks = seq(0.1, 0.9, 0.2),
##                             space="rgb"))

#an attempt at a high saturation colormap. Answer: need to think
#about saturation as a perceptual variable.

colorful.colors <- chain(
  c("blue", "red", "orange", "yellow"
    )
  , col2rgb, t, ./255, RGB, as("XYZ") #get color names into colorspace
  , coords, apply(., 2, `/`, .[,"Y"]), XYZ #equalize luminance
  , as("RGB"), coords, sRGB, coords, ./max(.), sRGB, hex #conflate RGB as sRGB because colorRamp is linear in sRGB...
  , colorRamp, .(seq(0,1,length=100)) , RGB
  )

colorful.tints <- chain(
  colorful.colors, as("XYZ"), coords
  , apply(., 2, '-', .[,"Y"]), XYZ, as("RGB")
  )

colorful.levels = chain(
  c("black", "white"), colorRamp(space="Lab"),
  .(seq(0,1,length=chain(colorful.colors,coords,dim, `[`(1)))),
  ./255, sRGB)

colorful.gradient <- pushRGB(colorful.levels, colorful.tints, pointwise=TRUE)

## library(reshape2)
## (ggplot(melt(volcano), aes(x=Var1, y=Var2, fill=value)) + geom_raster()
##  +     scale_fill_gradientn(colours=hex(colorful.gradient),
##                             space="rgb"))
