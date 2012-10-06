mu <- 3
beta <- 1
sig <- pi/beta/sqrt(3)
sig_height <- dlogis(sig, 0, beta)
old <- theme_set(theme_bw() + theme(panel.grid=element_blank()))
# a bell curve, labeled with estimated delta-x, and sensitivity.
(ggplot(data.frame(x=c(-8, 8)), aes(x))
 + stat_function(fun = function(x) dlogis(x, mu, beta))
 + scale_x_continuous('Envelope motion(\\Delta x)', breaks = c(0, mu), labels=c("0", '\\bar{\\Delta x}'))
 + scale_y_continuous('\\Pr\\left(\\Delta x\\right)', breaks = c())
 + geom_segment(x=mu, xend=mu+sig, y=sig_height, yend=sig_height)
 + geom_text(label="\\sim\\frac{1}{\\beta}", x=mu+sig/2, y=sig_height, hjust=0.5, vjust=1)
 )

# The sensitivity function
cs <- 3
sens = function(s) (2 - 2./(1+exp(-cs/s)))
(ggplot(data.frame(x=c(0,8), y=c(0,1)), aes(x,y))
 + stat_function(fun=sens)
 + scale_y_continuous('Sensitivity ($\\beta$)', breaks=c(0, 1), labels=c(0, "$\\beta_0$"))
 + coord_cartesian(ylim=c(0,1), xlim=c(0,8))
 + scale_x_continuous('Spacing ($S$)', breaks=c(0,cs), labels=c(0, "$S_C$"))
 + geom_point(x=cs, y=0.5)
 )
