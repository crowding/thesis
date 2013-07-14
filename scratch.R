#scratchpad for testing/debugging code I'm refining elsewhere

##########
 ## something about GNM is fucked and sensitive to different ways of contrasting...

  #a fit that ain't working for whatever reason
  model1 <- match_df(many.fits, cbind(bad.fits[3,], subject="nj", type="spacing"))$model[[1]]

#how you run through a mexican hat....

summary(gnm(cbind(n_cw, n_ccw)
            ~ displacementTerm(spacing, displacement, start=c(cs=6, beta_dx=10))
            + content
            + I(content*abs(content))
            + content_global
            + content_global:full_circle
            ,
            data=model1$data,
            family=binomial(link=logit.2asym(0.025, 0.025)),
            ))

 summary(gnm(cbind(n_cw, n_ccw)
            ~ displacementTerm(spacing, displacement, start=c(cs=6, beta_dx=10))
            + surroundTerm(content_local, extent, start=c(width0.5=2, strength=6))
            + content
            + I(content*abs(content))
            ,
            data=model1$data,
            family=binomial(link=logit.2asym(0.025, 0.025)),
            ))

 summary(gnm(cbind(n_cw, n_ccw)
            ~ displacementTerm(spacing, displacement, start=c(cs=6, beta_dx=10))
            + softMinSurroundTerm(content, extent, start=c(width=10, strength=6))
            + content
            + I(content*abs(content))
            ,
            data=model1$data,
            family=binomial(link=logit.2asym(0.025, 0.025)),
            ))

summary(gnm(cbind(n_cw, n_ccw)
            ~ displacementTerm(spacing, displacement, start=c(cs=6, beta_dx=10))
            + content
            + I(content*abs(content))
            + content_local
            + content_local:factor( extent)
            + content_local:full_circle
            ,
            data=model1$data,
            family=binomial(link=logit.2asym(0.025, 0.025)),
            ))

s <- sample(length(model1$residuals), 10)
model.matrix(~ side - 1, ndata,
             contrasts.arg=list(side="contr.treatment "))[s,]
options()
umodel <- update(model0, .~.+addition, data=ndata, family=model0$family)

model.matrix(model1)[s,]
model.matrix(model2)[s,]
ndata$addition[s,]


summarise2 <- function (.data, ...) {
  env <- list2env(.data, parent = parent.frame())
 
  cols <- eval(substitute(alist(...)))
   for (col in names(cols)) {
    env[[col]] <- eval(cols[[col]], env)
  }
  quickdf(mget(names(cols), env))
}

summarize3 <- macro(function(.data, ...) {
  template(quickdf(eval(quote({
          ...(Map(list(...), names(list(...)),
                  f = function(x,n)
                  template( .(as.name(n)) <- .(x))))
          mget(.(names(list(...))), environment())
        }),
        list2env(.(.data)))))
})

library(ptools)
summarize4 <- macro(function(.data, ...) {
  e <- list(...)
  n <- names(e)
  template(.(quickdf)(.(do.call)(
    `{`,
    .(c(Map(e=e, n=n, function(e, n) template( .(`<-`)(.(n),.(e)))),
        template(.(mget)(.(n), .(environment)())))),
    envir=.(list2env)(.(.data)))))
})

head(summarize4(baseball, g2 = g * 2, g4 = g2 * 2))

library(microbenchmark)

microbenchmark(
#  summarise(baseball, g2 = g * 2, g4 = g2 * 2),
  summarise2(baseball, g2 = g * 2, g4 = g2 * 2),
  summarize3(baseball, g2 = g * 2, g4 = g2 * 2),
  summarize4(baseball, g2 = g * 2, g4 = g2 * 2)
  )


makePromiseEnv <- function(expressions, parent=parent.frame()) {
   f <- function() environment()
   formals(f) <- as.pairlist(expressions)
   environment(f) <- parent
   f()
}

e <- makePromiseEnv(alist(a = {print("hello"); 4}, b = {print("again"); 6}))

makePromiseEnv2 <- function(expressions, envir=parent.frame()) {
  f <- eval(substitute(function() EE()), list(EE=environment))
  arglist <- expressions
  arglist[] <- list(quote(expr=))       #delete defaults, keep names
  formals(f) <- as.pairlist(arglist)
  do.call(f, expressions, envir=envir)
}

e <- makePromiseEnv2(alist(a=ls()))

d_env <- environment()
d_env$log <- function(x) s( (1/x) * d(x) )
d_env$`^` <- function(x, y) s( y*x^(y-1) * d(y) * d(x) )
d_env$`+` <- function(x, y) if(nargs()==1) s( d(x) ) else s( d(x)  + d(y) )
d_env$`*` <- function(x, y) s( y*d(x)  + x*d(y) )
d_env*`-` <- function(x, y) if(nargs()==1) s( -d(x) ) else -d(x) - d(y) 
d_env$`exp` <- function(x) s( exp(x) * d(x) )
d_env$`sin` <- function(x) s( cos(x) * d(x) )
d_env$`cos` <- function(x) s( -sin(x) * d(x) )

d <- function(a, b) {}

(ggplot() + aes(with(circle.models$model[[1]]$data, content/spacing),
                with(recast_data(circle.models$model[[1]]$data), content_global),
                color=circle.models$model[[1]]$data$target_number_shown)
 + geom_point())

chain(
  segment,
  extract_segment(fold=TRUE, spindle=TRUE, collapse=FALSE),
#  models$pbm$data,
#  refold(fold=TRUE),
  recast_data(carrier.factor=3, envelope.factor=3),
  subset(content==1),
  ggplot,
  (.+aes(x=spacing, y=number_shown_as_spacing,
         color=target_number_shown, group=target_number_shown)
   + geom_line()))

chain(
  quad.recast.models
  , match_df(data.frame(subject="pbm",
                        content.local=FALSE, envelope.local=TRUE))
  , .$model[[1]]$data
  , ggplot
  , +aes(content_local, content_global, color=target_number_shown)
  , +geom_point())

chain(
  new.data,
  ggplot
  , +aes(x=extent, y=abs(content_global), color=factor(spacing)
         , group=factor(spacing))
  , +geom_point(position="dodge")
  , +geom_line())


# show a contour plot of spacing versus number for a particular idea.


apply(HairEyeColor, c("Hair", "Eye"), sum)
melt(acast(melt(HairEyeColor), Hair ~ Eye, sum))

modelA <- subset(informed.models, subject=="pbm")$model[[1]]
modelB <- subset(adj.models, carrier.local==TRUE & envelope.local==TRUE & subject=="pbm")$model[[1]]

names(cbind(modelA$data, modelB$data))

(ggplot(modelB$data)
 + aes(spacing, content_global, color=target_number_shown)
 + facet_grid(content~.)
 + geom_point())

function(.) {
    . <- {
        browser()
        .
    }
    . <- mutate(., envelope.factor = envelope.factor, carrier.factor = carrier.factor)
    . <- mutate_when_missing(., eccentricity = 20/3, target_number_shown = round(2 * 
        pi * eccentricity/spacing), target_number_all = target_number_shown, 
        content_cw = (content + 1/4), content_ccw = (1 - content)/4, 
        side = factor("all", levels = c("all", "bottom", "left", 
            "right", "top")))
    . <- mutate(., content_local = content/spacing, content_global = content/2/pi/eccentricity * 
        abs.pmin(target_number_all/carrier.factor, target_number_shown), 
        full_circle = target_number_shown == target_number_all, 
        extent = spacing * target_number_shown, number_shown_as_spacing = 2 * 
            pi * eccentricity/envelope.factor/abs.pmin(target_number_shown, 
            target_number_all/envelope.factor))
    .
}


logliks <- mkchain(.$fits, put(.$fit, NULL), mutate(optimized=lapply(optimized, `[[`, "lp__")))
merge(logliks(a), logliks(b), "subject")


a$stan_predict <- stan_predict
environment(a$stan_predict) <- a
test <- predict(a, a$data[1,])
test <- predict(a, a$data[1:5,], summary=colwise_se_frame)
