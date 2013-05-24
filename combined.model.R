## model <- circle.models[["nj", "model"]]
library(ptools)
library(plyr)
library(ggplot2)
library(gnm)
library(psyphy)
source("library.R")
source("slopeModel.R")
source("density.modeling.R")
setup_theme()

use_folding <- TRUE

options(width=130)
infile <- "density.modeling.RData"
outfile <- "combined.model.RData"
plotfile <- "combined.model.pdf"

join <- function(...) suppressMessages(plyr::join(...))

 #We're going to try a nonlinear term for "surround size."
#Define a "surround" term and a surround function

#Relative summation strength as a function of width at narrow spacing.
#i.e. at narrow widths the surround is saturated, at wide it drops off
#(relatively). range 0 to 1
#Relative summation strength as a function of extent
#i.e. smaller fields stimulate a wider thing range 0 to 1
term_and_fun <- function(params, vars, form) {
    qe(list(
      term = nonlinearTerm(...(params))(...(vars))(.(form)),
      fun = (function(.a=...(quote_args %()% params), .b=...(quote_args %()% vars))
             .(form))))
}

bind[term=surroundTerm, fun=surroundFun] <-
  term_and_fun(
    params = alist(width, strength),
    vars = alist(density, extent),
    form = quote(
      strength * density
      * (2*exp(extent/sqrt(width^2) * pi^2/3)
         / (1+exp(extent/sqrt(width^2) * pi^2/3)) - 1)))

surroundDemo <- exploreFun %<<% dots(
  surroundFun, c(width=10, strength=1, density=1, extent=10))

bind[term=softMinSurroundTerm, fun=softMinSurroundFun] <-
  term_and_fun(
      params = alist(width, strength),
      vars = alist(carrier, number, spacing),
      form = quote(
        carrier*number*strength #strength density
        *(-log( #times soft-min of
            exp(-1 * 10) #unity
          + exp(-width/number/spacing * 10)) #< 1 if extent > field
         /10)))

softMinSurroundDemo <- exploreFun %<<% dots(
  softMinSurroundFun, c(width=10, strength=1, carrier=1, number=5, spacing=2))

bind[term=spacingDisplacementTerm, fun=spacingDisplacementFun] <-
  term_and_fun(
    params=alist(cs, beta_dx),
    vars=alist(spacing, displacement),
    form=quote(beta_dx * displacement * (2 - 2/(1+exp(-cs/spacing)))))

spacingDisplacementDemo <- exploreFun %<<% dots(
  spacingDisplacementFun, c(cs=2, beta_dx=1, spacing=2, displacement=1))

bind[term=softMinDisplacementTerm, fun=softMinDisplacementFun] <-
  term_and_fun(
    params=alist(beta_dx, cs, dsharpness),
    vars=alist(displacement, spacing),
    form=quote(
      -displacement * beta_dx #maximum sensitivity times didsplacement, times
      * log( #soft min of
        exp(-1*dsharpness^2) #sensitivity 1 (uncrowded or)
            + exp(-spacing/cs*dsharpness^2) #sensitivity limited by crowding
            )/dsharpness^2
      ))

softMinDisplacementDemo <- exploreFun %<<% dots(
  softMinDisplacementFun,
  c(beta_dx=1, cs=5, dsharpness=3, displacement=1, spacing=5))

#This is actually just a combination of soft-min displacement and number....
bind[term=softNumberDisplacementTerm, fun=softNumberDisplacementFun] <-
  term_and_fun(
    params=alist(beta_dx, cs, field, csharpness, fsharpness),
    vars=alist(displacement, number, spacing),
    form=quote(
      displacement * beta_dx  #maximum sensitivity times displacement,
      *-(log(                 #times soft min of
             exp(-1 * csharpness^2)          #sensitivity 1 (uncrowded), or
             + exp(-log(      #soft MAX of
                        exp(field/cs/number * fsharpness^2) # less than 1 if number > "critical"
                        + exp(spacing/cs * fsharpness^2) # Less than 1 if spacing > "critical"
                        ) / fsharpness^2 * csharpness^2)
             ))/csharpness^2))

softNumberDisplacementDemo <- exploreFun %<<% dots(
  softNumberDisplacementFun, c(beta_dx=1, cs=2, field=10, csharpness=2,
                                fsharpness=2, displacement=1, number=8, spacing=2)
  )

softNumberDisplacement2Fun <- softNumberDisplacementFun
softNumberDisplacement2Demo <- exploreFun %<<% dots(
  softNumberDisplacementFun, c(beta_dx=1, cs=2, field=10, csharpness=2,
                               fsharpness=2, displacement=1, number=100, spacing=2)
  )

termDemos <- chain(l=ls(),
                   str_match("(.*)Fun$"), .[,2], na.exclude,
                   paste0("Demo"), intersect(l),
                   mget(globalenv()))

showSurrounds <- function(new=FALSE) {
  dev.new(pointsize=6)
  par(mar=c(2.1, 2.1, 3.1, 0), mgp=c(1,0,0), tcl=0.2)
  par(bg="transparent")
  s = split.screen(c(length(termDemos), 1), erase=FALSE)
  tryCatch({Map(s, termDemos, f=function(s, d) {
    screen(s)
    plot.new()
    par(oma=c(0,0,1,0))
    d()
  })})
  close.screen(s)
}

#things to try adding to the spacingish model to account for the
#spacing experiment.
model.types <- chain(
  declare_data(
    list(type = "spacing", dconstrain=list(c()),
         fmla = (cbind(n_cw, n_ccw) ~
                 displacementTerm(spacing, displacement,
                                  start=c(cs=6, beta_dx=10))
                 + content + I(content * abs(content)))),
    list("softSpacing", list(c(dsharpness=3)),
         (cbind(n_cw, n_ccw) ~
          softMinDisplacementTerm(spacing, displacement,
                                 start=c(cs=5, beta_dx=14))
          + content + I(content * abs(content)))
         ),
    list("softNumberSpacing", list(c(csharpness=sqrt(3), fsharpness=sqrt(3))),
         (cbind(n_cw, n_ccw) ~
          softNumberDisplacementTerm(displacement, target_number_shown, spacing,
                           start=c(cs=4, beta_dx=14, field=pi*20/3))
          + content + I(content * abs(content)))
         ),
    list("null", list(c()),
         cbind(n_cw, n_ccw) ~ 1)
    ),
  colwise(factor_in_order)()
  )

#here are all the different models we will consider.
model.additions <- chain(
  merge(data.frame(type=c("softSpacing", "spacing", "softNumberSpacing"), stringsAsFactors=FALSE),
        declare_data(
          list(addition = "none", constrain = list(c()),
               fmla = . ~ .),
          list("global", list(c()),
               . ~ . + content_global),
          ## list("local", list(c()),
          ##      . ~ . + content_local),
          list("local_extent", list(c()),
               . ~ . + content_local:extent),
          ## list("local_fullcircle", list(c()),
          ##      . ~ . + content_local + content_local:full_circle),
          list("global_fullcircle", list(c()),
               . ~ . + content_global + content_global:full_circle),
          ## list("global_fullcircle_displacement", list(c()),
          ##      . ~ . + displacement + content_global + content_global:full_circle),
          ## list("local_global_fullcircle", list(c()),
          ##      . ~ . + content_global
          ##      + content_local + content_global:full_circle),
          #should actually be equivalent
          ## list("local_global_fullcircle_2", list(c()),
          ##      . ~ . + content_global
          ##      + content_local + content_local:full_circle),
          ## list("global_extent_fullcircle", list(c()),
          ##      . ~ . + content_global + extent + content_global:full_circle),
          ## list("local_extent_fullcircle", list(c()),
          ##      . ~ . + content_local + extent + content_local:full_circle),
          list("logit_surround", list(c()),
               . ~ . + surroundTerm(content_local, extent,
                                    start = c(width=pi*20/3, strength=6))),
          list("soft_min_surround", c(sharpness=3),
               . ~ . + softMinSurroundTerm(content, target_number_shown, spacing,
                                   start = c(width=pi*20/3, strength=6))),
          list("soft_min_surround_sharp", list(c()),
               . ~ . + softMinSurroundTerm(content, target_number_shown, spacing,
                                   start = c(width=pi*20/3, strength=6, sharpness=3)))
          #content:side not identifiable by GNM for some reason...
          ## list("global_fullcircle_contentside",
          ##      . ~ . + content_global
          ##      + content_global:full_circle + content:factor(side)),
          ## list("globalside",
          ##      . ~ . + content_global:side),
          ## list("global_fullcircle_contentside",
          ##      . ~ . + content_global + I(content_global*full_circle)
          ##      - content + content:factor(side)),
          )),
  rbind(
    declare_data(
      list(type="null", addition="none", constrain=list(c()), fmla=. ~ .)
      , .
      )),
  colwise(factor_in_order)()
  )
#with that setup, here are the fields which uniquely identify a "model type"
model.identifiers <- (names(model.types)
                      %v% names(model.additions)
                      %-% c("model", "fmla", "dconstrain",
                            "constrain", "start"))

#we toy with weights to tease the fitting for number style. Not really helping...
segment.weight <- 1
#circle.weight <- 2

fit_many_models <- function(model.df, combined.data,
                            starts = model.types,
                            scenarios = model.scenarios)
  {
  #for each given model....
  results <- adply(model.df, 1, function(row) {
    bind[model=bind[model], ...=model.group] <- row
    print(model.group)
    # skip if we don't have numdensity for this observer
    if (all(merge(combined.data, model.group)$full_circle))
      return (quickdf(list(model=list())))
    #for all types of starting points (spacing/number)...
    adply(starts, 1, function(row) {
      bind[fmla=bind[starting.fmla], dconstrain=bind[dconstrain],
           ...=start.group] <- row
      print(cbind(model.group, start.group))
      new.data <- merge(combined.data, model.group, names(model.group))
      new.data$side <- factor(new.data$side)
      # Then try the model update, using all the data.
      # Then refit under each of the scenarios
      adply(
        merge(model.additions, start.group)[names(model.additions)], 1,
        function(row) {
          bind[fmla=bind[fmla], constrain=bind[constrain], ...=addition.group] <- row
          print(cbind(model.group, start.group, addition.group))
          new.fmla <- update(starting.fmla, fmla)
          bind[new.model, warnings=warnings, error=error] <- captureWarnings({
            makeAModel <- function(fit.fmla, fit.data, fit.family, constrain, dconstrain) {
              #holy crap model objects are annoyingly obtuse about
              #scope, so fucking keep the scope I guess I'll have to capture the args
              m <- gnm(fit.fmla, data=fit.data, verbose=FALSE,
                     family=fit.family,
                     constrain=(c(names(constrain) %||% character(0),
                                  names(dconstrain) %||% character(0))),
                          constrainTo=c(constrain, dconstrain))
              if(!is.null(m)) m$original.env <- environment()
              m
            }
            environment(makeAModel) <- globalenv()
            makeAModel(new.fmla, new.data, model$family, constrain, dconstrain)
          })
          #return the model with some statistics....
          if(!is.null(new.model)) {
            statify(data.frame(
              #todo this could be better factored out as running
              #diagnostics on the completed dataset, or just DGAF
              coef=I(list(coef(new.model))),
              na.coef= sum(is.na(coef(new.model))),
              model=I(list(new.model)),
              warnings = I(list(warnings)),
              n.warnings = length(warnings),
              error=I(list(error)),
              is.error=!is.null(error),
              made.fit = TRUE
              ))
          }
          else data.frame(
              warnings = I(list(warnings)),
              n.warnings = length(warnings),
              error=I(list(error)),
              is.error=!is.null(error),
              made.fit = FALSE
            )
        })
    })
  })
  asisify(results)
}

#all of the stats we can compute
stat.funs <- functions.of(model)(
  deviance.all=deviance(model),
  deviance.circle=subset_deviance(model, full_circle),
  deviance.segment=subset_deviance(model, !full_circle),
  edf=length(residuals(model)) - model$df,
  aic.penalty=2 * (length(residuals(model)) - model$df.residual),
  bic.penalty=(log(length(residuals(model))) *
               (length(residuals(model)) - model$df.residual)),
  aic=extractAIC(model)[2],
  bic=extractAIC(model,
    k=log(length(residuals(model))))[2],
  hosmer.all=hosmerlem(model)$stat,
  hosmer.circle=hosmerlem(model,
    subset(model$data, full_circle))$stat,
  hosmer.segment=hosmerlem(model,
    subset(model$data, !full_circle))$stat
  )

statify <- function(model.frame, functions=stat.funs) {
  stats <- ldply(
    model.frame$model,
    function(model) quickdf(lapply(functions, function(y) y(model))))
  model.frame[names(stats)] <- stats
  model.frame
}

coef_frame <- function(model.frame) {
  ddply(model.frame, names(model.frame) %-% "model", function(x)
        quickdf(as.list(coef(x$model[[1]]))))
}

infile <- "density.modeling.RData"
infile2 <- "slopeModel.RData"
infile3 <- "motion_energy.csv"
outfile <- "combined.model.RData"
plotfile <- "combined.model.pdf"

main <- function(infile="density.modeling.RData",
                 infile2="slopeModel.RData",
                 infile3="motion_energy.csv",
                 outfile="combined.model.RData",
                 plotfile="combined.model.pdf") {

  old <- gcinfo(TRUE)
  on.exit(gcinfo(old))

  load(infile2)
  load(infile)

  motion.energy <- add_energies(read.csv(infile3))

  if (interactive()) {
    while(length(dev.list()) < 1) dev.new()
    plot.dev <- dev.list()[1]
  } else {
    cairo_pdf(plotfile, onefile=TRUE)
    plot.dev <- dev.cur()
    on.exit(dev.off(plot.dev), add=TRUE)
  }

  dev.set(plot.dev)
  showSurrounds()

  model.df <<- model.df

  #all the data we have. We are casting the field size as 1 i.e. "full" since
  #there are hopefulle momlineary terms that will fit it explicitly.
  splits <<- splits <-
    c(segment.config.vars, segment.experiment.vars, "side") %-% "exp_type"
  combined.data <<- combined.data <- chain(
    data,
    subset(exp_type %in% c("content", "spacing", "numdensity")),
    mutate(side=factor(ifelse(exp_type %in% "numdensity", side, "all"))),
    recast_data(envelope.factor=2, carrier.factor=2),
    rename(c(content_global="content_hemi",
              number_shown_as_spacing="number_as_spacing_hemi")),
    recast_data(envelope.factor=1),
    cut_extents,
    mkrates(splits))

  prediction.dataset <<- prediction.dataset <- chain(
    c(segment.config.vars, segment.experiment.vars,
      "bias", "full_circle", "side"),
    subset(combined.data, side!="all", select=.),
    unique,
    recast_data)

  #collect a bunch of fits to compare.
  {
    many.fits <<- many.fits <- fit_many_models(model.df, combined.data)
    gc()
  }

  #attach the null deviance to each model
  many.fits <- chain(
    fits=many.fits,
    subset(type=="null"),
    idata.frame,
    ddply(names(model.df) %-% "model", summarize,
          null.deviance.all=deviance.all,
          null.deviance.segment=deviance.segment,
          null.deviance.circle = {gc(); deviance.circle}),
    merge(fits, by=names(model.df) %-% "model"),
    idata.frame)
  if(interactive()) many.fits <<- many.fits

  #force some compaction
  save(file=outfile, list=ls())

  # Summarize the quality of fit across models...
  chain(many.fits,
        drop_recursive,
#        subset(na.coef == 0),
#       subset(n.warnings == 0),
        merge(data.frame(subset=c("bic", "deviance", "segment", "circle"))),
        ddply(c("type","addition"), mutate,
              grpmean = mean(bic[type=="spacing" & is.finite(bic)])),
        arrange(desc(grpmean)),
        mutate(addition=factor_in_order(addition)),
        mutate(., stat = cbind(
                    bic=bic,
                    all=deviance.all- null.deviance.all,
                    segment=deviance.segment - null.deviance.segment,
                    circle=deviance.circle - null.deviance.circle
                    )[cbind(1:nrow(.), subset)]),
        ddply(c("subset", "subject"),
              mutate, stat=stat - min(stat[is.finite(stat)]) ),
        mutate(stat=ifelse(is.finite(stat), stat, Inf)),
        mutate(stat=stat/quantile(stat[is.finite(stat)],0.9)),
#        ddply("subset", mutate, stat=stat/quantile(stat[is.finite(stat)], 0.9)),
        force
        ) -> plotdata

  (ggplot(plotdata)
   + aes(x=addition, y=stat, color=subject, # shape=subject,
         #linetype=type, shape=type,
         group=interaction(subject, type))
   + facet_grid(subset ~ type)
   + coord_cartesian(ylim=c(-0.2,1.5))
   + geom_line(alpha=0.5)
   + geom_point()
   + theme(axis.text.x=element_text(angle=-90, hjust=0, size=rel(0.5)),
           legend.position="bottom"))

  bind[incomplete.model.types, complete.model.types] <- chain(
    many.fits, drop_recursive,
    ddply(model.identifiers, summarize, complete=all(made.fit)),
    alply(data.frame(complete=c(FALSE, TRUE)), 1, merge, .))
  print.if.nonempty(incomplete.model.types)

  bind[unproblematic.model.types, problematic.model.types] <- chain(
    many.fits, drop_recursive,
    ddply(model.identifiers, summarize,
          problematic=!all(made.fit & n.warnings==0 & na.coef==0)),
    alply(data.frame(problematic=c(FALSE, TRUE)), 1, merge, .))
  print.if.nonempty(problematic.model.types)

  #okay, well, pick the best overall model type...
  best.model.type <- chain(
    many.fits, drop_recursive,
    merge(unproblematic.model.types),
    ddply(model.identifiers, summarize, total.bic=sum(bic)),
    arrange(total.bic), .[1,]
    )

  # and let's plot its fits across the segment data
  selected.model <- best.model.type

  warned.fits <-
    chain(many.fits, drop_recursive, as.data.frame,
          subset(made.fit & n.warnings > 0),
          melt(c("type", "addition","subject")),
#          acast(type ~ addition ~ subject, length, margins=TRUE),
          acast(type + addition ~ subject, margins=TRUE, length))

  selected.model <-
    data.frame(type="spacing", addition="soft_min_surround")
  selected.subject <- data.frame(subject=="jb")

  # and lets' plot its fits across coefficient data

  #plot the predictions of some chosen model versus the data...
  #I really want softMinSurround to work...
  selected.model <- data.frame(type="spacing",
                                     addition="soft_min_surround")
  selected.subject <- data.frame(subject="nj")

  theM <- chain(many.fits, merge(cbind(selected.model, selected.subject)), .$model[[1]])

  wat$terms
  wat$constrain
  coef(wat)
  wat$coefficients
  

  #also as diagnostic, plot the predictions of the model on
  #a subject's full circle data.
  selected.model <- data.frame(type="spacing", addition="")
  plot_combined_fit(selected.model, selected.subject, prediction.dataset, combined.data)

  #Wow!! That's terrible, in a good way.


  #and finally, construct a frameof all the subjects' coefficients in
  #the preferred model (or the model left over from circles...)

  save(file=outfile, list=ls())
}

if(FALSE) {  evalq(update(theM,
                    formula=(
                             cbind(n_cw, n_ccw)
                             ~ displacement + content
                             + softMinSurroundTerm( trace=TRUE,
                                                          content, target_number_shown, spacing,
                                                          start=c(width=pi*20/3, strength=6))),
                    constrain=c( "width"),
                    constrainTo=c(20),
                    trace=TRUE,
                    iterStart=2,
                    iterMax=400),
                   theM$original.env) -> wat
                   }


cut_extents <- function(dataset) {
  mutate(dataset,
         cut_extent=cut(extent, c(0, 10, 11.5, 12.4, 14.5, 20, Inf)),
         put(contrasts(.$cut_extent), "helmert"))
  ## chain(d=dataset, count("extent", "n_obs"), arrange(extent),
  ##       mutate(cut_extent=cut(extent, c(0, 10, 11.5, 12.4, 14.5, 20, 100))),
  ##       ddply("cut_extent", mutate, cutsize=sum(freq)))
}

run_as_command()

# let's make sure I'm computing "local" and "global" correctly.
FALSE && {
  chain(combined.data,
        subset(select=c("content_local", "content_global",
                 "target_number_all", "target_number_shown", "full_circle")),
        unique, ggplot,
        .+aes(content_local, content_global,
              color=factor(target_number_shown),
              shape=full_circle) +
        facet_wrap(~target_number_all) + geom_point())

  #this would be a better way to test the "number versus spacing" deal
  #but it's not giving any reasonable results.
  #the whole tanh thing is an attempt to bound the relative field size between zero and 1
  numberDisplacementTerm <-
    (nonlinearTerm(load, beta_dx, fieldsize)(number, displacement, fullcircle)
     ( (2 - 2/(1+exp(-number/load
                     * ((tanh(fieldsize)*1)/2*(1-fullcircle) + (0+fullcircle)) )))
      * beta_dx * displacement))
  test.data <- subset(combined.data, subject=="pbm")
  number.factor.model <- gnm(
    (cbind(n_cw, n_ccw) ~
     numberDisplacementTerm(target_number_shown, displacement, full_circle,
                            start=c(cs=4, beta_dx=14, fieldsize=0.5))
     + content_global:full_circle
     + content + I(content * abs(content))
     + bias - 1),
    data=test.data,
    family=binomial(link=logit.2asym(g=0.025, lam=0.025))
    )

  number_factor_model <- gnm(
    (cbind(n_cw, n_ccw) ~
     spacingDisplacementTerm(number_shown_as_spacing, displacement,
                             start=c(cs=4, beta_dx=14, fieldsize=0.5))
     + content_global
     + content_local
     + content + I(content * abs(content))
     + bias - 1),
    data=test.data,
    family=binomial(link=logit.2asym(g=0.025, lam=0.025))
    )
}

plot_combined_fit <- function(selected.model,
                           selected.subject = data.frame(),
                           pdata=prediction.dataset,
                           cdata=combined.data) {
  selected.model.predictions <- predict_from_model_frame(
    merge(many.fits, selected.model)
    , newdata=pdata
    , fold=TRUE, spindle=TRUE, collapse=TRUE)
  #plot.spacing et al come from density.modeling.RData

  selected.model.predictions <-
    predict_from_model_frame(
      merge(many.fits,
            cbind(selected.model)),
      pdata,
      fold=TRUE, spindle=TRUE, collapse=TRUE)

  while(length(dev.list()) < 3) dev.new()

  dev.set(dev.list()[[1]])
  print((plot.spacing %+% segment.folded.spindled.mutilated
         + prediction_layers(selected.model.predictions, connect="number")
         ))

  splits <<- splits <- splits %-% "exp_type"
  mmm <- chain(selected.model,
               cbind(selected.subject),
               merge(many.fits),
               .$model[[1]])
  ddd <- chain(selected.model,
               cbind(selected.subject),
               cbind(data.frame(full_circle=TRUE)),
               merge(combined.data),
               mutate(exp_type="Full circle"))

  dev.set(dev.list()[[2]])
  plot_fit(mmm, data=ddd, style="bubble", splits=splits, fold=TRUE)

  #show me also profiles of the model fit?
  dev.set(dev.list()[[3]])

  local({
    #fucking model objects, got to be breaking lexical scope all the time
    assign("new.data", model$data, globalenv())
    assign("model", mmm, globalenv())
    prof <- profile(model) #fucking model objects
    plot(prof)
  })

}
