## model <- circle.models[["nj", "model"]]
library(vadr)
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

#%should be the same as above but different parameterization%
bind[term=logitSurroundTerm, fun=logitSurroundFun] <-
  term_and_fun(
    params = alist(width, strength),
    vars = alist(global, spacing, number),
    form = quote(
      global*strength *
      (2 * (exp(width/number/spacing * pi^2/3)
            / (1+exp(width/number/spacing * pi^2/3))) - 1)))

logitSurroundDemo <- exploreFun %<<% dots(
  logitSurroundFun, c(width=10, strength=1, global=10, spacing=1, number=10))

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
  if(interactive()) {
    dev.off()
    quartz(pointsize=8, dpi=120, width=10, height=10)
  }
  par(mar=c(2.6, 2.1, 2.6, 0.1), mgp=c(0,0,0), tcl=0.2)
  par(bg="transparent")
  s = split.screen(c(length(termDemos), 1), erase=FALSE)
  tryCatch({Map(s, termDemos, f=function(s, d) {
    screen(s)
    par(oma=c(0,0,0,0), mgp=c(1,0,0))
    plot.new()
    d()
  })})
  close.screen(s)
}
##showSurrounds()

# Here are the baseline model types (the spacing/displacement terms)
model.types <- chain(
  declare_data(
    list(type = "spacing", dconstrain=list(c()),
         fmla = (cbind(n_cw, n_ccw) ~
                 displacementTerm(spacing, displacement,
                                  start=c(cs=6, beta_dx=10))
                 + content + I(content * abs(content)))),
    ## list(type = "hemifield_number", dconstrain=list(c()),
    ##      fmla = (cbind(n_cw, n_ccw) ~
    ##              displacementTerm(number_as_spacing_hemi, displacement,
    ##                               start=c(cs=6, beta_dx=10))
    ##              + content + I(content * abs(content)))),
    ## list("softSpacing", list(c(dsharpness=3)),
    ##      (cbind(n_cw, n_ccw) ~
    ##       softMinDisplacementTerm(spacing, displacement,
    ##                              start=c(cs=5, beta_dx=14))
    ##       + content + I(content * abs(content)))
    ##      ),
    ## list("softNumberSpacing", list(c(csharpness=sqrt(3), fsharpness=sqrt(3))),
    ##      (cbind(n_cw, n_ccw) ~
    ##       softNumberDisplacementTerm(displacement, target_number_shown, spacing,
    ##                        start=c(cs=4, beta_dx=14, field=pi*20/3))
    ##       + content + I(content * abs(content)))
    ##      ),
    list("null", list(c()),
         cbind(n_cw, n_ccw) ~ 1)
    ),
  colwise(factor_in_order)()
  )

#here are all the different models we will consider.
model.additions <- chain(
  merge(data.frame(type=c("spacing", "hemifield_number"), stringsAsFactors=FALSE),
        declare_data(
          list(addition = "none", constrain = list(c()),
               fmla = . ~ .),
          list("global", list(c()),
               . ~ . + content_global),
          ## list("local", list(c()),
          ##      . ~ . + content_local),
          ## list("local_fullcircle", list(c()),
          ##      . ~ . + content_local + I(content_local*full_circle)),
          list("global_fullcircle", list(c()),
               . ~ . + content_global + I(content_global*full_circle)),
          list("hemi_fullcircle", list(c()),
               . ~ . + content_hemi + I(content_global*full_circle)),
          ## list("global_fullcircle_displacement", list(c()),
          ##      . ~ . + displacement + content_global + content_global:full_circle),
          ## list("local_global_fullcircle", list(c()),
          ##       . ~ . + content_local
          ##      + content_global + I(content_global*full_circle)),
          ## list("global_extent_fullcircle", list(c()),
          ##      . ~ . + content_global + extent + content_global:full_circle),
          ## list("local_extent_fullcircle", list(c()),
          ##      . ~ . + content_local + extent + content_local:full_circle),
          ## list("surround", list(c()),
          ##      . ~ . + surroundTerm(content_local, extent,
          ##                                start = c(width=10, strength=6))),
          list("surround_global", list(c()),
               . ~ . + content_global +
               surroundTerm(content_local, extent,
                            start = c(width=10, strength=6)))
          ## list("logit_surround", list(c()),
          ##      . ~ . + logitSurroundTerm(content_global, spacing, target_number_shown,
          ##                                start = c(width=pi*20/3, strength=6))),
          ## list("logit_surround_global", list(c()),
          ##      . ~ . + content_global +
          ##      logitSurroundTerm(content_global, spacing, target_number_shown,
          ##                   start = c(width=10, strength=6)))
          ## list("soft_min_surround", c(sharpness=3),
          ##      . ~ . + softMinSurroundTerm(content, target_number_shown, spacing,
          ##                          start = c(width=pi*20/3, strength=6)))
          ##
          ## ,
          ## list("soft_min_surround_sharp", list(c()),
          ##      . ~ . + softMinSurroundTerm(content, target_number_shown, spacing,
          ##                          start = c(width=pi*20/3, strength=6, sharpness=3)))
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
      ),
    . ),
  colwise(factor_in_order)()
  )
#with that setup, here are the fields which uniquely identify a "model type"
model.identifiers <- (names(model.types)
                      %v% names(model.additions)
                      %-% c("model", "fmla", "dconstrain",
                            "constrain", "start"))

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
              #scope, so to fucking keep the scope I guess I'll have
              #to capture the args
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
    function(model) {
      #print(model)
      quickdf(lapply(functions, function(y) {
        #print(y)
        y(model)
      }))
    })
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

main <- function(infile="slopeModel.RData",
                 infile2="density.modeling.RData",
                 infile3="motion_energy.csv",
                 outfile="combined.model.RData",
                 aplotfile="combined.model.pdf") {

  load(infile2)
  load(infile)
  plot.spacing <<- plot.spacing
  motion.energy <- add_energies(read.csv(infile3))

  #grumble
  for (i in ls()) assign(i, get(i), envir=globalenv())

  if (interactive()) {
    while(length(dev.list()) < 1) dev.new()
    plot.dev <- dev.list()[1]
  } else {
    cat("Writing plots to ", plotfile, "\n")
    cairo_pdf(aplotfile, onefile=TRUE)
    plot.dev <- dev.cur()
    on.exit(dev.off(plot.dev), add=TRUE)
  }

  dev.set(plot.dev)
  ##showSurrounds()

  model.df <<- model.df

  old <- gcinfo(TRUE)
  on.exit(gcinfo(old))

  #all the data we have. We are casting the field size as 1 i.e. "full" since
  #there are hopefully nonlinear terms that will fit it explicitly.
  splits <<- splits <-
    c(segment.config.vars, segment.experiment.vars, "side") %-% "exp_type"
  combined.data <<- combined.data <- chain(
    data,
    subset(exp_type %in% c("content", "spacing", "numdensity")),
    combined_model_recast,
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
    mutate(spec=interaction(addition, type, drop=TRUE)),
    subset(type=="null"),
    ddply(names(model.df) %-% "model", summarize,
          null.deviance.all=deviance.all,
          null.deviance.segment=deviance.segment,
          null.deviance.circle = {gc(); deviance.circle}),
    merge(fits, by=names(model.df) %-% "model"))
  if(interactive()) many.fits <<- many.fits

  # Summarize the quality of fit across models...
  # using a regression model to sort the incomplete-but-promising ones

  chain(
    many.fits,
    drop_recursive,
    subset(made.fit),
    model.frame=mutate(spec=interaction(addition, type, drop=TRUE)),
    lm(formula=bic ~ spec + subject - 1),
    coef,
    .[grep("^spec",names(.))],
    put(names(.), str_replace(names(.), "spec", "")),
    sort,
    data.frame(spec=names(.), model.score=.),
    merge(model.frame),
    .[(names(.)
       %^% c(names(model.types), names(model.additions)))
      %v% "model.score"],
    unique, arrange(model.score)
    ) -> model.quality

  #save here to be sure
  save(file=outfile, list=ls())

  chain(many.fits,
        drop_recursive,
        subset(made.fit),
        merge(model.quality),
        merge(data.frame(subset=c("bic", "deviance", "segment", "circle"))),
        ddply(c("type","addition"), mutate,
              grpmean = mean(bic[type=="spacing" & is.finite(bic)])),
        arrange(desc(model.score)),
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

  print(ggplot(plotdata)
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
  print.if.nonempty(best.model.type)

  # and let's plot its fits across the segment data

  warned.fits <-
    chain(many.fits, drop_recursive, as.data.frame,
          subset(made.fit & n.warnings > 0),
          melt(c("type", "addition","subject")),
          acast(type + addition ~ subject, margins=TRUE, length))
  print.if.nonempty(warned.fits)

  selected.model <- best.model.type
  selected.model <- data.frame(type="spacing", addition="surround_global")

  selected.coefs <- chain(many.fits, merge(selected.model), subset(made.fit),
                          .[c("subject", "coef")], mutate(coef=rbind%()%coef),
                          cbind(subject=.$subject, as.data.frame(.$coef)),
                          print)

  #selected.subject <- data.frame(subject="pbm")
  #theM <- chain(many.fits, merge(cbind(selected.model, selected.subject)), .$model[[1]])

  #and lets' plot its fits extent data.

  #plot the predictions of some chosen model versus the data...
  #I really want softMinSurround to work...

  #also as diagnostic, plot the predictions of the model on
  #a subject's full circle data.
  #selected.model <- data.frame(type="spacing", addition="")
  selected.fits <- subset(merge(many.fits, selected.model), made.fit)

  print(plot_segment_fit(selected.fits, prediction.dataset))

  #and can we profile the extents?
  print(make_extent_plots(selected.fits, combined.data))

  save(file=outfile, list=ls())
  for (i in ls()) assign(i, get(i), envir=globalenv())
}

make_extent_plots <- function(selected.fits, combined.data) {

  #we calculate a "global content sensitivity" i.e. the change in
  #response to an addition of a certaion abmount of blobal carrier
  #strength, as a function of the carrier sensitivity.

  #now let's calculate sensitivity
  #coefficients and standard errors for each bin
  measures <- chain(
    m=selected.fits, vapply(is.recursive, FALSE), .[!.],
    names, . %^% names(combined.data) %v% "model", m[.],
    adply(., 1, function(row) {
      bind[model=bind[model], ...=group] <- row
      chain(combined.data, merge(group), d=subset(content != 0),
            mutate(., model.fit=folding_predict(model, ., type="link", fold=TRUE)),
            d=mutate(., model.csens=content_sensitivity(., model)),
            glm(data=., formula=cbind(n_cw, n_ccw)
                ~ offset(model.fit) + content +
                cut_extent:content_global,
                family=model$family),
            folding_predict(type="terms", se.fit=TRUE, fold=TRUE),
            data.frame %()% .,
            cbind(d),
            mutate(model=0),
            ##.[c("extent", "model.csens", "cut_extent.content_global")],
            ddply("cut_extent", summarize,
                  extent=mean(extent),
                  model.csens=mean(model.csens),
                  adjustment = mean(fit.cut_extent.content_global/content_global),
                  content.sensitivity=mean(
                    model.csens +
                    fit.cut_extent.content_global/content_global),
                  content.sensitivity.se = mean(
                    se.fit.cut_extent.content_global/abs(content_global))),
            mutate(model=0),
            force)
    })
    , mutate(label=paste("Observer", toupper(subject)))
    )

  nominal.ecc <- 20/3
  fake.dataset <- chain(
    expand.grid(subject=unique(data$subject),
                exp_type="numdensity",
                displacement=0,
                content=1,
                side="left",
                spacing = seq(1, 2*pi*nominal.ecc/7, length.out=100),
                target_number_shown=7
                ),
    mutate(target_number_all = 2*pi*nominal.ecc/spacing),
    combined_model_recast)

  predictions <-
    chain(
      m=selected.fits, vapply(is.recursive, FALSE), .[!.],
      names, . %^% names(fake.dataset) %v% "model", m[.],
      adply(., 1, function(row) {
        bind[model=bind[model], ...=group] <- row
        chain(fake.dataset, d=merge(group),
              content_sensitivity(model),
              cbind(content.sensitivity=., d, model=0))
      }),
      mutate(label=paste("Observer", toupper(subject))))

  (ggplot(predictions)
   + aes(x=extent, y=content.sensitivity)
   + geom_line()
   + facet_wrap(~label)
   + geom_hline(y=0, color="gray50", linetype="11")
   + geom_pointrange(data=measures,
                     aes(ymin=content.sensitivity - content.sensitivity.se,
                         ymax = content.sensitivity + content.sensitivity.se))
   + labs(x="Stimulus extent", y="Carrier sensitivity")
   )

}

take_content_columns <- mkchain(
  x=.,
  colnames,
  str_match_matching("^fit\\..*content_.*"),
  .[!grepl(".*I\\.content.*", .)],
  .[!grepl("^fit\\.content$",.)],
  x[,.])

content_sensitivity <- function(dataset, model) {
  ## cc <- as.list(coef(model))
  chain(d=dataset,
          folding_predict(model, ., type="terms", fold=TRUE),
          list, data.frame(fit=.),
          take_content_columns,
          ## cbind(dataset[c("extent", "content_local", "content_global")]),
          ## dlply("extent", force), .[[5]],
          ## mutate(surr=surroundFun(cc$width, cc$strength, content_local, extent)
          ##             +cc$content_global*content_global)
          rowSums, ./dataset$content_global)
}



if(FALSE) {

  for (x in ls(model$original.env)) {
    assign(x, model$original.env[[x]], globalenv())
  }

  #Yuck! This is why we can't have nice things, look how it either
  #jumps off into the weeds!  I think that settles it, can't seriously
  #use gnm here, will have to switch to Stan or BUGS etc.

  ldply (c(1:20, 500), function(iter) {

      m2 <- gnm(
        formula=(
          cbind(n_cw, n_ccw) ~ displacementTerm(spacing, displacement,
                                                start = c(cs = 6, beta_dx = 10))
          + content + I(content * abs(content)) + content_global
          + surroundTerm(content_local, extent,
                         start = c(width = 10, strength = 6))),
        family = fit.family,
        data = fit.data,
        constrain = (c(names(constrain) %||% character(0)))
#        iterStart=2, iterMax=500
        )

    coef(m2)
  }) -> theIter

  print(theIter)
  plot_fit(m2, data=subset(ddd, exp_type != "numdensity"), style="bubble", splits=splits, fold=TRUE)

}

cut_extents <- mkchain(
  mutate(cut_extent=cut(extent, c(0, 10, 11.5, 12.4, 14.5, 20, Inf))),
  put(contrasts(.$cut_extent), "contr.helmert"))

combined_model_recast <- mkchain(
  mutate(side=factor(ifelse(exp_type %in% "numdensity", side, "all"))),
  recast_data(envelope.factor=2, carrier.factor=2),
  rename(c(content_global="content_hemi",
           number_shown_as_spacing="number_as_spacing_hemi")),
  recast_data(envelope.factor=1, carrier.factor=1),
  cut_extents)

plot_segment_fit <- function(selected.models,
                             pdata=prediction.dataset,
                             plot.dataset=segment.folded.spindled.mutilated) {

  selected.model.predictions <-
    predict_from_model_frame(
      selected.models,
      pdata,
      fold=TRUE, spindle=TRUE, collapse=TRUE)

  (plot.spacing %+% plot.dataset
        + density_prediction_layers(selected.model.predictions, connect="number")
        )
}

plot_combined_fit <- function(selected.model,
                           selected.subject = data.frame(),
                           pdata=prediction.dataset,
                           cdata=combined.data) {

  if(interactive()) while(length(dev.list()) < 3) dev.new()

  if(interactive()) dev.set(dev.list()[[1]])
  plot_segment_fit(selected.model, pdata)

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

  if(interactive()) dev.set(dev.list()[[2]])
  plot_fit(mmm, data=ddd, style="bubble", splits=splits, fold=TRUE)

  #show me also profiles of the model fit?
  if(interactive()) dev.set(dev.list()[[3]])

  local({
    #fucking model objects, got to be breaking lexical scope all the time
    e <- mmm$original.env
    for (n in ls(e)) {
      assign(n, e[[n]], globalenv())
    }
    prof <- profile(mmm) #fucking model objects
    plot(prof)
  })

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
    family=binomial(link=logit.2asym(g=0.025, lam=0.025)))

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
