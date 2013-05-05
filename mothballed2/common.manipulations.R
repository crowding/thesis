append <- function(filename, what) {
  f <- file(filename, 'a')
  tryCatch(finally=close(f), {
    writeLines(what, f)
  })
}

print.to.pdf <- function(plot,
                         file=outfilename(paste(".",deparse(substitute(plot)),".pdf",sep="")), ...) {
  pdf(file=file, ...)
  tryCatch(finally=dev.off(), {
    append(prodfile, file)
    print(plot)
  })
}

print.and.pdf <- function(plot,
                          file=outfilename(paste(".",deparse(substitute(plot)),".pdf",sep="")),
                          ...) {
  print(plot)
  print.to.pdf(plot, file=file, ...)
}


outfilename <- function(what) {
  sub('\\..*$', what,prodfile)
}

clearfile <- function(filename, what) {
  f <- file(filename, 'w')
  tryCatch(finally=close(f), {
    seek(f, 0)
    truncate(f)
  })
}

mean.and.sd <- function(x)data.frame(  mean=numcolwise(mean)(x)
                                     , sd=numcolwise(sd)(x)
                                     , n=numcolwise(length)(x))

what.varies <- function(runs) {
  zip <- function(l) { ##like python zip
    lapply(apply(do.call(rbind, l), 2, list), unlist, recursive=FALSE)
  }

  whatChanges <- function(randomizer) {
    whichParamsChange <- which(lapply(lapply(
      zip(randomizer$values),
         unique), length) > 1)
    changingFields <- sapply(lapply(
                                    randomizer$subs[whichQuestParamsChange],
                                    unlist), function(x) do.call(paste, as.list(c("trial", x, sep=""))))
  }
  
  changing.conditions <-
    lapply(runs$beforeRun.trials.randomizers,
           function(x)whatChanges(x[[1]]))

  changing.conditions <- lapply(changing.conditions, function(x) {
    if (any(matches <- grep("ccluders", x))) {
      #A thing I do not like about R: this is the easiest way to drop elements from a vector.
      #Think about a fix for the sanity package
      x <- x[-matches]
      x <- c(x,"visibilityCondition")
    }
    #nTargets always varies... ignore it.
    if (any(matches <- grep("nTargets", x))) {
      x <- x[-matches]
    }
    x
  })

  #"conditions" determines how we organize the lattice plot, below.
  conditions <- Reduce(union, (changing.conditions))
  #not sure what the double dot is doing there but eh
  conditions <- setdiff(conditions, c("trial..extra.r", "trial.extra.r"))
}

common.manipulations <- function(envir=parent.env(environment())) {
  ##do some basic calculations we always need and
  ##fix up some changes in the trial data structures.
  ##Note that this does it by destructively modifying your workspace!
  ##There is a way around that, though, by enclosing this in a local() block.
  with(envir, {
    library(plyr)
    ##Define the visibility condition for occluder experiments, or
    ##leave as "full" elsewhere.

    ##Ugh, drop any trials columns that begin with "params."
    ##We can look in "runs" for those if we ever need them
    trials <- trials[grep("^params", colnames(trials), invert=TRUE)]

    ##We only look at the concentric trials.
    trials <- subset(trials, trial.version__.function == "ConcentricTrial")
    
    if (!is.null(trials$trial.occluders)) {
      mutate(trials, 
             trial.extra.visibilityCondition =
               factor(ifelse(is.na(trial.useOccluders)
                             , "full"
                             , ifelse(trial.useOccluders
                                      , ifelse(sapply(trial.occluders,function(x)if(is.list(x))x$startAngle else NA) > pi
                                               , "left"
                                               , "right")
                                      , "full"
                                      )
                             ),  levels=c("left", "right", "full"))
             ) -> trials
    } else {
      trials$visibilityCondition <- factor("full", levels=c("left", "right", "full"))
    }
    
    ##determine for each trial if the "correct" response was
    ##given. Somewhat irrelevant when we know about hte induced motion effect
    mutate(trials,
           correct=(result.response
                    ==( - trial.extra.globalDirection
                       - trial.extra.localDirection * !trial.extra.globalDirection))) -> trials
    
    ##mark the motion condition in each trial based on whether local and
    ##global run the same direction.
    trials <- transform(trials,
                        motionCondition=
                        factor(array(c("congruent","counterphase","incongruent",
                                       "local",NA,"local",
                                       "incongruent","counterphase","congruent"),c(3,3))
                               [cbind(trial.extra.localDirection+2,trial.extra.globalDirection+2)],
                               levels=c("local","counterphase","congruent","incongruent")))

    ##compute the response time as follows.

    ##Find the 'cw' or 'ccw' trigger for each trial.
    t.responseTimestamp <- triggers[grep("cw",triggers$name),c("trials.i", "knobTime")]
    colnames(t.responseTimestamp) <- c("i", "responseTimestamp")
    trials <- merge(trials, t.responseTimestamp, all.x=TRUE)

    ##idk wth this happens?
    colnames(triggers)[colnames(triggers) == "next"] <- "next."
    
    ##find the timestamp of motion onset (this is before the actual first
    ##element, as given by motionFirstElementDelay)
    t.motionBegun <- triggers[triggers$name == "ConcentricTrial/run/startMotion",
                              c("trials.i", "next.")]
    colnames(t.motionBegun) <- c("i", "motionBegun")
    trials <- merge(trials, t.motionBegun, all.x=TRUE)

    trials <- transform(trials, responseTime = (responseTimestamp - motionBegun - trial.motion.process.t))
    if (is.null(trials$trial.maxResponseLatency)) trials$trial.maxResponseLatency <- NA

    ##Determine whether the response time was in the allowable range for
    ##each trial
    trials <- transform(trials, minResponseTime = trial.awaitInput - trial.motion.process.t)
    trials <- transform(trials, maxResponseTime = ifelse(is.na(trial.maxResponseLatency),
                                  Inf, trial.maxResponseLatency)
                        + minResponseTime)
    trials <- transform(trials, responseInWindow = (responseTime<maxResponseTime) & (responseTime > minResponseTime))
    if ("trial.occluders" %in% colnames(trials)) {
      trials <- transform(trials, n.occluders = sapply(trial.occluders,
                          function(x) if(is.null(x) || identical(x, NA)) 0 else length(x)))
    }

    ##extract the "directioncal content" from the trials. some trials
    ##have all targets all same direciton, other trials have targets
    ##in opposite directions overlaid in pairs.

    if (!"trial.extra.useFlankers" %in% colnames(trials)) {
      trials$trial.extra.useFlankers <- NA
      trials$trial.extra.flankerAngle <- list(c())
    }
    
    direction.content <- with(trials,
                              mapply(  trial.motion.process.velocity
                                     , trial.motion.process.color
                                     , trial.extra.useFlankers
                                     , trial.extra.flankerAngle
                                     , FUN=function(vel, col, fl, fa) {
                                       
                                       ##Note that if the trials had
                                       ##flankers then the flankers
                                       ##were listed first...
                                       if ( fl && (length(fa) > 0)) {
                                         if (dim(col)[2] > 1) {
                                           col <- col[,-(1:4)]
                                         }
                                         vel <- vel[-(1:4)]
                                       }
                                       ##
                                       if ( diff(range(vel)) > 0 ) { #mixture of CW and CCW
                                         if (dim(col)[2] < 2) { #pure counterphase
                                           c(col[1,1], col[1,1])
                                         } else { #direcitonal content
                                           if (sign(vel[1]) > 0) {
                                             c(col[1,1], col[1,2])
                                           } else {
                                             c(col[1,2], col[1,1])
                                           }
                                         }
                                       } else { #pure cw/ccw
                                         if (sign(vel[1]) > 0) {
                                           c(col[1,1], 0)
                                         } else {
                                           c(0, col[1,1])
                                         }
                                       }
                                     }
                                     )
                              )

    trials <- mutate(  trials
                     , trial.extra.content.cw = direction.content[1,]
                     , trial.extra.content.ccw = direction.content[2,]
                     )
                         
    ##assign a "direction contrast" to trials that were done before
    ##direction contrast was a real parameter Express contrast, displacement and response in
    ##common absolute directions (positive=CCW, negative=CW)
    if (is.null(trials$trial.extra.directionContrast)) {
      trials$trial.extra.directionContrast <- NA
    }
    if (any(is.na(trials$trial.extra.directionContrast))) {
      trials <-
        mutate( trials
               , trial.extra.directionContrast=
                  ( (trial.extra.content.cw - trial.extra.content.ccw)
                    / (trial.extra.content.cw + trial.extra.content.ccw)
                    * sign(trial.extra.localDirection)
                   )
               )
    }

    ##compute absolute direction contrasts and displacements.
    trials <-
      mutate(  trials
             , abs.localDirectionContrast =
                   sign(trial.extra.localDirection) * trial.extra.directionContrast,
             , abs.displacement =
                 sign(trial.extra.globalDirection) * trial.extra.globalVScalar
                 * trial.extra.r * trial.motion.process.dt
             , abs.response = -result.response
             , folded.content.with =
                 ifelse(  abs.localDirectionContrast > 0
                        , trial.extra.content.cw, trial.extra.content.ccw)
             , folded.content.against =
                 ifelse(  abs.localDirectionContrast > 0
                        , trial.extra.content.ccw, trial.extra.content.cw)             
             , folded.localDirectionContrast = ((folded.content.with - folded.content.against)
                                                / (folded.content.with + folded.content.against))
             , folded.displacement =
                 ifelse(abs.localDirectionContrast != 0,
                        abs.displacement * sign(abs.localDirectionContrast),
                        abs.displacement * sign(abs.displacement)),
             , folded.response =
                 ifelse(abs.localDirectionContrast != 0,
                        abs.response * sign(abs.localDirectionContrast),
                        abs.response * sign(abs.displacement))
             )

    trials <-
      within(trials,
             target.spacing <- 2 * pi * trial.motion.process.radius / trial.extra.nTargets)
    trials <- subset(trials, !is.nan(responseTime))
    trials <- transform(trials, log.target.spacing=log(target.spacing))

    ##pull out the subject ID into each trial, as well as the actual experiment file name.
    if (!"subject" %in% colnames(runs)) {
      if ("beforeRun.params" %in% colnames(runs)) {
        runs$subject <- unlist(lapply(runs$beforeRun.params,
                                      function(x) ifelse(length(dim(x))==3, x[['subject',1,1]], NA)), recursive=FALSE)
      } else {
        runs$subject = runs$beforeRun.params.subject
      }
    }
    if (! "source.file" %in% colnames(runs)) {
      if ("beforeRun.params" %in% colnames(runs)) {
        runs$source.file <- sapply(runs$beforeRun.params, function(x)x[['logfile',1,1]][[1]])
      } else {
        runs$source.file <- runs$beforeRun.params.logfile
      }
    }
    trials <- merge(trials, runs[,c("i" ,"subject", "source.file")],
                    by.x='runs.i', by.y='i',
                    all.x=TRUE)

    #we have to remove any var we used that we aren't going to save....
    rm(t.responseTimestamp)
    rm(t.motionBegun)
    rm(direction.content)
    
    if ("frame.skips" %in% ls()) rm(frame.skips)
  })
}

common.manipulations.variations <- function(env) {
  with(env, {
      ##oops...
    trials <- within(trials, subject[subject=="dtdt"] <- "dt")
  
    ## work out what we vere varying during these experiments...
    condition.columns <- what.varies(runs)
    
    ##here we put special case handling of things.
    condition.exprs <- parse(text=condition.columns)

    trials <-
      within(trials,
             target.spacing <- 2 * pi * trial.motion.process.radius / trial.extra.nTargets)
    trials <- subset(trials, trial.version__.function == "ConcentricTrial")
    trials <- transform(trials, log.target.spacing=log(target.spacing), target.spacing=NULL)

    ##strip away data/columns we aren't planning on using, for the time being.
    trials <- subset(trials, select=c(condition.columns,
                               'trial.motion.process.radius', 'result.success',
                               'log.target.spacing', 'motionCondition', 'subject',
                               'responseInWindow', 'responseTime', 'minResponseTime',
                               'maxResponseTime', 'correct', 'runs.i', 'trial.extra.nTargets'))
    rm(runs)
  })
}

common.manipulations.displacement <- function(env=parent.env(environment())) {
  library(plyr)
  with(env, {
      mutate(trials,
             abs.localDirectionContrast =
             sign(trial.extra.localDirection) * trial.extra.directionContrast,
             abs.displacement =
             sign(trial.extra.globalDirection) * trial.extra.globalVScalar
             * trial.extra.r * trial.motion.process.dt,
             abs.response = -result.response
             ) -> trials

      mutate(trials,
             folded.localDirectionContrast = abs(abs.localDirectionContrast),
             folded.displacement =
             ifelse(abs.localDirectionContrast != 0,
                    abs.displacement * sign(abs.localDirectionContrast),
                    abs.displacement * sign(abs.displacement)),
             folded.response =
             ifelse(abs.localDirectionContrast != 0,
                    abs.response * sign(abs.localDirectionContrast),
                    abs.response * sign(abs.displacement))
             ) -> trials
  })
}

