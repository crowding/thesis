library(plyr)
cdatafile <- "calibrationData.RData"
if (file.exists(cdatafile)) {
  load(cdatafile)
} else {
  calib <- local({
    load("../unpacked/pbm-2012-03-23__15-04-11-ConcentricAdjustmentPeriodic.RData")
    with(runs, list(
                  resolution=beforeRun.params.cal.rect[[1]]
                 , size=beforeRun.params.cal.spacing[[1]] * beforeRun.params.cal.rect[[1]][c(3,4)]
                 , distance=beforeRun.params.cal.distance
                 , interval=beforeRun.params.cal.interval
                 , gamma=beforeRun.params.cal.gamma
                 , luminance=beforeRun.params.cal.calibration.stage2.readings[[1]]
                 , gamma=beforeRun.params.cal.gamma[[1]]
                 ))
  })
  save(calib, file=cdatafile)
}

spatial <- mutate(  subset(trials, trial.motion.process.radius-6.67 < 0.1)
                 , velocity=abs(sapply(trial.motion.process.velocity, `[`, 1))
                 , temporal.freq=abs(velocity/wavelength)
                 )


temporal.freq =
  with(trials, mapply(function(x,y) abs(x[1]/y), trial.motion.process.velocity,
                      trial.motion.process.wavelength))

#tell me about example trial data....

example.trial <- tail(subset(trials, trial.version...function == "ConcentricTrial")) & subject == "jt" & result.success),n=1)

save("calibrationData.R", data)

fieldNames <- with.db.connection(SQLite(), "../discrimination.sqlite", fn=function(conn) {
  #Parameters used for example trials
  dbListFields(conn, "trials")
select_from_db(
})
