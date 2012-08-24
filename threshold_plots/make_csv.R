#import files from the 
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  source("db_functions.R")
  library(plyr)
  library(stringr)
})

renamings <- c(
trial_motion_process_radius   = "eccentricity",
folded_localDirectionContrast = "folded_direction_content",
abs_localDirectionContrast    = "abs_direction_content",
folded_displacement           = "folded_displacement",
abs_displacement              = "abs_displacement",
target_spacing                = "target_spacing",
trial_extra_nTargets          = "target_number",
subject                       = "subject",
abs_response                  = "drop",
folded_response               = "drop",
folded_response_with_carrier  = "folded_response_with_carrier",
abs_response_cw               = "abs_response_cw",
responseInWindow              = "drop"
)

main <- function(flist, dbfile, outfile) {
  files <- str_trim(readLines(flist))
  trials <- pull.from.sqlite(dbfile, data.frame(loaded_from = files))
  trials <- mutate(trials
                   , abs_response_cw = abs_response > 0
                   , folded_response_with_carrier= folded_response > 0)
  trials <- subset(trials
                   , responseInWindow==TRUE
                   , select=names(renamings[names(renamings) %in% names(trials)
                                           & renamings != "drop"]))
  trials <- rename(trials, renamings[names(renamings) %in% colnames(trials)])
  write.csv(trials, outfile, row.names=FALSE)
}

if ("--slave" %in% commandArgs()) {
  args <- commandArgs(trailingOnly=TRUE)
  do.call("main", as.list(args))
}
