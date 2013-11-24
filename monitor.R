#!/usr/bin/env Rscript
library(stringr)
library(plyr)
library(tools)
library(vadr)

getprocs <- function() {
  columns <- c("uid","user", "pid", "vsz", "tty", "command")
  system(paste("ps -ax -o ", paste0(columns, collapse=",")),
         intern=TRUE) -> processes
  pattern <- paste0(" *", str_dup("([^ ]+) +", length(columns)-1),
                   "(.*)", collapse="")
  data <- str_match(processes, pattern)
  data <- data[-1, -1]
  colnames(data) <- columns
  data <- as.data.frame(data, stringsAsFactors=FALSE)
  data <- mutate(data, uid=as.numeric(uid),
                 pid=as.numeric(pid), vsz=as.numeric(vsz))
  data
}

main <- function(thresh=6e6, interval=5) {
  thresh <- as.numeric(thresh)
  repeat {
    procs <- getprocs()
    m_ply(subset(procs, tty != "??" & vsz > thresh),
          function(pid, command, ...) {
            cat("KILLING PROCESS ", pid, "\n",
                command, "\n", "=====\n")
            pskill(pid)
          })
    Sys.sleep(interval)
  }
}

run_as_command()
