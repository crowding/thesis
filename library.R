

##' @title Interpret command line arguments and invoke a main
##' function with them.
##'
##' The idea is that to write a command line utility with R, you just
##' write a main() function use Rscript as your #! interpreter,
##' and at the end of your R script call run_as_command.
##'
##' @param func Which function to invoke. Defaults to whatever "main"
##' function is defined in the calling scope.
##' @param arguments The command line arguments to parse. By default,
##' uses the ones provided on the R/Rscript command line
##' @param require_toplevel Only run if invoked from the top level, as
##' from Rscript.
##' @param require_noninteractive Only run if in a non-interactive R
##' session.
##' @return Nothing. Things printed will naturally go out stdout and
##' errors during execution will naturally result in a nonzero exit
##' code.
##' @author Peter Meilstrup
run_as_command <- function(  func=parent.frame()$main
                , arguments=commandArgs(trailingOnly=TRUE)
                , require_toplevel=TRUE, require_noninteractive
                , parse_args=TRUE, verify_args=TRUE) {
  if (     (length(sys.frames()) == 1 || !require_toplevel )
        && (!interactive() || !require_noninteractive) ) {
    do.call(func, as.list(arguments))
  }
}

