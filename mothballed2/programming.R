###Useful R functions for programming, (encouraging terseness.) "Ptools"
library(arm)

index.data.frame <- function(df, row = 1:nrow(df), col=NULL, value) {
  ##A replacement for the inadequate
  ix <- index.data.frame.core(df, row, col)

  length(value) <- max(0, vapply(ix$indices, max, 0))
  mapply(function(i, c) {
    value[i] <<- df[i, c]
  }, indices, columns)
  output
}

index.data.frame.core <- function(df, row, col){
  if (is.null(col) && length(dim(row)) >= 2) {
    col <- row[,nrow(row)]
    row[,length(row)] <- NULL
  }

  row <- mkmatch(row, rownames(df))
  col <- mkmatch(col, colnames(df))
  indices <- tapply(1:length(col), row, identity)
  columns <- vapply(indices, function(i) col[i[1]], 0)
  list(indices=indices, columns=columns)
}

mkmatch <- function(index, names) {
 if (is.character(index) || is.factor(index)) {
   matched <- pmatch(index, names, duplicates.ok = TRUE)
   if (any(is.na(matched) & !is.na(index))) {
     stop("Unknown indices")
   }
 } else {
   matched <- index
 }
   matched
}

pad.missing.cells <- function(data, factors) {
  ##add a row containing NA for each combination of factors that is not represented.
  ##This is a workaround for a bug in current reshape / ggplot.
  chain( factors
       , combinations <- llply(., function(x)unique(data[,x,drop=FALSE]))
       , llply(nrow)
       , llply(seq)
       , do.call(expand.grid, .)
       , mapply(function(x, y) x[y,,drop=FALSE], combinations, .)
       , do.call(cbind,.)
       , merge(data, all.x=TRUE)
       )
}

mutate.where <- function(x, subset, ...) {
  ##a combination of mutate and subset.
  ##mutate those rows where subset evaluates to true, returning the entire modified data frame.
  e <- substitute(subset)
  r <- eval(e, x, parent.frame())
  if (!is.logical(r))
    stop("'subset' must evaluate to logical")
  r <- r & !is.na(r)
  cols <- as.list(substitute(list(...)))[-1]
  cols <- cols[names(cols) != ""]
  .data <- x[r,]
  for (col in names(cols)) {
    .data[[col]] <- eval(cols[[col]], .data, parent.frame())
  }
  for (col in names(cols)) {
    x[r,col] <- .data[,col]
  }
  x
}

keep.if <- function(x, expr, enclos=parent.frame()) {
  ##keep a subset if the expression evaluates to true. Use with ddply.
  if (eval(substitute(expr),x, enclos))
    x
  else
    x[c(),, drop=FALSE]
}

#fake out command line scripts for debugging...
run.command <- function(command) {
  blargs <- strsplit(command, ' ')[[1]]
  trace(commandArgs, at=3, tracer=substitute(args <-  c("--slave","--args", blargs[-1]), list(blargs=blargs)))
  on.exit(untrace(commandArgs))
  source(blargs[[1]])
}

#`subset<-(test,symbol=`.`)
##A shortcut for various assignments of the form:
##complicated.subset[with.long,names=TRUE] <- some.function.of(complicated.subset[with.long,names=TRUE],with.other.options=TRUE)
`%<-%` <- function(target, val) {
  val <- eval(substitute(substitute(val, list(.=quote(target)))))
  eval.parent(substitute(target <- val))
}

prefixing.assign <- function(prefix='', l=list(), env=parent.frame()) {
  for (n in names(l)) {
    assign(paste(prefix,n,sep=""),eval(substitute(l$n,list(n=n))),envir=env)
  }
}

almost.unique <- function(values, thresh = 0.0001) {
  values <- sort(values, na.last=TRUE)
  index <- pipe(values, diff, . > thresh, cumsum, c(0,.))
  tapply(values, index, mean)
}

cluster.to.unique <- function(values, thresh=0.0001) {
  pipe(values,
       .[ord <- order(.)],
       diff, .>thresh, cumsum, c(0,.),
       tapply(values[ord], ., function(x) {x[] <- mean(x); x}),
       unlist,
       .[inverse.permutation(ord)])
}

ammoc <- function(...) {
  ##evaluate all arguments in order, then return the first argument.
  ##Named because it is the left-to-right inverse of C's comma operator.
  ## the main use case for this is to remove a temp variable while keeping its value
  ## i.e. ammoc(.tmp, rm(.tmp))
  list(...)[[1]]
}

inverse.permutation <- function(perm) {
  ##if X is a vector expressing a permutation, for example the output
  ##of ORDER(), returns the inverse of that permutation.
  perm[perm] <- 1:length(perm)
  return(perm)
}

#this doesn't really handle "dots" arguments...?
curry <- function(..FUNCTION, ..., .currying.env=parent.frame()) {
  force(.currying.env)
  enclosing.env <- environment(..FUNCTION)
  
  defaults <- formals(..FUNCTION)
  curried.args <- as.list(match.call(..FUNCTION, sys.call(sys.parent())[-1], expand.dots=FALSE))[-1]
  
  for (n in names(curried.args))
    if(n == "...")
      stop("don't know how to curry varargs")
    else
      defaults[[n]] <- NULL

  out <- function(...) {
    calling.env <- parent.frame()
    eval.env <- new.env(parent=enclosing.env)
    calling.args <- as.list(match.call(expand.dots=FALSE))[-1]

    for (i in seq(len=length(defaults)))
      eval(substitute(delayedAssign(names(defaults)[[i]], .,
                    eval.env=eval.env,
                    assign.env=eval.env), list(.=defaults[[i]])))
    for (i in seq(len=length(curried.args)))
      eval(substitute(delayedAssign(names(curried.args)[[i]], ,
                    eval.env=.currying.env,
                    assign.env=eval.env), list(.=curried.args[[i]])))
    for (i in seq(len=length(calling.args)))
      if(names(calling.args)[[i]] == "...")
        stop("don't know how to curry varargs")
      else
        eval(substitute(delayedAssign(names(calling.args)[[i]], .,
                      eval.env=calling.env,
                      assign.env=eval.env), list(.=calling.args[[i]])))
    eval(body(..FUNCTION), eval.env)
  }
  formals(out) <- defaults
  out

  ##write the source attribute to be more human-readable
}

test.curry <- function() {
  require(testthat)
  adder <- function(x1, x2) x1+x2

  test_that("simple currying", {
    add3 <- curry(adder, 3)
    expact_that(add3(4) == 7)
  })

  test_that("curry named arguments", {
    divider <- function(x1, x2) x1/x2
    div2 <- curry(divider, x2=2)
    expect_that(divider(5) == 2.5)
  })

  test_that("unused argument detected at curry time", {
    expect_error(curry(adder, dsfargeg=3), throws_error("unused argument"))
  })

  test_that("unused argument detected at call time", {
    add3 <- curry(x1, 3)
    expect_error(curry(adder, dsfargeg=3), throws_error("unused argument"))    
  })
            
  test_that("curried argument evaluation is lazy", {
    3 <- x
    add3 <- curry(adder, x)
    x <- 4
    expect_that(add3(5) == 9)
  })

  test_that("curried arguments are evaluated in the currying context", {
    add3 <- NULL
    x <- 2
    local({
      3 <- x
      add3 <<- curry(adder, x)
    })
    x <- 4
    expect_that(add3(4) == 7)
  })

  test_that("curried arguments are evaluated anew each time", {
    evals <- 0
    evaluator <- function() {
      evals <<- evals + 1
      3
    }
    
    curry(adder, evaluator())
    expect_that(adder(3) == 6)
    expect_that(evals == 1)
    expect_that(adder(3) == 4)
    expect_that(evals == 2)
  })

  test_that("can curry quoted arguments", {
    sumX <- curry(with, expr=sum(x))
    expect_that(sumX(data.frame(x=c(1,2,3))) == 6)
  })

  test_that("can call with quoted arguments", {
    x123 <- curry(with, data=data.frame(x=c(1,2,3)))
    expect_that(with, data=data.frame(x=c(1,2,3)))
  })
  
  test_that("parent.frame accesses the final caller", {
  })
}

## As we can never remember how to use "substitute" on a non-quoted expression.
## I don't think this use of do.call is officially supported but it seems to work.
substitute.nq <- function(expr,...) {
  do.call(substitute, list(expr,...), envir=parent.frame())
}

## safer "data-like" version of pipe
chain <- function(...,dwim=TRUE) {
   arg.list <- as.list(match.call())
   arg.list[1] <- NULL
   if(dwim) arg.list[-1] <- pipe.dwim(arg.list[-1])
   sym <- gensym(envir=parent.frame())
   arg.list <- lapply(arg.list, substitute.nq, list(.=sym))
   arg.list <- lapply(arg.list, fn(substitute(a <- b, list(a=sym,b=.))))
   pipe.call <- do.call(call, c(list("{", enquote(substitute(. <- NULL, list(.=sym)))),
                                lapply(arg.list, enquote),
                                list(enquote(substitute(remove.returning(.), list(.=sym))))))
   eval.parent(pipe.call)
}

pipe <- chain

## a version that creates a function to be called later. useful in ddply etc.
mkchain <- function(..., dwim=TRUE) {
  arg.list <- as.list(match.call())
  arg.list[1] <- NULL
  if(dwim) arg.list <- pipe.dwim(arg.list)
  arg.list <- lapply(arg.list, substitute.nq, list(.=quote(..DATA)))
  arg.list <- lapply(arg.list, fn(substitute(a<-b, list(a=quote(..DATA),b=.))))
  pipe.call <- do.call(call, c(list("{"), lapply(arg.list, enquote), quote(quote(..DATA))))
  pipe.fn <- substitute(function(..DATA) ., list(.=pipe.call))
  out <- eval.parent(pipe.fn)
  attr(out, "source") <- NULL
  out
}

mkpipe <- mkchain

##this slightly ugly hack is required by the current way that pipe works.
remove.returning <- function(sym) {
  q <- substitute(sym)
  val <- force(sym)
  rm(list=as.character(q), envir=parent.frame())
  val
}

load.as.list <- function(...) {
  a = environment()
  load(envir=a, ...)
  as.list(a)
}

## a "macro-like" version of pipe. Note that in principle I could do
## the substitution once and then remember it, leading to a speed
## increase? Modifying the calling function wherever it is called from?
## replace.parent.call.with()? I think I can guarantee 
pipe.macro <- function(..., dwim=TRUE) {
  arg.list <- match.call()
  if(dwim) arg.list[-1:-2] <- pipe.dwim(arg.list[-1:-2])
  eval.parent(Reduce(function(x, y) substitute.nq(y, list(.=x)), arg.list[-1]))
}

## helper function that "guesses" the correct form of the arguments to
## pipe if they do not contain dots.
pipe.dwim <- function(arg.list) {
  lapply(arg.list,
         function(expr) if("." %in% all.names(expr)) {
           expr
         } else {
           switch(mode(expr),
                  name=call(as.character(expr),quote(.)), 
                  call=as.call(c(expr[[1]], quote(.), as.list(expr[-1]))),
                  stop("don't know how to pipe into a '", mode(expr), "'"))
         })
}

## "fn" is a slight shorthand for defining functions of one argument,
## useful in the arguments to higher-order functions.  Rather than
## function(.) .+1 we have fn(.+1) and otherwise the same.
fn <- function(expr) {
  ff <- eval.parent(substitute(function(.) e, list(e = substitute(expr))))
  attr(ff, "source") <- NULL
  ff
  ##need to do something with printing and use.source?
  ##mm <- match.call()
  ##mm$expr <- NULL
  ##mm[[1]] <- as.name("fn")
  ##attr(ff, "source") <- deparse(mm)
}

###below this fold, things get unreasonable.

## ##assign. overloading dot may produce some confusion with the plyr package.
## ##For example, 
## `.<-` <- function(x, expr) {
  
## }

## `.macro<-` <- function(x, expr) {
  
## }




             
## how to cache the results of macro substitution?
## modify the enclosing function in place?

## this comes from Gary Sabot and Thomas Lumley. via this post to s-news:
## http://www.biostat.wustl.edu/archives/html/s-news/2002-10/msg00064.html
gensym <- function(base=".v.",envir=parent.frame()){
  repeat{
    nm<-paste(base,paste(sample(letters,7,replace=TRUE),
                         collapse=""),sep=".")
    if (!exists(nm,where=envir))
      break
  }
  as.name(nm)
}

##future improvement: make gensym also scan all functions in the
##lexical environment, in preparation for a self-rewriting "pipe"

## here's a "defmacro" with locals
defmacro <-function(...,expr,local=NULL){
  expr<-substitute(expr)
  a<-substitute(list(...))[-1]
  nn<-names(a)
  if (is.null(nn)) nn<-rep("",length(a))
  for(i in seq(length=length(a))){
    if (nn[i]=="") {
      nn[i]<-paste(a[[i]])
      msg<-paste(a[[i]],"not supplied")
      a[[i]]<-substitute(stop(foo),list(foo=msg))
    }
  }
  names(a)<-nn
  a<-as.list(a)

  if(is.null(local)) {
     ff<-eval(substitute(function(){
        tmp<-substitute(body)
        eval(tmp,parent.frame())
     },list(body=expr)))
  } else {
     ff<-eval(substitute(function(){
        tmp<-substitute(body)
        locals<-lapply(local,gensym,envir=parent.frame())
        names(locals)<-local
        tmp<-do.call("substitute",list(tmp,locals))
        rval<-eval(tmp,parent.frame())
        rm(list=as.character(locals),envir=parent.frame())
        rval
     },list(body=expr)))
  }

  formals(ff)<-a
  mm<-match.call()
  mm$expr<-NULL
  mm[[1]]<-as.name("macro")
  attr(ff,"source")<-c(deparse(mm),deparse(expr))
  ff
}

## things this package ought to contain:

## a version of "substitute" that quotes all its arguments (thus is
## only a lexical transform, which can be replaced), and erases itself
## out of functions for speed. It might be used like this
##
## .(. = somefunction(.), .=complicated[[expr(argument)]])
##
## . <- function(...) {
##   args <- match.call()
##   eval.parent.rewriting.call(substitute.na(args[[2]], as.list(args[seq(3, length(args),1])))
## }

## the function that ties all my other functions' self erasure together.... will be a trick
## eval.parent.rewriting.call <- function(expr) {
##   eval.parent(expr)
## }

## an extended version of "gensym" that inspects the definitions of
## calling functions, and functions defining enclosing frames, for a
## list of all applicable symbols.
## gensym.lexical <- function(base=".v.",


## what are the pythagorean quadruples under 20?
test.pipe <- function() {
  ##quads <- subset(expand.grid(x=1:20, y=1:20, z=1:20, w=1:20), x^2+y^2+z^2 == w^2)

  ##each step's length is 5, 3, 9, 11, 7, 3
  path <- matrix(c(0, 0, 0,
                   0, 3, 4,
                   1, 1, 2,
                   0,-3,-6,
                   2, 3, 3,
                   0, 0,-3,
                   0, 0, 0), ncol=3, byrow=TRUE)
  basic.distance <- cumsum(sqrt(rowSums(apply(path, 2, diff)^2)))
  pipe.macro.distance <- pipe.macro(path, apply(.,2,diff), .^2, rowSums(.), sqrt(.), cumsum(.))
  pipe.macro.dwim.distance <- pipe.macro(path, apply(2,diff), .^2, rowSums, sqrt, cumsum)
  pipe.distance <- pipe(path, apply(.,2,diff), .^2, rowSums(.), sqrt(.), cumsum(.))
  pipe.dwim.distance <- pipe(path, apply(2,diff), .^2, rowSums, sqrt, cumsum)
}

## question: if you use eval.parent() to eval something in the parent,
## what do parent frames look like?
