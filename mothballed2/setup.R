

##install packages we use
install.packages("ggplot2","cacheSweave","xtable",)

## install Yihui Xie's patch improving pgfSweave
## ref: 
local({
  setwd("/Applications/Lyx.app/contents/resources/")
  patch <- url("http://www.lyx.org/trac/raw-attachment/ticket/7555/sweave-patch.diff")
  on.exit(close(patch),add=TRUE)
  patcher <- base::pipe("patch -N -p1", 'w')
  on.exit(close(patcher),add=TRUE)
  writeLines(readLines(patch),patcher)
})
