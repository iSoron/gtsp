#! /usr/bin/env Rscript
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

d<-scan("stdin", quiet=TRUE)
cat(sprintf("%12s %12.3f\n", "min", min(d)))
cat(sprintf("%12s %12.3f\n", "max", max(d)))
cat(sprintf("%12s %12.3f\n", "median", median(d)))
cat(sprintf("%12s %12.3f\n", "geom mean", gm_mean(d)))
