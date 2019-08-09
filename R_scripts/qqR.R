#Function which makes qq-plots for counts data (expected p value calculated acording to null distribution)
#Input: p-values obtained from results
#Output: qq-plot 

#Load required library
library(ggplot2)

qq = function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
  
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )

  plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + geom_abline(intercept=0,slope=1, col="red")
  #plot=plot+opts(title=title)
  plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
  plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
  if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
  plot
}