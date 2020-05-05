#!/usr/bin/env Rscript

library(ggplot2) # import the ggplot2 library
args <- commandArgs(TRUE) # all arguments are character types

coverage = read.table(args[1], sep = "\t", h=F) # read the dataframe

#m = paste("Coverage - median_coverage =", med, sep = " ")
maxi = max(coverage$V3)
mini = min(coverage$V3)
med = median(coverage$V3)
m = paste("Coverage - median_coverage =", med, "[", mini, ":", maxi, "]", sep = " ")

p = ggplot(coverage, aes(x = coverage$V2, y = coverage$V3)) + geom_line(color = "black", size = 0.5) # trace the curve
p1 = p + labs(list(title = m, x = "Base Position", y = "Number of Reads")) # add titles to the graph and axes
#p2 = p1 + geom_text(x = 1000, y = -10, label = m) + ylim(-10, maxi + 10) # add a tag to the plot giving the median_coverage value
#p3 = p2 + geom_line( aes( x = coverage$V2, y = med, color = 2 ) ) # to trace the median line on the plot


ggsave(args[2], plot = last_plot(), width = 20, dpi = 600) # save the plot
