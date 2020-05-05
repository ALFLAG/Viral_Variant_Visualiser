#!/usr/bin/env Rscript

library(ggplot2) # import the ggplot2 library
args <- commandArgs(TRUE) # all arguments are character types

density = read.table(args[1], h=T, sep = "\t") # read the dataframe
threshold = as.numeric(args[2]) # define threshold

t = paste("Threshold", args[2], sep = ' ') # define the threshold text to add on the plot

p = ggplot(density) # plot initialization
p1 = p + geom_point(aes(x = position, y = 0 , size=1, color = gene_id), show.legend = F, shape = 15) # add the consensus points, and remove the legend "size"
p2 = p1 + geom_point(aes(x = position, y = variant_percent, color = gene_id)) # add the variant points
p3 = p2 + geom_line(aes(x=position, y=threshold)) # add the threshold line
p4 = p3 + labs(list(title = "Nucleotide Variation", x = "Base Position", y = "Variant Percent")) # add axes and graph titles
p5 = p4 + theme(legend.position = "bottom") # modify the legend position
p6 = p5 + ylim(0,0.5) # modify the scale
p7 = p6 + geom_text(aes(x = position, y = variant_percent + 0.01, label = add_text, angle = 0)) # add ref and alt nucleotide to the grapĥ
p8 = p7 + geom_text(x = 0, y = threshold + 0.01, label = t) # add the threshold text

ggsave(args[3], plot = last_plot(), width = 20, dpi = 600) # save the graph