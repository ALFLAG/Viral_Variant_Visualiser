#!/usr/bin/env Rscript
#
# AF; last modification October 3rd 2018
#
# Description:  This script has been written in order to generate a graph
#               representing the viral genome and the possible variants inside
#               the viral population.
#               This script has been written in order to be executed as the last
#               step of a pipeline named Me.
#               It uses ggplot2 in order to draw this graph.
#
################################################################################
# ~ start of script ~

library(ggplot2) # import the ggplot2 library
args <- commandArgs(TRUE) # all arguments are character types
density = read.table(args[1], h=T, sep = "\t") # read the dataframe
threshold = as.numeric(args[2]) # define threshold
#density = read.table("test_res", h=T, sep = "\t") # read the dataframe
#threshold = 0.07
t1 = as.character(threshold) # define threshold as a character
t = paste("Threshold", t1, sep = ' ')

p = ggplot(density) # plot initialization
p1 = p + geom_point(aes(x = position, y = 0 , size = 1, color = gene_id), show.legend = F, shape = 15) # add the consensus points, and remove the legend "size"
p2 = p1 + geom_point(aes(x = position, y = as.numeric(variant_percent), color = gene_id)) # add the variant points
p3 = p2 + geom_line(aes(x=position, y=threshold)) # add the threshold line
p4 = p3 + labs(list(title = "Nucleotide Variation", x = "Base Position", y = "Variant Frequency")) # add axes and graph titles
p5 = p4 + theme(legend.position = "bottom") # modify the legend position
p6 = p5 + ylim(-0.06,1.2) # modify the scale
p7 = p6 + geom_text(aes(x = position, y = variant_percent + 0.01, label = add_alt, angle = 0)) # add alt nucleotide to the grapĥ
p8 = p7 + geom_text(aes(x = position, y = -0.03, label = add_ref, angle = 0)) # add ref nucleotide to the graph
p9 = p8 + geom_text(x = 0, y = threshold + 0.01, label = t) # add the threshold text
#p10 = p9 + geom_text(x = 500, y = -0.01, label = "consensus nucleotide")
#p10

#ggsave("graph_test.pdf", plot = last_plot(), width = 20, dpi = 600) # save the graph
ggsave(args[3], plot = last_plot(), width = 20, dpi = 600) # save the graph

# ~ end of script ~
