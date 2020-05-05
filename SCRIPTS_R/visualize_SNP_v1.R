#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(TRUE) # les arguments sont de types characters
density = read.table(args[1], h=T, sep = "\t")


# define vectors coordonates, flag and annotation vector
coordx = c() ; coordy = c() ; ref = "" ; pool = c() ; flag = 0

# this loop allow to create vector coordonates, and etiquettes
for (i in 1:length(density$position)) {
  if (density$SNP[i] == 1) {
    x = density$position[i]
    a = density$ref[i]
    b = density$alt[i]
    ref = paste(as.character(x),a,b,sep="_")
    pool = c(pool,ref)
    coordx = c(coordx, density$position[i])
    if (flag == 0) {
      coordy = c(coordy,0.4)
      flag = 1
    } else if (flag == 1) {
        coordy = c(coordy,0.6)
        flag = 0
      }
    }
}

p = ggplot(density, aes(x = position, y = SNP, color = gene_id)) + geom_point(shape = 15, aes(size = size_point))
p1 = p + labs(list(title = "Nucleotide variation", x = "base position", y = "Presence of SNP"))
p2 = p1 + annotate(geom="text", x=coordx, y=coordy, label= pool, color="black", angle = 90)
p2

# record the final graph
ggsave(args[2], plot = p2, width=10)
