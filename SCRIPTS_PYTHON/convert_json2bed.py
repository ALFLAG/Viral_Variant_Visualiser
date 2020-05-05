#! usr/bin/env python3
# -*- coding:utf-8 -*-
# AF, last modification August 9th 2018
# Name: convert_json2bed.py
#
# Description:  take a json file coordonates and return a bedfile.
#
################################################################################

# ~ start of script ~

#############
## Modules ##
#############

import argparse
import sys
import json

#################
## Main Script ##
#################

## entry variables definition 
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--injson", type = str,\
help = "the json file with genomique coordonate you want to convert in bed file")

parser.add_argument("-o","--outbed", type = str,\
help = "the bed file you want to create")

parser.add_argument("--header", type = str, help = "add header line or not")

arg = parser.parse_args()

## check for the rpesence of variables
if not arg.injson or not arg.outbed:
    parser.print_help()
    sys.exit("\nEither you forgot, or you misfilled the different options.\n\
Please, read the help section above.\n")

## read the entry file
filin = open(arg.injson, "r", encoding = "utf-8")
dico = json.load(filin)
filin.close()

## write the output file
output = arg.outbed

line_list = []
for key in dico:
    if key != "genomesize" and key != "virus_id":
        new_line= str(dico[key][0]) +"\t"+ str(dico[key][1])+"\t"+key+"\n"
        line_list.append(new_line)
    elif key == "virus_id":
        virus_id = dico[key]

## lambda function to sort the list of gene according to their start positions
line_list.sort(key = lambda line: int(line.split("\t")[1]))

with open(output,"w") as filout:
    if arg.header:
        header = "chr\tstart\tend\tgene\n"
        filout.write(header)
    else:
        for line in line_list:
            filout.write(virus_id[:-2]+"\t")
            filout.write(line)

# ~ end of script ~
