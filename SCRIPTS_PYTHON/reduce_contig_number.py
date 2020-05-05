#! /usr/bin/env python3
# -*- coding : utf8 -*-
# AF, last modification : December, 12th 2018
#
# Script name: reduce_contig_number.py
#
# Description:  this script is taking a contig.fasta file and determine if the
#                size of the contig worth to be used in a megablast task, in the
#                next rule of the snakefile.
#
################################################################################

# ~ start script ~ #

#############
## Modules ##
#############

import argparse
import sys
import re

###############
## Functions ##
###############


#################
## Main Script ##
#################

## call for command line argument
parser = argparse.ArgumentParser()

parser.add_argument('-c', '--contigs', type = str, help = 'the contig file you want to reduce.')
parser.add_argument('--outfile', '-o', type = str, help = 'the output file name')
parser.add_argument('--size_limit', '-s', type = str, default = "1000", help = '''enter the minimun size of the contigs you consider to be "good"
 in order to be filtered in.''')

args = parser.parse_args()
 
## check the presence of all arguments

if not args.contigs\
or not args.outfile\
or not args.size_limit:
    parser.print_help()
    sys.exit("""At least one argument is missing. All arguments are requiered.
    Please, read the help section above.""")
    
## generate a dictionnary to record
dico = {}

# open and parse the file
with open(args.contigs, "r") as filin:
    for line in filin:
        # if the line is a sequence description
        # then a new key is generated for the dictionnary
        if ">" in line:
            name = line
            dico[name] = ""
        # if the line is a sequence, then it is added to the dictionnary as
        # value for the corresponding key
        else:
            dico[name] = dico[name] + line.rstrip("\n")
filin.close()

# N.B. for each line, the last character (matching the '\n') is not removed,
# just to facilitate the writing of the output file

limit_size = int(args.size_limit)

with open(args.outfile, "w") as filout:
    none = True
    for key in dico:
        size = len(dico[key])
        if size >= limit_size:
            none = False
            filout.write(key)
            filout.write(dico[key])
            filout.write("\n")
    if none == True:
        for key in dico:
            size = len(dico[key])
            if size >= limit_size/2:
                none = False
                filout.write(key)
                filout.write(dico[key])
                filout.write("\n")
       
# ~ end of script ~
