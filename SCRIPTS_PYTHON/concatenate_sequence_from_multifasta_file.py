#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# AF, December 18th 2018
# Name : turn_multi2single_fasta.py
#
# Description : This script take one multifasta file as input and return one
#               file containing only a single concatenated sequence in fasta
#               format.
#
################################################################################

# ~ start of script ~


#############
## Modules ##
#############

import sys
import argparse

###############
## Functions ##
###############



#################
## Main Script ##
#################

## call for command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', type = str, help = "the input file in fasta format you want to convert")
parser.add_argument('--output', type = str, help = "the output file you want to create")

args = parser.parse_args()

## check for the presence of all arguments
if not args.input or not args.output:
    parser.print_help()
    sys.exit("Error while entering argument, at least one of them is missing.")

## define the new sequence you want to create as a string
seq = ""

## read and record the data
with open(args.input, 'r') as filin:
    for line in filin:
        if ">" not in line and line != "\n":
            seq = seq + line[:-1]

print(len(seq))

## write the new sequence in the output file
with open(args.output, "w") as filout:
    filout.write(">concatenated_sequence_from_"+args.input+"\n")
    filout.write(seq)

# ~ end of script ~
