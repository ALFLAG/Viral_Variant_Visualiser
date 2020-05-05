#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# AF, September 18th 2018
#
# Description:  Take one sequence and a viral acronym and rename the sequence
#
################################################################################

# ~ start of script ~

#############
## Modules ##
#############

import argparse
import re
import sys


#################
## Main Script ##
#################

# call for arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", type = str, help = "the consensus sequence to rename")
parser.add_argument("-o", type = str, help = "the sequence with the new name")
parser.add_argument("-a", type = str, help = "the accession number")
args = parser.parse_args()

# check for the presence of all arguments
if not args.i or not args.o or not args.a:
    parser.print_help()
    sys.exit("\nOne or more arguments are missing.\n"\
"Please, read the help section above.\n")

# define accession number variable
AN = str(args.a)[:-2]

# open files and rename the sequence
filin = open(args.i,"r")
with open(args.o,"w") as filout:
    for line in filin:
        if ">" in line: filout.write(">consensus-from_reference_" + AN + "\n")
        else: filout.write(line)

filin.close()

# ~ end of script ~
