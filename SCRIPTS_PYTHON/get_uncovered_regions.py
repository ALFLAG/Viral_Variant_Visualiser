#! /usr/bin/env python3
#  -*- coding:utf-8 -*-
#
# AF, created February 10th 2019
# Name: get_uncovered_region.py
#
################################################################################
# ~ start of script ~

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

## call for command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("--input", type = str, help = "the file with coverage depth")
parser.add_argument("--output", type = str, help = "the TAB-delimited output file")

args = parser.parse_args()

## check for the presence of both arguments
if not args.input or not args.output:
    parser.print_help()
    sys.exit()

## define the regular expression matching only uncovered region
regex1 = '([A-Za-z0-9,-_\.;]+\t[0-9]+)\t[0-1]\n'
REGEX1 = re.compile(regex1)

regex2 = '(gi\|[0-9]+\|gb\|[A-Za-z0-9,-_\.;]+\|\t[0-9]+)\t[0-1]\n'
REGEX2 = re.compile(regex2)

## define a list to store output position matching 0 coverage
uncovered_position_list = []

## read the input file
with open(args.input, "r") as filin:
    for line in filin:
        if REGEX1.match(line):
            uncovered_position_list.append(REGEX1.search(line).group(1))
        elif REGEX2.match(line):
            uncovered_position_list.append(REGEX2.search(line).group(1))

## write the output file
with open(args.output, "w") as filout:
    for position in uncovered_position_list:
        filout.write(position + "\n")

# ~ end of script ~
################################################################################
