#! /usr/bin/env python3
#  -*- coding:utf-8 -*-
# AF - July, 2nd 2018
#
# Script name = recup_best_non_viral_ref_from_get_org.py
#
################################################################################

# ~ start script ~

#############
## Modules ##
#############

import re
import argparse
import sys

#################
## Main Script ##
#################

# variable importation
parser = argparse.ArgumentParser()
parser.add_argument("-i", help = "the get_org result file we want to parse\
 in order to determine the most present acession number")
parser.add_argument("-o", help = "choose the name of your output file")
arg = parser.parse_args()

# check for the presence of the input_file
if not arg.i or not arg.o:
    parser.print_help()
    sys.exit("An error occured while entering the options.\
Please, read the help section above")

else:
    # read the file, record the accession numbers, and their number of occurence
    dico = {}
    with open(arg.i) as filin:
        regex = re.compile("gb:([A-Z0-9]+)\.[0-9]?\t")
        line = filin.readline()
        stop = False
        while stop == False:
            if "virus" in line or "VIRUS" in line or "Virus" in line:
                matching_line = regex.search(line)
                if matching_line:
                    AN = matching_line.group(1)
                    stop = True
            line = filin.readline()
    # retrieve the five most represented AN and write them in a separate file
    with open(arg.o,"w") as filout:
        filout.write(AN)

# ~ end of script ~
