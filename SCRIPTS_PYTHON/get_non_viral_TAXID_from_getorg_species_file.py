#! /usr/bin/env python3
#  -*- coding:utf-8 -*-
# AF - September, 10th 2018
#
# Script name = get_non_viral_AN_from_getorg_species_file.py
#
################################################################################
 
# ~ start script ~ 
 
#############
## Modules ##
#############

import re
import argparse
import sys
import operator

#################
## Main Script ##
#################

# variable importation
parser = argparse.ArgumentParser()
parser.add_argument("-i", help = "the get_org species file we want to parse\
 in order to retrieve AN of the most present non viral seq")
parser.add_argument("-o", help = "choose the name of your output file")
arg = parser.parse_args()

# check for the presence of the input_file
if not arg.i or not arg.o:
    parser.print_help()
    sys.exit("An error occured while entering the options.\
Please, read the help section above")

else: # read the file, record the accession numbers
    with open(arg.i) as filin:
        TAXID_list = []
        exp_reg = "\ttax_id:([0-9]+)\t"
        regex = re.compile(exp_reg)
        a = 0 # counter to limit the number of non viral reference
        while a < 5: # limited to five references
            line = filin.readline()
            if "virus" not in line  and "Virus" not in line and "VIRUS" not in line:
                matching_line = regex.search(line)
                TAXID = matching_line.group(1)
                if TAXID not in TAXID_list:
                    TAXID_list.append(TAXID)
                a += 1

with open(arg.o,"w") as filout:
    for TAXID in TAXID_list:
        filout.write(TAXID+"\n")

# ~ end of script ~
