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
        AN_list = []
        exp_reg = "gb:([A-Z0-9_]+)\ttax_id:"
        regex = re.compile(exp_reg)

        for line in filin:
            if "virus" not in line  and "Virus" not in line and "VIRUS" not in line: # if not virus in line get AN and level up non_virus_flag
                matching_line = regex.search(line)                           
                AN = matching_line.group(1)                                  
                AN_list.append(AN)

#        virus_flag = 0 # counter to limit the number of non viral reference
#        non_virus_flag = 0
#        # we get the none viral reference until the second viral reference is reach.
#        # therefore, for each sample, the number of none viral sequence will be different
#        while virus_flag < 3:# and non_virus_flag < 50:
#            line = filin.readline()
#            print(virus_flag, non_virus_flag)
#            print(line)
#            if line != "": # if the end of the file has not been reached, then continue
#                if "virus" not in line  and "Virus" not in line and "VIRUS" not in line: # if not virus in line get AN and level up non_virus_flag
#                    matching_line = regex.search(line)
#                    AN = matching_line.group(1)
#                    AN_list.append(AN)
#                    non_virus_flag += 1
#                else: # if virus in line, level up virus_flag
#                    virus_flag +=1
#            else: # if the end of the file is reached, then automatically set virus_flag to 3 and non_virus_flag to 21, in order to stop the while loop
#                virus_flag = 3
#                non_virus_flag = 29

with open(arg.o,"w") as filout:
    for AN in AN_list:
        filout.write(AN+"\n")

# ~ end of script ~
