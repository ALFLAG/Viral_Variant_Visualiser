#! /usr/bin/env python3
#  -*- coding:utf-8 -*-
# AF - July, 2nd 2018
#
# Script name = count_AN_presence_from_getorg.py
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
		exp_reg = "[A-Z0-9]+\:[0-9]+\:[0-9]+\tgi\|[0-9\-]+\|[a-z]+\|([A-Z0-9_]+)"
		regex = re.compile(exp_reg)
		for line in filin:
			if "virus" not in line or "Virus" not in line or "VIRUS" not in line:
				matching_line = regex.search(line)
				AN = matching_line.group(1)
				if AN in dico:
					dico[AN] += 1
				else:
					dico[AN] = 1

	# turn the dictionnary into a list of lists
	AN_occ_list = []
	for key in dico:
		AN_occ_list.append([key,dico[key]])

	# then sort the list (M to m) according to the number of occurence
	sorted_list = sorted(AN_occ_list, key=operator.itemgetter(1), reverse=True)
	
	# retrieve the five most represented AN and write them in a separate file
	with open(arg.o,"w") as filout:
		for i in range(0,5):
			filout.write(sorted_list[i][0] + "\n")

# ~ end of script ~
