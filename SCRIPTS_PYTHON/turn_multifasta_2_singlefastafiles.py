#!/usr/bin/env/ python3
# -*- coding: utf-8 -*-
# AF, last modification: September 4th 2018
#
# Name of script: turn_multifastafile_2_singlefastafiles.py
#
################################################################################

# ~ start script ~

#############
## Modules ##
#############


import sys
import argparse
import re

#################
## Main Script ##
#################

parser = argparse.ArgumentParser()

parser.add_argument("-i", type = str,\
help = "the multifasta file you want to split in multiple single fasta files")
parser.add_argument("-o", type = str,\
help = "the last outputfile, where a list of sequence id will be stored")
parser.add_argument("--dir", type = str,
help = "enter the path where to store sequences, add the slash at the end")

args = parser.parse_args()

if not args.i or not args.o or not args.dir:
	parser.print_help()
	sys.exit("""You forget to specify at least one of the options.
Please, read the help section above.""")

seq_name_list = [] # to record sequence id
with open(args.i, "r") as filin: # open the multifasta file
	regexseq = re.compile('[ATCGN]+\n')
	for line in filin:
		if ">" in line:
			line_list = []
			seq_name = line
			# on enregistre le nom de la sequence
		elif regexseq.search(line):
			line_list.append(line[:-1])
		elif line == "\n":
			seqname = ""
			seq = "".join(line_list)
			for letter in seq_name:
				if letter == " ":
					seqname = seqname + "_"
				elif letter == ":":
					seqname = seqname + "_"
				elif letter == ",":
					seqname = seqname + "_"
				else:
					seqname = seqname + letter
			seq_name_list.append(args.dir+"split_"+seqname[1:-1]+".fasta")
			with open(args.dir+"split_" + seqname[1:-1] + ".fasta", "w") as filout:
				filout.write(seqname+"\n")
				filout.write(seq)

with open(args.o, "w") as filout:
	filout.write(",".join(seq_name_list))

# ~ end of script ~
