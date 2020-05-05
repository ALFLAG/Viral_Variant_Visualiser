#! /usr/bin/env python3
#  -*- coding:utf-8 -*-
#
# AF, created in February 1st 2019
#
# Name: check_duplicates_ReadName.py
#
# Description:  this script is used in the pipeline Me. It is used to check if
#               duplicates names of read exists inside the file.
#               Because the pipeline use programs in order to merge reads when
#               it is possible, and therefore, it is possible tat multiple reads
#               got the same name. And this is problematic. specially for Mira
#               assembler that reports an error.
#
################################################################################

# ~ start of script ~

#############
## Modules ##
#############

import gzip
import argparse
import sys

#################
## Main Script ##
#################

## call for command line arguments
parser = argparse.ArgumentParser(description = """
This program will take a read file as input and
will check if the read names are only present once. If duplicates are present,
the name of the read file will be modify by adding the suffix '_bis' or '_ter'.
""")

parser.add_argument("--read_file", type = str, help = "the name of the input read file you want to check")
parser.add_argument("--outfile", type = str, help = "the output file name")

args = parser.parse_args()

## check for the presence of arguments
if not args.read_file or not args.outfile:
    parser.print_help()
    sys.exit("""
At least one of the command line arguments is missing.
Please, read the help section above.
""")

## record the read sequence and modify the read name when duplicate
# opening, reading and turn the file into a list of lines

print("reading of the input file\n")

try:
    with open(args.read_file, "r") as filin:
        line_list = filin.readlines()
except:
    with gzip.open(args.read_file, "r") as filin:
        line_list = filin.readlines()

# define the length of the line_list and substract 3
# why? because we do not want to have an error list index out of range
size = len(line_list) - 3

print("record the sequence in a dictionnary\n")

# parse the line_list, record all lines, but only check the sequences name
dico = {} # define the object that will store the data, a dictionnary

for i in range(0, size, 4):
    # define what's a fastq sequence
    A = str(line_list[i]).rstrip("\n").lstrip("b'").lstrip('b"').rstrip("\n'").rstrip('\n"') # seq name, because in case of paired data, there is a second part in the name.
    B = str(line_list[i+1]).rstrip("\n").lstrip("b'").lstrip('b"').rstrip("\n'").rstrip('\n"') # sequence
    C = str(line_list[i+2]).rstrip("\n").lstrip("b'").lstrip('b"').rstrip("\n'").rstrip('\n"') # +
    D = str(line_list[i+3]).rstrip("\n").lstrip("b'").lstrip('b"').rstrip("\n'").rstrip('\n"') # quality
    Aa = A.split(" ")
    # record the sequence in dico, if the sequence name is not already in it
    # values of the dctionnary are tuples, because they use less memory space.
    if len(Aa) == 1:
        if Aa[0] not in dico:
            dico[Aa[0]] = (B, C, D)
        elif Aa[0]+"_bis" not in dico:
            dico[Aa[0]+"_bis"] = (B, C, D)
        else:
            dico[Aa[0]+"_ter"] = (B, C, D)
    else:
        if Aa[0] not in dico:
            dico[Aa[0]] = (Aa[1], B, C, D)
        elif Aa[0]+"_bis" not in dico:
            dico[Aa[0]+"_bis"] = (Aa[1], B, C, D)
        else:
            dico[Aa[0]+"_ter"] = (Aa[1], B, C, D)

print("write the output file\n")
## write the output file
with open(args.outfile, "w") as filout:
    if len(Aa) > 1:
        for key in dico:
            filout.write(key+" ")
            filout.write(dico[key][0]+"\n")
            filout.write(dico[key][1]+"\n")
            filout.write(dico[key][2]+"\n")
            filout.write(dico[key][3]+"\n")
    else:
        for key in dico:
            filout.write(key+"\n")                                                    
            filout.write(dico[key][0]+"\n")                                          
            filout.write(dico[key][1]+"\n")                                          
            filout.write(dico[key][2]+"\n")                                          

# ~ end of script ~
