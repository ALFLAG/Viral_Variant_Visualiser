#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# AF, January 22nd 2019
#
# Name: FilterIn_for_mira.python3
# description:  this script is designed for the selection of raw reads that have
#               not be filtered out after the "remove host" step.
#               The reads that will be selected here will be used for mira
#               assembling step.
#
################################################################################

# ~ start of script ~

#############
## Modules ##
#############

import sys
import argparse
import gzip

###############
## Functions ##
###############

def record_raw_reads(infile):
    dico = {}
    with gzip.open(infile, 'r') as filin:
        raw_list = filin.readlines()
        for i in range(0,len(raw_list)-3,4):
            key = str(raw_list[i]).rstrip("\\n'").lstrip("b'").split(" ")
            A = str(raw_list[i+1]).rstrip("\\n'").lstrip("b'")
            B = str(raw_list[i+2]).rstrip("\\n'").lstrip("b'")
            C = str(raw_list[i+3]).rstrip("\\n'").lstrip("b'")
            C_prime = str(C.rstrip('\\n"').lstrip('"')) # necessary, do not know why, but some binary string are not b'' but b"". Do not understand.

            dico[key[0]] = (A,B,C_prime)

    return dico

#################
## Main script ##
#################

## call for the command line argument
parser = argparse.ArgumentParser()
parser.add_argument("--filt", type = str, help = """
Please, enter the name of the file with the cleaned and filtered viral reads in a fastq format""")
parser.add_argument("--data", type = str, help = """
Please, select the type of data you have, single-ends (se) or pair-end (pe)""")
parser.add_argument("--raw1", type = str, help = """
Please, enter the name of the file with the raw reads1 in a fastq.gz format.
If you use single-ends data, use only this argument.
If you use pair-ends data, use both raw arguments.""")
parser.add_argument("--raw2", type = str, help = """
Please, enter the name of the file with the raw reads2 in a fastq.gz format.
If you use single-ends data, do not use this argument.
And if you use pair-ends data, use both raw arguments.""")
parser.add_argument("--out", type = str, help = """
Please, enter the name of the output file""")

args = parser.parse_args()

## check for the presence of all arguments

if not args.filt or not args.out or not args.data:
    parser.print_help()
    sys.exit("""\nAt least one of the three following arguments is missing (data, out, or filt).
    Please, read the help section above and try again.\n""")

if args.data == "se":
    if not args.raw1:
        parser.print_help()
        sys.exit("""
You want to use se data, but you do not specified the raw1 argument. Please, correct this.
""")
    elif args.raw2:
        parser.print_help()
        sys.exit("""
You want to use single-ends data, therefore you should not use the raw2 argument, only use the raw1 argument.
Please, read the help section above.
""")
elif args.data == "pe":
    if not args.raw1 or not args.raw2:
        parser.print_help()
        sys.exit("""
You want to use pair-ends data, but at least one of the raw1 or raw2 argument is missing.
Please, read the help secion above.
""")
    
# reading the raw_reads files
# and record it as a dictionnary
if args.data == "se":
    print("\nsingle-ends data will be treated\n")
    print("raw file")
    dico = record_raw_reads(args.raw1)
 
elif args.data == "pe":
    print("\npair-ends data will be treated\n")
    print("raw file1")
    dico1 = record_raw_reads(args.raw1)

    print("raw file2")
    dico2 = record_raw_reads(args.raw2)


# use the filtered_read file
# we only record the names of the different sequences
print("filtered_file")
filt_name_list = []
with open(args.filt, 'r') as filin:
    filt_list = filin.readlines()
    for i in range(0, len(filt_list)-3, 4):
        name = filt_list[i].rstrip("\n").split(" ") # separate the read name and the extension 1:N:0:1, that inform if the read comes from the file 1 or the file 2
        filt_name_list.append(name[0])

print("comparing")
# search the names of the filtered read inside the raw_list
# if found, then filtered in the raw reads as output
with open(args.out,'w') as filout:
    if args.data == "se":
        for key in filt_name_list:
            fastq_seq = key +"\n"+ dico[key][0] +"\n"+ dico[key][1] +"\n"+ dico[key][2] +"\n"
            filout.write(fastq_seq)
    else:
        for key in filt_name_list:
            fastq_seq1 = key +"\n"+ dico1[key][0] +"\n"+ dico1[key][1] +"\n"+ dico1[key][2] +"\n"
            fastq_seq2 = key +"\n"+ dico2[key][0] +"\n"+ dico2[key][1] +"\n"+ dico2[key][2] +"\n"
            filout.write(fastq_seq1)
            filout.write(fastq_seq2)

# ~ end of script ~
