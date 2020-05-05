#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# modified from e.hirchaud by Alexandre FLageul
# July 3rd, 2018
#
# Description:  Get sequences from NCBI by a file contain list of id (gi or gb).
#               Check on nuclotide or protein db and output a multifasta file or
#               a genbank file.
#
################################################################################

# ~ start of script ~

#############
## Modules ##
#############

import os
import sys
import re
import argparse

from Bio import Entrez
from Bio import SeqIO

#################
## Main Script ##
#################

## variables importation ##
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type = str,\
    help="File with id for retrieve sequences, one ID per line.")
parser.add_argument('--db', type = str,\
    choices=['nucleotide', 'protein', 'Taxonomy'], default="nucleotide",\
    help="database seach : nucleotide (default) or protein or Taxonomy")
parser.add_argument("-f",'--fmt', type = str,\
    choices=['fasta', 'gb'], default="fasta",\
    help="output format 'fasta' (default) or 'genbank'.")
parser.add_argument('-o','--output', type = str, help="Name of output file")

args = parser.parse_args()

if not args.input:
	parser.print_help()
	sys.exit("""
You forgot to put the input file.
Please read the help section above""")

elif not args.output:
	parser.print_help()
	sys.exit("""
You forgot to put the output file.
Please read the help section above""")


## always tells NCBI who you are
Entrez.email = "alexandre.flageul@anses.fr"

## reading and recording of all accession numbers
id_list = []
with open(args.input) as filin:
	for line in filin:
		line = line.strip()
		id_list.append(line)

## get the sequence from online database
handle = Entrez.efetch(db=args.db, id=",".join(id_list), retmode="text",\
    rettype=args.fmt)

## write sequences in the output file
with open(args.output, "w") as filout:
	filout.write(handle.read())

# ~ end of script ~
