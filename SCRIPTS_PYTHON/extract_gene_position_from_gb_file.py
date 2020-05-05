#! /usr/bin/env python3
#  -*- coding:utf-8 -*-
# AF - last modification: August 10th 2018
# Script name : gene_position_gb_extraction.py
#
# Description : this script take a genbank file as input file and extract the
#                virus accession number, the genome size and the start_N_end
#                positions of all genes and store these data in a json file.
#
# IMPORTANT_TO_NOTE : * The second to last part of this script has specifically
#                        been written to work with Coronaviruses genome. You may
#                        want to comment this part in order use the script for
#                        other virus genome
#
################################################################################

# ~ start script ~ #

##############
### Modules ##
##############

import re
import json
import argparse
import sys

###############
## Functions ##
###############

def get_min(dico):
    min_list = []
    for key in dico:
        if key != "genomesize" and key != "virus_id":
            min_list.append(dico[key][0])
    return min(min_list)

def get_max(dico):
    max_list = []
    for key in dico:
        if key != "genomesize" and key != "virus_id":
            max_list.append(dico[key][1])
    return max(max_list)

#################
## Main Script ##
#################


## call for the command line arguments ##
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "the gb file you want to parse",\
                    type = str)
parser.add_argument("-o", "--output", help = "the output file to create",\
                    type = str)
parser.add_argument("-a", "--virus_acronym", help = "write the virus acronym",\
                    type = str)

arg = parser.parse_args()

## check for the presence of all command line arguements ##
if not arg.input or not arg.output or not arg.virus_acronym:
    parser.print_help()
    sys.exit("\nAn error occured while entering the arguments.\n\
All three arguments are necessary, please read the help section above.\n")


## call for a dictionnary to store data ##
dico = {}

## opening and reading the genbank file ##
with open(arg.input,"r") as filein:
    line = filein.readline()
    regex1 = re.compile('^LOCUS {7}[A-Z0-9_]+ +([0-9]+)')
    regex1bis = re.compile('^VERSION {5}([A-Z0-9_\.]+)')
        # match à la première ligne du fichier
    r = regex1.search(line)
    dico["genomesize"] = int(r.group(1))

    # get the accession number of the sequence associated with its version
    for line_bis in filein:
        if "VERSION" in line_bis:
            r_bis = regex1bis.search(line_bis)
            dico["virus_id"] = arg.virus_acronym + "_" + r_bis.group(1)
            break

    # because we are not interested by the head section anymore, we escape it
    # using the following while loop
    while "FEATURES" not in line:
        line = filein.readline()

    regex2 = \
re.compile('^\ {5}CDS\ {13}join\(([0-9]+)\.{2}[0-9]+,[0-9]+\.{2}([0-9]+)\)')
        # match lines with 'CDS' and 'join', for the fragmented genes

    regex3 = re.compile('^\ {5}CDS\ {13}([0-9]+)\.{2}([0-9]+)')
        # match lines with 'CDS' for the unfragmented genes

    regex4 = re.compile('^\ {21}/product=[\"\']([A-Za-z0-9 -_]+)[\"\']')
        # match lines with '/product' to recover the product name of the gene

    flag = False # set a flag to make my life easier

    while "ORIGIN" not in line:

        if flag == True:# this is set "ON" when a CDS has previously been found

            if regex4.match(line):# look for the product name of the gene
                r = regex4.search(line)
                product_id = r.group(1)
                dico[product_id] = geneposition
                flag = False # set "OFF", in order to look for another gene
                line = filein.readline()

            else:
                line = filein.readline()

        elif regex2.match(line): # look fragmented genes
            r = regex2.search(line)
            geneposition = [ int(r.group(1)) , int(r.group(2)) ]
            flag = True
            line = filein.readline()

        elif regex3.match(line): # look for unfragmented genes
            r = regex3.search(line)
            geneposition = [ int(r.group(1)) , int(r.group(2)) ]
            flag = True
            line = filein.readline()

        else:
            line = filein.readline()


## add UTR's region to dictionnary ##
dico["UTR5"] = [1,get_min(dico)-1]
dico["UTR3"] = [get_max(dico)+1, dico["genomesize"]]

## writting the data into a json file ##

filout = open(arg.output,"w",encoding = 'utf-8')
json.dump(dico, filout, ensure_ascii = False)
filout.close()

# ~ end of script ~ #
