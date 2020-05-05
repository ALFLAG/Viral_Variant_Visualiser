#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# AF version, September 4th, 2018
# Script name: add_position_to_consensus.py
#
# Description:  Take one Needleman-Wunsch result file and a json file containing
#               the gene names and positions of the reference.
#               The result operation is the creation of a json file containing
#               name and positions of the genes in the consensus sequence.
#
#               How does it work ?
#               in a nutshell, the program will extract both reference and
#               consensus sequences. Then it will attribute numbers/position to
#               the reference sequence, taking in account the presence of gap
#               inside the reference sequence. Gene names will be added to each
#               position, and the gene list will be transposed to the consensus
#               sequence. Gene names matching gaps in the consensus sequence
#               will be removed.
#
################################################################################

## ~ start of script ~ ##

#############
## Modules ##
#############

import re
import json
import argparse
import sys
import operator

#################
## Main script ##
#################

# call for arguments
parser = argparse.ArgumentParser()
parser.add_argument("-j", "--json", type = str,\
                help = "the json file containing the gene names and positions")
parser.add_argument("-i", "--input", type = str,\
                help = "the Needleman-Wunsch alignement file")
parser.add_argument("-o","--output", type = str,\
                help = "the output file in json format")
arg = parser.parse_args()

## check for the presence of all arguments
if not arg.json or not arg.input or not arg.output:
    parser.print_help()
    sys.exit("\nAn error occured when entering the arguments."\
" Please, read the help section above\n")

## search for the sequence portion of line in the needle file

a = 1 # flag initialization, 1 = ref (-asequence), 2 = cons (-bsequence)
consensus = [] # to store the reference sequence
reference = [] # to store the consensus sequence

# regex mathching the lines with sequences with nucleotide ambiguities
regex1 = re.compile('[A-Za-z0-9\.\_\-]+\ +[0-9]+\ {1}([ATCGKMSWRYBDHVNUatcgkmswrybdhvnu\-]+)\ +[0-9]+\n')         

# open and read the file
with open(arg.input,"r") as filin:
    for line in filin:
        blank = re.compile('^[ |.]+\n$')
        if "#" in line or line == "\n" or "|" in line or blank.match(line): # skip the header section
            continue
        else:
            if a == 1: # reference sequence
                reference.append((regex1.search(line)).group(1))
                a = 2 # sequence switching, the next sequence wiil be consensus
            elif a == 2:# consensus sequence
                consensus.append((regex1.search(line)).group(1))
                a = 1 # sequence switching, the next sequence will be reference


## generate sequences and not lists of partial sequence
cons , ref = "".join(consensus) , "".join(reference)


## load the data contained in the json file, it's the ref gene coordonates
filin = open(arg.json,"r")
dico = json.load(filin)
filin.close()

## record and sort data in a list of tuples
gene_list = []
for key in dico:
    if key != "virus_id"  and key != "genomesize":
        # key = gene_id ; position 0 = start_of_gene, position 1 = end_of_gene
        gene_list.append((key, dico[key][0], dico[key][1]))
    elif key == "virus_id":
        virus_id = dico[key]

## sort the list according to the order of gene (the start position of genes)
sorted_genelist = sorted(gene_list, key=operator.itemgetter(1), reverse = False)


## add gene positions (numbers) to each base of the reference
# everytime the position match a nucleotide, a new position is attached
# everytime the position match a gap, the previous position is attached
a = 0
position_list = []
for i in range(0,len(ref)):
    if ref[i] != "-":
        a += 1
        position_list.append(a)
    else:
        position_list.append(a)

## add gene_id to the list of position
genelist = [] # to store list of gene to each position
for i in range(0,len(position_list)): # run all over the list of numbers

    x = position_list[i] # each numbers is a position
    name = "" # initialization of the name for this specific position

    # flag initialization to get (or not) 'INTERGEN' to this position
    is_in_gene = False

    for gene in gene_list: # run all over the list of gene
        b_inf , b_sup = int(gene[1]) , int(gene[2]) # define borders of genes

        if b_inf <= x <= b_sup: # loof if the position belongs to the gene
            # if yes...
            # then the flag become green, and the name will not be INTERGEN
            is_in_gene = True
            name = name + "//" + gene[0] # name definition
    # do it again, in case of superimposed genes

    # if the flag is red, then the position is not present in the list of gene
    if is_in_gene == False:
        name = "//INTERGEN" # so we give the name INTERGEN to the position

    genelist.append(name[2:]) # add the name of the positon to the genelist

## give the name of each position to the consensus.
cons_gene = [ genelist[i] for i in range(0,len(genelist)) if cons[i] != "-" ]

## create the new dictionnary in to store date in a json format

dico = {} # dictionnary initialization
intergen = [] # to store intergene positions
for i in range(0, len(cons_gene)): # run all over the genelist
    gene = cons_gene[i] # and define the gene

    # check the presence of gene (that is not intergen) in the dictionnary
    if gene not in dico and gene != "INTERGEN":
        # if not present, the gene is added tothe dictionnary
        dico[gene] = (i+1, i + cons_gene.count(gene))

    # if the position match an INTERGEN
    elif gene == "INTERGEN":
        # then we give intergen a number, matching the position of the intergene
        # e.g. 1,2,3,4,5,6,102,103,104,105,200,201,202,203)
        intergen.append(i+1)

# add the virus name and the genome size to the dictionnary
id_virus = virus_id.split("_")[1]
dico["virus_id"] = "consensus-from_reference_" + id_virus
dico["genomesize"] = len(cons_gene)

# the dictionnary is not ready yet, we also need to add the INTERGENE regions
# we define lists, to store id and positions...
new_intergenA = [] # store tuples (pos, id)
new_intergenB = [] # store id

a = 1 # initialization of the intergenc region number
for i in range(0,len(intergen)-1):

    # check if the value at position i is directly below the value at position i+1
    # e.g. we have a 5 at position i, and a 6 at position i+1
    if intergen[i] + 1 == intergen[i+1]:
        # if yes, the positions i and i+1 belongs to the same intergenic region
        # therefore, they will have the same name
        # add the position and the name of this intergenic region to the list A
        # add the name of this intergenic region to the list B
        new_intergenA.append((intergen[i],"INTERGEN"+str(a)))
        new_intergenB.append("INTERGEN"+str(a))

    # if  value in i is not directly below  value in i+1
    # e.g. we have a 4 at position i and 19 at position i+1
    # then we will change of intergenique region during the next iteration /
    # at the next position
    else:
        a += 1

# add the intergenic region to the dictionnary
# using the same logic than for the gene regions
for i in range(0,len(new_intergenA)):
    gene = new_intergenA[i][1]
    if gene not in dico:
        pos_init = new_intergenA[i][0]
        dico[gene] = (pos_init, pos_init + new_intergenB.count(gene))
# the dictionnary is now ready!

## write the final json file containing the consensus gene coordonates
filout = open(arg.output,"w")
json.dump(dico, filout, ensure_ascii = False)
filout.close()

## ~ end of script ~ ##
