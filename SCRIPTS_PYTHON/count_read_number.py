#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# AF, last modification September 7th 2018
# Name: count_read_number.py
#
# Description:  Script to count the number of reads in a gzipped (or not) file,
#               ans dtermine automatically the reads format (fasta or fastq).
#
################################################################################

# ~ start of script ~

############
## Module ##
############

import argparse
import gzip

##############
## Function ##
##############

def count(filin, b=0):
    a = 1
    line = filin.readline()

    if b == 1: # the file is in gzipped format
        if line[2] == "@":
            fastq = True
        elif line[2] == ">":
            fasta = True
        while line:
            line=filin.readline()
            a += 1

    else: # the file is not in gzipped format
        if line[0] == "@":
            fastq = True
        elif line[0] == ">":
            fasta = True
        while line:
            line=filin.readline()
            a += 1

    if fastq == True:
        print("fastq_" + str(int(a/4)))
    elif fasta == True:
        print("fasta_"+str(int(a/2)))

#################
## Main Script ##
#################

parser = argparse.ArgumentParser()
parser.add_argument("-f", type = str, help = "the fastq file you want to parse")

arg = parser.parse_args()

if arg.f[-3:] == ".gz":
    with gzip.open(arg.f, "rb") as filin:
        count(filin,b=1)
else:
    with open(arg.f, "r") as filin:
        count(filin)
 
# ~ end of script ~
