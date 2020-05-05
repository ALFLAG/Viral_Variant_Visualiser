#! usr/bin/env/ python3                                                          
# -*- coding: utf-8 -*- 
# AF, November 26th 2018
#
# Description:  part of Me. calculate the requiered number of read to have a 80
#               mean coverage depth.
#               Display the result in the terminal, do not write output files.
#
################################################################################

# ~ start of script ~

#############
## Modules ##
#############

import argparse
import sys
import gzip
import re
import statistics

###############
## Functions ##
###############

#################
## Main Script ##
#################

## call for arguments
parser = argparse.ArgumentParser()
parser.add_argument('--reads', type = str, help = 'the read file fastq.gz fmt')
parser.add_argument('--cov', type = str, help = 'the coverage file')

args = parser.parse_args()

## check for the presence of argument
if not args.reads or not args.cov:
    parser.print_help()
    sys.exit('\nOne or more argument are missing.\n\
Please, read the help section above\n')


## Step1 - count number of read

# open the read file, and turn it into a list
# with gzip.open(args.reads, 'rb') as filin:
with open(args.reads, 'r') as filin:
    lines = filin.readlines()

# determine the length of the list and divide it by 4 to obtain
# the number of reads in the fastq.gz file
total_read_count = len(lines) / 4
#print('count', total_read_count)

# remove the list
del(lines)

## Step2 - determine the mean coverage

# define the line of the coverage file
regex = re.compile("[A-Z0-9\._]+\t[0-9]+\t([0-9]+)\n")

# initiate the storage list
cov_list = []
with open(args.cov,"r") as filin:
    for line in filin:
        # catch the coverage for the position and add it to the list
        coverage = int((regex.search(line)).group(1))
        cov_list.append(coverage)

mean_coverage = statistics.mean(cov_list)
#print('mean', mean_coverage)



# Step3 - set the wanted coverage mean
wanted_cov_depth = 80

# Step4 - calculate the requiered number of read for the wanted coverage
rd_nbre = ( total_read_count * wanted_cov_depth ) / mean_coverage

print(int(rd_nbre))

# ~ end of script ~
