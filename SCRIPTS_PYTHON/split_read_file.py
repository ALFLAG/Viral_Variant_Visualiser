#! /usr/bin/env python3

#############
## Modules ##
#############

import sys
import argparse
import math
import os

###############
## Functions ##
###############


def record_all_reads(infile):
    """This function takes the read file and record reads as a dictionnary.
    {[read_id]: (seq, +, quality)}
    """
    with open(infile, "r") as filin:
        i = 0
        dico = {}
        for line in filin:
            i += 1
            if i == 1: read_id = line[:-1]
            elif i == 2:
                seq, i = line, 0
                dico[read_id] = ( "".join( (read_id + "\n", seq) ) )
    print("There is {} reads to process.".format (str(len(dico))))
    return dico

def main(read_file, output):
    """Main script that split a fasta file in several file
       that contain 500000 fasta sequences top.
    """

    # determine output directory in order to give names to partial read files
    output_directory = "".join(output.split("/")[:-1])

    # check if output files exist. If yes, they are removed before keep going
    if os.path.exists(output):
        print("Removing existing files.")
        os.system("rm {}".format(output))
        os.system("rm {}".format(output_directory + "partial_reads_file_*"))

   # open the read file in fastq format
    all_reads = record_all_reads(read_file)

    # split the read file
    # two meters: a = for the partial read file number ; i to count number of reads
    a, i = 0, 0
    for read in all_reads:
        i += 1
        with open(output_directory + "partial_reads_file_" + str(a) + ".fasta", "a") as filout:
            filout.write(all_reads[read])
    
        # every time i reach 500.000,
        # then it is set to 0 and a is goes up in order to switch the file name
        if i == 100000:
            a += 1
            i = 0

    # finally record the file that give the number of file that were created
    with open(output, "w") as filout: filout.write(str(a + 1))



#################
## Main Script ##
#################


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--read_file", type = str,
            help = "the read file in fastq format")
    parser.add_argument("--output", type = str, help = "the output file name")
    args = parser.parse_args()
    
    if not args.read_file or not args.output: 
        parser.print_help()
        sys.exit()
        
    else: main(args.read_file, args.output)
