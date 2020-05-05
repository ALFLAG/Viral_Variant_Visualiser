#! /usr/bin/env python3

#############
## Modules ##
#############

import sys
import argparse

###############
## Functions ##
###############

def get_read_id(infile):
    """
    This function takes the output file from megablast (with -outfmt 6 and -max_target_seqs options)
    and record as a list all reads id that have matched the vdb database.
    """

    with open(infile, "r") as filin:
        rid = [ "@"+ (line.rstrip("\n") ).split("\t")[0] for line in filin ]

    return tuple(rid)

def record_all_reads(infile):
    """
    This function takes the read file and record reads as a dictionnary.
    {[read_id]: (seq, +, quality)}
    """

    with open(infile, "r") as filin:
        i = 0
        dico = {}
        for line in filin:
            i += 1
            if i == 1: read_id = line.rstrip("\n").split(" ")[0]
            elif i == 2: seq = line
            elif i == 3: plus = line
            elif i == 4:
                qual, i,= line, 0
                dico[read_id] = ( "".join( (read_id + "\n", seq, plus, qual) ) )

    return dico

def main(id_file, read_file, output):
    """
    """
    # open the read id file
    reads_id = get_read_id(id_file)
    
    # open the read file in fastq format
    all_reads = record_all_reads(read_file)

#    print(reads_id[0])
#    for key in all_reads:
#        print(key)

    # write output file
    with open(output, "w") as filout:
        for ID in reads_id: filout.write(all_reads[ID])


#################
## Main Script ##
#################


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--read_id", type = str,
            help = "the ouput file from Megablast results with -outfmt 6 and -max_target_seqs 1 options")
    parser.add_argument("--read_file", type = str,
            help = "the read file in fastq format")
    parser.add_argument("--output", type = str, help = "the output file name")
    args = parser.parse_args()
    
    if not args.read_id or not args.read_file or not args.output: 
        parser.print_help()
        sys.exit()
        
    else: main(args.read_id, args.read_file, args.output)
