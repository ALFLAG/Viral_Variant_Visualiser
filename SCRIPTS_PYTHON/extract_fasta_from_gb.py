#! /usr/bin/env python3
#
# Name: extract_fasta_from_gb.py
#
###############################################################################

#############
## Modules ##
#############

from Bio import SeqIO
import sys


##############
## Function ##
##############

def get_fasta(infile, outfile):
    """
To extract fasta sequence from a genbank file and record it in the output file.
    """

    for seq_record in SeqIO.parse(infile, "genbank"):
        AN = seq_record.id
        fasta = str(seq_record.seq)
        full = "\n".join( (">"+AN, fasta,"") )
        with open(outfile, "w") as filout:
            filout.write(full)

#################
## Main script ##
#################

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("""
Incorrect argument number.
Usage: python3 script input output
""")
    else:
        get_fasta(sys.argv[1], sys.argv[2])

# ~ end of script ~
