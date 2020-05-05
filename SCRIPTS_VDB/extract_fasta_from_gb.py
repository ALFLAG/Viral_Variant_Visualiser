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
To extract fasta sequence from a genbank file
and record it in a file in append mode.
    """

    for seqRecord in SeqIO.parse(infile, "genbank"):
        AN, DEF, fasta = seqRecord.id, seqRecord.description, str(seqRecord.seq)
        name = " ".join( (">"+AN, DEF) )
        tax = ";".join(seqRecord.annotations["taxonomy"])
        full = "\n".join( (name + "#" + tax, fasta,"") )
        with open(outfile, "a") as filout:
            filout.write(full)


def main_script(infile, outfile="all_viral_seq.fasta"):
    """
To parse a file containing all the gbfile name.
    """
    print("Open the input file.\n")

    with open(infile, "r") as filin:
        l = [ line[:-1] for line in filin ]
        lines = tuple(l)
    maxi = len(lines)

    print("Recording the fasta sequence.\n")
    i = 0    
    for name in lines:
        get_fasta(name, outfile)
        i += 1
        if i < maxi and i % 50000 == 0:
            print("{0} sequences have been treated.".format(i))
        elif i == maxi:
            print("All sequences have been treated.")


#################
## Main script ##
#################

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("""
Input file is missing.
Usage: python3 script input output
""")
    else:
        main_script(sys.argv[1], sys.argv[2])
#        get_fasta(sys.argv[1], sys.argv[2])
        sys.exit()

# ~ end of script ~
