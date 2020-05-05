#! /usr/bin/env python3
# -*- coding : utf8 -*-
# AF, last modification July, 20th 2018
#
# Script name: bestRefName-extraction_from_megaB.py
#
# Description:  retrieve the accession number from the best reference from a
#               Megablast file result. If multiple Megablast files are present
#               (contigs from spades and mira) the percentage of identity
#               of both file is compared and the best reference is kept.
#
################################################################################

# ~ start script ~ #

#############
## Modules ##
#############

import argparse
import sys
import re
import random

###############
## Functions ##
###############

def best_ref(fichier):
    """ This function  harvest a sequence name with the highest percentage of
identity and the percentage of identity.
    """
    # flag and regex definition
#    stop_flag = False
    regex = re.compile('[a-z]+\|([A-Z0-9_.]+)\|.+[\ ]?([0-9.e+]+)[\ ]{3}')

    # opening and recording the lines in a list
    with open(fichier,"r") as file_in:
        filin = file_in.readlines()

    # check two conditions:
    # if the previous line IS NOT MATCHING the regular expression
    # if the current line IS MATCHING the regular expression
    # to do it, start parsing the list at position 2 (index 1) and finish it
    # by respecting these conditions, we are sure to parse the entire megablast
    # result file and to extract data matching 1)a virus description and 2) that
    # this result is the first (the highest score) for the considered contig
    file_size = len(filin)
    none = True
    for i in range(1, file_size):
        if not regex.match(filin[i-1]) and regex.match(filin[i]):
            if "virus" in filin[i] \
            or "Virus" in filin[i] \
            or "VIRUS" in filin[i]:
            # if conditions are good, then extract the data
                none = False
                match = regex.search(filin[i])                               
                AN = match.group(1)                                      
                score = match.group(2)                                   
                break
    if none == True:
        return "NO_MATCH"
    else:
        return (AN, score)

def write_AN(out, AN):
    with open(out,"w") as filout:
        filout.write(AN)
    sys.exit()


#################
## Main Script ##
#################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = \
"for the extraction of the best reference's name")
    parser.add_argument("-s","--spades", type = str,\
 help = "for a contig generated from SPAdes")
    parser.add_argument("-m","--mira", type = str,\
 help = "for a contig generated from Mira")
    parser.add_argument("-o", type = str,\
 help = "the output file where to store the result")

    args = parser.parse_args()


    # in case of NO_MATCH
    out_text = """

The program stopped because it appears that no previously generated contig is
matching a virus after the megablast step.
Something probably went wrong. Sorry.

"""

    # check for the presence of all arguments
    if not args.spades and not args.mira:
        parser.print_help()
        sys.exit("\nYou have to specify at least one option -s or -m. \
Please, read the help section above.")
    elif not args.o:
        parser.print_help()
        sys.exit("\nYour forget to specify the output file.")
    else:
        if args.spades and args.mira:
            mira_tuple = best_ref(args.mira)
            spades_tuple = best_ref(args.spades)
            if mira_tuple == "NO_MATCH" and spades_tuple == "NO_MATCH":
                sys.exit(out_text)
            elif mira_tuple != "NO_MATCH" and spades_tuple == "NO_MATCH":
                write_AN(args.o, mira_tuple[0])                                
            elif spades_tuple != "NO_MATCH" and mira_tuple == "NO_MATCH":
                write_AN(args.o, spades_tuple[0])
            else:                                
                if mira_tuple[0] == spades_tuple[0]:
                    write_AN(args.o, spades_tuple[0])
                else:
                    # compare identity percent, we harvest the accession number with
                    # the highest percentage of identity ; in case of equality
                    # a random choice is done
                    score_m = float(mira_tuple[1])
                    score_s = float(spades_tuple[1])
                    if score_m == score_s:
                        ID_ = random.choice([score_m, score_s])
                        if ID_ == score_s:
                            write_AN(args.o, spades_tuple[0])
                        else:
                            write_AN(args.o, mira_tuple[0])
                    elif score_m > score_s:
                        write_AN(args.o, mira_tuple[0])
                    else:
                        write_AN(args.o, spades_tuple[0])
    
        elif args.spades and not args.mira:
            spades_tuple = best_ref(args.spades)
            if spades_tuple == "NO_MATCH":
                sys.exit(out_text)
            else:
                write_AN(args.o, spades_tuple[0])

        elif args.mira and not args.spades:
            mira_tuple = best_ref(args.mira)
            if mira_tuple == "NO_MATCH":
                sys.exit(out_text)
            else:
                write_AN(args.o, mira_tuple[0])
    
# ~ end of script ~ #
