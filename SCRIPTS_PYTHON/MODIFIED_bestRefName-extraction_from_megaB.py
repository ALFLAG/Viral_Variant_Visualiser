#! /usr/bin/env python3
# -*- coding : utf8 -*-
# AF, last modification July, 20th 2018
#
# Script name: bestRefName-extraction_from_megaB.py
#
# Description:  retrieve the accession number from the best reference from a
#               Megablast file result. If multiple Megablast files are present
#               (from spades contigs and mira contigs) te percentage of identity
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
    """ This function harvest a sequence name with the highest percentage of
identity and the percentage of identity.
    """
    # flag and regex definition
    stop_flag = False
    regex = re.compile('[a-z]+\|([A-Z0-9_.]+)\|.+[\ ]?([0-9.e+]+)[\ ]{3}')

    # opening the megablast file result and saving the lines in a list
    with open(fichier,"r") as file_in:
        filin = file_in.readlines()

    # parcours la liste et recupere seulement
    # Sequences producing significant alignments

    # harvest only the best match for each contig
    # (i.e. the first match for each contig)
    best_match = []
    top  = len(filin)
    for i in range(0, top):
        if 'Sequences producing significant alignments' in filin[i]:
            best_match.append(filin[i+2])
    del(filin)

    tuple_list = []
    for line in best_match:
        AN = (regex.search(line)).group(1)
        score = (regex.search(line)).group(2)
        tuple_list.append((AN,score))

    dico = {}
    for value in tuple_list:
        if value[0] not in dico: # if AN not in the dictionnary, then we add it
            dico[value[0]] = [] # then we add it and
            dico[value[0]].append(value[1]) # the associated score is recorded
        else:
            dico[value[0]].append(value[1])]


    # conmpte le nombre de fois que la référence a été trouvé
    # et celle qui est renvoyé le plus de fois est retournée
    # peut être n'est ce pas la meileur méthode, car beaucoup de petit score seront favorisé contre un seul très gros score
    maxtemp = 0
    for key in dico:
        occ_nbre = len(dico[key])
        if occ_nbre >= maxtemp:
            max_temp = occ_nbre
            best_AN = dico[key]

    return best_AN


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
    parser.add_argument("-s","--spades",\
 help = "for a contig generated from SPAdes", type = str)
    parser.add_argument("-m","--mira",\
 help = "for a contig generated from Mira", type = str)
    parser.add_argument("-o", type = str,\
 help = "the output file where to store the result")

    args = parser.parse_args()

    # check for the presence of all arguments
    if not args.spades and not args.mira:
        parser.print_help()
        sys.exit("\nYou have to specify at least one option -s or -m. \
Please, read the help section above.\n")
    elif not args.o:
        parser.print_help()
        sys.exit("\nYour forget to specify the output file.\n")
    else:
        if args.spades and args.mira:
            mira_ = best_ref(args.mira)
            spades_ = best_ref(args.spades)

            if mira_ == spades_:
                write_AN(args.o, spades_)

            else:
                # compare identity percent, we harvest the accession number with
                # the highest percentage of identity ; in case of equality
                # a random choice is done
                score_m = float(mira_tuple[1])
                score_s = float(spades_tuple[1])
                if score_m == score_s:
                    ID_ = random.choice([score_m, score_s])
                    if ID_ == ID_s:
                        write_AN(args.o, spades_tuple[0])
                    else:
                        write_AN(args.o, mira_tuple[0])
                elif ID_m > ID_s:
                    write_AN(args.o, mira_tuple[0])
                else:
                    write_AN(args.o, spades_tuple[0])

        elif args.spades and not args.mira:
            spades_tuple = best_ref(args.spades)
            write_AN(args.o, spades_tuple[0])

        elif args.mira and not args.spades:
            mira_tuple = best_ref(args.mira)
            write_AN(args.o, mira_tuple[0])

# ~ end of script ~ #
