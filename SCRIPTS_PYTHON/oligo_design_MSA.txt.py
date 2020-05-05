#! /usr/bin/env python
# -*- coding: utf8 -*-

# Generateur d'amorce d'hybridation pour des Acides Nucleiques specifiques

###############################################################################
###############################################################################
#########################                             #########################
#########################           MODULES           #########################
#########################                             #########################
###############################################################################
###############################################################################

import argparse
import sys

###############################################################################
###############################################################################
#########################                             #########################
#########################          FUNCTIONS          #########################
#########################                             #########################
###############################################################################
###############################################################################

###############################################################################
###############################################################################
#########################                             #########################
#########################         MAIN SCRIPT         #########################
#########################                             #########################
###############################################################################
###############################################################################

##################### Definitions des variables importées #####################
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",\
        help="this is the input file with sequence in fasta format from a MSA",\
        type = str)

parser.add_argument("-o", "--ouput",\
        help="this is the output file with oligo sequence in fasta format",\
        type = str, action = "store_true",\
        default = "primer_from_MSA.txt")

parser.add_argument("-g", "--gc_percent",\
        help="this is the minimum and maximum GC% you want for your targets as \
a list [min, max], type list of float",\
        type = list, action = "store_true",\
        default = [0.2,0.8])

parser.add_argument("-m", "--molar_conc",\
        help="what is the molar concentration of salt in the mix",\
        type = float, action = "store_true",\
        default = 0)

parser.add_argument("-f", "--formal",\
        help="what is the percentage of formaldehyde in the mix, type float \
between 0 and 1",\
        type = float, action = "store_true",\
        default = 0)

parser.add_argument("-M", "--mismatch",\
        help="how many mismatch are you agree to incorpore in your oligo ?",\
        type = int, action = "store_true",\
        default = 0)

parser.add_argument("-l", "--oligo_len",\
        help="What size do you want your oligos ?",\
        type = int, action = "store_true",\
        default = 20)

args = parser.parse_args()


f_in = args.input()

if args.output == TRUE:
    f_out = args.output()

if args.oligo_len == TRUE:
    spec_primer_lgth_input = args.oligo_len()

if args.gc_percent == TRUE:
    GC_min = args.gc_percent()[0] ; GC_max = args.gc_percent()[1]

if args.molar_conc == TRUE:
    molar_conc = args.molar_conc()

if args.formal == TRUE:
    forma_percent = args.formal()

if args.mismatch == TRUE:
    mis_m = args.mismatch()

###################### Lecture et Stockage des Sequences ######################
# On ouvre le fichier d'alignement multiple
file_in = open(f_in,"r")
temp_line_lgth = 0
# On créer une liste qui stockera toutes les séquences du fichier.
line_list = []
# On va lire l'intégralité du fichier.
for line in file_in:
    # On verifie que la ligne correspond à la sequence et pas au nom.
    if line[0] != ">":
        line_ = line.replace("\n","").replace("\r","")
        line_lgth = len(line_)
        # On regarde si les séquences ont toutes la même taille.
        if temp_line_lgth == 0:
            temp_line_lgth = line_lgth
            line_list.append(line_)
        elif line_lgth != temp_line_lgth:
            sys.exit('''
Sequences have not the same length.
Please, look again your multiple alignment to generate sequence of same length.
For example, add gaps to the short sequence ends or remove the part of the
longest sequences that are not usefull for the rest of the analysis.
''')
        else:
            line_list.append(line_)
# On termine par fermer le fichier d'entrée
file_in.close()
# On termine par conserver la liste en tant que tuple,
# cela prend moins de place dans la mémoire.
line_t = tuple(line_list) ; del(line_list)

####################       Fin de lecture du fichier       ####################



########################### Determination des amorces ##########################

###########################       Kmer methode        ##########################
# On va ici decomposer l'integralite des sequences en kmers de longueur
# specifiee par l'utilisateur.
# Pour une sequence donnee, une liste de kmers est generee.

# Pour l'attributiOn de numero a chaque sequence.
seq_number = 0

# On creer une liste pour ranger les dictionnaires seq_number/kmer_list associes
dico_list = []

# Cette boucle va permettre de generer la liste de tous les kmer pour toutes
# les sequences du fichier MSA.
# On commence par parcourir la liste de sequence.
for seq in line_t:
    # On attribue un numero pour chaque sequence
    seq_number = seq_number + 1
    # On creer une liste pour ranger les dictionnaires de position/kmer
    kmer_list = []
    # On parcours la sequence en cours 
    for i in range(0:length(seq)-(spec_primer_lgth_input-1)):
        # Pour chaque positiOn dans la sequence, un kmer est genere
        kmer = seq[i:i+spec_primer_lgth_input-1]
        # On enregistre la position de départ du kmer
        start = i
        # On enregistre la position et le kmer associe dans un dictionnaire
        dico = {"start" = start ; "kmer" = kmer}
        # On ajoute ce dictionnaire à la liste de kmer
        kmer_list.append(dico)
    # Une fois toute la sequence parcourue, un dictionnaire est creer pour
    # contenir le numero de la sequence et la liste de kmer associee.
    dico_seq = {"seq_number" = seq_number ; "kmer_tuple" = tuple(kmer_list)}
    # Ce dictionnaire est alors stocké dans une liste.
    dico_list.append(dico_seq)
# On supprime la liste de sequence de la RAM afin de liberer de l'espace
dico_t = tuple(dico_list)
del(dico_list)
del(kmer_list)

# Une fois toutes les sequences parcourues, une liste de kmer a ete genere pour
# chaque sequence.
# Il est alors temps de comparer les tuples de kmers entre eux.


# On definit la premiere sequence, qui correspond au premier dictionnaire de
# tuple ([dict(seq1, kmer_tuple)],      )
first_seq = dico_t[0]
# On creer une liste ou seront stockes tous les dico de kmer utilisable
raw_kmer_list = []
# Cette boucle va permettre de comparer les kmer entre eux a une meme position
# donnee.
# On parcours la liste de dictionnaire.
for kmer_dico in first_seq{"kmer_tuple"}:# kmer_dico est un dict -start, kmer-
    # dans le dictionnaire kmer_dico on veut acceder à la sequence kmer
    kmer_seq = kmer_dico{"kmer"}# On accede a la premiere sequence kmer
    green_flag = TRUE
    for i in range(1,len(dico_t)):
        next_dico = dico_t[i]
        next_kmer_tuple = next_dico{"kmer_tuple"}
        next_kmer_seq = next_kmer_tuple{"kmer"}
        if kmer_seq != next_kmer_seq:
            green_flag = FALSE
            break
    # Une fois tous les kmers compares pour une position donnee, si tout est
    # identique, alors on enregistre le kmer pour la prochaine étape
    if green_flag = TRUE:
        raw_kmer_list.append(kmer_dico)
raw_kmer_tuple = tuple(raw_kmer_list)
del(raw_kmer_list)

# Une fois que les kmers identique dans toutes les sequences ont ete
# selectionnes, il faut calculer le pourcentage en GC.
# Pour cela, on utilise une boucle qui va parcourir la liste de dictionnaire.
# On creer une liste qui stockera tous les dictionnaires utilisables.
raw_kmer_list_2 = []
for dico in raw_kmer_list:
    # On recupere la sequence du kmer
    seq = dico{"kmer"}
    # On compte le nombre de G et de C
    GC_number = seq.count("G") + seq.count("C")
    # On calcule le pourcentage en GC
    GC_percent = (100*GC_number)/spec_primer_lgth_input
    # On calcule le Tm du kmer
    AT_number = spec_primer_lgth_input - GC_number
    if spec_primer_lgth_input <= 20:
        Tm = AT_number*2 + GC_number*4
    elif:
        M_ = sys.argv[5]
        F_ = sys.argv[6]
        m_ = sys.argv[7]
        Tm_ = 81,5 + 16,6*math.log10(M_)+(0,41*GC_percent)
        Tm = Tm_–(500/spec_primer_lgth_input) – 0.6*F_ – 1,5*m_
        # M : concentration molaire en ions monovalents,
        # L : longueur de l’amorce (en base),
        # f : pourcentage de formamide,
        # m : pourcentage de mauvais appariement.
    # On regarde si le pourcentage en GC est conforme a la limite specifiee
    if GC_percent >= spec_GCpercent_limit[0]:
        if GC_percent <= spec_GCpercent_limit[1]:
            # On ajoute le pourcentage en GC au dictionnaire
            dico{"GC%"} = GC_percent
            dico{"Tm"} = melt_temp
            # On ajoute le dictionnaire a la liste
            raw_kmer_list_2.append(dico)
# On supprime la liste precedente de la RAM
del(raw_kmer_list)

file_out = open("result_probes.txt","w")
for dico in raw_kmer_list_2:
    start_position = dico{"start"}
    seq = dico{"kmer"}
    Tm = dico{"Tm"}
    GC_p = dico{"GC%"}
    file_out.write(
    ">position_{} \
    | GC_percent = {:.2f} \
    | Tm = {:.1f}°C\n\
    | {}".format(start_position, GC, Tm, seq))
file_out.close()
