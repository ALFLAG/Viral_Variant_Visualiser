#! /usr/bin/env python3

#############
## Modules ##
#############

import re
import argparse
import sys


##############
## Function ##
##############



def get_consensus_base(seq, position, infile):
    """
    Accept one string of ucleotides and count the occurence of each.
    
    """

    # new_seq pour transformer la sequence brute an sequence exploitable
    # i le compteur qui permet de parcourir la sequence
    # a la taille de la sequence, afin d'arrter la boucle
    # b la variable qui va servir à incrémenter la variable i
    new_seq, i, a = [], 0, len(seq)
    SEQ = seq.upper()
    element_set = set()
    while i < a:
        if SEQ[i] == "+":
            previous = SEQ[i-1]
            pos, caract, nb = i+1, "0", ""
            while caract in ("0","1","2","3","4","5","6","7","8","9"):

                # on enregistre les chiffres qui suivent le signe "+",
                # l'enregistrement complet se fait lorsque tous les tours de
                # boucles seront effectués
                nb = nb + SEQ[pos]

                # on enregistre le caractère qui suit le premier chiffre
                # si le caractère n'est pas un chiffre (en format string) alors la boucle ne redémarrera pas
                # si le caractère est un chiffre (en format string) alors la boucle redémarrera
                caract = SEQ[pos+1]

                # on incrémente la position pour redémarrer la boucle, si besoin
                pos += 1

            # nb correspond au nombre de nt qui suivent le nombre indiqué après le signe '+'
            # 1 pour ajouter l'espace occupé par le symbole + (soit 1 position)
            # len(nb) pour ajouter l'es pace occupé par la taille de la chaine de caractère enregistré dans la variable nb
            followed = SEQ[ i + len(nb) + 1 : i + len(nb) + int(nb) + 1]
            new_seq.append(previous + followed)
            element_set.add(previous + followed)


            b = int(nb) + 1 + len(nb)

        # dans la prochaine condition, on rencontre le signe négatif,
        # indiquant qu'il y aura une/des délétion(s) dans la/les colonnes qui
        # suivent. Le nombre de ces délétions est donné par le chiffre/nombre
        # qui suit le signe négatif
        # du coup, on retire le motif de la séquence
        elif SEQ[i] == "-":
            pos, caract, nb = i+1, "", ""

            # la boucle sert à récupérer le nombre de N qui suivent le symbone +
            while caract != "N":
                nb = nb + SEQ[pos] # on enregistre le premier chiffre qui suit le "-"
                caract = SEQ[pos+1] # on enregistre le caractère qui suit le premier chiffre, si N alors la boucle ne redémarrera pas
                pos += 1 # on incrémente la position pour passer à la boucle suivante, si il y a besoin

            # nb correspond au nombre de nt qui suivent le nombre indiqué après le signe '-'
            # 1 pour ajouter l'espace occupé par le symbole + (soit 1 position)
            # len(nb) pour ajouter l'es pace occupé par la taille de la chaine de caractère enregistré dans la variable nb
            b = int(nb) + 1 + len(nb)

        # Dans la prochaine condition, le caractère marque le début d'un read.
        # Le caractère suivant correspond qualité de la base qui suit.
        # il ne faut donc pas prendre en compte ce cartère, ni le suivant,
        # et aller directement au prochain caractère qui matche un nucleotide.
        elif SEQ[i] == "^": b = 2

        # Dans la prochaine condition, le caractère marque la fin du read.
        # il ne faut donc pas le prendre en compte et passer au prochain nt.
        elif SEQ[i] == "$": b = 1

        else:
            new_seq.append(SEQ[i])
            element_set.add(SEQ[i])
            b = 1

        # pour finir on passe à la prochaine base.
        i = i + b

    # maintenant qu'il y a une séquence exploitable,
    # il faut compter le nombre d'occurence de chaque base ou caractère '*',
    # et renvoyer celle qui est la plus présente
    base_list = [ (new_seq.count(elmt), elmt) for elmt in element_set ]
    
    # avec la fonction sort,
    # le tri d'une liste de liste, ou d'une liste de tuple,
    # se fait en fonction du premier element de la sous-liste (ou du tuple)..

    with open(infile+"_all_pos_nt", "a") as filout:
        filout.write(str(position) + "\t" + str(sorted(base_list)) + "\n")


    return sorted(base_list)[-1][1] 


### Choose what nucleotide to return
#
#quality_score = { '!':'0',  '"':'1',  '#':'2',  '$':'3',  '%':'4',  '&':'5',
#                  "'":'6',  '(':'7',  ')':'8',  '*':'9',  '+':'10', ',':'11',
#                  '-':'12', '.':'13', '/':'14', '0':'15', '1':'16', '2':'17',
#                  '3':'18', '4':'19', '5':'20', '6':'21', '7':'22', '8':'23',
#                  '9':'24', ':':'25', ';':'26', '<':'27', '=':'28', '>':'29',
#                  '?':'30', '@':'31', 'A':'32', 'B':'33', 'C':'34', 'D':'35',
#                  'E':'36', 'F':'37', 'G':'38', 'H':'39', 'I':'40'}
#
#
##utiliser les scores de qualité afin de pouvoir donner la qualité moyenne pour un nt donné à cette position
##    deux choses doivent être prise en compte pour retrouner le NT consensus à cette position
##    la quantité et la qualité
##    si il y a 1000 read avec * et 990 read avec T, on doit retourner T a condition que la qualité moyenne des T soit suffisante
##        pourquoi ? car sur un chromatogramme, on ne verrait que le pic du T, et on ne verrait pas du tout le pic du * puisqu'il n'y aurait rien à voir
##
##ATTENTION a ne pas tomber dans le probleme de bcftools, a savoir, ignorer les * dans certaines conditions (condition floue), ce qui aboutirait a avoir des résultats faux en terme de compositions de population virale
#
#    # if only one nt, then return it
#    if len(sorted(base_list)) == 1: return sorted(base_list)[-1][1]
#
#    # if several are possible
#    # check if the most present character is *, means deletion
#    elif sorted(base_list)[-1][1] == '*':
#        # if it's a star, count the occurence of the star and the second most present character
#
#        stars_count      = sorted(base_list)[-1][0]
#        other_base       = sorted(base_list)[-2][1]
#        other_base_count = sorted(base_list)[-2][0]
#
#        #if the stars_count is less than 200, then return it
#        if stars_count < 200: return sorted(base_list)[-1][1]
#
#        # if the stars_count is more than 200 and the other base in more than 20% of te stars count (mean at least 40 occurences)
#        # but if in this case 40 occurences is good enough, why in the first case it would not be good ? True, it does have any sense
#        elif stars_count > 200 and other_base_count > stars_count*0.20:
#            return sorted(base_list)[-2][1]
#        else:
#            return sorted(base_list)[-1][1] 
#    
#    # if the most present character is not *, then return the most present
#    else: return sorted(base_list)[-1][1] 


def main(infile):
    """
    Function that parse the output file from samtools mpileup -aa command.
    
    """

    #test="KR822424.1\t1\tN\t0\t*\t*\n"
    empty_regex = re.compile("[A-Za-z0-9\_\-\.]+\t[0-9]+\tN\t[0-9]+\t\*\t\*\n")

    base = []
    with open(infile, "r") as filin:
        for line in filin:
            AN = line.split("\t")[0]
            if empty_regex.match(line): base.append("N")
            else: base.append( get_consensus_base(line.split("\t")[4], line.split("\t")[1], infile ))

    seq_name = ">consensus-from_reference_{0}\n".format(AN.split(".")[0])
    last_seq = [ nt for nt in base if nt != "*"]

    print(seq_name + "".join(last_seq))

#################
## Main Script ##
#################

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pileup",
        help = "the pileup file from samtools mpileup command")

arg = parser.parse_args()

if not arg.pileup : parser.print_help()
else: main(arg.pileup)

