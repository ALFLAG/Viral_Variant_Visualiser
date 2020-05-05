#! /usr/bin/env python3
# non fini
#
# Name: compare_cds.py
#
################################################################################

#############
## Modules ##
#############

import argparse
import json
import sys
import copy
from Bio import Seq
from Bio import Alphabet

###############
## Functions ##
###############

def read_consensus_sequence(sequence):
    """
    Open the consensus sequence file and record only the sequence as a string.
    """
    print("Read the consensus sequence")

    seq = []
    with open(sequence, "r", encoding = "utf-8") as seqin:
        for line in seqin:
            if ">" not in line:
                seq.append(line[:-1])

    return "".join(seq)


def read_jsonfile(jsonfile):
    """
    Open the gene position file in json format.
    """
    print("Open the Json file, a.k.a the geneposition file")

    gene_file = open(jsonfile, "r", encoding = "utf-8")
    gepo = json.load(gene_file) # dictionnary containing gene position
    gene_file.close()
    return gepo


def read_variant_file(variant_file):
    """
    Open the variant file and extract the genes names.
    """
    print("Open and extract gene names from the variant file")

    with open(variant_file, "r", encoding = "utf-8") as filin:
        inlist = filin.readlines()[1:] # remove the header section

    var = [] # using a set, we avoid to add duplicate genes
    geneset = set()
    for elmt in inlist:
        line = elmt.split("\t")
        var.append(line)
        if "//" in line[5]:
            for g in line[5].split("//"):
                geneset.add(g)
        else:
            geneset.add(line[5])

    return tuple(var), tuple(geneset)


def get_consensus_genes(geneposition, consensus, genes):
    """
    Extract the genes consensus sequence and return a tuple containing:
        1) tuple of "start" and "end" position in the genome ;
        2) the gene sequence.
    """
    print("Extract the consensus genes from consensus sequence.")

    consensus_genes = {}
    for tpl in geneposition:
        if tpl[0] in genes:
            gene, start, end = tpl[0], tpl[1]-1, tpl[2]
            consensus_genes[gene] = ( (start, end), consensus[start:end] )
            print(gene, consensus[start:end][0:4])
    return consensus_genes


def get_mutated_genes(consensus_genes, variant_tpl):
    """
    Convert the genes consensus sequence as variant sequence.
    """
    print("Conversion of consensus genes to variant genes")

    # parse the list of variant and for each variant record wanted data
    var_list = []
    for var in variant_tpl:
        indice, gene = var[0], var[5]
        pos, ref, alt = int(var[1])-1, var[2], var[3]
        # get consensus gene and loop over it to generate the mutated sequence
        if gene in ["Intergen", "INTERGEN", "intergen"] or "UTR" in gene:
            isInterGen, mut_type, refseq, varseq = True, None, None, None
        else:
            isInterGen = False
            if "//" in gene:
                for g in gene.split("//"):
                   var_list.append(get_gene(indice, consensus_genes[g],
                                            pos, ref, alt, isInterGen) )
            else:
               var_list.append(get_gene(indice, consensus_genes[gene],
                                        pos, ref, alt, isInterGen) )

    return tuple(var_list)


def get_gene(indice, csgn, pos, ref, alt, isInterGen):
    """
    This function will turn the consensus gene inot the variant gene.
    Basically, it copy/paste the consensus gene except for the position(s)
    where the variant is
    """

    new_seq, start, end, gn = [], csgn[0][0], csgn[0][1], csgn[1]
    i = start
    while i != end:

        # if the position is not matching the variant position
        # then add the nucleotide in the mutated gene
        if i != pos:
            new_seq.append(gn[i-start])
            i += 1 # go to the next position
            # i = position in the full length genome
            # start = position of the beginning of the gene in the full
            # length genome
            # the difference is the position in the gene itself

        # if the position is matching the variant position
        # then compare the length of the ref and alt
            # if equal, just add the variant and go to the position
            # i + len(var)
            # if different, add the variant and calculate the difference,
                # if the diff is negative,
                    # the ref length is shorter than the var length
                # if the diff is positive,
                    # the var length is shorter than the ref length
                # finally go to the i + diff position
            # this dicotomy has incidence on the position of the next nt
            # basically, we always add the variant and respect the
            # continuity of the ref
        else:
            new_seq.append(alt)
            i += len(ref)
            diff = len(ref) - len(alt)
            if diff == 0: mut_type = "sub"
            elif diff < 0 and diff % 3 == 0: mut_type = "in_multiple_3"
            elif diff < 0 and diff % 3 != 0: mut_type = "in_shift"
            elif diff > 0 and diff % 3 == 0: mut_type = "del_multiple_3"
            elif diff > 0 and diff % 3 != 0: mut_type = "del_shift"

    if mut_type not in ("sub", "del_multiple_3", "in_multiple_3"):
        # in case of indel, perhaps the new seq length is not a multiple
        # of three, in that scenario, N must be added to the end of the
        # nuceotide sequence.
        N_to_add = 3 - ( len(new_seq) % 3 ) # give the number of N to add
        new_seq.append(N_to_add*"N")

    # record sequences in a Seq object from Biopython package 
    varseq = Seq.Seq( "".join(new_seq), Alphabet.generic_dna)
    refseq = Seq.Seq(gn, Alphabet.generic_dna)
    
    return (indice, mut_type, isInterGen, refseq, varseq)


def launch_compare(way, ref, var):
    """
    alignment A           alignment B    
    GFDKLIMPNCVMYYA       GFDKLIMPNCVMYYA
    |||||     |||||       ||||| |
    GFDKL-----VMYYA       GFDKLVMYYA
   
    The alignment A is correct, the B is uncorrect.
    The alignment B, yes I said "B", is the one we will use.
    Using index of both string, we record AA that are identical.
    As soon as the first difference is observed, the rest of the sequence is
    considered different. Result of the function is presented after.
    
    GFDKLIMPNCVMYYA                         GFDKLIMPNCVMYYA
    ||||| |         => GFDKL//////////                ||||| => //////////VMYYA
    GFDKLVMYYA                                   GFDKLVMYYA
    
    This is the purpose of this function. The two generated "slashed" sequences
    will be used to extract the region that is matching "/" in both sequences.
    """
    com = []
    mini, maxi = min( len(ref), len(var) ), max( len(ref), len(var) )
    
    if way is "F": start, stop, a = 0, mini, 1 # F = forward
    elif way is "R": start, stop, a = -1, -(mini + 1), -1 # R = reverse
    
    for i in range(start, stop, a):
        if ref[i] == var[i]: com.append( ref[i] )
        else:
            b = i
            while b != ((maxi * a) + start) :
                com.append("/")
                b += (1*a)
            break
    return tuple(com)

def get_differences(com1, com2):
    """
Raisonnement
dans le sens 1, seul les premiers AA sont identiques, les autres sont considérés différents
dans le sens 2, seul les derniers AA sont identiques, les autres sont considérés différents
l'idée est de récupéré à la fois l'index des positions dans le sens1 ET dans le sens2
    dans le sens1, index positif de 0 à len
    dans le sens2, index négatif de -1 à -(len+1)

all_ = ((A, /, /, /, /, /, /),
        (/, /, /, /, /, B, C))
   """
    print("Get region of difference.")
    diff_position = []

    for i in range(0, len(com1)):
        # in the next line, we record START and END positions of the
        # discordance, using the two lists index numerotation,
        # in the forward way (from 0 to len),
        # in the reverse way (from -(len+1) to -1)
        if com1[i] == com2[i]: diff_position.append( (i, i-len(com1)) )
            # pour chaque poition différente, on enregistre 1 tuple
                # ce tuple contient 2 valeurs:
                    # position via l'index forward
                    # position via l'index reverse

    start, end = diff_position[0][0], diff_position[-1][-1]
    del diff_position

    c1, c2 = list(copy.deepcopy(com1)), list(copy.deepcopy(com2))
    while "/" in c1: c1.remove("/")
    while "/" in c2: c2.remove("/")
    full_common = "".join( (c1 + c2) )
    
    return full_common, start, end


def translate_cons_var_sequences(all_tpl):
    """
    """
    print("Translation of consensus and variant sequences")
    all_output_lines = [(0,"mutation_type")] # record and store all protein variations

    for tpl in all_tpl:
        indice = int( tpl[0] )
        if tpl[2] == False: # check if the sequence is not intergenic sequence
            mut_type = tpl[1]
            ref_prot = tpl[3].translate()
            var_prot = tpl[4].translate()
            if mut_type is "sub":
                # because substitution can impact group of nt, it is necessary
                # to check them all.
                # in the final file, we want a mark that match the following
                # examples:
                    # if 1 substitution => P43L
                    # if 2 substitutions => PL65FE
                        # In this case, where the size of the substitution is
                        # greater than 1, the di aminoacid sequence is replaced
                        # by another di aminoacid sequence.
                        # The number between the two diaminoacids is the
                        # position of the first AA in the protein sequence
                # so, starting from the list of all AA substitutions,
                # we generate the final mark
                # by default, we suppose that mutation are synonymous
                synonymous = True
                final_mark = "synonymous"
                ref, var, pos = [], [], []
                for i in range(0, len(ref_prot)):
                    if ref_prot[i] != var_prot[i]:
                        synonymous = False
                        ref.append(ref_prot[i])
                        var.append(var_prot[i])
                        pos.append(i+1)

                if synonymous is False:
                    final_mark = "{}_{}_{}".format(
                                             "".join(ref), pos[0], "".join(var))
                all_output_lines.append( (indice, final_mark) )

            elif mut_type in ("del_shift", "in_shift"):
                # write what is different and stop writting when a '*' is met
                ref, var, pos = [], [], []
                maxi = max( len(ref_prot), len(var_prot) )
                for i in range(0, maxi):
                    if ref_prot[i] != var_prot[i]:
                        ref.append(ref_prot[i])
                        var.append(var_prot[i])
                        pos.append(i+1)
                        if "*" in (ref_prot[i], var_prot[i]):
                            break
                final_mark = "{}_{}_{}".format(
                                             "".join(ref), pos[0], "".join(var))
                all_output_lines.append( (indice, final_mark) )

            elif mut_type in ("del_multiple_3", "in_multiple_3"):
                final_mark = launch_particular(ref_prot, var_prot)
                all_output_lines.append( (indice, final_mark) )

        else: # if the sequence is intergenic sequence
            all_output_lines.append( (indice, "INTERGEN") )

    return all_output_lines

def launch_particular(ref_prot, var_prot):

    # forward way will out something like (A,Z,E,/,/,/,/,/,/,/,/)
    com1 = launch_compare("F", ref_prot, var_prot)

    # reverse way: the function will out the same thing that forward way,
        # but the fact that we parse the tuple reversly at the end will reverse
        # the order of the tuple
    com2 = launch_compare("R", ref_prot, var_prot)[::-1]

    full_common, start, end = get_differences(com1, com2)
    # full common must match the shortest sequence, so:
    # if it matches the var_protein, then it is a deletion,
    # if it matches the ref_protein, then it is an insertion
    # if there is no match, see the else section
    if full_common == var_prot:
        # in this case we have a deletion
        # we take the deleted part of the reference sequence AND
        # AA that are directly located upstream and downstream
        # the deleted part
        if end + 1 == 0: AA_ref = ref_prot[start - 1 :]
        else: AA_ref = ref_prot[start-1 : end+1]

        AA_var = "".join( (AA_ref[0], AA_ref[-1]) )
        final_mark = "{}_{}_{}".format(AA_ref, start, AA_var)

    elif full_common == ref_prot:
        # in this case we have an insertion
        # we take the inserted part of the variant sequence AND
        # AA that are directly located upstream and downstream
        # the inserted part
        if end + 1 == 0: AA_var = var_prot[start - 1 : ]
        else: AA_var = var_prot[start - 1 : end + 1]

        AA_ref = "".join( (AA_var[0], AA_var[-1] ))
        final_mark = "{}_{}_{}".format(AA_ref, start, AA_var)

    else:
        # in this case, the insertion or deletion is a multiple of
        # 3 nt, but not in phase with the codons of the reference
        # sequence and therefore, after translation and alignment of the
        # proteins we will obtain an indel associated with a new AA
        #
        # WARNING
        # if the new AA is matching one AA present at the
        # extremity of the deleted part, then the else condition
        # will not be launched because the scenario will match a
        # scenario with indel in phase with the codons of the
        # reference sequence i.e.the if and elif sections

        # IMPORTANT
        # if the else section is True, then var_prot and full_common are
        # to launch the present function, and in this second round of the
        # function the IF or ELIF condition will be OBLIGATORY True

        final_mark = launch_particular(var_prot, full_common)
        # FUCKING REALLY IMPORTANT
        # if the var_prot and the full_common object are exchanged by each other,
        # the programme will NOT RAISE ANY ERROR, but the result WILL BE FALSE
        # I have done the test, trust ME, this is really important

    return final_mark

def write_output(variant, all_mutant):
    """
    """
    print("writing the output file")
    with open(variant, "r") as filin:
        all_lines = filin.readlines()

    new_dico = {}
    for tpl in all_mutant: new_dico[ str(tpl[0]) ] = str(tpl[1])

    new_line = []
    for i in range(0, len(all_lines)):
        lineA = all_lines[i][:-1]
        if str(i) in new_dico: lineB = new_dico[str(i)] + "\n"
        else: lineB = "parceque le gbfile est mal annoté\n"

        new_line.append( "\t".join( (lineA, lineB) ))

    new_file_name = variant + "_with_mutation_type"
    with open(new_file_name, "w") as filout:
        for line in new_line: filout.write(line)

    return True

def convert_json_file(dico):
    """
The first loop allow to parse all records in the dictionary containing the genes
names and locations. During this loop, each item is added to a new dictionnary,
gene name as a key and values are the location of this gene.
If "//" is present in the gene name, this means the region is shared by multiple
genes. In that case, the "gene" name is splited according to "//" and every
subname is recorded as a key of the new dictionnary, and the associated values
are start and end position of the shared region. If the subname is already
present in the new dictionnary, the positions are still recored.

Once all region are individually recorded, for all keys in the new dictionnary,
the start and end positions of the gene (i.e. the min and max of the value list
associated to each key) are the only value to be kept and match the full length
gene. 
    """
    print("Convert the position file.")

    new_dico = {}
    for gene in dico:
        if gene == "genomesize" or gene == "virus_id": pass
        elif "//" in gene:
            for g in gene.split("//"):
                if g not in new_dico:
                    new_dico[g] = []
                    new_dico[g].append(dico[gene][0])
                    new_dico[g].append(dico[gene][1])
                else:
                    new_dico[g].append(dico[gene][0])
                    new_dico[g].append(dico[gene][1])
        else:
            if gene not in new_dico:
                new_dico[gene] = []
                new_dico[gene].append(dico[gene][0])
                new_dico[gene].append(dico[gene][1])
            else:
                new_dico[gene].append(dico[gene][0])
                new_dico[gene].append(dico[gene][1])
 
    gepolist = []
    for gene in new_dico:
        gepolist.append( (gene, min(new_dico[gene]), max(new_dico[gene])) )

    return tuple(gepolist)

def main_script(consensus, variant, json):
    """
    """
    # first open the consensus sequence
    cons = read_consensus_sequence(consensus)

    # then open the variant file
    # the var object will contain each line of the variant file as a tuple
    # the genes object will contain a tuple of all genes concerned by variation
    var, genes = read_variant_file(variant)

    # then open the json file containing the gene position
    dico = read_jsonfile(json)

    # get the full length gene locations, all genes present in the genome
        # NB this function can be used directly in order to get the full length
        # genes for a genbank pblication
    gepo = convert_json_file(dico)

    # then extract consensus genes and convert them as variant genes
    consensus_genes = get_consensus_genes(gepo, cons, genes)

    all_data = get_mutated_genes(consensus_genes, var) 
        # return a tuple containing:
            # indice matching the variant,
            # mutation_type
            # Trus/False about the intergenic region
            # the consensus gene,
            # the variant gene.

    all_mutations = translate_cons_var_sequences(all_data)
    all_mutation_sorted = sorted(all_mutations)

    # write the final output
    # reopen the variant file and add the new data at the end of each line

    flag = write_output(variant, all_mutation_sorted)

    if flag is True: sys.exit("""Analysis of variants is over !""")
    else: sys.exit("""A problem appears while writing the output file.""")

#################
## Main Script ##
#################

if __name__ == "__main__":

    # check arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--geneposition", type = str, help = "")
    parser.add_argument("--sequence", type = str, help = "")
    parser.add_argument("--variant", type = str, help = "")
    args = parser.parse_args()

    if not args.geneposition or not args.sequence or not args.variant:
        parser.print_help()
        sys.exit()
    else:
        main_script(args.sequence, args.variant, args.geneposition)

# ~ end of script ~
