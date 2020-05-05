#! /usr/bin/env python
# -*- coding: utf8 -*-

################################################################################
################################################################################
##########                                                            ##########
##########                       My Functions                         ##########
##########                                                            ##########
################################################################################
################################################################################


def rev_comp(seq):
    """ This function takes a DNA or RNA ambigous or not sequence and return the
complementary and reverse sequence.
    """
    base_dict = {"A":"T","U":"A","C":"G","T":"A","G":"C","K":"M","M":"K","S":"S","W":"W","R":"Y","Y":"R","B":"V","D":"H","H":"D","V":"B","N":"N"}
    seq = seq.upper()
    comp_list = []
    for base in seq:
        comp_list.append(base_dict[base])
        comp = "".join(comp_list)
    return comp[::-1]


def translate_DNA(seq, mode = 1):
    """This function takes a DNA or a RNA sequence and return the corresponding
protein sequence. Only the first open reading frame is translated.
When a STOP codon is met, the translation does not stop. You can specify a mode,
i.e. the way to show the protein sequence (e.g. 'mode = 1' for a 1 letter
aminoacid writing ; 'mode = 2' for 3 letter aminoacid writing).
    """
    import sys
    if mode == 1:
        Codon_AA_Dict = {"GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",\
"UUA":"L", "UUG":"L", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",\
"CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",\
"AAA":"K", "AAG":"K",\
"AAU":"N", "AAC":"N",\
"AUG":"M",\
"GAU":"D", "GAC":"D",\
"UUU":"F", "UUC":"F",\
"UGU":"C", "UGC":"C",\
"CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",\
"CAA":"Q", "CAG":"Q",\
"UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "AGU":"S", "AGC":"S",\
"GAA":"E", "GAG":"E",\
"ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",\
"GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",\
"CAU":"H", "CAC":"H",\
"UAU":"Y", "UAC":"Y",\
"AUU":"I", "AUC":"I", "AUA":"I",\
"GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",\
"UAA":"*", "UGA":"*", "UAG":"*",\
"UGG":"W",\
"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",\
"TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",\
"CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",\
"AAA":"K", "AAG":"K",\
"AAT":"N", "AAC":"N",\
"ATG":"M",\
"GAT":"D", "GAC":"D",\
"TTT":"F", "TTC":"F",\
"TGT":"C", "TGC":"C",\
"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",\
"CAA":"Q", "CAG":"Q",\
"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",\
"GAA":"E", "GAG":"E",\
"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",\
"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",\
"TGG":"W",\
"CAT":"H", "CAC":"H",\
"TAT":"Y", "TAC":"Y",\
"ATT":"I", "ATC":"I", "ATA":"I",\
"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",\
"TAA":"*", "TGA":"*", "TAG":"*"}
    elif mode == 2:
        Codon_AA_Dict = { "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",\
"UUA":"Leu", "UUG":"Leu", "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu",\
"CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg", "AGA":"Arg", "AGG":"Arg",\
"AAA":"Lys", "AAG":"Lys",\
"AAU":"Asn", "AAC":"Asn",\
"AUG":"Met",\
"GAU":"Asp", "GAC":"Asp",\
"UUU":"Phe", "UUC":"Phe",\
"UGU":"Cys", "UGC":"Cys",\
"CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",\
"CAA":"Gln", "CAG":"Gln",\
"UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser", "AGU":"Ser", "AGC":"Ser",\
"GAA":"Glu", "GAG":"Glu",\
"ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",\
"GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly",\
"UGG":"Trp",\
"CAU":"His", "CAC":"His",\
"UAU":"Tyr", "UAC":"Tyr",\
"AUU":"Ile", "AUC":"Ile", "AUA":"Ile",\
"GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",\
"UAA":"*", "UGA":"*", "UAG":"*",\
"GCT":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",\
"TTA":"Leu", "TTG":"Leu", "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu",\
"CGT":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg", "AGA":"Arg", "AGG":"Arg",\
"AAA":"Lys", "AAG":"Lys",\
"AAT":"Asn", "AAC":"Asn",\
"ATG":"Met",\
"GAT":"Asp", "GAC":"Asp",\
"TTT":"Phe", "TTC":"Phe",\
"TGT":"Cys", "TGC":"Cys",\
"CCT":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",\
"CAA":"Gln", "CAG":"Gln",\
"TCT":"Ser", "TCC":"Ser", "TCA":"Ser", "TCG":"Ser", "AGT":"Ser", "AGC":"Ser",\
"GAA":"GlT", "GAG":"GlT",\
"ACT":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",\
"GGT":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly",\
"TGG":"Trp",\
"CAT":"His", "CAC":"His",\
"TAT":"Tyr", "TAC":"Tyr",\
"ATT":"Ile", "ATC":"Ile", "ATA":"Ile",\
"GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val",\
"TAA":"*", "TGA":"*", "TAG":"*"}
    else:
        sys.exit("""\nThe mode is not correct. Please modify it. 
'mode = 1' for a single aminoacid code, or
'mode = 2' for the 3 letters aminoacid code.
""")

    prot_seq = []
    resting_nt = ""
    size_seq = len(seq)
    if size_seq%3 == 0:
        protein = translate_prot(seq, size_seq, mode)
   else:
        rest = size_seq//3
        protein = translate_prot(seq, size_seq-rest, mode)
        
def translate_prot(sequence, size, mode):
    prot_seq = []
    for i in range(0, size, 3)
            codon = seq[i:i+3].upper()
            AA = Codon_AA_Dict[codon]
            prot_seq.append(AA)
    return prot_seq
