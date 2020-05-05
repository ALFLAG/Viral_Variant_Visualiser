#! usr/bin/env/python3
# -*- coding:utf-8 -*-
# AF - last modification: August 21st 2018
# Script name : convert_vcffile_to_readablefile.py
#
# Description:  using a json file and a vcf file, return a tab delimited file
#               in order to feed a R script for SNP visualization
#               and another tab-delimited file which resume the reference and
#               alternative nucleotides at one position.
#
################################################################################

## ~ start script ~ ##

#############
## Modules ##
#############

import re                                                                        
import vcf
import numpy
import json
import argparse
import sys

###############
## Functions ##
###############


def find_key(dico_json, genomeposition):
    """ This function retrieves the gene name and the start and end position of
this gene.
    """
    gene = "intergene"
    for key in dico_json:
        if key != "virus_id" and key != "genomesize":
            base_inf , base_sup = dico_json[key][0] , dico_json[key][1]
            if base_inf <= genomeposition <= base_sup:
                gene = key
                break
            else:
                continue
    return gene, base_inf, base_sup

def write_line(temp_count, pos, gene, dico):
    """ This function write the line that will be added in the final output file.
    'position\tSNP\tref\talt\tvariant_percent\tadd_ref\tadd_alt\tgene_id\tsize_point\tlseq\trseq\thomopolymer\n'
    """
    if 'INTERGEN' in gene: # to replace the names INTERGEN1 -2 -3 etc...
        gene = 'intergen'
    global threshold
    global A

    if temp_count == 0:
        line = str(pos) + "\tNA\tNA\tNA\tNA\tNA\tNA\t"+gene+"\t1\tNA\tNA\tNA\tNA\n"
    elif temp_count == 1:
        ref, alt = dico[str(pos)][0], dico[str(pos)][1]
        freq, snp = dico[str(pos)][2], dico[str(pos)][3]
        lseq, rseq = dico[str(pos)][4], dico[str(pos)][5]

        if float(freq) >= threshold:
            A+=1
            is_homo = is_homopolymer(lseq, rseq, ref, alt)
            if is_homo is True:
                line = "\t".join( (str(pos), str(snp), ref, alt, freq, ref, alt, gene, "1", str(A), lseq, rseq, "yes") ) + "\n"
#                line = str(pos) + "\t"+ str(snp) + "\t" + ref +"\t"+ alt +"\t"+ freq +"\t"+ ref+ "\t" +alt + "\t" +gene+"\t1\t"+str(A)+"\t"+lseq+"\t"+rseq+"\tyes\n"
            else:
#                line = str(pos) + "\t"+ str(snp) + "\t" + ref +"\t"+ alt +"\t"+ freq +"\t"+ ref+ "\t" +alt + "\t" +gene+"\t1\t"+str(A)+"\t"+lseq+"\t"+rseq+"\tno\n"
                line = "\t".join( (str(pos), str(snp), ref, alt, freq, ref, alt, gene, "1", str(A), lseq, rseq, "no") ) + "\n"
        else:
#            line = str(pos) + "\t"+ str(snp) + "\t" + ref +"\t"+ alt +"\t"+ freq +"\tNA\tNA\t" +gene+"\t1\tNA\tNA\tNA\tNA\n"
            line = "\t".join( (str(pos), str(snp), ref, alt, freq, "NA\tNA", gene, "1\tNA\tNA\tNA\tNA\n") )
    return line

def is_homopolymer(lseq, rseq, ref, alt):
    """This function test if the region is homopolymeric.
    Return True if the region is considerd to be homopolymeric,
    otherwise, it will return false.
    The region is considered to homopolymeric if there is at least 3 identical
    nucleotides.
    """
    if ">" not in alt: # check if alt allele is not <DEL>
        if len(ref) < len(alt):# e.g. A vs AT
            if alt[-1] == rseq[0] and alt[-1] == rseq[1]: return True
        elif len(ref) > len(alt):# e.g. AT vs A
            if ref[-1] == rseq[0] and ref[-1] == rseq[1]: return True
        else:
            if ref[0] == rseq[0] and ref[0] == lseq[-1]: return True
            else: return False

    elif len(ref) == 1:
        if lseq[-1] == ref and rseq[0] == ref:
            return True
        elif lseq[-1] == ref and lseq[-2] == ref:
            return True
        elif rseq[0] == ref and rseq[1] == ref:
            return True
        else:
            return False

##########
## Main ##
##########

## call for the command line arguments ##
parser = argparse.ArgumentParser()
parser.add_argument('--vcfs', type = str, help='the vcf file to count')
parser.add_argument('--json', type = str, help = "json gene position file")
parser.add_argument('--out', type = str, help = "the output file")
parser.add_argument('--threshold', type = str, help='the threshold you use to cut the data', default = 0.1)

args = parser.parse_args()

## check for the presence of all command line arguments ##
if not args.vcfs or not args.json or not args.out or not args.threshold:           
    parser.print_help()
    sys.exit("""\nAn error occured while entering the arguments.
Please, read the help section above.\n""")
else:
    # recover genomesize data
    threshold = float(args.threshold)
    filin = open(args.json,"r")
    dico_json = json.load(filin)       
    genomesize = dico_json["genomesize"]
    filin.close()

    # lists initializing
    genomeposition = []
    temp_counts = []

    # dataframe
    snp_positions = numpy.zeros(genomesize+1, dtype=bool)
    
    # recover the snp position inside the vcf file
    file_vcf = open(args.vcfs, "r")
    for record in vcf.Reader(file_vcf):
        snp_positions[record.POS-1] = True
        # for each position written inside the vcf file,
        # the value in the dataframe is changed from False to True
    file_vcf.close()

    # turn the False or True values from the dataframe into 0 or 1
    for n in range(0, genomesize):
        temp_counts.append(numpy.count_nonzero(snp_positions[n]))
        genomeposition.append(n+1)
    
    # look for 3 groups: SNP_position , REF , ALT
    regex1 = '[a-zA-Z0-9\._]+[\t]([0-9]+)[\t][a-zA-Z0-9\._-]+[\t]([ATCGKMSWRYBDHVN\.]+)'\
'\t(([ATCGKMSWRYBDHVN\,\.]+)||<DUP>||<DEL>||<INV>||<INS>||<FUS>)\t[0-9]+\t[A-Za-z0-9\.\;\,]+\tSAMPLE='
    # look for 1 group: Variant frequency
    regex2 = '[0-9\/\,\.]+:[0-9\/\,\.]+:[0-9\/\,\.]+:[0-9\/\,\.]+:'\
'([0-9\/\,\.]+):[0-9\/\,\.]+:[0-9\.\/\,\.%]+'

    ## regex compilation
    reg1 , reg2 = re.compile(regex1) , re.compile(regex2)
    # get the left sequence and the right sequence from the reference base
    lseq, rseq = re.compile(";LSEQ=([A-Z0]+);"), re.compile(";RSEQ=([A-Z0]+);")

    # dictionnary initialization to store data
    dico = {} # forward in the script key = snp_location

    # parse vcf file and harvest snp data
    with open(args.vcfs, "r") as file_vcf:
        for line in file_vcf:
            if "#" in line: # skip the header section
                continue
            else: # data section
                new_list = [] # add in this order = ref, alt, freq, snp_number, lseq, rseq
                print(line)
                # search for the motif and capture the wanted group
                SNP_position = (reg1.search(line)).group(1)
                ref = (reg1.search(line)).group(2) ; new_list.append(ref)
                alt = (reg1.search(line)).group(3) ; new_list.append(alt)
                freq = (reg2.search(line)).group(1) ; new_list.append(freq)
                snp = len(alt.split(",")) ; new_list.append(snp)
                leftseq = (lseq.search(line)).group(1); new_list.append(leftseq)
                rightseq = (rseq.search(line)).group(1); new_list.append(rightseq)
                dico[SNP_position] = new_list

    # Once harvested, data need to be recorded in a file


    A = 0 # this flag is used in the write_line function in order to add indices to line where variants are upper than threshold
    with open(args.out, "w") as filout:
        summary_list = []
        regex = '([0-9]+)\t1\t([A-Z\>\<]+)\t([A-Z\>\<]+)\t([0-9\.]+)\t[A-Z\>\<]+\t[A-Z\>\<]+\t([A-Z0-9a-z\-_ /]+)\t1\t([0-9]+)\t([A-Z]+\t[A-Z]+\t(no|yes))\n'
        REGEX = re.compile(regex)    
        filout.write("position\tSNP\tref\talt\tvariant_percent\tadd_ref\tadd_alt\tgene_id\tsize_point\tindice\tlseq\trseq\tisHomo\n")
        a = 0 # flag to find the first genomic region after initialization of i

        # this part is requiered to write the correct name of the genomic region
        for i in range(0,len(genomeposition)):
            pos = genomeposition[i]
            if a != 0:
                if base_inf <= genomeposition[i] <= base_sup:
                    line = write_line(temp_counts[i], pos, gene_id, dico)
                else:
                    gene_id, base_inf, base_sup = find_key(dico_json, genomeposition[i])
                    line = write_line(temp_counts[i], pos, gene_id, dico)
                filout.write(line)
            else:
                gene_id, base_inf, base_sup = find_key(dico_json, genomeposition[i])
                line = write_line(temp_counts[i], pos, gene_id, dico)
                filout.write(line)
                a = 1
            
            # generate another file, which resume the variation
            #print(line)
            if REGEX.match(line):
                indice = (REGEX.search(line)).group(6)
                gene = (REGEX.search(line)).group(5)
                ref = (REGEX.search(line)).group(2)
                alt = (REGEX.search(line)).group(3)
                position = (REGEX.search(line)).group(1)
                freq = (REGEX.search(line)).group(4)
                homo = (REGEX.search(line)).group(7)
                new_line = str("\t".join( (indice, position, ref, alt, freq, gene, homo) ))
                summary_list.append(new_line + "\n")

    with open(args.out + "_summary", "w") as filout:
        filout.write("indice\tposition\tref\talt\tfreq\tgene\tlseq\trseq\tisHomo*\n")
        for line in summary_list:
            filout.write(line)
        filout.write("""
*NB: an homopolymer region is set to 'yes' if there is a succession of at least 3 identical nucleotides.
     it looks like a restrictive measure, but Ion Torrent sequencing is very bad on such region, so make sure you verify these variants.""")

## ~ end of script ~ ##
