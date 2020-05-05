#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# AF, August 13th 2018
# Status:
#    step:
#
# Description:    the present script is the first part of the pipeline edited with
#                Snakemake (SnakefileConsensus+SNP).
#                To lauch this script, you first need to complete the
#                "configuration_file.txt" you will find in the same directory
#                than this script.
#
################################################################################

# ~ start of script ~

#############
## Modules ##
#############


import sys
import argparse
import datetime
import os
import json
import subprocess


#################
## Main Script ##
#################


## arguments definition ##
parser = argparse.ArgumentParser("\n\nThis is the first part of Me.\n"\
"It is supposed to be used to detect genomic variants in a global viral\n"\
"population. Please, only use it to treat one project at a time,\n"\
"but inside this project, you can treat an unlimited number of samples\n"\
"(actually, the only limit is the computing capacity).\n")

parser.add_argument("-i", type = str,\
    help = "Enter the readfile name (no spaces). If you do not know it, "\
"please ask to the sequencing plateform. If you have multiple samples, "\
"then you have multiple readfiles. When you enter the readfiles names in "\
"a comma separated way, please make sure it is in the same order than "\
"the samples names. Otherwise it will be a mess for you.")

parser.add_argument("--generate_consensus", type = str,\
    help = "Do you want to perform consensus ? [yes or no]")
parser.add_argument("--consensus_file", type = str,\
    help = "Path to the consensus file")

parser.add_argument("--remove_duplicate", type = str,\
    help = "Remove duplicate reads for variant detection ? ['yes' or 'no']")
parser.add_argument("--mail", type = str,\
    help = "Enter your email address")
parser.add_argument("--threshold", type = str,\
    help = "set the threshold vlue as a float [0:1]")
parser.add_argument("-s", type = str,\
    help = "Enter the sample names (no spaces), comma separated")
parser.add_argument("--project", type = str,\
    help = "The NGS_project_code the sequencing platform communicate to you.")
parser.add_argument("--data", type = str,\
    help = 'Enter the type of data you have.\
    Is it Paired-end (pe) or Single-end (se) data?')
parser.add_argument("--db", type = str,\
    help = 'Choose the database to search viral sequence: \
        vipr = ViPR database ;\
        nt = NCBI nucleotide database ;\
        urvdb = Reference Viral DataBase ;\
        perso = personnal database, if you choose this option, \
            please read the documentation to know where to place and \
            how to build your database')


args = parser.parse_args()


## checking the different arguments, and variable attribution ##
if     not args.s \
    or not args.i \
    or not args.mail \
    or not args.project \
    or not args.generate_consensus \
    or not args.remove_duplicate \
    or not args.data \
    or not args.threshold:
    parser.print_help()
    sys.exit(    "\nOne or more arguments are missing. "\
                "Please, read the help section above\n")

if args.generate_consensus == "no" and not args.consensus_file:
    parser.print_help()
    sys.exit("You forget to provide the consensus file.")


## turn the samples name string into a list ##
if "_" in args.s:
    sys.exit("WARNING: the underscore '_' character is forbidden in sample name")

samples = args.s
sample_list = [ech for ech in samples.split(",")]

## turn the readfiles names string into a list ##
readfiles = args.i
rf_list = [name for name in readfiles.split(",")]


## check for the number of read file, is it correct compared to
# the number of samples ?
data = args.data
if data.upper() == "SE":
    if len(rf_list) != len(sample_list):
        parser.print_help()
        sys.exit("""
You selected single ends data, but the number of sample is different from the 
number of read file.
If you have multiple SE read file for one sample, please merge them and restart.
""")

else:
    if len(rf_list)/2 != len(sample_list):
        parser.print_help()
        sys.exit("""You want to use paired end data, which mean you are 
supposed to have 2 read file per sample, but it seems that is not correct.
Please, read the help section above""")

## generate the files complete path
line_list = []
if data.upper() == "SE":
    for i in range(0, len(rf_list)):
        line = "/home/aflageul/raw_data/"+args.project+"/"+str(rf_list[i])
        line_list.append(line)
elif data.upper() == "PE":
    for i in range(0,len(rf_list)-1,2):
        sub_list = []
        line1 = "/home/aflageul/raw_data/"+args.project+"/"+str(rf_list[i])
        line2 = "/home/aflageul/raw_data/"+args.project+"/"+str(rf_list[i+1])
        sub_list.append(line1)
        sub_list.append(line2)
        line_list.append(sub_list)


## create the dictionnary for the config file ##
if args.generate_consensus == "yes":
    json_dico = {'PROJECT_CODE'       :args.project,
                 'DUPLICATES'         :args.remove_duplicate,
                 'Echantillons'       :sample_list,
                 'READ_FILE'          :line_list,
                 'THRESHOLD'          :args.threshold,
                 'DATA'               :args.data,
                 'DATABASE'           :args.db,
                 'GENERATE_CONSENSUS' :"yes"
    }

elif args.generate_consensus == "no":
    consensus = "/home/aflageul/raw_data/"+args.project+"/"+args.consensus_file
    json_dico = {'DUPLICATES'         :args.remove_duplicate,
                 'PROJECT_CODE'       :args.project,
                 'Echantillons'       :sample_list,
                 'READ_FILE'          :line_list,
                 'THRESHOLD'          :args.threshold,
                 'DATA'               :args.data,
                 'DATABASE'           :args.db,
                 'CONSENSUS_FILE'     :consensus,
                 'GENERATE_CONSENSUS' :"no"
    }
# generatethe consensus full path



## create the config file to use ##
config_file = "/home/aflageul/GALAXY/ME/configuration_file.json"

filout = open(config_file, "w", encoding = "utf-8")
json.dump(json_dico, filout, ensure_ascii = False)
filout.close()


## create the config_file record ##
date = str(datetime.datetime.now())[0:10]
record_name = "{0}_{1}_{2}_configuration_file.json".format(
    ",".join(sample_list), args.project, date)
path = "/home/aflageul/GALAXY/ME/records_config_files/"

record = open(path + record_name, "w", encoding = "utf-8")
json.dump(json_dico, record, ensure_ascii = False)
record.close()


## path to report ##

user_name = ""
for letter in args.mail:
    if letter == "@": break
    elif letter == ".": user_name = user_name + "-"
    else: user_name = user_name + letter

report_path = "/home/aflageul/GALAXY/ME/REPORTS/{0}_{1}_{2}.html".format(
    user_name, ",".join(sample_list), date)


## write emails

email_success = """
mail -s 'Analysis complete' {0} << EOF
Hello, it\'s Me again.

I\'m glad to announce you that your analysis for the following sample {1} is complete.
Have a good day!
Regards,

Me
EOF""".format(args.mail, ",".join(sample_list))

email_fail = """
mail -s 'Analysis failed to launch' {0} << EOF
Hello, it\'s Me again.

I\'m afraid your analysis {1} has failed. I do not know the error, so please, read the log file and try again.
If the error persists, please refer to my god damn creator!
Nevertheless, try to have a good day!

Regards,

Me
EOF""".format(args.mail, ",".join(sample_list)) 

## launch the second part of Me, the pipeline ##

email = email_success
try:
# without the "check" argument of the subprocess.run command, in case of error,
# the sended email will be the email_success
    subprocess.run("bash -c 'source /home/aflageul/miniconda3/bin/activate "+\
"/home/aflageul/miniconda3/envs/snakes_update && "+\

"snakemake -s /home/aflageul/GALAXY/ME/snakefile_ViralVariantVisualizer_v2 "+\
"--cluster \"qsub -q long.q -V\" --cores --reason --use-conda && "+\
#"--cores --reason --use-conda && "+\

"conda deactivate'", shell = True, check = True)

except:
    email = email_fail

# send email
os.system(email)

## ~ end of script ~
