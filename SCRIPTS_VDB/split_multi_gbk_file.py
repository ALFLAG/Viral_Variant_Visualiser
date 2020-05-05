#! /usr/bin/env python3


## Modules
import sys
import re

## Functions

def write_gbfile(AN, line):
    """This function write the output file.
    """
    with open(AN+".gbk", "w") as filout:
        filout.write(line)

def main(infile):
    """The main script will parse the multi-gbfile
    """
    with open(infile, "r") as filin:
        gbfile = []
        start = False
        count = 0
        regex = re.compile("VERSION +([A-Z0-9_\.]+)\n")

        for line in filin:
            if "LOCUS" in line:
                count += 1
                start = True
                gbfile.append(line)
            elif "VERSION" in line:
                gbfile.append(line)
                AN = (regex.search(line)).group(1)
            elif "//\n" in line:
                gbfile.append(line)
                full = "".join(gbfile)
                write_gbfile(AN, full)
                gbfile = []
            elif start is True:
                gbfile.append(line)
        
    print( "the file {0} has been splited in {1} file(s)".format(infile, count) )

#################
## Main script ##
#################

if __name__ == "__main__":
    infile = sys.argv[1]
    main(infile)
