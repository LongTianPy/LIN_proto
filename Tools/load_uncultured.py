#!/usr/bin/python
"""
"""

# IMPORT
import os
from os import listdir
from os.path import isfile, isdir, join
import sys
import shutil

# FUNCTIONS


# MAIN
if __name__ == "__main__":
    folder = sys.argv[1]
    label_file = join(folder,"labels.txt")
    target_folder = "/home/linproject/Workspace/Psy_166/"
    f = open(label_file,"r")
    lines = [i.strip().split("\t") for i in f.readlines()]
    f.close()
    for i in lines:
        fasta_file = i[0] + ".fna"
        names = i[1].split(" ")
        if names[0] == "C.":
            Genus = "_".join(names[:2])
            Species = names[2]
            Strain = names[3]
        else:
            Genus = names[0]
            Species = names[1]
            Strain = names[2]
        shutil.copy(join(folder,fasta_file),target_folder)
        Type = "No"
        GPS = "NA"
        Date = "NA"
        Accession = ""
        Country = "NA"
        Paper = "NA"
        Source = "NA"
        attributes = "^^".join(
            [Genus, Species, Strain, Type, GPS, str(Date), str(Accession), str(Country), str(Paper), str(Source)])
        cmd = "python /home/linproject/Projects/LIN_proto/workflow.py -i {0} -u 2 -s 4 -t {1}".format(fasta_file,attributes)
        os.system(cmd)
