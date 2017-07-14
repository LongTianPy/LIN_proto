#!/usr/bin/python
"""
"""

# IMPORT
import os
from os import listdir
from os.path import isfile,isdir,join
import sys

# FUNCTIONS
#working_dir = "/home/linproject/Workspace_playaround/bbmap_sketches"
#genome_dir = "/home/linproject/Workspace/Psy_166"

def create_sketches(working_dir, genome_dir):
    cmd = "sh sketch.sh in={0} out={1}"
    genomes = [join(genome_dir,f) for f in listdir(genome_dir) if isfile(join(genome_dir,f)) and f.endswith(".fasta")]
    genomes_nopath = [f for f in listdir(genome_dir) if isfile(join(genome_dir,f))]
    for i in range(len(genomes)):
        os.system(cmd.format(genomes[i], join(working_dir,genomes_nopath[i]+".sketch")))

def run_comparesketch(working_dir,result_dir):
    files = [file for file in listdir(working_dir) if isfile(join(working_dir, file)) and file.endswith(".fasta.sketch")]
    cmd = "sh comparesketch.sh in={0} *.sketch minid=0 records=4000 > {1}"
    for i in files:
        os.system(cmd.format(i,result_dir+i.split(".")[0]+".txt"))


# MAIN
if __name__ == '__main__':
    working_dir = sys.argv[1]
    genome_dir = sys.argv[2]
    result_dir = sys.argv[3]
    create_sketches(working_dir=working_dir, genome_dir=genome_dir)
    os.chdir(working_dir)
    run_comparesketch(working_dir=working_dir,result_dir=result_dir)

