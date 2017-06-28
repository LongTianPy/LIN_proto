#!/usr/bin/python
"""
"""

# IMPORT
import os
from os import listdir
from os.path import isfile,isdir,join

# FUNCTIONS
working_dir = "/home/linproject/Workspace_playaround/bbmap_sketches"
genome_dir = "/home/linproject/Workspace/Psy_166"

def create_sketches():
    cmd = "sh /home/linproject/Projects/bbmap/sketch.sh in={0} out={1}"
    genomes = [join(genome_dir,f) for f in listdir(genome_dir) if isfile(join(genome_dir,f))]
    genomes_nopath = [f for f in listdir(genome_dir) if isfile(join(genome_dir,f))]
    for i in range(len(genomes)):
        os.system(cmd.format(genomes[i], genomes_nopath[i]+".sketch"))


# MAIN
if __name__ == '__main__':
    create_sketches()