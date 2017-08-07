#!/usr/bin/python
"""
"""

# IMPORT
import os
from Bio import SeqIO
from MySQLdb import Connect

# VARIABLES
all_sketches = "/home/linproject/Workspace/bbmap/all_sketches/"
tmp_sketches = "/home/linproject/Workspace/bbmap/tmp_sketches/"
bbmap_bin = "/home/linproject/Projects/bbmap/"

# FUNCTIONS
def change_seq_id(c,genome_id):
    c.execute("select FilePath from Genome where Genome_ID={0}".format(genome_id))
    filepath = c.fetchone()[0]
    f = open(filepath,"r")
    seqs = list(SeqIO.parse(f,"fasta"))
    f.close()
    for seq in seqs:
        seq.id = genome_id
    with open(filepath,"w") as f:
        SeqIO.write(seqs,"fasta")
    return filepath

def make_sketch(c,genome_id,filepath):
    cmd = "sh {0}sketch.sh in={0} out="
# MAIN