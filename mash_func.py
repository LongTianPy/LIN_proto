#!/usr/bin/python
"""A bunch of functions related to mash or sourmash
"""

# IMPORT
import os
from Bio import SeqIO
from os.path import join

# FUNCTIONS
def write_both_strand(Genome_ID,cursor,sourmash_dir):
    cursor.execute("select FilePath from Genome where Genome_ID={0}".format(Genome_ID))
    tmp = cursor.fetchone()
    filepath = tmp[0]
    f = open(filepathm,"r")
    records = list(SeqIO.parse(f,"fasta"))
    f.close()
    f = open("{0}{1}.fasta".format(sourmash_dir,Genome_ID),"w")
    for record in records:
        f.write(">" + str(record.id) + "\n")
        f.write(str(record.seq) + "\n")
        f.write(">" + str(record.id) + "_rc\n")
        f.write(str(record.seq.reverse_complement()) + "\n")
    f.close()

def create_signature(Genome_ID,sourmash_dir):
    sourmash_cmd = "sourmash compute {0}.fasta {0}".format(join(sourmash_dir,Genome_ID))
    os.system(sourmash_cmd)



# MAIN