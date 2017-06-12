#!/usr/bin/python
"""
This scripts takes a few gene sequences and BLAST the database (Pseudomonas syringae) for the homologues sequences in
the genomes in the database.
"""

# IMPORT
from Bio import SeqIO
import sys
import os
from os import listdir
import export_db
from export_db import connect_to_db, export_table
import pandas as pd
# import multiprocessing as mp

db = "/var/www/html/blast/db/genome_new.fasta"
db_seq = SeqIO.index(db,"fasta")

# FUNCTIONS
def blast_gene(gene_seq_file):
    cmd = "blastall -p blastn -i {0} -o {0}.blast_out -d {1} -a 8 -m 8 -b 500 -e 1e-10".format(gene_seq_file,db)
    os.system(cmd)
    return gene_seq_file+".blast_out"

def parse_blast_out(gene_seq_file, blast_out):
    f = open(blast_out,"r")
    lines = [i.strip().split("\t") for i in f.readlines()]
    f.close()
    f = open(gene_seq_file)
    record = list(SeqIO.parse(f,"fasta"))
    f.close()
    seq_len = len(record[0].seq)
    df = {}
    for line in lines:
        if line[1] not in df and int(line[3])==seq_len:
            df[line[1]] = [int(line[8]), int(line[9])]
        else:
            continue
    return df

def extract_seq(df):
    out_dict = {}
    for key in df.keys():
        if df[key][0] < df[key][1]:
            out_dict[key] = str(db_seq[key].seq)[df[key][0]-1:df[key][1]]
        else:
            out_dict[key] = str(db_seq[key].seq[df[key][0]-1:df[key][1]].reverse_complement())
    return out_dict

def get_meta(db="LINdb"):
    c = connect_to_db(db)
    df = export_table(c)
    return df


def write_file(gene_seq_file,out_dict):
    meta = get_meta(db="LINdb")
    f = open(gene_seq_file+".homologues","w")
    for i in out_dict.keys():
        genus = meta.get_value(int(i),"Genus")
        species = meta.get_value(int(i),"Species")
        strain = meta.get_value(int(i),"Strain")
        f.write(">{0} {1} {2}\n".format(genus,species,strain))
        f.write(out_dict[i]+"\n")
    f.close()




# MAIN
if __name__ == '__main__':
    gene_seq_folder = sys.argv[1]
    gene_seq_files = [i for i in listdir(gene_seq_folder) if i.endswith("fas") or i.endswith("fasta")]
    for gene_seq_file in gene_seq_files:
        blast_out = blast_gene(gene_seq_file)
        parsed_blast_out = parse_blast_out(gene_seq_file,blast_out)
        out_dict = extract_seq(parsed_blast_out)
        write_file(gene_seq_file,out_dict)
