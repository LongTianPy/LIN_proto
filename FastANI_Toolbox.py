#!/usr/bin/python

# IMPORT
from Bio import SeqIO
import os
from os import listdir
from os.path import isfile, isdir, join

# FUNCTION
def split_query_seq(query_seq, frag_size=1020):
    """
    input: object from concatenated query genome string
    """


    
def concate_reference_files(top10, cursor):
    concat_ref_file = open('concat_ref.fasta','w')
    for each_genome_ID in top10:
        concat_ref_file.write(">{0}\n".format(each_genome_ID))
        cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(each_genome_ID))
        tmp = c.fetchone()[0]
        each_ref = open(tmp,"r")
        records = list(SeqIO.parse(each_ref,"fasta"))
        each_ref.close()
        seq = ""
        for record in records:
            seq_modified = str(record.seq).replace("N","")
            seq = seq + seq_modified
        concat_ref_file.write(seq + "\n")
    concat_ref_file.close()

def makeblastdb():
    cmd = "makeblastdb -dbtype nucl -in concat_ref.fasta -title ref_genome -out ref_genome_blastdb"
    os.system(cmd)
    
def run_blastn(fragment):

    return

def FastANI(new_Geonme_ID, top10_list, cursor):

