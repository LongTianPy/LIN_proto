#!/usr/bin/python
"""
I pulled bacteria genomes from NCBI refseq and got 6000+ genomes.
"""

# IMPORT
import filecmp
from os import listdir
from os.path import isfile, join

# FUNCTIONS
def connect_to_db():
    conn = connect_to_db("localhost","root")
    c = conn.cursor()
    c.execute("use LINdb_RefSeq")
    return conn, c

def check_identical(new_file,file_paths,c):
    for i in range(len(file_paths)):
        if not filecmp.cmp(new_file,file_paths[i]):
            duplicated_genome = i
            break
        else:
            continue
    if duplicated_genome:
        return duplicated_genome
    else:
        return None


def load_refseq():
    refseq_path = "/home/linproject/NCBI_refseq/"
    meta_data = "/home/linproject/NCBI_refseq/assembly_summary.txt"
    genome_path = "/home/linproject/NCBI_refseq/GbBac/"
    conn, c = connect_to_db()
    c.execute("select Genome_ID, FilePath from Genome")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    FilePath = [i[1] for i in tmp]
    refseq_genomes = [join(genome_path,f) for f in listdir(genome_path) if isfile(join(genome_path,f)) and f.endswith("fna")]


# MAIN