#!/usr/bin/python
"""This script creates the blast database.
"""

# IMPORT
from Bio import SeqIO
from MySQLdb import Connect
import os
from os.path import isfile, join, isdir

# FUNCTIONS
def concat_genomes():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("USE LINdb_RefSeq")
    c.execute("SELECT Genome_ID, FilePath FROM Genome")
    tmp = c.fetchall()
    Genome_ID = [str(i[0]) for i in tmp]
    FilePath = [i[1] for i in tmp]
    del tmp
    c.close()
    db_dir = "/var/www/html/blast/db/"
    db_fasta_file = "genome.fasta"
    if not isfile(join(db_dir,db_fasta_file)):
        out_handler = open(join(db_dir,db_fasta_file),"w")
        for i in range(len(Genome_ID)):
            f = open(FilePath[i],"r")
            records = list(SeqIO.parse(f,"fasta"))
            f.close()
            out_handler.write(">{0}\n".format(Genome_ID[i]))
            for record in records:
                out_handler.write(str(record.seq))
            out_handler.write("\n")
        out_handler.close()
    else:
        out_handler = open(join(db_dir,db_fasta_file),"a")
        db_records = list(SeqIO.parse(out_handler,"fasta"))
        num_of_records = len(db_records)
        if num_of_records < len(Genome_ID):
            for i in range(num_of_records+1,len(Genome_ID)):
                f = open(FilePath[i],"r")
                records = list(SeqIO.parse(f,"fasta"))
                f.close()
                out_handler.write(">{0}\n".format(Genome_ID[i]))
                for record in records:
                    out_handler.write(str(record.seq))
                out_handler.write("\n")
        out_handler.close()

def makeblastdb():
    db_dir = "/var/www/html/blast/db/"
    db_fasta_file = "genome.fasta"
    makeblastdb_cmd = "makeblastdb -in {0} -dbtype nucl -hash_index -logfile {1}error_log".format(join(db_dir,db_fasta_file),db_dir)
    os.system(makeblastdb_cmd)

# MAIN
if __name__ == "__main__":
    concat_genomes()
    makeblastdb()