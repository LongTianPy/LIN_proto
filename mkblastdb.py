#!/usr/bin/python
"""This script creates the blast database.
"""

# IMPORT
from Bio import SeqIO
from MySQLdb import Connect
import os
from os.path import isfile, join, isdir
import sys

# FUNCTIONS
def concat_genomes(db):
    conn = Connect("127.0.0.1","LINbase","Latham@537")
    c = conn.cursor()
    c.execute("USE {0}".format(db))
    c.execute("SELECT Genome.Genome_ID, Genome.FilePath, LIN.LIN FROM Genome,LIN WHERE Genome.Genome_ID=LIN.Genome_ID AND LIN.Scheme_ID=4")
    tmp = c.fetchall()
    Genome_ID = [str(i[0]) for i in tmp]
    FilePath = [i[1] for i in tmp]
    del tmp
    c.close()
    db_dir = "/home/linproject/Workspace/BLAST/db"
    db_fasta_file = "{0}.fasta".format(db)
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
        out_handler = open(join(db_dir,db_fasta_file),"a+")
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
    return db_fasta_file

def makeblastdb(db_fasta_file,db):
    db_dir = "/home/linproject/Workspace/BLAST/db/"
    makeblastdb_cmd = "makeblastdb -in {0} -dbtype nucl -logfile {1}error_log -out {2}".format(join(db_dir,db_fasta_file),db_dir,db_dir+db)
    os.system(makeblastdb_cmd)

# MAIN
if __name__ == "__main__":
    db = sys.argv[1]
    db_fasta_file = concat_genomes(db=db)
    makeblastdb(db_fasta_file=db_fasta_file,db=db)