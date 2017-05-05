#!/usr/bin/python
"""This script used Mash to determine the subset of genomes to perform ANI calculation
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os



# FUNCTIONS
def mash_sketch(mash_exec, working_dir):
    cmd = "{0} sketch -k 27 -s 10000 -o {1}sig {1}*.fasta".format(mash_exec, working_dir)
    os.system(cmd)

def mash_dist(mash_exec, working_dir, new_Genome_ID):
    cmd = "{0} dist -t {1}sig.msh {1}/{2}.fasta > ".format(mash_exec, working_dir, str(new_Genome_ID))
    os.system(cmd)

def db_connect():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("USE LINdb")
    return c

# MAIN
if __name__ == "__main__":
    mash_exec = "/home/linproject/Scripts/mash-Linux64-v1.1.1/mash"
    working_dir = "/home/linproject/Workspace/mash/"

    c = db_connect()

    c.execute("SELECT Genome_ID, FilePath from Genome")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    FilePath = [i[1] for i in tmp]
    Genome_table = pd.DataFrame()
    Genome_table["FilePath"] = FilePath
    Genome_table.index = Genome_ID

    for i in range(1,len(Genome_ID)):
        subset_Genome_ID = Genome_ID[:i+1]
        for each_subset_Genome_ID in subset_Genome_ID:
            f = open(Genome_table.get_value(each_subset_Genome_ID,"FilePath"),"r")
            records = list(SeqIO.parse(f,"fasta"))
            f.close()
            records_to_write = []
            for record in records:
                rc = SeqRecord(seq=record.seq.reverse_complement(), id=record.id + "_rc")
                records_to_write.append(record)
                records_to_write.append(rc)
            with open(working_dir+"{0}.fasta".format(str(each_subset_Genome_ID)),"w") as f:
                SeqIO.write(records_to_write,f,"fasta")
        mash_sketch(mash_exec=mash_exec, working_dir=working_dir)

