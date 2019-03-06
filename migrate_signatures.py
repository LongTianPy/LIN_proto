#!/usr/bin/python
"""
"""

# IMPORT
from MySQLdb import Connect
import os
import shutil
from os.path import join, isfile, isdir
import pandas as pd

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2.0/all_sketches/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2.0/rep_bac/"

# FUNCTIONS
def connect_to_db():
    conn = Connect("127.0.0.1", "LINbase","Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def create_signature(genome,dest):
    cmd = "sourmash compute -o {0} {1} -k 21,31,51 -n 2000 -q > /dev/null 2>&1".format(dest, genome)
    os.system(cmd)
    return dest

def migrate(c):
    c.execute("SELECT LIN.Genome_ID, Genome.FilePath, LIN.LIN FROM Genome,LIN WHERE Genome.Genome_ID=LIN.Genome_ID")
    tmp = c.fetchall()
    genome_id = [int(i[0]) for i in tmp]
    filepath = [i[1] for i in tmp]
    lin = [i[2] for i in tmp]
    df = pd.DataFrame()
    pd["FilePath"] = filepath
    pd["LINgroup"] = [",".join(i.split(",")[:6]) for i in lin]
    df.index = genome_id
    for i in genome_id:
        this_sig = create_signature(df.loc[i,'FilePath'], sourmash_dir+'{0}.sig'.format(i))
        if not isdir(join(sourmash_dir,df.loc[i,'LINgroup'])):
            os.mkdir(join(sourmash_dir,df.loc[i,'LINgroup']))
            shutil.copy(this_sig,rep_bac_dir)
        shutil.copy(this_sig, join(sourmash_dir,df.loc[i,'LINgroup']))

if __name__ == '__main__':
    conn,c = connect_to_db()
    migrate(c)



# MAIN