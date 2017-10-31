#!/usr/bin/python
"""
"""

# IMPORT
import pandas as pd
from MySQLdb import Connect
import os
import shutil
from os.path import isfile,isdir,join

sourmash_dir = "/home/linproject/Workspace/Sourmash2/all_sketches/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2/rep_bac/"

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq")
    return conn, c

def extract_metadata(c):
    metadata = pd.DataFrame()
    c.execute("SELECT Genome.Genome_ID, Genome.FilePath,LIN.LIN FROM Genome,LIN WHERE Genome.Genome_ID=LIN.Genome_ID AND LIN.Scheme_ID=4")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    FilePath = [i[1] for i in tmp]
    LIN = [i[2] for i in tmp]
    metadata["FilePath"] = FilePath
    metadata["LIN"] = LIN
    metadata.index = Genome_ID
    return metadata

conn,c = connect_to_db()
metadata = extract_metadata(c)

def extract_LINgroups(threshold):
    c.execute("select Cutoff from Scheme where Scheme_ID=4")
    tmp = c.fetchone()[0]
    cutoff = tmp.split(",")
    position = cutoff.index(threshold) + 1
    LINgroups = {}
    rep_bac = []
    for i in metadata.index:
        LIN = metadata.get_value(i,"LIN")
        LINgroup = ",".join(LIN.split(",")[:position])
        if LINgroup not in LINgroups:
            LINgroups[LINgroup] = [i]
            rep_bac.append(i)
        else:
            LINgroups[LINgroup].append(i)
    return LINgroups, rep_bac

def manage_signatures(LINgroups,rep_bac):
    for i in rep_bac:
        sig_file = "{0}.sig".format(str(i))
        shutil.copyfile(sourmash_dir+sig_file,rep_bac_dir+sig_file)
    for LINgroup in LINgroups.keys():
        LINgroup_path = sourmash_dir+LINgroup + "/"
        if not isdir(LINgroup_path):
            os.mkdir(LINgroup_path)
        for each in LINgroups[LINgroup]:
            sig_file = "{0}.sig".format(str(each))
            shutil.copyfile(sourmash_dir+sig_file,LINgroup_path+sig_file)

if __name__ == '__main__':
    LINgroups,rep_bac = extract_LINgroups("95")


# MAIN