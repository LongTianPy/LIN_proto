#!/usr/bin/python
""" Workflow for genome identification of LINbase
    No metadata needed.

    Input is a uuid that points to the folder with the same name fasta file inside.
"""

# IMPORT
import os
from os import listdir
import sys
from os.path import join, isfile
from MySQLdb import Connect
import pandas as pd
import uuid

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2/all_sketches/"
# working_dir = "/home/linproject/Workspace/Genome_identification/"
# genome_dir = "/home/linproject/Workspace/Genome_identification/uploaded_genome"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2/rep_bac/"
FastANI_cmd = "fastANI -q {0} -r {1} -o /home/linproject/Workspace/Genome_identification/{2}"
scheme = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.9925, 0.995, 0.9975, 0.9990000000000001,
          0.99925, 0.9995, 0.9997499999999999, 0.9998999999999999, 0.9999899999999999]



# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost", "LINbase","Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def create_sketch(filepath,output):
    dest = output
    cmd = "sourmash compute -o {0} {1} -k 31 -n 1000".format(dest,filepath)
    os.system(cmd)
    return dest

def compare_sketch(tmp_sig,LINgroup,output):
    if LINgroup == "rep_bac":
        dest = rep_bac_dir
    else:
        dest = sourmash_dir + LINgroup + "/"
    folder_size = len([file for file in os.listdir(dest) if isfile(join(dest,file))])
    cmd = "sourmash search {0} {1}*.sig -n {2} > {3}"
    cmd = cmd.format(tmp_sig, dest, folder_size, output)
    os.system(cmd)
    return output

def parse_result(result_file):
    f = open(result_file,"r")
    lines = [i.strip().split(" \t ") for i in f.readlines()[3:]]
    f.close()
    df = pd.DataFrame()
    if len(lines) == 0:
        return df
    else:
        mash_d = [float(i[1]) for i in lines]
        df["Jaccard_similarity"] = mash_d
        df.index = [i[2].split("/")[-1].split(".")[0] for i in lines]
        return df

# MAIN
def genome_identification(dir):
    conn, c = connect_to_db()
    working_dir = dir
    input_genome = [join(dir,f) for f in listdir(dir) if f.endswith(".fasta")][0]
    output_stamp = str(uuid.uuid4())
    tmp_sig = create_sketch(input_genome,join(working_dir,output_stamp)+".sig")
    rep_bac_result = compare_sketch(tmp_sig,"rep_bac",join(working_dir,output_stamp+".mash_out"))
    df = parse_result(rep_bac_result)
    if df.empty: # Then there's no 95% level LINgroups matched
        # Check if there is any LINgroup matched to it.
        c.execute("select distinct(LINgroup),(length(LINgroup)-length(replace(LINgroup,',',''))) as level from Description where Description_Item_ID<4 and (length(LINgroup)-length(replace(LINgroup,',','')))<=5;")
        tmp = c.fetchall()
        lingroups = {str(i[0]):int(i[1]) for i in tmp}
        current_lingroup = ''
        current_max_value = 0
        current_max_genome_id = ''
        for lingroup in lingroups:
            representative_bacterium_LIN = lingroup + "," + ",".join(['0']*(20-lingroups[lingroup]))
            c.execute("select LIN.Genome_ID,Genome.FilePath from LIN, Genome where LIN.Genome_ID=Genome.Genome_ID and LIN='{0}'".format(representative_bacterium_LIN))
            tmp = c.fetchall()
            representative_bacterium_Genome_ID = tmp[0][0]
            representative_bacterium_FilePath = tmp[0][1]
            run_FastANI = FastANI_cmd.format(input_genome,representative_bacterium_FilePath,output_stamp+'_'+str(representative_bacterium_Genome_ID))
            os.system(run_FastANI)
            with open(working_dir+output_stamp+'_'+str(representative_bacterium_Genome_ID),"r") as f:
                try:
                    line = f.readlines()[0].strip().split("\t")
                    ani = float(line[2])/100
                except:
                    ani = 0
            if ani > current_max_value:
                current_max_value = ani
                current_max_genome_id = representative_bacterium_Genome_ID
                current_lingroup = lingroup
            else:
                current_max_value = current_max_value
                current_max_genome_id = current_max_genome_id
                current_lingroup = current_lingroup
        if current_max_value != 0:
            level = len(current_lingroup.split(","))
            if scheme[level-1] <= current_max_value:
                c.execute("select LIN from LIN where Genome_ID={0}".format(current_max_genome_id))
                best_LIN = c.fetchone()[0]
                result = {"LINgroup":current_lingroup,"LIN":best_LIN,"FastANI":current_max_value}
                # print("{0}\t{1}\t{2}".format(current_lingroup,current_max_genome_id,current_max_value))
            else:
                result = {}
        else:
            result = {}
    else:
        rep_bac_Genome_ID = int(df.index[0])
        c.execute("select LIN from LIN where Genome_ID={0} and Scheme_ID=4".format(rep_bac_Genome_ID))
        rep_bac_LIN = c.fetchone()[0]
        rep_bac_LINgroup = ",".join(rep_bac_LIN.split(",")[:6])
        rep_bac_result = compare_sketch(rep_bac_LINgroup,output_stamp)
        df = parse_result(rep_bac_result)
        current_max_genome_id= int(df.index[0])
        c.execute("select FilePath from Genome where Genome_ID={0}".format(current_max_genome_id))
        current_max_filepath = c.fetchone()[0]
        run_FastANI = FastANI_cmd.format(input_genome,current_max_filepath,output_stamp+"_"+str(current_max_genome_id))
        os.system(run_FastANI)
        with open(working_dir+output_stamp+'_'+str(current_max_genome_id),"r") as f:
            line = f.readlines()[0].strip().split("\t")
        ani = float(line[2])/100
        current_max_value = ani
        c.execute("select LIN from LIN where Genome_ID={0}".format(current_max_genome_id))
        best_LIN = c.fetchone()[0]
        result = {"LINgroup":rep_bac_LINgroup,"LIN":best_LIN,"FastANI":current_max_value}
        # print("{0}\t{1}\t{2}".format(rep_bac_LINgroup,current_max_genome_id,current_max_value))
    return result