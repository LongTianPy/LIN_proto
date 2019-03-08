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
import filecmp

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2.0/all_sketches/"
# working_dir = "/home/linproject/Workspace/Genome_identification/"
# genome_dir = "/home/linproject/Workspace/Genome_identification/uploaded_genome"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2.0/rep_bac/"

scheme = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.9925, 0.995, 0.9975, 0.9990000000000001,
          0.99925, 0.9995, 0.9997499999999999, 0.9998999999999999, 0.9999899999999999]



# FUNCTIONS
def connect_to_db():
    conn = Connect("127.0.0.1", "LINbase","Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def create_sketch(filepath,output):
    dest = output
    cmd = "sourmash compute -o {0} {1} -k 21,31,51 -n 2000 > /dev/null 2>&1".format(dest,filepath)
    os.system(cmd)
    return dest

def compare_sketch(tmp_sig,LINgroup,output):
    if LINgroup == "rep_bac":
        dest = rep_bac_dir
    else:
        dest = sourmash_dir + LINgroup + "/"
    folder_size = len([file for file in os.listdir(dest) if isfile(join(dest,file))])
    cmd = "sourmash search {0} {1}*.sig -n {2} -o {3} -k 31 --threshold 0.0001 -q 2> /dev/null 2>&1"
    cmd = cmd.format(tmp_sig, dest, folder_size, output)
    os.system(cmd)
    return output

def parse_result(result_file):
    df = pd.read_csv(result_file, sep=",", header=0)
    if df.empty:
        return df
    else:
        ids = []
        for each in df['filename']:
            id = int(each.split('/')[-1].split('.')[0])
            ids.append(id)
        df.index = ids
        return df

def check_belonged_LINgroups(conservevd_LIN,c):
    c.execute("select LINgroup_ID,LINgroup from LINgroup")
    tmp = c.fetchall()
    LINgroup_ID = [int(i[0]) for i in tmp]
    LINgroup = [i[1] for i in tmp]
    belongs_to = []
    for i in range(len(LINgroup_ID)):
        if conservevd_LIN.startswith(LINgroup[i]):
            belongs_to.append(LINgroup_ID[i])
    return belongs_to

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

# MAIN
def Genome_Identification(dir):
    conn, c = connect_to_db()
    metadata = extract_metadata(c)
    working_dir = dir
    FastANI_cmd = "fastANI -q {0} -r {1} -o {2} > /dev/null 2>&1"
    input_genome = [join(dir,f) for f in listdir(dir) if f.endswith(".fasta")][0]
    file_duplication = 0
    for i in metadata.index:
        if filecmp.cmp(input_genome,metadata.get_value(i,"FilePath")):
            file_duplication = 1
            SubjectGenome = int(i)
            break
    if file_duplication == 0:
        output_stamp = str(uuid.uuid4())
        tmp_sig = create_sketch(input_genome, join(working_dir, output_stamp) + ".sig")
        rep_bac_result = compare_sketch(tmp_sig, "rep_bac", join(working_dir, output_stamp + ".mash_out"))
        df = parse_result(rep_bac_result)
        if df.empty:  # Then there's no 95% level LINgroups matched
            # Check if there is any LINgroup matched to it.
            c.execute("select distinct(LINgroup),(1+length(LINgroup)-length(replace(LINgroup,',',''))) as level from "
                      "LINgroup where DescriptionType_ID<4 and (length(LINgroup)-length(replace(LINgroup,',','')))<=5 "
                      "order by LINgroup ASC;")
            tmp = c.fetchall()
            lingroups = {str(i[0]): int(i[1]) for i in tmp}
            current_lingroup = ''
            current_max_value = 0
            current_max_genome_id = ''
            for lingroup in lingroups:
                representative_bacterium_LIN = lingroup + "," + ",".join(['0'] * (20 - lingroups[lingroup]))
                c.execute(
                    "select LIN.Genome_ID,Genome.FilePath from LIN, Genome where LIN.Genome_ID=Genome.Genome_ID and LIN='{0}'".format(
                        representative_bacterium_LIN))
                tmp = c.fetchall()
                representative_bacterium_Genome_ID = tmp[0][0]
                representative_bacterium_FilePath = tmp[0][1]
                run_FastANI = FastANI_cmd.format(input_genome, representative_bacterium_FilePath,
                                                 join(working_dir,output_stamp + '_' + str(representative_bacterium_Genome_ID)))
                os.system(run_FastANI)
                with open(join(working_dir, output_stamp + '_' + str(representative_bacterium_Genome_ID)), "r") as f:
                    try:
                        line = f.readlines()[0].strip().split("\t")
                        ani = float(line[2]) / 100
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
                if scheme[level - 1] <= current_max_value:
                    c.execute("select LIN from LIN where Genome_ID={0}".format(current_max_genome_id))
                    best_LIN = c.fetchone()[0]
                    belongs_to = check_belonged_LINgroups(current_lingroup, c)
                    result = {"LINgroup"    : current_lingroup, "best LIN": best_LIN, "FastANI": current_max_value,
                              "LINgroup_IDs": belongs_to}
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
            rep_bac_result = compare_sketch(tmp_sig, rep_bac_LINgroup, output_stamp)
            df = parse_result(rep_bac_result)
            current_max_genome_id = int(df.index[0])
            c.execute("select FilePath from Genome where Genome_ID={0}".format(current_max_genome_id))
            current_max_filepath = c.fetchone()[0]
            run_FastANI = FastANI_cmd.format(input_genome, current_max_filepath,
                                             join(working_dir, output_stamp + '_' + str(current_max_genome_id)))
            os.system(run_FastANI)
            with open(join(working_dir, output_stamp + '_' + str(current_max_genome_id)), "r") as f:
                line = f.readlines()[0].strip().split("\t")
            ani = float(line[2]) / 100
            current_max_value = ani
            c.execute("select LIN from LIN where Genome_ID={0}".format(current_max_genome_id))
            best_LIN = c.fetchone()[0]
            for i in range(len(scheme)):
                if i < 19:
                    if ani > scheme[i] and ani < scheme[i + 1]:
                        LINgroup = best_LIN[:i + 1]
                    else:
                        i += 1
                else:
                    LINgroup = best_LIN
            belongs_to = check_belonged_LINgroups(LINgroup, c)
            result = {"LINgroup"    : LINgroup, "best LIN": best_LIN, "FastANI": current_max_value,
                      "LINgroup_IDs": belongs_to}
        # print("{0}\t{1}\t{2}".format(rep_bac_LINgroup,current_max_genome_id,current_max_value))
    else:
        c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
        best_LIN = c.fetchone()[0]
        belongs_to = check_belonged_LINgroups(best_LIN, c)
        result = {"best LIN": best_LIN, "LINgroup":best_LIN, "FastANI":1,"LINgroup_IDs":belongs_to}
    return result