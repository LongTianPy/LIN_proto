#!/usr/bin/python
"""Had a confusion of whether or not we should look at the next position using LINgroup indexing, here I'm writing
a script to compare them.
"""

# IMPORT
import pandas as pd
from MySQLdb import Connect
import os
from os.path import isfile,isdir,join
import shutil

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use LINdb_Psy")
    return c

def fetch_genomes(cursor):
    cursor.execute("select Genome.Genome_ID, Genome.FilePath, LIN.LIN from Genome, LIN where Genome.Genome_ID=LIN.Genome_ID")
    tmp = cursor.fetchall()
    Genome_ID=[int(i[0]) for i in tmp]
    FilePath = [i[1] for i in tmp]
    LIN = [i[2] for i in tmp]
    df = pd.DataFrame()
    df["FilePath"] = FilePath
    df["LIN"] = LIN
    df.index = Genome_ID
    return df, Genome_ID

def fetch_scheme(cursor):
    cursor.execute("select Cutoff from Scheme where Scheme_ID=3")
    tmp = cursor.fetchone()
    cutoff = tmp[0].split(",")
    cutoff = [float(i) / 100 for i in cutoff]
    return cutoff

def fetch_current(df, Genome_ID, idx):
    current_genome = Genome_ID[idx]
    current_db = df.loc[Genome_ID[:idx],]
    return current_genome, current_db

def old_indexing(previous_lin,current_level,working_dir,cursor,similarity_pool_old,cutoff,current_genome,current_genome_filepath,subject_genomes,reverse_LIN_dict):
    print previous_lin
    subject_genomes = ",".join([str(each) for each in subject_genomes])
    cursor.execute(
        "select Genome_ID, LIN from LIN where LIN LIKE '{1}%' and Genome_id in ({0})".format(subject_genomes,previous_lin))
    tmp = cursor.fetchall()
    df_piece = pd.DataFrame()
    genomes_piece = [int(i[0]) for i in tmp]
    df_piece["LIN"] = [i[1].split(",") for i in tmp]
    df_piece.index = genomes_piece
    LIN_dictionary = {}
    for each_genome in genomes_piece:
        each_LIN = df_piece.get_value(each_genome,"LIN")
        each_leading_part = ",".join(each_LIN[:current_level+1])
        if each_leading_part not in LIN_dictionary:
            LIN_dictionary[each_leading_part] = {each_LIN[current_level+1]:[each_genome]}
        else:
            if each_LIN[current_level+1] not in LIN_dictionary[each_leading_part]:
                LIN_dictionary[each_leading_part][each_LIN[current_level+1]] = [each_genome]
            else:
                LIN_dictionary[each_leading_part][each_LIN[current_level+1]].append(each_genome)
    if len(set(LIN_dictionary.keys())) > 1:
        LIN_ANI_storage = {}
        LIN_ANI_max_storage = {}
        for each_LIN_dictionary_key in LIN_dictionary.keys():
            LIN_ANI_storage[each_LIN_dictionary_key] = []
            for each_next_number in LIN_dictionary[each_LIN_dictionary_key].keys():
                subject_LIN = each_LIN_dictionary_key + "," + each_next_number + "".join([",0"]*(20-1-current_level-1))
                subject_genome_ID = reverse_LIN_dict[subject_LIN]
                if str(subject_genome_ID) in similarity_pool_old:
                    similarity = similarity_pool_old[str(subject_genome_ID)]
                else:
                    cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(subject_genome_ID))
                    subject_genome_filepath = cursor.fetchone()[0]
                    shutil.copyfile(current_genome_filepath,working_dir+"{0}.fasta".format(current_genome))
                    shutil.copyfile(subject_genome_filepath,working_dir+"{0}.fasta".format(subject_genome_ID))
                    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}output/ -m ANIblastall --nocompress".format(working_dir)
                    os.system(pyani_cmd)
                    ANIb_result = pd.read_table(working_dir+"output/ANIblastall_percentage_identity.tab", header=0, index_col=0).get_value(int(current_genome),str(subject_genome_ID))
                    os.system("rm -rf {0}*".format(working_dir))
                    similarity = ANIb_result
                    similarity_pool_old[str(subject_genome_ID)] = similarity
                LIN_ANI_storage[each_LIN_dictionary_key].append(similarity)
            LIN_ANI_max_storage[each_LIN_dictionary_key] = max(LIN_ANI_storage[each_LIN_dictionary_key])
        if max(LIN_ANI_max_storage.values()) > cutoff[current_level]:
            leading_part_w_max_ANI = max(LIN_ANI_max_storage, key=LIN_ANI_max_storage.get)# The best current route
            return leading_part_w_max_ANI, current_level+1
        else:
            leading_part_w_max_ANI = ",".join(previous_lin.split(",") + ["0"] * (19 - current_level))
            current_level = 19
            return  leading_part_w_max_ANI, current_level
    else:
        return LIN_dictionary.keys()[0], current_level+1

def new_indexing(previous_lin,current_level,cursor,similarity_pool_old,similarity_pool_new,cutoff,current_genome_filepath,subject_genomes):
    subject_genomes = ",".join([str(each) for each in subject_genomes])
    cursor.execute(
        "select Genome_ID, LIN from LIN where Genome_ID in ({0}) and LIN LIKE '{1}%'".format(subject_genomes, previous_lin))
    tmp = cursor.fetchall()
    df_piece = pd.DataFrame()
    genomes_piece = [int(i[0]) for i in tmp]
    df_piece["LIN"] = [i[1].split(",") for i in tmp]
    df_piece.index = genomes_piece
    LIN_dictionary = {}
    for each_genome in genomes_piece:
        each_LIN = df_piece.get_value(each_genome,"LIN")
        each_leading_part = ",".join(each_LIN[:current_level+1])
        if each_leading_part not in LIN_dictionary:
            LIN_dictionary[each_leading_part] = each_genome
        else:
            continue
    if len(set(LIN_dictionary.keys())) > 1:
        LIN_ANI_storage = {}
        # LIN_ANI_max_storage = {}
        for each_LIN_dictionary_key in LIN_dictionary.keys():
            each_subject_genome = LIN_dictionary[each_LIN_dictionary_key]
            assert str(each_subject_genome) in similarity_pool_old
            similarity_pool_new[str(each_subject_genome)] = similarity_pool_old[each_subject_genome]
            LIN_ANI_storage[each_LIN_dictionary_key].append(similarity_pool_new[str(each_subject_genome)])
        if max(LIN_ANI_storage.values()) > cutoff[current_level]:
            leading_part_w_max_ANI = max(LIN_ANI_storage, key=LIN_ANI_storage.get)
            return leading_part_w_max_ANI, current_level + 1
        else:
            leading_part_w_max_ANI = ",".join(previous_lin.split(",") + ["0"] * (19 - current_level))
            current_level = 19
            return leading_part_w_max_ANI, current_level
    else:
        return LIN_dictionary.keys()[0], current_level+1

def main():
    times_old = []
    times_new = []
    identical = []
    genome = []
    cursor = connect_to_db()
    full_df, Genome_ID = fetch_genomes(cursor)
    cutoff = fetch_scheme(cursor)
    reverse_LIN_dict = {full_df.get_value(each_genome, "LIN"): each_genome for each_genome in Genome_ID}
    if not isdir("/home/linproject/Workspace/LINgroup_indexing_test/old"):
        os.mkdir("/home/linproject/Workspace/LINgroup_indexing_test/old")
    working_dir = "/home/linproject/Workspace/LINgroup_indexing_test/old/"
    for genome_idx in range(1, len(Genome_ID)):
        current_genome, current_db = fetch_current(full_df, Genome_ID, genome_idx)
        genome.append(current_genome)
        current_genome_filepath = full_df.get_value(current_genome, "FilePath")
        subject_genomes = current_db.index
        similarity_pool_old = {}
        similarity_pool_new = {}
        current_level_old = 0
        current_level_new = 0
        previous_lin_old = ""
        previous_lin_new = ""
        while current_level_old < 19:
            previous_lin_old, current_level_old = old_indexing(previous_lin=previous_lin_old,
                                                               current_level=current_level_old,
                                                               working_dir=working_dir,
                                                               cursor=cursor,similarity_pool_old=similarity_pool_old,
                                                               cutoff=cutoff,current_genome=current_genome,
                                                               current_genome_filepath=current_genome_filepath,
                                                               subject_genomes=subject_genomes,
                                                               reverse_LIN_dict=reverse_LIN_dict)
        while current_level_new < 19:
            previous_lin_new, current_level_new = new_indexing(previous_lin=previous_lin_new,
                                                               current_level=current_level_new,
                                                               cursor=cursor,similarity_pool_new=similarity_pool_new,
                                                               similarity_pool_old=similarity_pool_old,
                                                               cutoff=cutoff,
                                                               current_genome_filepath=current_genome_filepath,
                                                               subject_genomes=subject_genomes)
        if previous_lin_new == previous_lin_old:
            identical.append("Y")
        else:
            identical.append("N")
        times_new.append(len(similarity_pool_new))
        times_old.append(len(similarity_pool_old))
    result_df = pd.DataFrame()
    result_df["Times_old"] = times_old
    result_df["Times_new"] = times_new
    result_df["Identical"] = identical
    result_df.index = genome
    result_df.to_csv("/home/linproject/Workspace/LINgroup_indexing_test/result.csv")

# MAIN
if __name__ == "__main__":
    main()