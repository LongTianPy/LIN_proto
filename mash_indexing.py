#!/usr/bin/python
"""This method combines LINgroup indexing and mash together.
    Reasons:
        For mash, it's not able to detect mash distance when ANI is lower than 90%.
        For LINgroup indexing, the efficiency largely depends on how similar is the new submission
            to the genomes in the database.
        In light of this, we need to change our current workflow.
"""

# IMPORT
import pandas as pd
import mash_func
from MySQLdb import Connect
import shutil
import os
from os.path import isdir
import LIN_Assign

# Pre-set variables
sourmash_dir = "/home/linproject/Workspace/Sourmash/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash/rep_bac/"
workspace_dir = '/home/linproject/Workspace/New/workspace/'
working_dir = '/home/linproject/Workspace/New/workspace/'

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use LINdb")
    return conn, c

def prepare_input(c):
    c.execute("select Genome_ID, FilePath from Genome")
    tmp = c.fetchall()
    genome_id = [int(i[0]) for i in tmp]
    filepath = [i[1] for i in tmp]
    df_Genome = pd.DataFrame()
    df_Genome["FilePath"] = filepath
    df_Genome.index = genome_id
    return df_Genome, genome_id

def go_through_LIN_table(previous_route, current_level,LIN_table,cursor,reverse_LIN_dict,New_Genome_filepath,working_dir,New_Genome_ID,User_ID,similarity_pool,cutoff):
    """

    :param previous_route: a comma separated string of the leading part of LIN decided by last step, if it's the
    first position, the previous route would be ""
    :param current_level: current position to look at to determine the extension of the route
    :param LIN_table:
    :return: extended route
    """

    this_threshold = cutoff[current_level]
    cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE Scheme_ID=4 and LIN LIKE '{0}%'".format(previous_route))
    tmp = cursor.fetchall()
    LIN_table_piece = pd.DataFrame()
    genomes_piece = [int(i[0]) for i in tmp]
    LIN_table_piece["LIN"] = [i[1].split(",") for i in tmp]
    LIN_table_piece.index = genomes_piece
    LIN_dictionary = {}
    for each_genome in genomes_piece:
        each_LIN = LIN_table_piece.get_value(each_genome,"LIN")
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
                subject_LIN = each_LIN_dictionary_key + "," + each_next_number + "".join([",0"]*(19-1-current_level-1))
                subject_genome_ID = reverse_LIN_dict[subject_LIN]
                if str(subject_genome_ID) in similarity_pool:
                    similarity = similarity_pool[str(subject_genome_ID)]
                else:
                    cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(subject_genome_ID))
                    subject_genome_filepath = cursor.fetchone()[0]
                    shutil.copyfile(New_Genome_filepath,working_dir+"{0}.fasta".format(New_Genome_ID))
                    shutil.copyfile(subject_genome_filepath,working_dir+"{0}.fasta".format(subject_genome_ID))
                    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}{1}_output/ -m ANIblastall --nocompress".format(working_dir,User_ID)
                    os.system(pyani_cmd)
                    ANIb_result = pd.read_table(working_dir+"{0}_output/ANIblastall_percentage_identity.tab".format(User_ID), header=0, index_col=0).get_value(int(New_Genome_ID),str(subject_genome_ID))
                    os.system("rm -rf {0}*".format(working_dir))
                    similarity = ANIb_result
                    similarity_pool[str(subject_genome_ID)] = similarity
                LIN_ANI_storage[each_LIN_dictionary_key].append(similarity)
            LIN_ANI_max_storage[each_LIN_dictionary_key] = max(LIN_ANI_storage[each_LIN_dictionary_key])
        if max(LIN_ANI_max_storage.values()) > this_threshold:
            leading_part_w_max_ANI = max(LIN_ANI_max_storage, key=LIN_ANI_max_storage.get) # The best current route
            return leading_part_w_max_ANI, current_level+1
        else:
            if previous_route != "":
                leading_part_w_max_ANI = ",".join(previous_route.split(",") + ["0"] * (19 - current_level))
            else:
                leading_part_w_max_ANI = ",".join(["0"] * (18 - current_level))
            current_level = 19
            return leading_part_w_max_ANI, current_level
    else:
        return LIN_dictionary.keys()[0], current_level+1


def LINgroup_indexing(cursor, New_Genome_ID, New_Genome_filepath , working_dir, User_ID):
    cursor.execute("SELECT Cutoff FROM Scheme WHERE Scheme_ID=4")
    tmp = cursor.fetchone()
    cutoff = tmp[0].split(",")
    cutoff = [float(i) / 100 for i in cutoff]
    cursor.execute("SELECT Genome_ID, LIN FROM LIN where Scheme_ID=4")
    tmp = cursor.fetchall()
    if len(tmp) == 0:
        new_LIN = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
        top1_Genome_ID = New_Genome_ID
        top1_similarity = 1
        conserved_LIN = ""
    else:
        LIN_table = pd.DataFrame()
        genomes = [int(i[0]) for i in tmp]
        LIN_table["LIN"] = [i[1].split(",") for i in tmp]
        LIN_table.index = genomes
        reverse_LIN_dict = {",".join(LIN_table.get_value(each_genome, "LIN")):each_genome for each_genome in genomes}
        if len(LIN_table.index) == 1:
            subject_genome_ID = LIN_table.index[0]
            cursor.execute("SELECT FilePath from Genome WHERE Genome_ID={0}".format(LIN_table.index[0]))
            Subject_Genome_filepath = cursor.fetchone()[0]
            shutil.copyfile(New_Genome_filepath, working_dir + "{0}.fasta".format(New_Genome_ID))
            shutil.copyfile(Subject_Genome_filepath, working_dir + "{0}.fasta".format(LIN_table.index[0]))
            pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}{" \
                        "1}_output/ -m ANIblastall --nocompress".format(
                working_dir, User_ID)
            os.system(pyani_cmd)
            ANIb_result = pd.read_table(working_dir + "{0}_output/ANIblastall_percentage_identity.tab".format(User_ID),
                                         header=0, index_col=0).get_value(int(New_Genome_ID),str(subject_genome_ID))
            os.system('rm -rf {0}*'.format(working_dir))
            similarity = ANIb_result
            top1_Genome_ID = LIN_table.index[0]
            top1_similarity = similarity
            new_LIN_object = LIN_Assign.getLIN(Genome_ID=top1_Genome_ID, Scheme_ID=4, similarity=top1_similarity,c=cursor,current_genome=New_Genome_ID)
            new_LIN = LIN_Assign.Assign_LIN(new_LIN_object, c=cursor).new_LIN
        else:
            similarity_pool = {}
            previous_route = "" # To initiate
            current_level = 0
            while current_level < 19:
                previous_route, current_level = go_through_LIN_table(previous_route=previous_route,
                                                                     current_level=current_level, LIN_table=LIN_table,
                                                                     cursor=cursor, reverse_LIN_dict=reverse_LIN_dict,
                                                                     New_Genome_filepath=New_Genome_filepath,
                                                                     working_dir=working_dir,
                                                                     New_Genome_ID=New_Genome_ID, User_ID=User_ID,
                                                                     similarity_pool=similarity_pool,
                                                                     cutoff=cutoff)
            # print previous_route
            cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE Scheme_ID=4 and LIN LIKE '{0}%'".format(previous_route))
            tmp = cursor.fetchall()
            final_candidate_LIN_table = pd.DataFrame()
            final_candidate_LIN_table["LIN"] = [i[1] for i in tmp]
            final_candidate_LIN_table.index = [int(i[0]) for i in tmp]
            LIN_ANI_storage = {}
            for each_final_candidate in final_candidate_LIN_table.index:
                if str(each_final_candidate) not in similarity_pool:
                    subject_genome_ID = each_final_candidate
                    cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(subject_genome_ID))
                    subject_genome_filepath = cursor.fetchone()[0]
                    shutil.copyfile(New_Genome_filepath, working_dir + "{0}.fasta".format(New_Genome_ID))
                    shutil.copyfile(subject_genome_filepath, working_dir + "{0}.fasta".format(subject_genome_ID))
                    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}{1}_output/ -m ANIblastall --nocompress".format(
                        working_dir, User_ID)
                    os.system(pyani_cmd)
                    ANIb_result = pd.read_table(
                        working_dir + "{0}_output/ANIblastall_percentage_identity.tab".format(User_ID), header=0,
                        index_col=0).get_value(int(New_Genome_ID), str(subject_genome_ID))
                    os.system("rm -rf {0}*".format(working_dir))
                    similarity = ANIb_result
                    similarity_pool[str(subject_genome_ID)] = similarity
                    LIN_ANI_storage[str(subject_genome_ID)] = similarity
                else:
                    LIN_ANI_storage[str(each_final_candidate)] = similarity_pool[str(each_final_candidate)]
            # LIN_ANI_storage = {str(each_final_candidate):similarity_pool[str(each_final_candidate)] for each_final_candidate in final_candidate_LIN_table.index}
            final_best_Genome_ID = str(max(LIN_ANI_storage,key=LIN_ANI_storage.get))
            # final_best_LIN = final_candidate_LIN_table.get_value(final_best_Genome_ID,"LIN")
            final_best_ANI = LIN_ANI_storage[final_best_Genome_ID]
            new_getLIN_object = LIN_Assign.getLIN(Genome_ID=int(final_best_Genome_ID), Scheme_ID=4, similarity=final_best_ANI,
                                              c=cursor)
            new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_getLIN_object,c=cursor,current_genome=New_Genome_ID).new_LIN
            SubjectGenome= int(final_best_Genome_ID)
            ANIb_result = final_best_ANI
        conserved_LIN = new_getLIN_object.conserved_LIN
    return new_LIN, SubjectGenome, ANIb_result,conserved_LIN

def mash_indexing(cursor, new_Genome_ID, User_ID,conn):
    """
    We will first look at if there is match from MinHash, if not, use LINgroup indexing until determine 
    the G level LINgroup
    :param cursor: 
    :param New_Genome: 
    :return: G level LINgroup
    """
    cursor.execute("select FilePath from Genome where Genome_ID={0}".format(new_Genome_ID))
    new_FilePath = cursor.fetchone()[0]
    shutil.copy(new_FilePath, sourmash_dir + "{0}.fasta".format(new_Genome_ID))
    new_SigPath = mash_func.create_signature(Genome_ID=new_Genome_ID, sourmash_dir=sourmash_dir, cursor=cursor)
    df_rep_bac = mash_func.sourmash_searching(sourmash_dir=sourmash_dir,LINgroup="rep_bac",
                                              current_sig_path=new_SigPath,current_genome=new_Genome_ID)
    if len(df_rep_bac) != 0:
        top_rep_bac = df_rep_bac.index[0]
        cursor.execute("select LIN from LIN where Genome_ID={0} and Scheme_ID=4".format(top_rep_bac))
        top_LIN_rep_bac = cursor.fetchone()[0]
        top_LINgroup_rep_bac = ",".join(top_LIN_rep_bac.split(",")[:6])
        df_LINgroup = mash_func.sourmash_searching(sourmash_dir=sourmash_dir,LINgroup=top_LINgroup_rep_bac,
                                                   current_sig_path=new_SigPath,current_genome=new_Genome_ID)
        SubjectGenomes = df_LINgroup.index
        similarity_pool = {}
        for SubjectGenome in SubjectGenomes:
            cursor.execute("select FilePath from Genome where Genome_ID={0}".format(SubjectGenome))
            SubjectGenome_FilePath = cursor.fetchone()[0]
            shutil.copy(new_FilePath,workspace_dir+"{0}.fasta".format(new_Genome_ID))
            shutil.copy(SubjectGenome_FilePath,workspace_dir+"{0}.fasta".format(SubjectGenome))
            pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}{1}_output/ " \
                        "-m ANIblastall --nocompress".format(workspace_dir, User_ID)
            os.system(pyani_cmd)
            ANIb_result = pd.read_table(workspace_dir + "{0}_output/ANIblastall_percentage_identity.tab".format(User_ID),
                                        header=0, index_col=0).get_value(int(new_Genome_ID), str(SubjectGenome))
            os.system("rm -rf {0}*".format(workspace_dir))
            similarity_pool[str(SubjectGenome)] = ANIb_result
        SubjectGenome = str(max(similarity_pool,key=similarity_pool.get))
        ANIb_result = similarity_pool[SubjectGenome]
        new_getLIN_object = LIN_Assign.getLIN(Genome_ID=SubjectGenome, Scheme_ID=4,
                                              similarity=ANIb_result,
                                              c=cursor)
        new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_getLIN_object, c=cursor, current_genome=new_Genome_ID).new_LIN
        conserved_LIN = new_getLIN_object.conserved_LIN
    else:
        new_LIN, SubjectGenome, ANIb_result,conserved_LIN = LINgroup_indexing(cursor=cursor,New_Genome_ID=new_Genome_ID,
                                                                New_Genome_filepath=new_FilePath,
                                                                working_dir=workspace_dir,User_ID=User_ID)

    cursor.execute("INSERT INTO LIN (Genome_ID, Scheme_ID, LIN, SubjectGenome, ANI) values ({0}, 4, '{1}', '{2}', {3})"
              .format(new_Genome_ID, new_LIN, SubjectGenome, ANIb_result))
    conn.commit()
    return new_LIN, SubjectGenome, ANIb_result,new_SigPath,conserved_LIN


# MAIN
if __name__ == "__main__":
    # sourmash_dir = "/home/linproject/Workspace/Sourmash/"
    conn, c = connect_to_db()
    df_Genome, Genome_ID = prepare_input(c)
    # c.execute("select Genome_ID from LIN")
    # tmp = [i[0] for i in c.fetchall()]
    # startpoint = len(tmp)
    # if startpoint == 0:
    #     new_Genome_ID = Genome_ID[0]
    #     new_FilePath = df_Genome.get_value(new_Genome_ID, "FilePath")
    #     shutil.copy(new_FilePath, sourmash_dir + "{0}.fasta".format(new_Genome_ID))
    #     new_SigPath = mash_func.create_signature(Genome_ID=new_Genome_ID, sourmash_dir=sourmash_dir, cursor=c,
    #                                              conn=conn)
    #     c.execute("insert into LIN (Genome_ID,Scheme_ID,SubjectGenome,ANI,LIN) VALUES"
    #               " ({0},3,{0},1,'0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')".format(Genome_ID[startpoint]))
    #     conn.commit()
    #     if not isdir(sourmash_dir+ "0,0,0,0,0,0,0"):
    #         os.mkdir(sourmash_dir+ "0,0,0,0,0,0,0")
    #     shutil.copy(new_SigPath,sourmash_dir+ "0,0,0,0,0,0,0/")
    #     shutil.copy(new_SigPath, sourmash_dir + "rep_bac/")
    #     startpoint += 1
    for i in range(370,len(Genome_ID)):
        new_Genome_ID = Genome_ID[i]
        # new_FilePath = df_Genome.get_value(new_Genome_ID,"FilePath")
        # shutil.copy(new_FilePath, sourmash_dir + "{0}.fasta".format(new_Genome_ID))
        # new_SigPath = mash_func.create_signature(Genome_ID=new_Genome_ID,sourmash_dir=sourmash_dir,cursor=c,conn=conn)
        new_LIN, SubjectGenome, ANIb_result,new_SigPath = mash_indexing(cursor=c,new_Genome_ID=new_Genome_ID,
                                                            User_ID=2,conn=conn)
        new_LINgroup = ",".join(new_LIN.split(",")[:6])
        if not isdir(sourmash_dir + new_LINgroup):
            os.mkdir(sourmash_dir + new_LINgroup)
            shutil.copy(new_SigPath,sourmash_dir+"rep_bac/")
        shutil.copy(new_SigPath,sourmash_dir+new_LINgroup + "/")



