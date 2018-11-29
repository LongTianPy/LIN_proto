#!/usr/bin/python
"""
This script contains fundamental functions to realize LINgroup indexing.
LINgroup indexing is a method to select subject genome to calculate similarity

"""

# IMPORT
import os
import shutil
import pandas as pd
import LIN_Assign
# OBJECT
# class LIN_decision_tree(object):
#     def __init__(self, LIN_table):
#         self.LIN_table = LIN_table
#         self.tree()
#     def tree(self,LIN_table=None):
#         if LIN_table == None:
#             LIN_table = self.LIN_table



# FUNCTION
def go_through_LIN_table(previous_route, current_level,LIN_table,cursor,reverse_LIN_dict,New_Genome_filepath,working_dir,New_Genome_ID,User_ID,similarity_pool,cutoff):
    """

    :param previous_route: a comma separated string of the leading part of LIN decided by last step, if it's the
    first position, the previous route would be ""
    :param current_level: current position to look at to determine the extension of the route
    :param LIN_table:
    :return: extended route
    """

    this_threshold = cutoff[current_level]
    cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE LIN LIKE '{0}%'".format(previous_route))
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
                subject_LIN = each_LIN_dictionary_key + "," + each_next_number + "".join([",0"]*(20-1-current_level-1))
                subject_genome_ID = reverse_LIN_dict[subject_LIN]
                if str(subject_genome_ID) in similarity_pool:
                    similarity = similarity_pool[str(subject_genome_ID)]
                else:
                    cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(subject_genome_ID))
                    subject_genome_filepath = cursor.fetchone()[0]
                    shutil.copyfile(New_Genome_filepath,working_dir+"{0}.fasta".format(New_Genome_ID))
                    shutil.copyfile(subject_genome_filepath,working_dir+"{0}.fasta".format(subject_genome_ID))
                    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}{1}_output/ -m ANIb --nocompress".format(working_dir,User_ID)
                    os.system(pyani_cmd)
                    ANIb_result = pd.read_table(working_dir+"{0}_output/ANIb_percentage_identity.tab".format(User_ID), header=0, index_col=0).get_value(int(New_Genome_ID),str(subject_genome_ID))
                    os.system("rm -rf {0}*".format(working_dir))
                    similarity = ANIb_result
                    similarity_pool[str(subject_genome_ID)] = similarity
                LIN_ANI_storage[each_LIN_dictionary_key].append(similarity)
            LIN_ANI_max_storage[each_LIN_dictionary_key] = max(LIN_ANI_storage[each_LIN_dictionary_key])
        if max(LIN_ANI_max_storage.values()) > this_threshold:
            leading_part_w_max_ANI = max(LIN_ANI_max_storage, key=LIN_ANI_max_storage.get) # The best current route
            return leading_part_w_max_ANI, current_level+1
        else:
            leading_part_w_max_ANI = ",".join(previous_route.split(",") + ["0"] * (19 - current_level))
            current_level = 19
            return leading_part_w_max_ANI, current_level
    else:
        return LIN_dictionary.keys()[0], current_level+1


def LINgroup_indexing(cursor, New_Genome_ID, working_dir, User_ID):
    cursor.execute("SELECT Cutoff FROM Scheme WHERE Scheme_ID=4")
    tmp = cursor.fetchone()
    cutoff = tmp[0].split(",")
    cutoff = [float(i) / 100 for i in cutoff]
    cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(New_Genome_ID))
    New_Genome_filepath = cursor.fetchone()[0]
    cursor.execute("SELECT Genome_ID, LIN FROM LIN")
    tmp = cursor.fetchall()
    if len(tmp) == 0:
        new_LIN = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
        top1_Genome_ID = New_Genome_ID
        top1_similarity = 1
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
                        "1}_output/ -m ANIb --nocompress".format(
                working_dir, User_ID)
            os.system(pyani_cmd)
            ANIb_result = pd.read_table(working_dir + "{0}_output/ANIb_percentage_identity.tab".format(User_ID),
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
                print(current_level)
                previous_route, current_level = go_through_LIN_table(previous_route=previous_route,
                                                                     current_level=current_level, LIN_table=LIN_table,
                                                                     cursor=cursor, reverse_LIN_dict=reverse_LIN_dict,
                                                                     New_Genome_filepath=New_Genome_filepath,
                                                                     working_dir=working_dir,
                                                                     New_Genome_ID=New_Genome_ID, User_ID=User_ID,
                                                                     similarity_pool=similarity_pool,
                                                                     cutoff=cutoff)
                print(current_level)
            # print previous_route
            cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE LIN LIKE '{0}%'".format(previous_route))
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
                    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}{1}_output/ -m ANIb --nocompress".format(
                        working_dir, User_ID)
                    os.system(pyani_cmd)
                    ANIb_result = pd.read_table(
                        working_dir + "{0}_output/ANIb_percentage_identity.tab".format(User_ID), header=0,
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
            top1_Genome_ID = int(final_best_Genome_ID)
            top1_similarity = final_best_ANI
            print("The size of current database is " + str(len(LIN_table.index)))
            print("The number of calculations done is " + str(len(similarity_pool.keys())) + "\n\n")
    return new_LIN, top1_Genome_ID, top1_similarity

