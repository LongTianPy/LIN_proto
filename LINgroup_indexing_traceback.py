#!/usr/bin/python
"""
This method traces back to the LINgroup indexing process to see how it reduced
the times of calculation
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd


# FUNCTIONS
def find_current_step(previous_route, current_level, similarity_pool,c, current_genome_index,subjectlin,Genome_ID):
    sub_genome_id = ",".join([str(each_genome_id) for each_genome_id in Genome_ID[:current_genome_index]])
    c.execute("SELECT Genome_ID, LIN FROM LIN WHERE LIN LIKE '{0}%' and Genome_ID in ({1})".format(previous_route,sub_genome_id))
    tmp = c.fetchall()
    LIN_table_piece = pd.DataFrame()
    genomes_piece = [int(i[0]) for i in tmp]
    LIN_table_piece["LIN"] = [i[1].split(",") for i in tmp]
    LIN_table_piece.index = genomes_piece
    LIN_dictionary = {}
    reverse_LIN_dict = {",".join(LIN_table_piece.get_value(each_genome, "LIN")):each_genome for each_genome in genomes_piece}
    for each_genome in genomes_piece:
        each_LIN = LIN_table_piece.get_value(each_genome,"LIN")
        each_leading_part = ",".join(each_LIN[:current_level+1])
        if each_leading_part not in LIN_dictionary:
            LIN_dictionary[each_leading_part] = each_genome
        else:
            continue
    if len(set(LIN_dictionary.keys())) > 1:
        LIN_ANI_storage = {}
        LIN_ANI_max_storage = {}
        for each_LIN_dictionary_key in LIN_dictionary.keys():
            each_subject_genome = LIN_dictionary[each_LIN_dictionary_key]
            if str(each_subject_genome) in similarity_pool:
                continue
            else:
                similarity_pool[str(each_subject_genome)] = []
            # LIN_ANI_storage[each_LIN_dictionary_key] = []
            # for each_next_number in LIN_dictionary[each_LIN_dictionary_key].keys():
            #     subject_LIN = each_LIN_dictionary_key + "," + each_next_number + "".join([",0"]*(20-1-current_level-1))
            #     subject_genome_ID = reverse_LIN_dict[subject_LIN]
            #     if str(subject_genome_ID) in similarity_pool:
            #         continue
            #     else:
            #         similarity_pool[str(subject_genome_ID)] = []
        new_route = subjectlin[current_level+1]
        return new_route, current_level+1
    else:
        return subjectlin[current_level+1], current_level+1



def LINgroup_indexing_traceback():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use LINdb_Psy")
    c.execute("SELECT Genome_ID, SubjectGenome, LIN FROM LIN")
    tmp = c.fetchall()
    Genome_ID = [i[0] for i in tmp]
    SubjectGenome = [i[1] for i in tmp]
    LIN = [i[2] for i in tmp]
    df = pd.DataFrame()
    df["SubjectGenome"] = SubjectGenome
    df["LIN"] = LIN
    df.index = Genome_ID
    for current_genome_index in range(1,len(Genome_ID)):
        sub_table = df.loc[Genome_ID[:current_genome_index],]
        subjectgenome = df.get_value(Genome_ID[current_genome_index],"SubjectGenome")
        subjectlin = df.get_value(subjectgenome,"LIN").split(",")
        previous_route = ""
        current_level = 0
        similarity_pool = {}
        while current_level < 19:
            previous_route, current_level = find_current_step(previous_route=previous_route, current_level=current_level,
                                                              similarity_pool=similarity_pool, c=c,
                                                              current_genome_index=current_genome_index,
                                                              subjectlin=subjectlin,Genome_ID=Genome_ID)
        print Genome_ID[current_genome_index]
        print subjectlin
        print len(similarity_pool.keys())


# MAIN
if __name__=="__main__":
    LINgroup_indexing_traceback()