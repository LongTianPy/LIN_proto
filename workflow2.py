#!/usr/bin/python
"""New workflow wrapper, make the code and thoughts cleaner and straightforward.

"""

# IMPORT
import LIN_Assign
import LINgroup_indexing
import mash_indexing
import MySQLdb
from MySQLdb import Connect
import pandas as pd
import os
from os.path import isdir, isfile, join
import sys
from LoadingExternalInfo import LoadInfo
import ExtractInfo
import IntermediateResult
import logging
import logging.handlers
import argparse
from datetime import datetime
from pytz import timezone
from Bio import SeqIO
import filecmp
import uuid
import sendEmail
import shutil
import multiprocessing as mp
from functools import partial
import filecmp

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2/all_sketches/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2/rep_bac/"
sourmash_tmp = "/home/linproject/Workspace/Sourmash2/tmp_2/"
sourmash_result = "/home/linproject/Workspace/Sourmash2/result/"
bbmap_bin = "/home/linproject/Projects/bbmap/"
bbmap_dir = "/home/linproject/Workspace/bbmap/"
bbmap_rep_bac_dir = "/home/linproject/Workspace/bbmap/rep_bac/"
bbmap_all_sketches = "/home/linproject/Workspace/bbmap/all_sketches/"
bbmap_tmp = "/home/linproject/Workspace/bbmap/tmp_sketches/"
bbmap_results = "/home/linproject/Workspace/bbmap/results/"
original_folder  = '/home/linproject/Workspace/LINdb/'
tmp_folder = '/home/linproject/Workspace/tmp_upload/'
workspace_dir = '/home/linproject/Workspace/New/workspace/'

class DecisionTree(object):
    def __init__(self,ANI,cov,wkid):
        self.decide(ANI=ANI,cov=cov,wkid=wkid)
    def decide(self,ANI,cov,wkid):
        if wkid == "N/A":
            if ANI > 0.7:
                if cov > 0.145:
                    self.same_family = True
                else:
                    self.same_family = False
            else:
                self.same_family = False
        elif wkid == 0:
            if cov > 0.046:
                self.same_family = True
            else:
                self.same_family = False
        elif wkid > 0:
            if ANI > 0.7:
                if cov > 0.229:
                    self.same_family = True
                else:
                    False
            else:
                self.same_family = False


# FUNCTIONS
### Parse arguments
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="LIN platform backend"
    )
    parser.add_argument("-i", dest="new_genome", help="xxxxxx.fasta")
    parser.add_argument("-u", dest="User_ID", help="An interger")
    parser.add_argument("-s", dest="Interest_ID", help="Interest ID")
    parser.add_argument("-t", dest="Attributes", help="Attributes")
    parser.add_argument("-p", dest="privacy", help="Is it private information")
    args = parser.parse_args()
    return args

### Connect to database
def connect_to_db():
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq")
    return conn, c

### Get metadata from database
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

def load_new_metadata(c,db,args):
    c.execute("INSERT INTO Submission (User_ID, Time) VALUES ({0},'{1}')".format(args.User_ID, standardtime))
    db.commit()
    c.execute("SELECT Submission_ID FROM Submission where User_ID={0} and Time='{1}'".format(User_ID,standardtime))
    Submission_ID = int(c.fetchone()[0])
    os.system("cp {0} {1}".format(tmp_folder+args.new_genome,original_folder+args.new_genome))
    c.execute("INSERT INTO Genome (Interest_ID, Submission_ID, FilePath) VALUES ({0}, {1}, '{2}')"
              .format(args.Interest_ID, Submission_ID, original_folder+args.new_genome))
    db.commit()
    c.execute("SELECT Genome_ID FROM Genome WHERE Submission_ID={0}".format(Submission_ID))
    new_Genome_ID = int(c.fetchone()[0])
    c.execute("SELECT Attribute_IDs FROM Interest WHERE Interest_ID={0}".format(args.Interest_ID))
    tmp = c.fetchone()[0].split(",")
    Attribute_ID_list = [int(id) for id in tmp]
    Attributes = args.Attributes.split("^^")
    for i in range(len(Attribute_ID_list)):
        attributevalue = Attributes[i].replace("_"," ")
        sql = "insert into AttributeValue (Genome_ID,Interest_ID,Attribute_ID,AttributeValue,User_ID,Private) VALUES ({0},{1},{2},'{3}',{4},{5})".format(new_Genome_ID,args.Interest_ID, Attribute_ID_list[i],attributevalue,args.User_ID,args.privacy)
        # print(sql)
        c.execute(sql)
        db.commit()
    return new_Genome_ID


def create_sketch(filepath):
    dest = sourmash_tmp+"tmp.sig"
    cmd = "sourmash compute -o {0} {1} -k 31 -n 1000".format(dest,filepath)
    os.system(cmd)

def compare_sketch(LINgroup):
    if LINgroup == "rep_bac":
        dest = rep_bac_dir
    else:
        dest = sourmash_dir + LINgroup + "/"
    folder_size = len([file for file in os.listdir(dest) if isfile(join(dest,file))])
    cmd = "sourmash search {0} {1}*.sig -n {2} > {3}"
    cmd = cmd.format(sourmash_tmp+"tmp.sig", dest, folder_size, sourmash_result+"tmp_result.txt")
    os.system(cmd)
### By expectation, this returns The MinHash top hit, estimated ANI, and the Jaccard similarity
### Load the new genome into database with metadata

def parse_result():
    f = open(sourmash_result+"tmp_result.txt","r")
    lines = [i.strip().split(" \t ") for i in f.readlines()[3:]]
    f.close()
    df = pd.DataFrame()
    if len(lines) == 0:
        return df
    else:
        mash_d = [float(i[1]) for i in lines]
        df["Jaccard_similarity"] = mash_d
        df.index = [i[2].split("/")[-1].split(".")[0] for i in lines]
        # df = df[df["Jaccard_similarity"] > (df.get_value(df.index[0], "Jaccard_similarity") - 0.05)]
        return df

### Assign LIN
def LINgroup_indexing(cursor,metadata,new_genome_filepath):
    cursor.execute("SELECT Cutoff FROM Scheme WHERE Scheme_ID=4")
    cutoff = cursor.fetchone()[0].split(",")
    cutoff = [float(i) / 100 for i in cutoff]
    if metadata.empty:
        new_LIN = ",".join(["0"]*len(cutoff))
        top1_Genome_ID = None
        top1_similarity = 1
        top1_coverage = 1
        conserved_LIN = ""
    else:
        reverse_LIN_dict = {metadata.get_value(each_genome, "LIN"):each_genome for each_genome in metadata.index}
        if len(metadata.index) == 1:
            subject_genome_ID = metadata.index[0]
            subject_genome_filepath = metadata.get_value(subject_genome_ID,"FilePath")
            sub_working_dir = workspace_dir + str(subject_genome_ID) + "/"
            if not isdir(sub_working_dir):
                os.mkdir(sub_working_dir)
            shutil.copyfile(new_genome_filepath, sub_working_dir + "tmp.fasta")
            shutil.copyfile(subject_genome_filepath, sub_working_dir + "{0}.fasta".format(subject_genome_ID))
            pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                        "-i {0} -o {0}output -m ANIblastall --nocompress -f".format(sub_working_dir)
            os.system(pyani_cmd)
            ANIb_result = pd.read_table(sub_working_dir+"output/ANIblastall_percentage_identity.tab",sep="\t",header=0,
                                        index_col=0).get_value('tmp',str(subject_genome_ID))
            cov_result = pd.read_table(sub_working_dir+"output/ANIblastall_alignment_coverage.tab",sep="\t",header=0,
                                        index_col=0).get_value('tmp',str(subject_genome_ID))
            os.system("rm -rf {0}".format(sub_working_dir))
            if isdir(sub_working_dir):
                os.system("rmdir {0}".format(sub_working_dir))
            predict = DecisionTree(ANI=ANIb_result, cov=cov_result, wkid="N/A")
            if predict.same_family:
                top1_similarity = ANIb_result
                top1_coverage = cov_result
                top1_Genome_ID = subject_genome_ID
                new_LIN_object = LIN_Assign.getLIN(Genome_ID=top1_Genome_ID,Scheme_ID=4,similarity=top1_similarity,
                                                   c=cursor)
                new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object,c=cursor).new_LIN
            else:
                top1_similarity = ANIb_result
                top1_coverage = cov_result
                top1_Genome_ID = subject_genome_ID
                new_LIN = "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
                conserved_LIN = ""
        else:
            similarity_pool = pd.DataFrame()
            previous_route = ""
            current_level = 0
            while current_level < 19:
                previous_route,current_level = go_through_LIN_table(previous_route,current_level,cursor,
                                                                    reverse_LIN_dict,new_genome_filepath,
                                                                    workspace_dir,similarity_pool,cutoff)
            if previous_route == ",".join(["0"]*20):
                top1_similarity = 0.1
                top1_coverage = 0.1
                top1_Genome_ID = 1
                conserved_LIN = ""
                new_LIN_object = LIN_Assign.getLIN(Genome_ID=top1_Genome_ID,Scheme_ID=4,similarity=top1_similarity,
                                                   c=cursor)
                new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object,c=cursor).new_LIN
            else:
                cursor.execute(
                    "SELECT Genome_ID, LIN FROM LIN WHERE Scheme_ID=4 and LIN LIKE '{0}%'".format(previous_route))
                tmp = cursor.fetchall()
                final_candidate_LIN_table = pd.DataFrame()
                final_candidate_LIN_table["LIN"] = [i[1] for i in tmp]
                final_candidate_LIN_table.index = [int(i[0]) for i in tmp]
                LIN_ANI_storage = pd.DataFrame()
                def calculate_ANI(each_final_candidate):
                    if each_final_candidate not in similarity_pool.index:
                        subject_genome_ID = each_final_candidate
                        cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(subject_genome_ID))
                        subject_genome_filepath = cursor.fetchone()[0]
                        sub_working_dir = workspace_dir + str(subject_genome_ID) + "/"
                        if not isdir(sub_working_dir):
                            os.mkdir(sub_working_dir)
                        shutil.copyfile(new_genome_filepath, sub_working_dir + "tmp.fasta")
                        shutil.copyfile(subject_genome_filepath,
                                        sub_working_dir + "{0}.fasta".format(subject_genome_ID))
                        pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                                    "-i {0} -o {0}output -m ANIblastall --nocompress -f".format(sub_working_dir)
                        os.system(pyani_cmd)
                        ANIb_result = pd.read_table(sub_working_dir + "output/ANIblastall_percentage_identity.tab",
                                                    sep="\t", header=0,
                                                    index_col=0).get_value('tmp', str(subject_genome_ID))
                        cov_result = pd.read_table(sub_working_dir + "output/ANIblastall_alignment_coverage.tab",
                                                   sep="\t", header=0,
                                                   index_col=0).get_value('tmp', str(subject_genome_ID))
                        os.system("rm -rf {0}".format(sub_working_dir))
                        if isdir(sub_working_dir):
                            os.system("rmdir {0}".format(sub_working_dir))
                        predict = DecisionTree(ANI=ANIb_result, cov=cov_result, wkid="N/A")
                        sub_df = pd.DataFrame(0, index=[subject_genome_ID], columns=["ANI", "Coverage", "Same_family"])
                        sub_df.loc[subject_genome_ID, "ANI"] = ANIb_result
                        sub_df.loc[subject_genome_ID, "Coverage"] = cov_result
                        if predict.same_family:
                            sub_df.loc[subject_genome_ID, "Same_family"] = 1
                        else:
                            sub_df.loc[subject_genome_ID, "Same_family"] = 0
                        return sub_df
                    else:
                        sub_df = similarity_pool.loc[each_final_candidate,]
                        return sub_df
                # pool = mp.Pool(4)
                # results = pool.map(calculate_ANI,final_candidate_LIN_table.index)
                # for i in results:
                #     LIN_ANI_storage = LIN_ANI_storage.append(i)
                for each_final_candidate in final_candidate_LIN_table.index:
                    sub_df = calculate_ANI(each_final_candidate=each_final_candidate)
                    LIN_ANI_storage = LIN_ANI_storage.append(sub_df)
                sorted_df = LIN_ANI_storage.sort_values("ANI",ascending=False)
                final_best_Genome_ID = sorted_df.index[0]
                final_best_ANI = sorted_df.get_value(final_best_Genome_ID,"ANI")
                final_best_coverage = sorted_df.get_value(final_best_Genome_ID,"Coverage")
                new_LIN_object = LIN_Assign.getLIN(Genome_ID=final_best_Genome_ID,Scheme_ID=4,similarity=final_best_ANI,c=cursor)
                new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object,c=cursor).new_LIN
                conserved_LIN = ",".join(new_LIN_object.conserved_LIN)
                top1_Genome_ID = final_best_Genome_ID
                top1_similarity = final_best_ANI
                top1_coverage = final_best_coverage
    return new_LIN, top1_Genome_ID, top1_similarity,top1_coverage,conserved_LIN

def go_through_LIN_table(previous_route, current_level,cursor,reverse_LIN_dict,new_genome_filepath,working_dir,similarity_pool,cutoff):
    this_threshold = cutoff[current_level]
    cursor.execute("SELECT Genome_ID,LIN FROM LIN WHERE Scheme_ID=4 and LIN LIKE '{0}%'".format(previous_route))
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
            # print(current_level)
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
            LIN_ANI_storage[each_LIN_dictionary_key] = pd.DataFrame()
            def parallel_each_position(each_next_number,similarity_pool):
                subject_LIN = each_LIN_dictionary_key + "," + each_next_number + "".join([",0"]*(20-1-current_level-1))
                subject_genome_ID = reverse_LIN_dict[subject_LIN]
                if subject_genome_ID in similarity_pool.index:
                    similarity = similarity_pool = similarity_pool.get_value(subject_genome_ID,"ANI")
                else:
                    cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(subject_genome_ID))
                    subject_genome_filepath = cursor.fetchone()[0]
                    sub_working_dir = workspace_dir + str(subject_genome_ID) + "/"
                    if not isdir(sub_working_dir):
                        os.mkdir(sub_working_dir)
                    shutil.copyfile(new_genome_filepath,sub_working_dir+"tmp.fasta")
                    shutil.copyfile(subject_genome_filepath, sub_working_dir + "{0}.fasta".format(subject_genome_ID))
                    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                            "-i {0} -o {0}output -m ANIblastall --nocompress -f".format(sub_working_dir)
                    os.system(pyani_cmd)
                    ANIb_result = pd.read_table(sub_working_dir + "output/ANIblastall_percentage_identity.tab", sep="\t",
                                                header=0,
                                                index_col=0).get_value('tmp', str(subject_genome_ID))
                    cov_result = pd.read_table(sub_working_dir + "output/ANIblastall_alignment_coverage.tab", sep="\t",
                                               header=0,
                                               index_col=0).get_value('tmp', str(subject_genome_ID))
                    os.system("rm -rf {0}".format(sub_working_dir))
                    if isdir(sub_working_dir):
                        os.system("rmdir {0}".format(sub_working_dir))
                    predict = DecisionTree(ANI=ANIb_result, cov=cov_result, wkid=0)
                    sub_df = pd.DataFrame(0,index=[subject_genome_ID],columns=["ANI","Coverage","Same_family"])
                    sub_df.loc[subject_genome_ID,"ANI"] = ANIb_result
                    sub_df.loc[subject_genome_ID,"Coverage"] = cov_result
                    if predict.same_family:
                        sub_df.loc[subject_genome_ID, "Same_family"] = 1
                    else:
                        sub_df.loc[subject_genome_ID, "Same_family"] = 0
                    return sub_df
            # partial_parallel_each_position = partial(parallel_each_position,similarity_pool=similarity_pool)
            # pool = mp.Pool(4)
            # results = pool.map(partial_parallel_each_position,LIN_dictionary[each_LIN_dictionary_key].keys())
            for each_next_number in LIN_dictionary[each_LIN_dictionary_key].keys():
                sub_df = parallel_each_position(each_next_number=each_next_number,similarity_pool=similarity_pool)
                LIN_ANI_storage[each_LIN_dictionary_key] = LIN_ANI_storage[each_LIN_dictionary_key].append(sub_df)
                similarity_pool = similarity_pool.append(sub_df)
            small_df = LIN_ANI_storage[each_LIN_dictionary_key]
            each_df = small_df[small_df["Same_family"]==1]
            if each_df.empty:
                LIN_ANI_max_storage[each_LIN_dictionary_key] = 0
            else:
                LIN_ANI_max_storage[each_LIN_dictionary_key] = max(each_df["ANI"])
        if max(LIN_ANI_max_storage.values()) == 0:
            next_level = 20
            return ",".join(["0"]*20), next_level
        else:
            if max(LIN_ANI_max_storage.values()) > this_threshold:
                leading_part_w_max_ANI = max(LIN_ANI_max_storage,key=LIN_ANI_max_storage.get)
                next_level = current_level + 1
                return leading_part_w_max_ANI, next_level
            else:
                if previous_route != "":
                    leading_part_w_max_ANI = ",".join(previous_route.split(",") + ["0"] * (19 - current_level))
                else:
                    leading_part_w_max_ANI = ",".join(["0"] * (18 - current_level))
                next_level = 20
                return leading_part_w_max_ANI, next_level
    else:
        return LIN_dictionary.keys()[0], current_level+1

def update_LINgroup(Genome_ID, c, new_LIN, conn):
    c.execute("select * from Description")
    tmp = c.fetchall()
    Description_ID = [int(i[0]) for i in tmp]
    LINgroup = [i[1] for i in tmp]
    members = [i[-1] for i in tmp]
    df = pd.DataFrame()
    for i in range(len(Description_ID)):
        if new_LIN.startswith(LINgroup[i]):
            this_member = ",".join(members[i].split(",").append(Genome_ID))
            c.execute("update Description set Genome_IDs={0} where Description_ID={1}".format(this_member,
                                                                                              Description_ID[i]))
            conn.commit()
### Email


# MAIN
if __name__ == '__main__':
    args = get_parsed_args()
    new_genome = args.new_genome
    new_genome_filepath = tmp_folder + new_genome
    User_ID = int(args.User_ID)
    Interest_ID_new_genome = int(args.Interest_ID)
    Attributes = args.Attributes
    privacy = args.privacy
    if privacy == "True":
        privacy = 1
        args.privacy = 1
    else:
        privacy = 0
        args.privacy = 0
    eastern = timezone("EST")
    currenttime = eastern.localize(datetime.now())
    fmt_time_display = '%Y-%m-%d %H:%M:%S %Z%z'
    standardtime = currenttime.strftime(fmt_time_display)
    db, c = connect_to_db()
    metadata = extract_metadata(c=c)
    file_duplication = 0
    for i in metadata.index:
        if filecmp.cmp(new_genome_filepath,metadata.get_value(i,"FilePath")):
            file_duplication = 1
            SubjectGenome = int(i)
            break
    if file_duplication == 0:
        create_sketch(tmp_folder+new_genome)
        if metadata.empty:
            new_LIN = ",".join(["0"] * 20)
            top1_Genome_ID = 1
            top1_similarity = 1
            top1_coverage = 1
            conserved_LIN = ""
            SubjectGenome = top1_Genome_ID
            ANIb_result = top1_similarity
            cov_result = top1_coverage
        else:
            compare_sketch(LINgroup="rep_bac")
            df = parse_result()
            if df.empty:
                print("###########################################################")
                print("System message:")
                print("No Jaccard similarity detected, will use LINgroup indexing.")
                print("###########################################################")
                ## LINgroup indexing
                new_LIN, SubjectGenome, ANIb_result, cov_result, conserved_LIN = LINgroup_indexing(cursor=c,metadata=metadata,new_genome_filepath=new_genome_filepath)
            else:
                rep_bac_Genome_ID = int(df.index[0])
                rep_bac_LIN = metadata.get_value(rep_bac_Genome_ID,"LIN")
                rep_bac_LINgroup = ",".join(rep_bac_LIN.split(",")[:6])
                compare_sketch(LINgroup=rep_bac_LINgroup)
                df = parse_result()
                if df.get_value(df.index[0],"Jaccard_similarity") == 1:
                    print("###########################################################")
                    print("System message:")
                    print("100% Jaccard similarity detected, checking duplication.")
                    print("LIN will be assigned if new genome.")
                    print("###########################################################")
                    # Same genome found
                    sub_df = df[df["Jaccard_similarity"]==1]
                    ANIb_result = 0
                    cov_result = 0
                    SubjectGenome = 0
                    # There is a table about same genome, better record it
                    #[new_LIN, ANIb_result,cov_result,conserved_LIN] = [None]*4
                    for each_subject_genome_ID in sub_df.index[:3]:
                        subject_genome_filepath = metadata.get_value(int(each_subject_genome_ID),"FilePath")
                        sub_working_dir = workspace_dir + str(each_subject_genome_ID) + "/"
                        if not isdir(sub_working_dir):
                            os.mkdir(sub_working_dir)
                        shutil.copyfile(new_genome_filepath, sub_working_dir+"tmp.fasta")
                        shutil.copyfile(subject_genome_filepath,sub_working_dir+"{0}.fasta".format(each_subject_genome_ID))
                        pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                                    "-i {0} -o {0}output -m ANIblastall --nocompress -f".format(sub_working_dir)
                        os.system(pyani_cmd)
                        this_ANIb_result = pd.read_table(sub_working_dir + "output/ANIblastall_percentage_identity.tab",
                                                    sep="\t",
                                                    header=0,
                                                    index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        this_cov_result = pd.read_table(sub_working_dir + "output/ANIblastall_alignment_coverage.tab", sep="\t",
                                                   header=0,
                                                   index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        os.system("rm -rf {0}".format(sub_working_dir))
                        if this_ANIb_result > 0.99999:
                            ANIb_result = this_ANIb_result
                            cov_result = this_cov_result
                            SubjectGenome = each_subject_genome_ID
                            break
                        else:
                            if this_ANIb_result > ANIb_result:
                                ANIb_result = this_ANIb_result
                                cov_result = this_cov_result
                                SubjectGenome = each_subject_genome_ID
                            else:
                                continue
                    new_LIN_object = LIN_Assign.getLIN(Genome_ID=SubjectGenome, Scheme_ID=4,
                                                       similarity=ANIb_result, c=c)
                    new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
                    conserved_LIN = ",".join(new_LIN_object.conserved_LIN)
                else:
                    print("###########################################################")
                    print("System message:")
                    print("Jaccard similarity detected, calculating ANIs.")
                    print("LIN will be assigned.")
                    print("###########################################################")
                    for each_subject_genome_ID in df.index[:3]:
                        sub_working_dir = workspace_dir + str(each_subject_genome_ID) + "/"
                        if not isdir(sub_working_dir):
                            os.mkdir(sub_working_dir)
                        subject_genome_filepath = metadata.get_value(int(each_subject_genome_ID), "FilePath")
                        shutil.copyfile(new_genome_filepath, sub_working_dir + "tmp.fasta")
                        shutil.copyfile(subject_genome_filepath, sub_working_dir + "{0}.fasta".format(each_subject_genome_ID))
                        pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                                    "-i {0} -o {0}output -m ANIblastall --nocompress -f".format(sub_working_dir)
                        os.system(pyani_cmd)
                        ANIb_result = pd.read_table(sub_working_dir + "output/ANIblastall_percentage_identity.tab", sep="\t",
                                                    header=0,
                                                    index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        cov_result = pd.read_table(sub_working_dir + "output/ANIblastall_alignment_coverage.tab", sep="\t",
                                                   header=0,
                                                   index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        os.system("rm -rf {0}".format(sub_working_dir))
                        if isdir(sub_working_dir):
                            os.system("rmdir {0}".format(sub_working_dir))
                        predict = DecisionTree(ANI=ANIb_result, cov=cov_result, wkid=df.get_value(each_subject_genome_ID,"Jaccard_similarity"))
                        if predict.same_family:
                            break
                        else:
                            continue
                    if predict.same_family:
                        new_LIN_object = LIN_Assign.getLIN(Genome_ID=each_subject_genome_ID,Scheme_ID=4,similarity=ANIb_result,c=c)
                        new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
                        conserved_LIN = ",".join(new_LIN_object.conserved_LIN)
                        SubjectGenome = each_subject_genome_ID
                    else:
                        new_LIN_object = LIN_Assign.getLIN(Genome_ID=each_subject_genome_ID, Scheme_ID=4, similarity=ANIb_result,
                                                           c=c)
                        new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=cursor).new_LIN
                        conserved_LIN = ""
                        SubjectGenome = each_subject_genome_ID
        c.execute("SELECT EXISTS(SELECT LIN FROM LIN WHERE LIN='{0}')".format(new_LIN))
        duplication = c.fetchone()[0] # 0 = no, 1 = yes
        if duplication == 0:
            print("###########################################################")
            print("System message:")
            print("New genome uploaded.")
            print("LIN will be assigned.")
            print("###########################################################")
            new_genome_ID = load_new_metadata(c=c,db=db,args=args)
            this_95_LINgroup = ",".join(new_LIN.split(",")[:6])
            this_95_LINgroup_path = sourmash_dir + this_95_LINgroup + "/"
            c.execute("SELECT EXISTS(SELECT LIN FROM LIN WHERE LIN LIKE '{0}%')".format(this_95_LINgroup))
            LINgroup_exists = c.fetchone()[0]
            c.execute("INSERT INTO LIN (Genome_ID, Scheme_ID,SubjectGenome,ANI,Coverage,LIN) values "
                      "({0},4,{1},{2},{3},'{4}')".format(new_genome_ID,SubjectGenome,ANIb_result,cov_result,new_LIN))
            db.commit()
            os.system("cp {0} {1}".format(sourmash_tmp+"tmp.sig",sourmash_dir+str(new_genome_ID)+".sig"))
            if LINgroup_exists == 0: # It's a new rep_bac
                os.system("cp {0} {1}".format(sourmash_tmp+"tmp.sig", rep_bac_dir + str(new_genome_ID)+".sig"))
                if not isdir(this_95_LINgroup_path):
                    os.mkdir(this_95_LINgroup_path)
                    os.system("cp {0} {1}".format(sourmash_tmp+"tmp.sig",this_95_LINgroup_path+str(new_genome_ID)+".sig"))
            else:
                os.system(
                    "cp {0} {1}".format(sourmash_tmp + "tmp.sig", this_95_LINgroup_path + str(new_genome_ID) + ".sig"))
            update_LINgroup(Genome_ID=new_genome_ID,c=c,new_LIN=new_LIN,conn=db)
            c.execute("SELECT LIN_ID FROM LIN WHERE Scheme_ID=4 AND Genome_ID={0}".format(new_genome_ID))
            LIN_ID = c.fetchone()[0]
            Job_uuid = str(uuid.uuid4())
            c.execute("INSERT INTO Job (LIN_ID,User_ID,Job_uuid,Conserved_LIN) VALUES ({0},{1},'{2}','{3}')".format(LIN_ID,
                                                                                                                    User_ID,
                                                                                                                    Job_uuid,
                                                                                                                    conserved_LIN))
            db.commit()
            c.execute("SELECT Email FROM User WHERE User_ID={0}".format(User_ID))
            user_email = c.fetchone()[0]
            email_cmd = "python /home/linproject/Projects/LIN_proto/sendEmail.py {0} Submission_result {1}".format(
                user_email, Job_uuid)
            # os.system(email_cmd)
        else:
            print("###########################################################")
            print("System message:")
            print("Duplicate submission found, recording.")
            print("###########################################################")
            c.execute("INSERT INTO Duplicated_upload (Reference_Genome_ID,Who_uploads_too) VALUES ({0},{1})".format(SubjectGenome,User_ID))
            db.commit()
            c.execute("SELECT Email FROM User WHERE User_ID={0}".format(User_ID))
            user_email = c.fetchone()[0]
            email_cmd = "python /home/linproject/Projects/LIN_proto/duplicated_upload.py {0} Submission_result {1}".format(
                    user_email,SubjectGenome)
            os.system(email_cmd)
    else:
        print("###########################################################")
        print("System message:")
        print("Duplicate submission found, recording.")
        print("###########################################################")
        c.execute(
            "INSERT INTO Duplicated_upload (Reference_Genome_ID,Who_uploads_too) VALUES ({0},{1})".format(SubjectGenome,
                                                                                                          User_ID))
        db.commit()
        c.execute("SELECT Email FROM User WHERE User_ID={0}".format(User_ID))
        user_email = c.fetchone()[0]
        email_cmd = "python /home/linproject/Projects/LIN_proto/duplicated_upload.py {0} Submission_result {1}".format(
                user_email, SubjectGenome)
        os.system(email_cmd)
    os.system("rm {0}".format(new_genome_filepath))
    c.close()
    db.close()