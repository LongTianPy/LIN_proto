#!/usr/bin/python
"""
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
from Bio import SeqIO, Entrez
import filecmp
import uuid
import sendEmail
import shutil
import multiprocessing as mp
from functools import partial
from uuid import uuid4
import workflow2

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2/all_sketches/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2/rep_bac/"
sourmash_batch_tmp = "/home/linproject/Workspace/Sourmash2/batch_tmp/"
sourmash_batch_tmp_rep_bac = "/home/linproject/Workspace/Sourmash2/batch_tmp/rep_bac/"
sourmash_batch_result = "/home/linproject/Workspace/Sourmash2/batch_result/"
original_folder  = '/home/linproject/Workspace/LINdb/'
tmp_folder = '/home/linproject/Workspace/tmp_upload/'
workspace_dir = '/home/linproject/Workspace/New/workspace/'
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']
Entrez.email = "aaa@bb.cc"


# FUNCTIONS
def current_time():
    fmt_time_display = '%H:%M:%S %m/%d/%Y'
    current_time = datetime.now().strftime(fmt_time_display)
    return current_time

def connect_to_db():
    conn = Connect("localhost", "LINbase","Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def parse_FastANI_output(fastani_res):
    with open(fastani_res,"r") as f:
        lines = [i.strip().split("\t") for i in f.readlines()]
    genomes = list(set([i[0] for i in lines]))
    df = pd.DataFrame(0,index=genomes,columns=genomes)
    for i in range(len(lines)):
        idx = lines[i][0]
        col = lines[i][1]
        ani = float(lines[i][2])/100
        df.loc[idx,col] = ani
    return df

def cluster_genomes(fastani_matrix, threshold):
    clusters = {}
    cluster_number = 0
    cluster_pool = []
    for i in fastani_matrix.columns:
        cluster = list(fastani_matrix[fastani_matrix[i]>threshold].index)
        if cluster not in cluster_pool:
            clusters[cluster_number] = cluster
            cluster_number += 1
    return clusters

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

def LIN_assignment_within_LINgroup(df,new_genome_filepath,metadata,c):
    """

    :param df: parsed minhash result
    :return:
    """
    if df.get_value(df.index[0], "Jaccard_similarity") == 1:
        print("###########################################################")
        print("System message:")
        print("100% Jaccard similarity detected, checking duplication.")
        print("LIN will be assigned if new genome.")
        print("###########################################################")
        # Same genome found
        sub_df = df[df["Jaccard_similarity"] == 1]
        ANIb_result = 0
        cov_result = 0
        SubjectGenome = 0
        # There is a table about same genome, better record it
        # [new_LIN, ANIb_result,cov_result,conserved_LIN] = [None]*4
        for each_subject_genome_ID in sub_df.index[:3]:
            subject_genome_filepath = metadata.get_value(int(each_subject_genome_ID), "FilePath")
            sub_working_dir = workspace_dir + uuid.uuid4() + "/"
            if not isdir(sub_working_dir):
                os.mkdir(sub_working_dir)
            shutil.copyfile(new_genome_filepath, sub_working_dir + "tmp.fasta")
            shutil.copyfile(subject_genome_filepath, sub_working_dir + "{0}.fasta".format(each_subject_genome_ID))
            pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                        "-i {0} -o {0}output -m ANIb --nocompress -f".format(sub_working_dir)
            os.system(pyani_cmd)
            this_ANIb_result = pd.read_table(sub_working_dir + "output/ANIb_percentage_identity.tab",
                                             sep="\t",
                                             header=0,
                                             index_col=0).get_value('tmp', str(each_subject_genome_ID))
            this_cov_result = pd.read_table(sub_working_dir + "output/ANIb_alignment_coverage.tab", sep="\t",
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
                        "-i {0} -o {0}output -m ANIb --nocompress -f".format(sub_working_dir)
            os.system(pyani_cmd)
            ANIb_result = pd.read_table(sub_working_dir + "output/ANIb_percentage_identity.tab", sep="\t",
                                        header=0,
                                        index_col=0).get_value('tmp', str(each_subject_genome_ID))
            cov_result = pd.read_table(sub_working_dir + "output/ANIb_alignment_coverage.tab", sep="\t",
                                       header=0,
                                       index_col=0).get_value('tmp', str(each_subject_genome_ID))
            os.system("rm -rf {0}".format(sub_working_dir))
            if isdir(sub_working_dir):
                os.system("rmdir {0}".format(sub_working_dir))
            predict = workflow2.DecisionTree(ANI=ANIb_result, cov=cov_result,
                                   wkid=df.get_value(each_subject_genome_ID, "Jaccard_similarity"))
            if predict.same_family:
                break
            else:
                continue
        if predict.same_family:
            new_LIN_object = LIN_Assign.getLIN(Genome_ID=each_subject_genome_ID, Scheme_ID=4, similarity=ANIb_result,
                                               c=c)
            new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
            conserved_LIN = ",".join(new_LIN_object.conserved_LIN)
            SubjectGenome = each_subject_genome_ID
        else:
            new_LIN_object = LIN_Assign.getLIN(Genome_ID=each_subject_genome_ID, Scheme_ID=4, similarity=ANIb_result,
                                               c=c)
            new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
            conserved_LIN = ""
            SubjectGenome = each_subject_genome_ID
    c.execute("SELECT EXISTS(SELECT LIN FROM LIN WHERE LIN='{0}')".format(new_LIN))
    duplication = c.fetchone()[0]  # 0 = no, 1 = yes
    if duplication == 0:
        print("###########################################################")
        print("System message:")
        print("New genome uploaded.")
        print("LIN will be assigned.")
        print("###########################################################")
        c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
        best_LIN = c.fetchone()[0]
        belongs_to = workflow2.check_belonged_LINgroups(conserved_LIN, c)
        result = {"new LIN"     : new_LIN, "best LIN": best_LIN, "ANI": ANIb_result, "LINgroup": conserved_LIN,
                  "LINgroup_IDs": belongs_to, "Coverage":cov_result}
    else:
        print("###########################################################")
        print("System message:")
        print("Duplicate submission found, recording.")
        print("###########################################################")
        c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
        best_LIN = c.fetchone()[0]
        result = {"best LIN": best_LIN}
    return result



def LIN_assignment_new_LINgroup(cluster,metadata,c,result_pool,tmp_sig_folder,rep_bac_folder,result_folder,id_table):
    def compare_local_sketch(tmp_sig, LINgroup, basename):
        if LINgroup == "rep_bac":
            dest = rep_bac_folder
        else:
            dest = join(result_folder,LINgroup)
        folder_size = len([file for file in os.listdir(dest) if isfile(join(dest, file))])
        cmd = "sourmash search {0} {1}/*.sig -n {2} > {3}"
        cmd = cmd.format(tmp_sig, dest, folder_size, join(tmp_sig_folder, basename)+".minhash_out")
        os.system(cmd)
        return join(tmp_sig_folder, basename)+".minhash_out"
    rep_genome = cluster[0]
    rep_bac_sig = uuid4() + ".sig"
    rep_bac_tmp_sig = create_sketch(rep_genome,join(tmp_sig_folder,rep_bac_sig))
    id_table[rep_genome] = rep_bac_sig
    folder_size_rep_bac_folder = len([file for file in os.listdir(rep_bac_folder)])
    if folder_size_rep_bac_folder > 0:
        df = parse_result(compare_local_sketch(rep_bac_tmp_sig,"rep_bac",rep_bac_sig))
        if df.empty:
            new_LIN, SubjectGenome, ANIb_result, cov_result, conserved_LIN = workflow2.LINgroup_indexing(cursor=c,
                                                                                                         metadata=metadata,
                                                                                                         new_genome_filepath=rep_genome)
            c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
            best_LIN = c.fetchone()[0]
            belongs_to = workflow2.check_belonged_LINgroups(conserved_LIN, c)
            result = {"new LIN"     : new_LIN, "best LIN": best_LIN, "ANI": ANIb_result, "LINgroup": conserved_LIN,
                      "LINgroup_IDs": belongs_to, "Coverage": cov_result}
        else:
            result = LIN_assignment_within_LINgroup(df,rep_genome,metadata,c)
    else:
        new_LIN, SubjectGenome, ANIb_result, cov_result, conserved_LIN = workflow2.LINgroup_indexing(cursor=c, metadata=metadata,
                                                                                       new_genome_filepath=rep_genome)
        c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
        best_LIN = c.fetchone()[0]
        belongs_to = workflow2.check_belonged_LINgroups(conserved_LIN, c)
        result = {"new LIN"     : new_LIN, "best LIN": best_LIN, "ANI": ANIb_result, "LINgroup": conserved_LIN,
                  "LINgroup_IDs": belongs_to, "Coverage": cov_result}
    result_pool[rep_genome] = result
    lingroup_this_cluster = ",".join(result["new LIN"].split(",")[:6])
    # rep_bac_sig = uuid4() + ".sig"
    # id_table[rep_genome] = rep_bac_sig
    lingroup_this_cluster_folder = join(result_folder,lingroup_this_cluster)
    if not isdir(lingroup_this_cluster_folder):
        os.mkdir(lingroup_this_cluster_folder)
        shutil.copy(rep_bac_tmp_sig,rep_bac_folder)
    shutil.copy(rep_bac_tmp_sig, lingroup_this_cluster_folder)
    if len(cluster) > 1:
        for genome in cluster[1:]:
            basename_sig = uuid4() + ".sig"
            id_table[genome] = basename_sig
            each_tmp_sig = create_sketch(genome,join(tmp_sig_folder,basename_sig))
            df = parse_result(compare_local_sketch(each_tmp_sig,lingroup_this_cluster,basename_sig))
            if df.empty:
                new_LIN, SubjectGenome, ANIb_result, cov_result, conserved_LIN = workflow2.LINgroup_indexing(cursor=c,
                                                                                                             metadata=metadata,
                                                                                                             new_genome_filepath=genome)
                c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
                best_LIN = c.fetchone()[0]
                belongs_to = workflow2.check_belonged_LINgroups(conserved_LIN, c)
                result = {"new LIN"     : new_LIN, "best LIN": best_LIN, "ANI": ANIb_result, "LINgroup": conserved_LIN,
                          "LINgroup_IDs": belongs_to, "Coverage": cov_result}
            else:
                result = LIN_assignment_within_LINgroup(df, genome, metadata, c)
            result_pool[genome] = result
            shutil.copy(each_tmp_sig,lingroup_this_cluster_folder)
    return result_pool

def process_each_genome_belonging_to_existing_LINgroup(each_genome,LINgroup,meta,id_table,tmp_sig_folder):
    prefix_uuid = uuid4()
    id_table[each_genome] = prefix_uuid + ".sig"
    tmp_sig = create_sketch(each_genome,join(tmp_sig_folder,prefix_uuid) + ".sig")
    each_minhash_out = compare_sketch(tmp_sig,LINgroup,tmp_sig + ".minhash_out")
    df = parse_result(each_minhash_out)
    result = LIN_assignment_within_LINgroup(df=df,new_genome_filepath=each_genome,metadata=meta)

    return result

def initial_clustering(input_dir,output_dir,fastani_threshold):
    print("#####################################################################################")
    print("Step 1: Clustering input genomes by FastANI, current time: {0}".format(current_time()))
    print("Retrieving input genomes. Current time: {0}".format(current_time()))
    step1_folder = join(output_dir, "Step_1")
    if not os.path.isdir(step1_folder):
        os.mkdir(step1_folder)
    ref_genomes = [join(input_dir,f) for f in os.listdir(input_dir)]
    refList_file = join(step1_folder,"refList.txt")
    with open(refList_file, "w") as f:
        for genome in ref_genomes:
            f.write(genome + "\n")
    step1_fastANI_result = join(step1_folder,"1_fastANI_result.txt")
    fastANI_cmd = "fastANI --rl {0} --ql {0} -o {1}".format(refList_file, step1_fastANI_result)
    print("Calculating pairwise FastANI, this may take a while. Current time: {0}".format(current_time()))
    os.system(fastANI_cmd)
    print("Calculation finished. Clustering.")
    fastani_matrix = parse_FastANI_output(step1_fastANI_result)
    initial_clusters = cluster_genomes(fastani_matrix, fastani_threshold)
    print("Returning {0} clusters at {1}% FastANI level.\n".format(len(initial_clusters), fastani_threshold*100))
    with open(join(step1_folder, "1_initial_clusters.txt"), "w") as f:
        for each_clsuter in initial_clusters.keys():
            f.write("{0}\t{1}\n".format(str(each_clsuter), "\t".join(initial_clusters[each_clsuter])))
    return initial_clusters

def compare_init_cluster_w_existing_LINgroups(init_clusters,output_dir,c):
    duplication_pool = {}
    print("#####################################################################################")
    print("Step 2: Checking if each cluster belongs to existing 95% level LINgroups.")
    step2_folder = join(output_dir,"Step_2")
    if not os.path.isdir(step2_folder):
        os.mkdir(step2_folder)
    print("Checking with MinHash. Current time: {0}".format(current_time()))
    LINgroup_membership = {} # rep_bac:[cluster(s)]
    cluster_belonging = {} # Cluster: rep_bac
    cluster_rep_sig_dictionary = {} # Cluster_rep: uuid.sig
    print("Check if there is any file duplication with current database. Current time: {0}".format(current_time()))
    c.execute("SELECT LIN.Genome_ID, LIN.LIN, Genome.FilePath FROM Genome,LIN WHERE Genome.Genome_ID=LIN.Genome_ID")
    tmp = c.fetchall()
    meta = pd.DataFrame()
    Genome_ID = [int(i[0]) for i in tmp]
    LIN = [i[1] for i in tmp]
    FilePath = [i[2] for i in tmp]
    meta["LIN"] = LIN
    meta["FilePath"] = FilePath
    meta.index = Genome_ID
    tmp_sig_folder = join(step2_folder,"tmp_sig")
    tmp_result_folder = join(step2_folder,"tmp_result")
    if not os.path.isdir(tmp_sig_folder):
        os.mkdir(tmp_sig_folder)
    if not os.path.isdir(tmp_result_folder):
        os.mkdir(tmp_result_folder)
    for each_cluster in init_clusters:
        cluster_rep = init_clusters[each_cluster][0]
        cluster_rep_uuid = uuid4()
        cluster_rep_tmp_sig = create_sketch(cluster_rep,join(tmp_sig_folder+".sig"))
        cluster_rep_minhash_out = compare_sketch(cluster_rep_tmp_sig, "rep_bac", join(tmp_result_folder,cluster_rep_uuid+".minhash_out"))
        df= parse_result(cluster_rep_minhash_out)
        if df.empty:
            if "new" not in LINgroup_membership:
                LINgroup_membership["new"] = [each_cluster]
            else:
                LINgroup_membership["new"].append(each_cluster)
            cluster_belonging[each_cluster] = "new"
        else:
            rep_bac_Genome_ID = int(df.index[0])
            c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(rep_bac_Genome_ID))
            rep_bac_LIN = c.fetchone()[0]
            rep_bac_LINgroup = ",".join(rep_bac_LIN.split(",")[:6])
            if rep_bac_LINgroup not in LINgroup_membership:
                LINgroup_membership[rep_bac_LINgroup] = [each_cluster]
            else:
                LINgroup_membership[rep_bac_LINgroup].append(each_cluster)
            cluster_belonging[each_cluster] = rep_bac_LINgroup
        cluster_rep_sig_dictionary[cluster_rep] = cluster_rep_tmp_sig
    print("Checking if there is file duplication in each cluster. Current time: {0}".format(current_time()))
    for each_LINgroup in LINgroup_membership:
        if each_LINgroup != "new":
            meta_this_lingroup = meta[meta["LIN"].str.startswith(each_LINgroup)]
            for each_cluster in LINgroup_membership[each_LINgroup]:
                for each_file in each_cluster:
                    for i in meta_this_lingroup.index:
                        if filecmp.cmp(each_file,meta_this_lingroup.get_value(i,"FilePath")):
                            LINgroup_membership[each_LINgroup][each_cluster].remove(each_file)
                            duplication_pool[each_file] = {"best LIN":meta_this_lingroup.get_value(i,"LIN")}
    print("Returning if any cluster belongs to any 95% level LINgroup. Current time: {0}".format(current_time()))
    with open(join(step2_folder, "LINgroup_membership.txt"), "w") as f:
        for i in LINgroup_membership.keys():
            f.write("{0}\t{1}\n".format(i, "\t".join(LINgroup_membership[i])))
    os.system("rm {0}/*".format(sourmash_batch_tmp))
    return LINgroup_membership, cluster_belonging, cluster_rep_sig_dictionary,meta

def assign_LIN_to_each_cluster(LINgroup_membership,meta,output_dir):
    print("Step3: Assign LINs by cluster. ")
    step3_folder = join(output_dir, "Step_3")
    if not os.path.isdir(step3_folder):
        os.mkdir(step3_folder)
    tmp_sig_folder = join(step3_folder, "tmp_sig")
    rep_bac_folder = join(step3_folder, "rep_bac")
    result_folder = join(step3_folder, "result")
    for i in [tmp_sig_folder, rep_bac_folder, result_folder]:
        if not os.path.isdir(i):
            os.mkdir(i)
    result_pool = {}
    # Sort LINgroups first, putting potential LINgroups to the back
    sorted_LINgroups = sorted(LINgroup_membership.keys())
    for each_LINgroup in sorted_LINgroups:
        if each_LINgroup != "new":
            print("Assigning LINs to the genomes of LINgroup {0}. Current time: {1}".format(each_LINgroup, current_time()))
            included_clusters = LINgroup_membership[each_LINgroup]
            for each_cluster in included_clusters:
                for each_genome in each_cluster:
                    each_result = process_each_genome_belonging_to_existing_LINgroup(each_genome,each_LINgroup,meta)
                    result_pool[each_genome] = each_result
        else:
            print("Assigning LINs to genomes belonging to new LINgroups, this will take a while. Current time: {0}".format(current_time()))
            for each_cluster in LINgroup_membership["new"]:
                result_pool = LIN_assignment_new_LINgroup(each_cluster,meta,c,result_pool,tmp_sig_folder,rep_bac_folder,result_folder)






def batch_submission(input_dir, output_dir, fastani_threshold,Username,InterestName):
    conn, c = connect_to_db()
    c.execute("SELECT User_ID FROM User WHERE Username='{0}'".format(Username))
    User_ID = int(c.fetchone()[0])
    c.execute("SELECT Interest_ID FROM Interest WHERE InterestName='{0}'".format(InterestName))
    Interest_ID_new_genome = int(c.fetchone()[0])
    # Starting step 1
    initial_clusters = initial_clustering(input_dir, output_dir, fastani_threshold)
    # Starting step 2
    LINgroup_membership, cluster_belonging, cluster_rep_bac_dictionary,meta = compare_init_cluster_w_existing_LINgroups(initial_clusters,output_dir)
    # Starting step 3


# MAIN