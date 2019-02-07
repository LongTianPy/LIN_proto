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
# import sendEmail
import shutil
import multiprocessing as mp
from functools import partial
from uuid import uuid4
import workflow2
from scipy.cluster import hierarchy
import json

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2.0/all_sketches/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2.0/rep_bac/"
sourmash_batch_tmp = "/home/linproject/Workspace/Sourmash2/batch_tmp/"
sourmash_batch_tmp_rep_bac = "/home/linproject/Workspace/Sourmash2/batch_tmp/rep_bac/"
sourmash_batch_result = "/home/linproject/Workspace/Sourmash2/batch_result/"
original_folder  = '/home/linproject/Workspace/LINdb/'
tmp_folder = '/home/linproject/Workspace/tmp_upload/'
workspace_dir = '/home/linproject/Workspace/New/workspace/'
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']
Entrez.email = "aaa@bb.cc"


# FUNCTIONS
def connect_to_db():
    conn = Connect("127.0.0.1", "LINbase", "Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c


def create_signature(genome_filepath, working_dir):
    print("Creating 21,31,51-mer signatures of given files")
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)
    if not genome_filepath.endswith('/'):
        genome_filepath = genome_filepath + "/"
    signatures_dir = join(working_dir, 'signatures')
    if not isdir(signatures_dir):
        os.mkdir(signatures_dir)
    all_sig = join(signatures_dir, 'all')
    if not isdir(all_sig):
        os.mkdir(all_sig)
    else:
        sig_files = [file for file in os.listdir(all_sig) if file.endswith("sig")]
    genome_files = [f for f in os.listdir(genome_filepath) if os.path.isfile(join(genome_filepath, f))]
    file_map = {}
    if len(sig_files) == len(genome_files):
        for each in genome_files:
            file_map[each] = "{0}/{1}.sig".format(all_sig, ".".join(each.split(".")[:-1]))
        print("Signatures found, skipping...")
        return all_sig, file_map
    else:
        cmd = "sourmash compute {0}{1} -k 21,31,51 -n 1000 -o {2}/{3}.sig > /dev/null 2>&1"
        for each in genome_files:
            os.system(cmd.format(genome_filepath, each, all_sig, ".".join(each.split(".")[:-1])))
            file_map[each] = "{0}/{1}.sig".format(all_sig, ".".join(each.split(".")[:-1]))
        print("Signatures created.")
    return all_sig, file_map

def compare_signatures_pairwisely(all_sig, working_dir,k):
    print("Pairwisely comparing signatures, k={0}".format(k))
    if not working_dir.endswith("/"):
        working_dir += "/"
    cmd = "sourmash compare {0}/*.sig -k {2} -o {1}all_sig_matrix.txt --csv {1}all_sig_k{2}.csv -q"
    os.system(cmd.format(all_sig, working_dir,k))
    df = pd.read_csv("{0}".format(join(working_dir,"all_sig_k{0}.csv".format(k))), sep=",", header=0,dtype='float')
    df.index = df.columns
    return df


def generate_coarse_distance_matrix(genome_filepath, working_dir, k):
    all_sig, file_map = create_signature(genome_filepath, working_dir)
    df = compare_signatures_pairwisely(all_sig, working_dir, k)
    return df, file_map


def import_existing_distance_matrix(existing_df):
    df = pd.read_csv(existing_df, sep=",", header=0,dtype='float')
    df.index = df.columns
    return df


def cluster_by_threshold(df, threshold):
    print("Clustering genomes at jaccard similarity threshold of {0}".format(threshold))
    samples = list(df.columns)
    sample_distances = {}
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            sample_distances[(samples[i], samples[j])] = 1 - df.loc[samples[i], samples[j]]
    keys = [sorted(k) for k in sample_distances.keys()]
    values = sample_distances.values()
    sorted_keys, distances = zip(*sorted(zip(keys, values)))
    Z = hierarchy.linkage(distances)
    labels = sorted(set([key[0] for key in sorted_keys] + [sorted_keys[-1][-1]]))
    # hierarchy.dendrogram(Z, labels=labels)
    cutree = hierarchy.cut_tree(Z, height=threshold)
    clusters = {}
    for i in range(len(cutree)):
        if str(cutree[i][0]) not in clusters:
            clusters[str(cutree[i][0])] = [labels[i]]
        else:
            clusters[str(cutree[i][0])].append(labels[i])
    with open(join(working_dir,"clusters_k{0}.json"),"w") as f:
        json.dump(clusters,f)
    return clusters


def pick_representative(clusters, working_dir,file_map):
    print("Selecting the first genome of each cluster as the representative genome")
    uploaded_rep_bac = {clusters[i][0]:i for i in clusters.keys()}
    uploaded_rep_bac_dir = join(working_dir, "signatures", "rep_bac")
    sig_to_fasta = {}
    if isdir(uploaded_rep_bac_dir):
        shutil.rmtree(uploaded_rep_bac_dir)
    os.mkdir(uploaded_rep_bac_dir)
    # cmd = "sourmash compute {0} -k 21,31,51 -n 1000 -o {1}/{2}.sig > /dev/null 2>&1"
    for each in uploaded_rep_bac.keys():
        shutil.copy(file_map[each], uploaded_rep_bac_dir)
        # os.system(cmd.format(each, uploaded_rep_bac_dir, ".".join(each.split("/")[-1].split(".")[:-1])))
        sig_to_fasta[".".join(each.split("/")[-1].split(".")[:-1])+".sig"] = each
    return uploaded_rep_bac, uploaded_rep_bac_dir, sig_to_fasta

def parse_sourmash_search(file):
    with open(file,"r") as f:
        lines = [i.strip().split(",") for i in f.readlines()]
    if len(lines)<2:
        return ''
    else:
        return lines[1][2]

def compare_uploaded_with_existing_rep_bac(k,uploaded_rep_bac, uploaded_rep_bac_dir, sig_to_fasta, rep_bac_dir,working_dir,c):
    print("Comparing each representative genome to existing 95% LINgroups in LINbase")
    uploaded_cluster_to_LINgroup = {}
    cmd = "sourmash search {0} {1}*.sig -o {2} -k {3} -q"
    # rep_bac = [f[:-4] for f in os.listdir(uploaded_rep_bac_dir) if isfile(join(uploaded_rep_bac_dir,f) and f.endswith(".sig"))]
    for i in sig_to_fasta.keys():
        os.system(cmd.format(join(uploaded_rep_bac_dir, i), rep_bac_dir, join(working_dir,"tmp.csv"),k))
        best_rep_bac = parse_sourmash_search(join(working_dir,"tmp.csv"))
        if best_rep_bac != '':
            genome_id = int(best_rep_bac.split('.')[0])
            c.execute('SELECT LIN FROM LIN WHERE Genome_ID={0} AND Scheme_ID=4'.format(genome_id))
            lin = c.fetchone()[0]
            lingroup = ','.join(lin.split(',')[:6])
        else:
            lingroup = ''
        uploaded_cluster_to_LINgroup[uploaded_rep_bac[sig_to_fasta[i]]] = lingroup
    return uploaded_cluster_to_LINgroup

def coarse_search(genome_filepath, working_dir, k, threshold = 0.08,precomputed_sim_matrix = None):
    conn, c = connect_to_db()
    if precomputed_sim_matrix is None:
        df, file_map = generate_coarse_distance_matrix(genome_filepath, working_dir, k)
    else:
        df = import_existing_distance_matrix(precomputed_sim_matrix)
    clusters = cluster_by_threshold(df,threshold)
    uploaded_rep_bac, uploaded_rep_bac_dir, sig_to_fasta = pick_representative(clusters, working_dir)
    uploaded_cluster_to_LINgroup = compare_uploaded_with_existing_rep_bac(k, uploaded_rep_bac, uploaded_rep_bac_dir,
                                                                          sig_to_fasta, rep_bac_dir, working_dir, c)
    print("Exporting the files")
    with open(join(working_dir,"cluster_to_LINgroup.json"),"w") as f:
        json.dump(uploaded_cluster_to_LINgroup,f)
    return uploaded_cluster_to_LINgroup




# MAIN
if __name__ == '__main__':
    genome_filepath = sys.argv[1]
    working_dir = sys.argv[2]
    threshold = float(sys.argv[3])
    k = int(sys.argv[4])
    if len(sys.argv) == 6:
        precomputed = sys.argv[5]
    else:
        precomputed = None
    uploaded_cluster_to_LINgroup = coarse_search(genome_filepath,working_dir, k,threshold,precomputed)

