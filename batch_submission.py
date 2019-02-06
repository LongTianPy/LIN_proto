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
from scipy.cluster import hierarchy

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
def connect_to_db():
    conn = Connect("localhost", "LINbase", "Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c


def create_signature(genome_filepath, working_dir):
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)
    if not genome_filepath.endswith('/'):
        genome_filepath = genome_filepath + "/"
    signatures_dir = join(working_dir, 'signatures')
    os.mkdir(signatures_dir)
    all_sig = join(signatures_dir, 'all')
    os.mkdir(all_sig)
    genome_files = [f for f in os.listdir(genome_filepath) if os.path.isfile(join(genome_filepath, f))]
    cmd = "sourmash compute {0}{1} -k 21,31,51 -n 1000 -o {2}/{3}.sig > /dev/null 2>&1"
    for each in genome_files:
        os.system(cmd.format(genome_filepath, each, all_sig, ".".join(each.split(".")[:-1])))
    return all_sig


def compare_signatures_pairwisely(all_sig, working_dir):
    if not working_dir.endswith("/"):
        working_dir += "/"
    cmd = "sourmash compare {0}/*.sig -k 31 -o {1}all_sig_matrix.txt --csv {1}all_sig.csv -q"
    os.system(cmd.format(all_sig, working_dir))
    df = pd.read_csv("{0}all_sig.csv".format(working_dir), sep=",", header=0)
    df.index = df.columns
    return df


def generate_coarse_distance_matrix(genome_filepath, working_dir):
    all_sig = create_signature(genome_filepath, working_dir)
    df = compare_signatures_pairwisely(all_sig, working_dir)
    return df


def import_existing_distance_matrix(existing_df):
    df = pd.read_csv(existing_df, sep="\t", header=0, index_col=0)
    return df


def cluster_by_threshold(df, threshold):
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
    return clusters


def pick_representative(clusters, working_dir):
    uploaded_rep_bac = {clusters[i][0]:clusters.keys()[i] for i in clusters.keys()}
    uploaded_rep_bac_dir = join(working_dir, "signatures", "rep_bac")
    sig_to_fasta = {}
    os.mkdir(uploaded_rep_bac_dir)
    cmd = "sourmash compute {0} -k 21,31,51 -n 1000 -o {2}/{3}.sig > /dev/null 2>&1"
    for each in uploaded_rep_bac.keys():
        os.system(cmd.format(each, uploaded_rep_bac_dir, ".".join(each.split("/")[-1].split(".")[:-1])))
        sig_to_fasta[".".join(each.split("/")[-1].split(".")[:-1])+".sig"] = each
    return uploaded_rep_bac, uploaded_rep_bac_dir, sig_to_fasta

def parse_sourmash_search(file):
    with open(file,"r") as f:
        lines = [i.strip().split(",") for i in f.readlines()]
    if len(lines)<2:
        return ''
    else:
        return lines[1][2]

def compare_uploaded_with_existing_rep_bac(uploaded_rep_bac, uploaded_rep_bac_dir, sig_to_fasta, rep_bac_dir,working_dir,c):
    uploaded_cluster_to_LINgroup = {}
    cmd = "sourmash search {0} {1}*.sig -o {2}"
    # rep_bac = [f[:-4] for f in os.listdir(uploaded_rep_bac_dir) if isfile(join(uploaded_rep_bac_dir,f) and f.endswith(".sig"))]
    for i in sig_to_fasta.keys():
        os.system(cmd.format(join(uploaded_rep_bac_dir, i), rep_bac_dir, join(working_dir,"tmp.csv")))
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

def coarse_search(genome_filepath, working_dir, threshold = 0.08,precomputed_sim_matrix = None):
    conn, c = connect_to_db()
    if precomputed_sim_matrix is None:
        df = generate_coarse_distance_matrix(genome_filepath, working_dir)
    else:
        df = import_existing_distance_matrix(precomputed_sim_matrix)
    clusters = cluster_by_threshold(df,threshold)
    uploaded_rep_bac, uploaded_rep_bac_dir, sig_to_fasta = pick_representative(clusters, working_dir)
    uploaded_cluster_to_LINgroup = compare_uploaded_with_existing_rep_bac(uploaded_rep_bac, uploaded_rep_bac_dir,
                                                                          sig_to_fasta, rep_bac_dir, working_dir, c)
    return uploaded_cluster_to_LINgroup



# MAIN
if __name__ == '__main__':
    genome_filepath = sys.argv[1]
    working_dir = sys.argv[2]
    threshold = sys.argv[3]
    uploaded_cluster_to_LINgroup = coarse_search()