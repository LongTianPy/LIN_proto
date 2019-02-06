#!/usr/bin/python
"""
Basically it's a LIN assignment process. With the following steps
1. Check if initializing procedures are needed. If not,
    1.1 create a database with needed tables and columns in them.
    1.2 create folders for k-mer signatures
2.
"""

# IMPORT
import os
import pandas as pd
from os.path import join
from scipy.cluster import hierarchy
from scipy.spatial import distance

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2/all_sketches/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2/rep_bac/"

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost", "LINbase","Latham@537")
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
    genome_files = [f for f in os.listdir(genome_filepath) if os.path.isfile(join(genome_filepath,f))]
    cmd = "sourmash compute {0}{1} -k 21,31,51 -n 1000 -o {2}/{3}.sig > /dev/null 2>&1"
    for each in genome_files:
        os.system(cmd.format(genome_filepath,each,all_sig,".".join(each.split(".")[:-1])))
    return all_sig

def compare_signatures_pairwisely(all_sig,working_dir):
    if not working_dir.endswith("/"):
        working_dir += "/"
    cmd = "sourmash compare {0}/*.sig -k 21 -o {1}all_sig_matrix.txt --csv {1}all_sig.csv -q"
    os.system(cmd.format(all_sig, working_dir))
    df = pd.read_csv("{0}all_sig.csv".format(working_dir), sep=",",header=0)
    df.index = df.columns
    return df

def generate_coarse_distance_matrix(genome_filepath, working_dir):
    all_sig = create_signature(genome_filepath,working_dir)
    df = compare_signatures_pairwisely(all_sig, working_dir)
    return df

def import_existing_distance_matrix(existing_df):
    df = pd.read_csv(existing_df,sep="\t",header=0,index_col=0)
    return df

def cluster_by_threshold(df,threshold):
    samples = list(df.columns)
    sample_distances = {}
    for i in range(len(samples)):
        for j in range(i+1,len(samples)):
            sample_distances[(samples[i],samples[j])] = 1- df.loc[samples[i],samples[j]]
    keys = [sorted(k) for k in sample_distances.keys()]
    values = sample_distances.values()
    sorted_keys, distances = zip(*sorted(zip(keys,values)))
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

def pick_representative(clusters,working_dir):
    uploaded_rep_bac = {clusters[i][0] for i in clusters.keys()}
    uploaded_rep_bac_dir = join(working_dir,"signatures","rep_bac")
    os.mkdir(uploaded_rep_bac_dir)
    cmd = "sourmash compute {0} -k 21,31,51 -n 1000 -o {2}/{3} > /dev/null 2>&1"
    for each in uploaded_rep_bac.keys():
        os.system(cmd.format(each,uploaded_rep_bac_dir,".".join(each.split("/")[-1].split(".")[:-1])))
    return uploaded_rep_bac, uploaded_rep_bac_dir

def compare_uploaded_with_existing_rep_bac(uploaded_rep_bac,uploaded_rep_bac_dir, sourmash_rep_bac):
    pass


# MAIN