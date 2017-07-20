#!/usr/bin/python
"""
"""

# IMPORT
import pandas as pd
import multiprocessing as mp
import sys
import os
from os import listdir
from os.path import isfile, isdir, join

# FUNCTIONS
def create_job_map(working_dir):
    all_dirs = [dir for dir in listdir(working_dir) if isdir(join(working_dir,dir))]
    dirs = []
    job_map = []
    idx_count = 1
    for i in all_dirs:
        if "output" in listdir(i):
            if "ANIblastall_percentage_identity.tab" in listdir(join(i,"output")):
                files = [file for file in listdir(i) if file.endswith(".fasta")]
                prefix = [str(file.split(".")[0]) for file in files]
                job_map.append([idx_count,i,prefix[0],prefix[1]])
                idx_count += 1
    return job_map

def create_empty_dfs(working_dir):
    files = [file for file in listdir(working_dir) if file.endswith(".fasta")]
    cols = [str(file.split(".")[0]) for file in files]
    idxs = [int(i) for i in cols]
    ani = pd.DataFrame(columns=cols,index=idxs)
    cov = pd.DataFrame(columns=cols,index=idxs)
    aln = pd.DataFrame(columns=cols,index=idxs)
    return ani, cov, aln

def fill_dfs(job_pair,ani,cov,aln):
    idx = job_pair[0]
    dir = job_pair[1]
    file1 = job_pair[2]
    file2 = job_pair[3]
    ani_df = pd.read_table(join(dir,"/output/ANIblastall_percentage_identity.tab"),header=0,index_col=0)
    cov_df = pd.read_table(join(dir,"/output/ANIblastall_alignment_coverage.tab"),header=0,index_col=0)
    aln_df = pd.read_table(join(dir,"/output/ANIblastall_alignment_lengths.tab"),header=0,index_col=0)
    ani.loc[int(file1),file2] = ani_df.get_value(int(file1),file2)
    ani.loc[int(file2),file1] = ani_df.get_value(int(file2),file1)
    cov.loc[int(file1), file2] = cov_df.get_value(int(file1), file2)
    cov.loc[int(file2), file1] = cov_df.get_value(int(file2), file1)
    aln.loc[int(file1), file2] = aln_df.get_value(int(file1), file2)
    aln.loc[int(file2), file1] = aln_df.get_value(int(file2), file1)
    return ani, cov, aln

# MAIN
if __name__ == '__main__':
    working_dir = sys.argv[1]
    job_map = create_job_map(working_dir=working_dir)
    ani, cov, aln = create_empty_dfs(working_dir=working_dir)
    for each_job in job_map:
        ani, cov, aln = fill_dfs(job_pair=each_job, ani=ani, cov=cov, aln=aln)
    ani.to_csv("../ANI.csv")
    cov.to_csv("../cov.csv")
    aln.to_csv("../aln.cov")

