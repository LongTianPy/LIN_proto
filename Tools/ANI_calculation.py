#!/usr/bin/python
"""
    The multi-threading or -processing function of pyani isn't working.
    This script manually set it to multi-processing using python's multiprocess module
"""

# IMPORT
import multiprocessing as mp
import os
from os import listdir
from os.path import isdir, isfile, join
import shutil
import pandas as pd
import sys
from functools import partial
import uuid

# FUNCTIONS
def create_job_map(working_dir):
    files = [file for file in listdir(working_dir) if isfile(join(working_dir,file))]
    job_pairs = []
    for i in range(len(files)):
        for other_file in files[i:]:
            job_pairs.append(files[i] + "+" +other_file)
    return job_pairs

def check_done_jobs(working_dir,job_pairs):
    done_jobs = []
    undones = []
    dirs = [join(working_dir,dir) for dir in listdir(working_dir) if isdir(join(working_dir,dir))]
    for dir in dirs:
        files = [file for file in listdir(dir) if isfile(join(dir,file)) and file.endswith(".fasta")]
        if isdir(join(dir,"output")) and "ANIblastall_percentage_identity.tab" in listdir(join(dir,"output")):
            done_jobs.append("+".join(files))
        else:
            # shutil.rmtree(dir)
            undones.append("+".join(files))
    print(len(done_jobs))

    tmp_undones = set(job_pairs) - set(done_jobs)
    undones = set(undones) | tmp_undones
    return undones


def use_pyani(pair_str,ANI,cov,aln):
    pair = pair_str.split("+")
    query = ".".join(pair[0].split(".")[:-1])
    subject = ".".join(pair[1].split(".")[:-1])
    workstation = str(uuid.uuid4())
    if not isdir(workstation):
        os.mkdir(workstation)
    else:
        workstation = str(uuid.uuid4())
        os.mkdir(workstation)
    shutil.copy(pair[0],workstation)
    shutil.copy(pair[1],workstation)
    cmd = "python3 /home/longtian/dragonstooth/python/bin/average_nucleotide_identity.py -i {0} -o {0}/output -m ANIblastall --nocompress".format(workstation)
    os.system(cmd)
    # ANI_df = pd.read_table("{0}/output/ANIblastall_percentage_identity.tab".format(workstation),header=0,index_col=0)
    # cov_df = pd.read_table("{0}/output/ANIblastall_alignment_coverage.tab".format(workstation),header=0,index_col=0)
    # aln_df = pd.read_table("{0}/output/ANIblastall_alignment_lengths.tab".format(workstation),header=0,index_col=0)
    # ANI.set_value(query,subject,ANI_df.get_value(query,subject))
    # ANI.set_value(subject,query,ANI_df.get_value(subject,query))
    # cov.set_value(query, subject, cov_df.get_value(query, subject))
    # cov.set_value(subject, query, cov_df.get_value(subject, query))
    # aln.set_value(query, subject, aln_df.get_value(query, subject))
    # aln.set_value(subject, query, aln_df.get_value(subject, query))

def fill_dfs(each_dir,ANI,cov,aln):
    files = [file for file in listdir(each_dir) if isfile(join(each_dir,file)) and file.endswith(".fasta")]
    prefix = [str(i.split(".")[0]) for i in files]
    ani_df = pd.read_table(each_dir+"/output/ANIblastall_percentage_identity.tab",header=0,index_col=0)
    cov_df = pd.read_table(each_dir+"/output/ANIblastall_alignment_coverage.tab",header=0,index_col=0)
    aln_df = pd.read_table(each_dir+"/output/ANIblastall_alignment_lengths.tab",header=0,index_col=0)
    ANI.set_value(int(prefix[0]), prefix[1], ani_df.get_value(int(prefix[0]),prefix[1]))
    ANI.set_value(int(prefix[1]), prefix[0], ani_df.get_value(int(prefix[1]), prefix[0]))
    cov.set_value(int(prefix[0]), prefix[1], cov_df.get_value(int(prefix[0]), prefix[1]))
    cov.set_value(int(prefix[1]), prefix[0], cov_df.get_value(int(prefix[1]), prefix[0]))
    aln.set_value(int(prefix[0]), prefix[1], aln_df.get_value(int(prefix[0]), prefix[1]))
    aln.set_value(int(prefix[1]), prefix[0], aln_df.get_value(int(prefix[1]), prefix[0]))

def mp_fill_dfs(working_dir,ANI,cov,aln):
    dirs = [join(working_dir,dir) for dir in listdir(working_dir) if isdir(join(working_dir,dir))]
    partial_fill_dfs = partial(fill_dfs,ANI=ANI,cov=cov,aln=aln)
    pool_size = 200
    pool = mp.Pool(processes=pool_size)
    pool.map(partial_fill_dfs,dirs)
    return ANI,cov,aln



# MAIN
if __name__ == '__main__':
    working_dir = sys.argv[1]
    # dirs = [dir for dir in listdir(working_dir) if isdir(dir)]
    # for dir in dirs:
    #     shutil.rmtree(join(working_dir,dir))
    # files = [str(file.split(".")[:-1]) for file in listdir(working_dir) if isfile(join(working_dir, file))]
    # job_pairs = create_job_map(working_dir=working_dir)
    # undone_job_pairs = check_done_jobs(working_dir=working_dir, job_pairs=job_pairs)
    # print(undone_job_pairs)
    #
    #
    ANI = pd.DataFrame(0,index=files,columns=files)
    cov = pd.DataFrame(0,index=files,columns=files)
    aln = pd.DataFrame(0,index=files,columns=files)
    # partial_use_pyani = partial(use_pyani,ANI=ANI,cov=cov,aln=aln)
    # pool_size = 200
    # pool = mp.Pool(processes=pool_size)
    # pool.map(partial_use_pyani,job_pairs)
    ANI_tb, cov_tb, aln_tb = mp_fill_dfs(working_dir=working_dir,ANI=ANI, cov=cov, aln=aln)
    if not isdir("/work/dragonstooth/longtian/Data/output"):
        os.mkdir("/work/dragonstooth/longtian/Data/output")
    ANI_tb.to_csv("/work/dragonstooth/longtian/Data/output/ANI.csv")
    cov_tb.to_csv("/work/dragonstooth/longtian/Data/output/coverage.csv")
    aln_tb.to_csv("/work/dragonstooth/longtian/Data/output/alignment_length.csv")

