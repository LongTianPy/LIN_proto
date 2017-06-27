#!/usr/bin/python
"""
    The multi-threading or -processing function of pyani isn't working.
    This script manually set it to multi-processing using python's multiprocess module
"""

# IMPORT
import multiprocessing as mp
import os
from os import listdir
from os.path import isdir
import shutil
import pandas as pd
import sys
from functools import partial
import uuid

# FUNCTIONS
def create_job_map(working_dir):
    files = [file for file in listdir(working_dir)]
    job_pairs = []
    for i in range(len(files)):
        for other_file in files[i:]:
            job_pairs.append([files[i],other_file])
    return  job_pairs

def use_pyani(pair,ANI,cov,aln):
    workstation = str(uuid.uuid4())
    if not isdir(workstation):
        os.mkdir(workstation)
    else:
        workstation = str(uuid.uuid4())
        os.mkdir(workstation)
    shutil.copy(pair[0],workstation)
    shutil.copy(pair[1],workstation)
    cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} -o {0}/output" \
          "-m ANIblastall --nocompress".format(workstation)
    os.system(cmd)
    ANI_df = pd.read_table("{0}/ANIblastall_percentage_identity.tab".format(workstation),header=0,index_col=0)
    cov_df = pd.read_table("{0}/ANIblastall_alignment_coverage.tab".format(workstation),header=0,index_col=0)
    aln_df = pd.read_table("{0}/ANIblastall_alignment_lengths.tab".format(workstation),header=0,index_col=0)
    ANI.set_value(pair[0],pair[1],ANI_df.get_value(pair[0],pair[1]))
    ANI.set_value(pair[1],pair[0],ANI_df.get_value(pair[1],pair[0]))
    cov.set_value(pair[0], pair[1], cov_df.get_value(pair[0], pair[1]))
    cov.set_value(pair[1], pair[0], cov_df.get_value(pair[1], pair[0]))
    aln.set_value(pair[0], pair[1], aln_df.get_value(pair[0], pair[1]))
    aln.set_value(pair[1], pair[0], aln_df.get_value(pair[1], pair[0]))

# MAIN
if __name__ == '__main__':
    working_dir = sys.argv[1]
    files = [file for file in listdir(working_dir)]
    job_pairs = create_job_map(working_dir=working_dir)
    ANI = pd.DataFrame(0,index=files,columns=files)
    cov = pd.DataFrame(0,index=files,columns=files)
    aln = pd.DataFrame(0,index=files,columns=files)
    partial_use_pyani = partial(use_pyani,ANI=ANI,cov=cov,aln=aln)
    pool_size = 288
    pool = mp.Pool(processes=pool_size)
    pool.map(partial_use_pyani,job_pairs)
    os.mkdir("output")
    ANI.to_csv("output/ANI.csv")
    cov.to_csv("output/coverage.csv")
    aln.to_csv("output/alignment_length.csv")

