#!/usr/bin/python
"""This is for the ARC mpi to perform ANI calculation
"""

# IMPORT
from mpi4py import MPI
import os
from os import listdir
from os.path import isfile, isdir ,join
import shutil
import sys
import random


# FUNCTIONS
def create_job_map(working_dir):
    files = [file for file in listdir(working_dir) if isfile(join(working_dir,file)) and file.endswith(".fasta")]
    job_pairs = []
    for i in range(len(files))[:-1]:
        for other_file in files[i+1:]:
            job_pairs.append(files[i] + "+" +other_file)
    return job_pairs

def use_pyani(job_pair):
    pair = job_pair.split("+")
    [query, subject] = pair
    workstation = job_pair
    if not isdir(workstation):
        try:
            os.mkdir(workstation)
            shutil.copy(query,workstation)
            shutil.copy(subject,workstation)
            cmd = "python3 /home/longtian/dragonstooth/python/bin/average_nucleotide_identity.py -i {0} -o {0}/output -m ANIblastall --nocompress"
            os.system(cmd.format(workstation))
        except FileExistsError:
            pass

# MAIN
if __name__ == '__main__':
    working_dir = sys.argv[1]
    files = [str(file.split(".")[0]) for file in listdir(working_dir) if isfile(join(working_dir, file))]
    job_pairs = create_job_map(working_dir)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    number = len(job_pairs)
    for i in range(number):
        idx = random.randint(0,number)
        job_pair = job_pairs[idx]
        use_pyani(job_pair=job_pair)
