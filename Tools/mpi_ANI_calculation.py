#!/usr/bin/python
"""This is for the ARC mpi to perform ANI calculation
"""

# IMPORT
from mpi4py import MPI
import os
from os import listdir
from os.path import isfile, isdir ,join
import shutil
import ANI_calculation
import sys


# FUNCTIONS


# MAIN
if __name__ == '__main__':
    working_dir = sys.argv[1]
    files = [file for file in listdir(working_dir) if file.endswith(".fasta") and isfile(join(working_dir,file))]
    job_pairs = ANI_calculation.create_job_map(working_dir)
