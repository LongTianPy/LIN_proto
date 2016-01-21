#!/usr/bin/python

# IMPORT
import os
import multiprocessing as mp
import sys
from os import listdir
from os.path import isfile, join

# FUNCTIONS
def DistanceCalc(left,right,pos): # Left is one genome from new genome set, right is one genome from the original set
    cmd = "kpal distance tmp_count original_count -l %s -r %s > tmp"%(left,right)
    os.system(cmd)
    with open("tmp","r") as f:
        line = f.readlines()[0].strip().split(' ')
    return [pos, line[2]]

def NewCount(filepath):
    filepath = sys.argv[1]
    if filepath.endswith('/'):
        filepath = filepath
    else:
        filepath = filepath + '/'
    cmd1 = "kpal count -k 12 %s*.fasta tmp_count"%filepath
    os.system(cmd1)
    NewFiles = [f for f in listdir(filepath) if isfile(join(filepath, f)) and f.endswith('.fasta')].sort()
    OrigFiles = [f for f in listdir('~/Data/kPALevaluation/Psy/') if isfile(join('~/Data/kPALevaluation/Psy/',f)) and f.endswith('fasta')].sort() # File path here should be the workspace where all the original genomes are in.
    Num_OrigFiles = len(OrigFiles)
    index_OrigFiles = [[i,OrigFiles[i]] for i in range(Num_OrigFiles)]
    pool = mp.Pool(processes = 8)
    results_matrix = [[]]*len(NewFiles)
    for i in NewFiles:
        results = [pool.apply(DistanceCalc,args=(i,j[1],j[0])) for j in index_OrigFiles]
        results_tmp = results.sort()
        results_tmp2 = [k[1] for k in results_tmp]
        results_matrix[i] = results_tmp2
    return results_matrix

if __name__ == "__main__":
    NewCount(filepath=sys.argv[1])