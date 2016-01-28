#!/usr/bin/python

# IMPORT
import os
import multiprocessing as mp
import sys
from os import listdir
from os.path import isfile, join
import time

# FUNCTIONS
def DistanceCalc(obj_list): # Left is one genome from new genome set, right is one genome from the original set
    """
    :param obj_list: A list with sublist where 0 is left, 1 is right and 2 is position
    :return: returns position and distance value
    """
    cmd = "kpal distance tmp_count /home/vinatzerlab/Data/kPALevaluation/Psy/countk12 -l %s -r %s > %s.txt"%(obj_list[0], obj_list[1], obj_list[2])
    os.system(cmd)
    with open("%s.txt"%obj_list[2],"r") as f:
        line = f.readlines()[0].strip().split(' ')
    print '\t'.join(obj_list)
    os.system('rm %s.txt'%obj_list[2])
    return [obj_list[2], line[2]]

def pseudoFunc(obj_list):
	"""This is a pseudo-function to test if multi-threading works"""
	return len(obj_list)

def NewCount(filepath):
    filepath = sys.argv[1]
    if filepath.endswith('/'):
        filepath = filepath
    else:
        filepath = filepath + '/'
    cmd1 = "kpal count -k 12 %s*.fasta tmp_count"%filepath
    os.system(cmd1)
    NewFiles = [f[:-6] for f in listdir(filepath) if isfile(join(filepath, f)) and f.endswith('.fasta')]
    NewFiles.sort()
    Num_NewFiles = len(NewFiles)
    OrigFiles = [f[:-6] for f in listdir('/home/vinatzerlab/Data/kPALevaluation/Psy/') if isfile(join('/home/vinatzerlab/Data/kPALevaluation/Psy/',f)) and f.endswith('fasta')] # File path here should be the workspace where all the original genomes are in.
    OrigFiles.sort()
    Num_OrigFiles = len(OrigFiles)
    index_OrigFiles = [[i,OrigFiles[i]] for i in range(Num_OrigFiles)]
    obj_list = []
    pos = 0
    for i in range(Num_NewFiles):
        for j in range(Num_OrigFiles):
            obj_list.append([NewFiles[i],OrigFiles[j],pos])
            pos += 1
    return obj_list
    # pool = mp.Pool(processes = 8)
    # results_matrix = [[]]*len(NewFiles)
    # for i in range(len(NewFiles)):
    #     results = [pool.apply(DistanceCalc,args=(NewFiles[i],j[1],j[0])) for j in index_OrigFiles]
    #     results.sort()
    #     results_tmp = [k[1] for k in results]
    #     results_matrix[i] = results_tmp
    # return results_matrix

if __name__ == "__main__":
    obj_list = NewCount(filepath=sys.argv[1])
    print (len(obj_list))
    pool_size = 8
    pool = mp.Pool(processes=pool_size)
    results = pool.map(DistanceCalc, obj_list)
    #results = pool.map(pseudoFunc,obj_list)
    time.sleep(20)
    results.sort()
    f = open('result.txt','r')
    for i in results:
	f.write(i)
	f.write('\n')
    f.close()
    print results
