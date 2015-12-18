#!/usr/bin/python

__author__ = 'longtian'

# IMPORT
import numpy as np
from scipy import stats
from scipy import spatial
import os
from os import path
import ParseGenome

# FUNCTIONS
def OriMatrix(givendir, KMers, k): # Should be a workspace directory, use for test, we use '/Users/longtian/Desktop/LIN/LIN_proto/KMer/test'
    files = [i for i in os.listdir(givendir)]
    fullfilepath = [givendir+str(i) for i in files]
    print files
    genomes = [ParseGenome.LoadGenome(i) for i in fullfilepath]
    print genomes[0][0]
    print genomes[0][1][:100]
    pi_ks = [ParseGenome.CalcKMer_SlidWindow_Simple(i, KMers, k) for i in genomes]
    pi_ks = np.array(pi_ks)
    pi_ks_zscores = stats.zscore(pi_ks,axis = 1)
    rownames = [i[0] for i in genomes]
    colnames = [i[0] for i in genomes]
    print 'pi_ks_zscores'+'\n'+str(len(pi_ks_zscores))
    print pi_ks_zscores
    return rownames, colnames, pi_ks_zscores

def DisMatrix(pi_ks_zscores):
    SimMatrix = spatial.distance.pdist(pi_ks_zscores, 'euclidean')
    print 'SimMatrix done\n'
    print SimMatrix
    try:
        f = open('/home/vinatzerlab/Desktop/DistMat_DisMatrix.txt','w')
        for i in SimMatrix:
            f.write(i)
            f.write('\n')
    except:
        print 'f.write doesnt work'
        return SimMatrix
    return SimMatrix

def writeMatrix(SimMatrix, colnames):
    np.savetxt('/home/vinatzerlab/Desktop/DistMat_lambda.txt', SimMatrix, delimiter='t',header='\t'+'\t'.join(colnames))






