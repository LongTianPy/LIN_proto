#!/usr/bin/python

# WHERE MAIN FUNCTION IS

__author__ = 'longtian'

# IMPORT
from Bio import SeqIO
from itertools import product
import re
import numpy as np
from scipy import stats
from scipy import spatial
import os
from os import path
import optparse
import GenerateMatrix
import KMerPrep
import ParseGenome

def main(k, givendir):
    k = int(k)
    KMers = KMerPrep.KMerLib(k)
    rownames, colnames, pi_ks_zscores = GenerateMatrix.OriMatrix(givendir=givendir, KMers=KMers, k=k)
    SimMatrix = GenerateMatrix.DisMatrix(pi_ks_zscores)
    GenerateMatrix.writeMatrix(SimMatrix,colnames=colnames)

if __name__ == '__main__':
    # parser object for managing input options
    parser = optparse.OptionParser()

    # essential data
    parser.add_option('-p', dest = 'givendir' ,   default='', help = 'Path of the workspace that contains all the files' )
    parser.add_option('-k', dest = 'k' ,  default = '' , help = 'Length of the substring, K-mer' )
    #parser.add_option('-o',  dest = 'output' , default = '' , help = 'Outputpath used to store results' ) #Not using it for now
    # load the inputs
    (options , args) = parser.parse_args()

    # process the inputs
    # note: all commandline inputs are str by default
    givendir = options.givendir
    k = options.k
    #out = options.output

    main(k=k, givendir=givendir)
