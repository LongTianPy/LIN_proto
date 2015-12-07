#!/usr/bin/python

# Test run on single file

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

def main(givenfile, k):
    k = int(k)
    KMers = KMerPrep.KMerLib(k)
    genome = ParseGenome.LoadGenome(givenfile)
    pi_ks = ParseGenome.CalcKMerProb(genome[1],KMers,k)
    pi_ks = np.array(pi_ks)
    pi_ks_zscores = stats.zscore(pi_ks,axis = 1)
    rownames = [i[0] for i in genomes]
    colnames = [i[0] for i in genomes]
