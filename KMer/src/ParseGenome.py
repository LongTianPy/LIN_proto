#!/usr/bin/python

__author__ = 'longtian'

# IMPORT
from Bio import SeqIO
from KMerPrep import KMerLib
import re
import sys


# FUNCTIONS
def LoadGenome(genomefilepath): # Read and load genomes, concatenate into one whole string
    f = open(genomefilepath,'r')
    records = list(SeqIO.parse(f,'fasta'))
    f.close()
    contigs = [str(i.seq) for i in records]
    genome = [genomefilepath]+[''.join(contigs)]
    return genome

def CountKMer(genomefilepath,Kmers):
    genome = LoadGenome(genomefilepath)
    KMer_counts = [0]*len(Kmers) # Used for store K-mer counting records
    for i in range(len(Kmers)):
        pattern = '(?=%s)'%Kmers[i]
        pattern = re.compile(pattern,re.IGNORECASE)
        KMer_counts[i] = len([m.start() for m in re.finditer(pattern,genome[1])])
    return genomefilepath, KMer_counts

def KMerFreq(KMer_counts):


if __name__ == '__main__':
    genomefilepath = sys.argv[1]
    k = sys.argv[2]
    KMers = KMerLib(k)

