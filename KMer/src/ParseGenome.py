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

def CalcKMerProb(genome,KMers):
    L = len(genome)
    KMer_counts_k = [0]*len(KMers) # Used for store K-mer counting records
    KMer_counts_1k_1 = [0]*len(KMers)
    KMer_counts_2k = [0]*len(KMers)
    KMer_counts_2k_1 = [0]*len(KMers)
    for i in range(len(KMers)):
        pattern_k = '(?=%s)'%KMers[i]
        pattern_1k_1 = '(?=%s)'%KMers[i][:-1]
        pattern_2k = '(?=%s)'%KMers[i][1:]
        pattern_2k_1 = '(?=%s)'%KMers[i][1:-1]
        pattern_k = re.compile(pattern_k,re.IGNORECASE)
        pattern_1k_1 = re.compile(pattern_1k_1,re.IGNORECASE)
        pattern_2k = re.compile(pattern_2k,re.IGNORECASE)
        pattern_2k_1 = re.compile(pattern_2k_1,re.IGNORECASE)
        KMer_counts_k[i] = len([m.start() for m in re.finditer(pattern_k,genome[1])])
        KMer_counts_1k_1[i] = len([m.start() for m in re.finditer(pattern_1k_1,genome[1])])
        KMer_counts_2k[i] = len([m.start() for m in re.finditer(pattern_2k,genome[1])])
        KMer_counts_2k_1[i] = len([m.start() for m in re.finditer(pattern_2k_1,genome[1])])
    Total_substrings = L-k+1
    Prob_appearance_k = [float(i)/Total_substrings for i in KMer_counts_k]
    Prob_appearance_1k_1 = [float(i)/Total_substrings for i in KMer_counts_1k_1]
    Prob_appearance_2k = [float(i)/Total_substrings for i in KMer_counts_2k]
    Prob_appearance_2k_1 = [float(i)/Total_substrings for i in KMer_counts_2k_1]
    expected_appearance_prob = [Prob_appearance_1k_1[i]*Prob_appearance_2k[i]/Prob_appearance_2k_1[i] for i in range(len(Prob_appearance_2k))]
    pi_k = [(Prob_appearance_k[i]-expected_appearance_prob[i])/expected_appearance_prob[i] if expected_appearance_prob[i] != 0 else 0 for i in range(len(Prob_appearance_k)) ]
    return pi_k


def KMerFreq(genomefilepath,k):
    genomename, KMer_counts, L = CountKMer(genomefilepath,k)
    Total_substrings = L-k+1
    Prob_appearance = [float(i)/Total_substrings for i in KMer_counts]



