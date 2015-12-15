#!/usr/bin/python

__author__ = 'longtian'

# IMPORT
from Bio import SeqIO
import re


# FUNCTIONS
def LoadGenome(genomefilepath): # Read and load genomes, concatenate into one whole string
    f = open(genomefilepath,'r')
    records = list(SeqIO.parse(f,'fasta'))
    f.close()
    contigs = [str(i.seq) for i in records]
    genome = [genomefilepath]+[''.join(contigs)]
    return genome

def CalcKMerProb(genome,KMers,k):
    L = len(genome[1])
    # KMer_counts_k = [0]*len(KMers) # Used for store K-mer counting records
    # KMer_counts_1k_1 = [0]*len(KMers)
    # KMer_counts_2k = [0]*len(KMers)
    # KMer_counts_2k_1 = [0]*len(KMers)
    # for i in range(len(KMers)):
    #     pattern_k = '(?=%s)'%KMers[i]
    #     pattern_1k_1 = '(?=%s)'%KMers[i][:-1]
    #     pattern_2k = '(?=%s)'%KMers[i][1:]
    #     pattern_2k_1 = '(?=%s)'%KMers[i][1:-1]

    pattern_k = ['(?=%s)'%i for i in KMers]
    pattern_1k_1 = ['(?=%s)'%i[:-1] for i in KMers]
    pattern_2k = ['(?=%s)'%i[1:] for i in KMers]
    pattern_2k_1 = ['(?=%s)'%i[1:-1] for i in KMers]


        # pattern_k = re.compile(pattern_k,re.IGNORECASE)
        # pattern_1k_1 = re.compile(pattern_1k_1,re.IGNORECASE)
        # pattern_2k = re.compile(pattern_2k,re.IGNORECASE)
        # pattern_2k_1 = re.compile(pattern_2k_1,re.IGNORECASE)
    KMer_counts_k = map(lambda x:len([m.start() for m in re.finditer(re.compile(x,re.IGNORECASE),genome[1])]), pattern_k)
    KMer_counts_1k_1 = map(lambda x:len([m.start() for m in re.finditer(re.compile(x,re.IGNORECASE),genome[1])]), pattern_1k_1)
    KMer_counts_2k = map(lambda x:len([m.start() for m in re.finditer(re.compile(x,re.IGNORECASE),genome[1])]), pattern_2k)
    KMer_counts_2k_1 = map(lambda x:len([m.start() for m in re.finditer(re.compile(x,re.IGNORECASE),genome[1])]), pattern_2k_1)






        # KMer_counts_k[i] = len([m.start() for m in re.finditer(pattern_k,genome[1])])
        # KMer_counts_1k_1[i] = len([m.start() for m in re.finditer(pattern_1k_1,genome[1])])
        # KMer_counts_2k[i] = len([m.start() for m in re.finditer(pattern_2k,genome[1])])
        # KMer_counts_2k_1[i] = len([m.start() for m in re.finditer(pattern_2k_1,genome[1])])
    Total_substrings = L-k+1
    Prob_appearance_k = [float(i)/Total_substrings for i in KMer_counts_k]
    Prob_appearance_1k_1 = [float(i)/Total_substrings for i in KMer_counts_1k_1]
    Prob_appearance_2k = [float(i)/Total_substrings for i in KMer_counts_2k]
    Prob_appearance_2k_1 = [float(i)/Total_substrings for i in KMer_counts_2k_1]
    expected_appearance_prob = [Prob_appearance_1k_1[i]*Prob_appearance_2k[i]/Prob_appearance_2k_1[i] if Prob_appearance_2k_1[i] != 0 else 0 for i in range(len(Prob_appearance_2k))]
    pi_k = [(Prob_appearance_k[i]-expected_appearance_prob[i])/expected_appearance_prob[i] if expected_appearance_prob[i] != 0 else 0 for i in range(len(Prob_appearance_k)) ]
    return pi_k

def CalcKMerProb_slidingwindow(genome,KMers,k):
    sequence = genome[1]
    L = length(sequence)
    # k
    KMer_counts_k = [0]*len(KMers)
    for i in range(L-k+1):
        substring = sequence[i:(i+k)]
        idx = KMers.index(substring)
        KMer_counts_k[idx] += 1
    # KMer from 1 to k-1
    KMers_1k1 = [i[:-1] for i in KMers]
    KMer_counts_1k_1 = [0]*len(K_1Mers)
    # KMer from 2 to k
    KMers_2k = [i[1:] for i in KMers]
    KMer_counts_2k = [0]*len(KMers_2k)
    for i in range(L-k-1+1):
        substring = sequence[i:(i+k-1)]
        idx1 = KMers_1k1.index(substring)
        idx2 = KMers_2k.index(substring)
        KMer_counts_1k_1[idx1] += 1
        KMer_counts_2k[idx2] += 1
    # KMer from 2 to k-1
    KMers_2k1 = [i[1:-1] for i in KMers]
    KMer_counts_2k_1 = [0]*len(KMers_2k1)
    for i in range(L-k-2+1):
        substring = sequence[i:(i+k-2)]
        idx = KMers_2k1.index(substring)
        KMer_counts_2k_1[idx] += 1
    Total_substrings = L-k+1
    Prob_appearance_k = [float(i)/Total_substrings for i in KMer_counts_k]
    Prob_appearance_1k_1 = [float(i)/Total_substrings for i in KMer_counts_1k_1]
    Prob_appearance_2k = [float(i)/Total_substrings for i in KMer_counts_2k]
    Prob_appearance_2k_1 = [float(i)/Total_substrings for i in KMer_counts_2k_1]
    expected_appearance_prob = [Prob_appearance_1k_1[i]*Prob_appearance_2k[i]/Prob_appearance_2k_1[i] if Prob_appearance_2k_1[i] != 0 else 0 for i in range(len(Prob_appearance_2k))]
    pi_k = [(Prob_appearance_k[i]-expected_appearance_prob[i])/expected_appearance_prob[i] if expected_appearance_prob[i] != 0 else 0 for i in range(len(Prob_appearance_k)) ]
    return pi_k



