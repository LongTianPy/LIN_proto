#!/usr/bin/python

"""
    This script will read the k-mer profile into array/matrix in the memory and calculate them with newly uploaded
    genomes.
"""

# IMPORT
import h5py
import numpy as np

# FUNCTIONS
def read_by_iteration(count_profile):
	f = h5py.File(count_profile,'r')
	data = f['profiles']
	genomes = data.keys()
	kmer_profiles = data[genomes[0]][:]
	i = 1
	while i < len(genomes):

		kmer_profiles = np.vstack((kmer_profiles,genomes[i]))
        i += 1
	return kmer_profiles

# Use os.system(cmd) count the k-mer profile of newly uploaded genome(s)
def KmerCountNew(filepath):
    """
    :param filepath:
    :return: A k-mer counting profile will be generated on the hardisk
    """
    if filepath.endswith('/'):
        filepath = filepath
    else:
        filepath = filepath + '/'
    cmd = "kpal count %s *.fasta -k 12 tmp_count"%filepath
    os.system(cmd)

# Use h5py, read k-mer profiles of both the original and new one
Original_kmer =
# And Read them into 2 arrays

# Concatenate these two arrays together, either concatenate/vstack will do, but better try which one is faster

# Probably need to transform to k-mer frequency matrix first

# Calculate the pairwise distance (Could be ALLvsALL or only calculate the distances between the original ones and new
# ones)

# Create Tree based on this distance matrix