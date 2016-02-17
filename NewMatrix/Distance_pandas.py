#!/usr/bin/python
"""
This script is going to use pandas, aka data.frame to store data and apply distance/similarity calculation directly on
the data.frame, rather than store data into arrays, as is in GenerateMatrixMemory.
"""
# IMPORT
import pandas as pd
import h5py
import os
import sys
import scipy.spatial
import numpy as np


# FUNCTIONS
def read_into_dataframe(countfile):
    print "Started to read %s"%countfile
    f = h5py.File(countfile,'r')
    data = f['profiles']
    data = dict(data)
    df = pd.DataFrame(data)
    print "Loading completed."
    return df

def isfilepath(filepath):
    if filepath.endswith('/'):
        filepath = filepath
    else:
        filepath = filepath + '/'
    return filepath

def KmerCountNew(filepath):
    """
    :param filepath:
    :return: A k-mer counting profile will be generated on the hardisk
    """
    cmd = "kpal count -k 12 %s tmp_count"%(filepath)
    os.system(cmd)

def generate_distance_test(subjectfilepath, queryfilepath):
    """
    :param subjectfilepath:
    :param queryfilepath:
    :return: Print the distance matrix in a csv file.
    """
    queryfilepath = isfilepath(queryfilepath)
    KmerCountNew(queryfilepath)
    new_kmer_profile_path = filepath+'tmp_count'
    original_kmer_profile = read_into_dataframe(subjectfilepath)
    new_kmer_profile = read_into_dataframe(new_kmer_profile_path)
    print "Concatenating two data frames by columns..."
    total_mker_profile = pd.concat([original_kmer_profile, new_kmer_profile], axis=1)
    print "... Done."
    del original_kmer_profile
    del new_kmer_profile
    frequency_transform = lambda column: column/np.sum(column)
    print "Transforming k-mer countings to frequencies..."
    total_frequency = total_mker_profile.apply(frequency_transform)
    print "... Done."
    del total_mker_profile
    cosine_similarity = lambda column1, column2: scipy.spatial.distance.cosine(column1,column2)
    print "Calculating distance..."
    result = total_frequency.apply(lambda col1: total_frequency.apply(lambda col2: cosine_similarity(col1, col2)))
    print "... Done."
    print "Writing distance matrix to %s"%queryfilepath
    result.to_csv('%sdistance.csv'%queryfilepath)

def generate_distance(subjectpath,queryfilepath):
    """
    In the real case, we hope we only have one newly uploaded genome at a time, which means we don't calculate the whole
    data frame -- just the new one with the old one. And concatenate the new result to the existing distance matrix of
    the original genomes.
    :param subjectpath:
    :param queryfilepath
    :return: Print the distance matrix of the whole.
    """
    # Here we only have one fasta file
    KmerCountNew(queryfilepath)
    original_kmer = read_into_dataframe(subjectpath)
    new_kmer = read_into_dataframe('tmp_count')
    new_kmer_name = new_kmer.keys()[0]
    print "Concatenating two data frames..."
    total_mker_profile = pd.concat([original_kmer, new_kmer], axis=1)
    del original_kmer
    del new_kmer
    print "... Done."
    frequency_transform = lambda column: column/np.sum(column)
    print "Transforming counts to frequencies..."
    total_frequency = total_mker_profile.apply(frequency_transform)
    print "... Done."
    del total_mker_profile
    cosine_similarity = lambda column1, column2: scipy.spatial.distance.cosine(column1,column2)
    new_kmer_column = total_frequency[new_kmer_name]
    print "Calculating cosine similarities..."
    result = total_frequency.apply(lambda new_kmer_column: total_frequency.apply(lambda col2: cosine_similarity(new_kmer_column, col2)))
    print "... Done.\n\n"
    print "Writing distance matrix."
    result.to_csv('distance.csv')


def main(subjectpath,queryfilepath):
    generate_distance(subjectpath=subjectpath, queryfilepath=queryfilepath)



if __name__ == '__main__':
    filepath = sys.argv[2]
    subjectpath = sys.argv[1]
    main(queryfilepath=filepath,subjectpath=subjectpath)