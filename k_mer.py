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

# def generate_distance_test(subjectfilepath, queryfilepath):
#     """
#     :param subjectfilepath:
#     :param queryfilepath:
#     :return: column_names: the genomes that the new genome compares with;
#              result_new2old: similarity values where small means similar.
#     """
#     queryfilepath = isfilepath(queryfilepath)
#     KmerCountNew(queryfilepath)
#     new_kmer_profile_path = filepath+'tmp_count'
#     original_kmer_profile = read_into_dataframe(subjectfilepath)
#     new_kmer_profile = read_into_dataframe(new_kmer_profile_path)
#     print "Concatenating two data frames by columns..."
#     total_mker_profile = pd.concat([original_kmer_profile, new_kmer_profile], axis=1)
#     print "... Done."
#     del original_kmer_profile
#     del new_kmer_profile
#     frequency_transform = lambda column: column/np.sum(column)
#     print "Transforming k-mer countings to frequencies..."
#     total_frequency = total_mker_profile.apply(frequency_transform)
#     print "... Done."
#     del total_mker_profile
#     cosine_similarity = lambda column1, column2: scipy.spatial.distance.cosine(column1,column2)
#     print "Calculating distance..."
#     result = total_frequency.apply(lambda col1: total_frequency.apply(lambda col2: cosine_similarity(col1, col2)))
#     print "... Done."
#     print "Writing distance matrix to %s"%queryfilepath
#     result.to_csv('%sdistance.csv'%queryfilepath)

def generate_distance(subjectpath,queryfilepath):
    """
    In the real case, we hope we only have one newly uploaded genome at a time, which means we don't calculate the whole
    data frame -- just the new one with the old one. And concatenate the new result to the existing distance matrix of
    the original genomes.
    :param subjectpath: The k-mer frequency table of original genomes
    :param queryfilepath A FASTA file
    :return: Print the distance matrix of the whole.
    """
    # Here we only have one fasta file
    KmerCountNew(queryfilepath)
    # original_kmer = read_into_dataframe(subjectpath)
    subjectpath=isfilepath(subjectpath)
    subject_frequency_file = subjectpath+'frequency'
    subject_distance_matrix_file = subjectpath+'distance.csv'
    original_frequency = pd.read_hdf(subject_frequency_file,'profiles').to_sparse()
    new_kmer = read_into_dataframe('tmp_count').to_sparse()
    new_kmer_name = str(new_kmer.keys()[0])
    frequency_transform = lambda column: column/np.sum(column)
    print "transforming new counts to frequencies"
    # original_frequency = original_kmer.apply(frequency_transform)
    new_frequency = new_kmer.apply(frequency_transform)
    print "... Done."
    del new_kmer
    cosine_similarity = lambda column1, column2: scipy.spatial.distance.cosine(column1,column2)
    new_kmer_column = new_frequency[new_kmer_name]
    print "Calculating cosine similarities..."
    # result = total_frequency.apply(lambda new_kmer_column: total_frequency.apply(lambda col2: cosine_similarity(new_kmer_column, col2)))
    # result_new2old = original_frequency.apply(lambda column :cosine_similarity(column, new_kmer_column)) # Usually it's one column
    result_new2old = [1-cosine_similarity(new_kmer_column,original_frequency[i]) for i in original_frequency.keys()]
    column_names = original_frequency.keys()
    return column_names, result_new2old
    # newrow = [new_kmer_name]+result_new2old
    # print newrow
    # # Starting to manipulate the final distance matrix
    # print "Manipulating distance matrix"
    # original_distance_matrix = pd.read_csv(subject_distance_matrix_file)
    # original_distance_matrix[new_kmer_name]=result_new2old
    # original_distance_matrix.loc[original_distance_matrix.shape[0]]=newrow+[0]
    # print "Writing distance matrix."
    # original_distance_matrix.to_csv('new_distance.csv')



if __name__ == '__main__':
    filepath = sys.argv[2]
    subjectpath = sys.argv[1]
    print generate_distance(queryfilepath=filepath,subjectpath=subjectpath)