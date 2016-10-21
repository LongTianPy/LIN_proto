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
from MySQLdb import Connect
import MySQLdb
import logging


# FUNCTIONS
def read_into_dataframe(countfile):
    f = h5py.File(countfile, 'r')
    data = f['profiles']
    data = dict(data)
    df = pd.DataFrame(data)
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
    cmd = "kpal count -k 12 %s /home/linproject/Workspace/kpal_count/tmp_count" % (filepath)
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

def generate_distance(queryfilepath, Genome_ID, User_ID):
    """
    In the real case, we hope we only have one newly uploaded genome at a time, which means we don't calculate the whole
    data frame -- just the new one with the old one. And concatenate the new result to the existing distance matrix of
    the original genomes.
    :param subjectpath: The k-mer frequency table of original genomes
    :param queryfilepath A FASTA file
    :return: The similarity dictionary of the new genome with each of the original genomes.
    """
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    # Here we only have one fasta file
    logging.info("Calculating k-mer profile of the new submission.")
    KmerCountNew(queryfilepath)

    subject_frequency_file = '/home/linproject/Workspace/Psy_166/init/frequency_sample'

    logging.info("Reading original frequency matrix.")
    original_frequency = pd.read_hdf(subject_frequency_file, 'profiles')

    logging.info("Reading k-mer profile into data frame.")
    new_kmer = read_into_dataframe('/home/linproject/Workspace/kpal_count/tmp_count')

    logging.info("Changing the name of this data frame column to the Genome ID")
    new_kmer.columns = [Genome_ID]

    logging.info("Cleaning the processed k-mer profile")
    os.system('rm /home/linproject/Workspace/kpal_count/tmp_count')
    new_kmer_name = Genome_ID
    frequency_transform = lambda column: column / np.sum(column)

    logging.info("Transforming k-mer profile of new submission to frequency")
    # original_frequency = original_kmer.apply(frequency_transform)
    new_frequency = new_kmer.apply(frequency_transform)

    del new_kmer
    cosine_similarity = lambda column1, column2: scipy.spatial.distance.cosine(column1, column2)
    new_kmer_column = new_frequency[new_kmer_name]
    logging.info("Calculating cosine similarities...")
    # result = total_frequency.apply(lambda new_kmer_column: total_frequency.apply(lambda col2: cosine_similarity(
    # new_kmer_column, col2)))
    # result_new2old = original_frequency.apply(lambda column :cosine_similarity(column, new_kmer_column)) # Usually
    # it's one column
    result_new2old = [1 - cosine_similarity(new_kmer_column, original_frequency[i]) for i in original_frequency.keys()]
    similarities = pd.DataFrame()
    similarities['Genome'] = original_frequency.keys()
    similarities[new_kmer_name] = result_new2old
    # for i in range(len(result_new2old)):
    #     similarities[original_frequency.keys()[i]] = [result_new2old[i]]
    # original_frequency[new_kmer_name] = new_kmer_column
    # original_frequency.to_csv('/home/vinatzerlab/Desktop/updated_frequency.csv')
    similarities_sorted = similarities.sort_index(by=[new_kmer_name], ascending=False)
    new_hdf = pd.concat([original_frequency, new_frequency], axis=1)
    del original_frequency
    del new_frequency

    logging.info("Writing new frequency matrix.")

    # Remember un-cooment here
    new_hdf.to_hdf(subject_frequency_file, key='profiles')

    return similarities_sorted
    # return similarities.sort(axis=0, ascending=False, kind="mergesort")
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
    filepath = sys.argv[1]
    print generate_distance(queryfilepath=filepath)
