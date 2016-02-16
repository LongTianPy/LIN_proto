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
    cmd = "kpal count -k 12 %s*.fasta %stmp_count"%(filepath, filepath)
    os.system(cmd)

# MAIN
def main(subjectfilepath, queryfilepath):
    queryfilepath = isfilepath(queryfilepath)
    KmerCountNew(queryfilepath)
    new_kmer_profile_path = filepath+'tmp_count'
    original_kmer_profile = read_into_dataframe(subjectfilepath)
    new_kmer_profile = read_into_dataframe(new_kmer_profile_path)
    print "Concatenating two data frames by columns"
    total_mker_profile = pd.concat([original_kmer_profile, new_kmer_profile], axis=1)
    print "... Done."
    del original_kmer_profile
    del new_kmer_profile
    euclidean_distance = lambda column1, column2: pd.np.linalg.norm(column1 - column2)
    print "Calculating euclidean distance"
    result = total_mker_profile.apply(lambda col1: total_mker_profile.apply(lambda col2:distance(col1, col2)))
    print "... Done."
    print result

if __name__ == '__main__':
    filepath = sys.argv[2]
    subjectpath = sys.argv[1]
    main(queryfilepath=filepath,subjectfilepath=subjectpath)