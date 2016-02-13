#!/usr/bin/python
"""
This script is going to use pandas, aka data.frame to store data and apply distance/similarity calculation directly on
the data.frame, rather than store data into arrays, as is in GenerateMatrixMemory.
"""
# IMPORT
import pandas
import h5py
import os
import sys

# FUNCTIONS
def read_into_dataframe(countfile):
    f = h5py.File(countfile,'r')
    data = f['profiles']

