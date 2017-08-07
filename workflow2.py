#!/usr/bin/python
"""New workflow wrapper, make the code and thoughts cleaner and straightforward.

"""

# IMPORT
import LIN_Assign
import LINgroup_indexing
import mash_indexing
import MySQLdb
from MySQLdb import Connect
import pandas as pd
import os
from os.path import isdir, isfile, join
import sys
from LoadingExternalInfo import LoadInfo
import ExtractInfo
import IntermediateResult
import logging
import logging.handlers
import argparse
from datetime import datetime
from pytz import timezone
from Bio import SeqIO
import filecmp
import uuid
import sendEmail
import shutil

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash/rep_bac/"
original_folder  = '/home/linproject/Workspace/Psy_166/'
workspace_dir = '/home/linproject/Workspace/New/workspace/'


# FUNCTIONS

### Parse arguments
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="LIN platform backend"
    )
    parser.add_argument("-i", dest="new_genome", help="xxxxxx.fasta")
    parser.add_argument("-u", dest="User_ID", help="An interger")
    parser.add_argument("-s", dest="Interest_ID", help="Interest ID")
    parser.add_argument("-t", dest="Attributes", help="Attributes")
    parser.add_argument("-p", dest="privacy", help="Is it private information")
    args = parser.parse_args()
    return args

### Pre-process to see if there's duplication
### Use either bbmap or sourmash to do this and return the result if it's a new genome


### Load the new genome into database with metadata

### Assign LIN

### Email


# MAIN