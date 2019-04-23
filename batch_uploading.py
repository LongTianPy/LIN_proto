#!/usr/bin/python
"""Given a table of NCBI accession numbers, create a command line list to call workflow2.genome_submission sequentially
"""

# IMPORT
import pandas as pd
import sys
from Bio import Entrez
import argparse

# VARIABLES
Entrez.email = "aaa@bbb.ccc"
wgs_base_url = 'ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/'


# FUNCTIONS
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="Download sequence by accession numbers or so, input is always a file with accession numbers, each per line.")
    parser.add_argument("-i", dest='file', help="CSV file containing NCBI accession numbers")
    parser.add_argument("-o", dest='out_dir', help="Target file folder.")
    args = parser.parse_args()
    return args

def file_meta(accession):
    handler = Entrez.esearch(db='nucleotide',term='NC_000964')
    handler = Entrez.efetch(db='nucleotide', id='255767013',retmode='xml')
    record = Entrez.read(handler)


# MAIN