#!/usr/bin/python
"""
This script implements functions of fast k-mer calculation to narrow search subjects, ANIb calculation for accurate
genomic percentage identity and LIN assignment to the new genome.
"""

# IMPORT
import ANI_Wrapper_2
import k_mer
import LIN_Assign
import MySQLdb
from MySQLdb import Connect
import pandas

# INPUT
#
# Sequentially, we will first process the genome with k-mer calculation. The input is the name of the genome "genome";
# the "Interest_ID", by which we know the part of genomes to compare.
# The output of this step is a data frame of similarities between this new genome with all those subject genomes.
# One condition is if the biggest similarity is below 60% which is the bottomline of LIN assignment cutoff, we suggest
# not using the following workflow but running the LIN assignment program Haitham development with JSpecies.
#
# The second step is ANIb calculation. When we find out the number of most similar genomes to the new genome by k-mer,
# we process them to ANIb calculation. In this step, we need first find out the those 2nd level subject genomes, by a
# side-function.
# The output of ANIb would be the name of the most similar genome and their percentage identity.
#
# Then the last step is to assign a new LIN to this genome. It takes the genome name "genome", Scheme_ID, and the
# similarity

# MAIN
def main():
    # There should be a script that uploads the genome sequence to a subfolder in the workspace where the name is
    # randomly and uniquely generated
    workspace_dir = '/home/vinatzerlab/Data/LIN_Workspace/'
    subfolder = 'Test_3-7-2016/'
    workspace_dir = workspace_dir + subfolder
    # And we have the file name of the genome
    new_genome = '3337_Psy-DSM10604.fasta' # Fetched from the front end
    # As well as its Interest_ID
    Interest_ID_new_genome = 1 # We hard-code it here, but it should be able to be read from the front end
    db = Connect('localhost','root')
    c = db.cursor()
    c.execute('use LINdb_test')
    # # Fetch the file paths of all the genomes from the database that have the same interest ID
    # c.execute('SELECT FilePath FROM Genome WHERE Interest_ID = {0}'.format(Interest_ID_new_genome))
    # FilePaths = c.fetchall()
    # FilePaths = [i[0] for i in FilePaths]
    # We need first to read their k-mer frequencies, which, are calculated beforehand and store in the server.
    # This reminds me of adding one more table for the location of those frequency files.
    ## And supposedly, the frequency of new genome shuold also be added to this file once it's calculated.
    c.execute('SELECT FrequencyFilePath FROM FrequencyFile WHERE Interest_ID = {0}'.format(Interest_ID_new_genome))
    FrequencyFilePath = c.fetchone()[0]
    similarity = k_mer.generate_distance(FrequencyFilePath, workspace_dir)
    # Then we need to sort them, or get the top n similar genomes





