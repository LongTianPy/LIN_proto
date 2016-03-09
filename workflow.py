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
import os
import sys

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
def main(new_genome):
    # There should be a script that uploads the genome sequence to a subfolder in the workspace where the name is
    # randomly and uniquely generated
    workspace_dir = '/home/linproject/Workspace/New/'
    subfolder = 'Test_3-7-2016/'
    workspace_dir = workspace_dir + subfolder
    # And we have the file name of the genome
    # Fetched from the front end
    new_genomeID = new_genome[:-6]
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
    similarity = k_mer.generate_distance(FrequencyFilePath, workspace_dir+new_genome) # Already sorted
    # Check the biggest value, if it is bigger than the bottomline of the cutoff being used
    if similarity[new_genomeID].max() < 0.6:
        print "No similar genome found, run ANIb calculation sequentially to all genome is recommended."
        sys.exit()
    else:
        print "Looking for 10 most similar genome from our database."
        top10 = similarity.head(5)['Genome'].values
        # Get their file paths and copy them to the workspace
        for i in top10:
            c.execute("SELECT FilePath FROM Genome WHERE FilePath like '%{0}%'".format(i))
            cmd = "cp {0} {1}".format(c.fetchone()[0], workspace_dir)
            os.system(cmd)
        # Now we have all of them in the workspace
        ANIb_result = ANI_Wrapper_2.unified_anib(workspace_dir)[new_genomeID]
        top1_genome = ANIb_result[[x for x in ANIb_result.index if x != new_genomeID]].idxmax()
        print top1_genome
        top1_similarity = ANIb_result[[x for x in ANIb_result if x != new_genomeID]].max()
        if top1_similarity >= 1:
            print "This is most likely " + top1_genome
        else:
            new_LIN_object = LIN_Assign.getLIN(genome=top1_genome, Scheme_ID=3, similarity=top1_similarity)
            new_LIN = LIN_Assign.Assign_LIN(new_LIN_object)
        db.commit()
        return new_LIN

if __name__ == '__main__':
    new_genome = sys.argv[1] # Actually fetched from front end
    print main(new_genome=new_genome)






