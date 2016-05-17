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
import pandas as pd
import os
import sys
from sklearn.cluster import KMeans
from LoadingExternalInfo import LoadInfo
import ExtractInfo
import IntermediateResult
# import ExtractInfo

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
    original_folder  = '/home/linproject/Workspace/Zika/'
    workspace_dir = '/home/linproject/Workspace/New/'
    subfolder = 'workspace/'
    workspace_dir = workspace_dir + subfolder
    InfoFile = "/home/linproject/Workspace/Zika/Attribute_full.csv"
    # And we have the file name of the genome
    # Fetched from the front end
    new_genomeID = new_genome.split('.')[0]
    # As well as its Interest_ID
    Interest_ID_new_genome = 2 # We hard-code it here, but it should be able to be read from the front end
    db = Connect('localhost','root')
    c = db.cursor()
    c.execute('use LINdb_zika')
    c.execute('INSERT INTO Genome (Interest_ID, Submission_ID, FilePath, GenomeName) values ({0}, 1, "{1}", "{2}")'
              .format(Interest_ID_new_genome ,original_folder+new_genome, new_genomeID))
    db.commit()
    ## For Zika virus case only, comment out when initialization is done.
    # LoadInfo(InfoFile,c,new_genomeID,Interest_ID_new_genome)
    # db.commit()
    # # Fetch the file paths of all the genomes from the database that have the same interest ID
    # c.execute('SELECT FilePath FROM Genome WHERE Interest_ID = {0}'.format(Interest_ID_new_genome))
    # FilePaths = c.fetchall()
    # FilePaths = [i[0] for i in FilePaths]
    # We need first to read their k-mer frequencies, which, are calculated beforehand and store in the server.
    # This reminds me of adding one more table for the location of those frequency files.
    ## And supposedly, the frequency of new genome shuold also be added to this file once it's calculated.
    print original_folder+new_genome;
    similarity = k_mer.generate_distance(original_folder+new_genome) # Already sorted
    # Check the biggest value, if it is bigger than the bottomline of the cutoff being used
    # if similarity[new_genomeID].max() < 0.6:
    #     print "No similar genome found, run ANIb calculation sequentially to all genome is recommended."
    #     sys.exit()
    # else:
    print "Looking for similars genome from our database."
    if len(similarity['Genome'])<=10:
        n_top = len(similarity['Genome'])
    else: # As for selecting numbers of clusters, 3 might not be very good since it is possible when we have a lot of
          # genomes to select while there might be way more than 10 genomes in the top cluster. So we are setting the
          # the number of clusters based on the total number of genomes to select from. Since 10 is a preferred number
          # of genomes to perform pairwise blasting, I guess we can do (M/10)+1
        n_clusters = 3
        km = KMeans(n_clusters=n_clusters)
        km.fit(similarity[new_genomeID].reshape(-1,1))
        centroid_idx = list(km.cluster_centers_).index(max(km.cluster_centers_))
        top_cluster_idx = [i for i,x in enumerate(km.labels_) if x==centroid_idx]
        n_top = len(top_cluster_idx)
        print "We are comparing your genome with {0} genomes in our database.".format(n_top)
    top10 = similarity.head(n_top)['Genome'].values
    print top10
    try:
        IntermediateResult.write_kmer_result(top10=top10,db_cursor=c)
        IntermediateResult.send_email("kmer",User_ID=2,db_cursor=c)
    except:
        pass
    # top10_LINs = [ExtractInfo.get_top10_LIN(i,c) for i in top10] # This can be used to send preliminary results
    # print top10_LINs
    # Get their file paths and copy them to the workspace
    similarities = pd.DataFrame()
    for i in top10:
        c.execute("SELECT FilePath FROM Genome WHERE FilePath like '%{0}%'".format(i))
        cmd = "cp {0} {1}".format(c.fetchone()[0], workspace_dir)
        os.system(cmd)
        os.system('cp {0} {1}'.format(original_folder+new_genome,workspace_dir))
        # Now we have all of them in the workspace
        ANIb_result = ANI_Wrapper_2.unified_anib(workspace_dir)[new_genomeID]
        os.system('rm -rf {0}*'.format(workspace_dir))
        similarity = ANIb_result.loc[i]
        similarities[i]=[similarity]
    top1_genome = similarities.idxmax(axis=1)[0]
    top1_similarity = similarities.max(axis=1)[0]
    print top1_genome
    print top1_similarity
    # if top1_similarity >= 1:
    #     print "This is most likely to be" + top1_genome
    # else:
    new_LIN_object = LIN_Assign.getLIN(genome=top1_genome, Scheme_ID=3, similarity=top1_similarity,c=c)
    print "The most similar record is " + top1_genome+ " , whose LIN is " +','.join(new_LIN_object.LIN) + '.'
    print "The similarity to it is " + str(top1_similarity*100) + "%."
    new_LIN = LIN_Assign.Assign_LIN(new_LIN_object,c=c).new_LIN
    print "The LIN assigned to your genome is " + new_LIN
    c.execute('SELECT Genome_ID FROM Genome where FilePath like "%{0}%"'.format(new_genome))
    Genome_ID = int(c.fetchone()[0])
    c.execute("INSERT INTO LIN (Genome_ID, Scheme_ID, LIN, SubjectGenome, ANI) values ({0}, 3, '{1}', '{2}', {3})"
              .format(Genome_ID, new_LIN, top1_genome, top1_similarity))
    db.commit()
    try:
        IntermediateResult.write_ANI_result(new_genomeID=new_genomeID,new_LIN_object=new_LIN_object,new_LIN=new_LIN,db_cursor=c)
        IntermediateResult.send_email(file_source="ANI",User_ID=1,db_cursor=c)
    except:
        pass




if __name__ == '__main__':
    new_genome = sys.argv[1] # Actually fetched from front end
    print main(new_genome=new_genome)
