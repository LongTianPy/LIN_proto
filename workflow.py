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
import logging
import logging.handlers
import argparse
from datetime import datetime
from pytz import timezone
from Bio import SeqIO

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

def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))

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

def get_contig_number(fastafile):
    f = open(fastafile,"r")
    records = list(SeqIO.parse(f,"fasta"))
    f.close()
    return len(records)

# MAIN
def main(argv=None): # The genome file name we are expecting for is a
    # There should be a script that uploads the genome sequence to a subfolder in the workspace where the name is
    # randomly and uniquely generated
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()
    new_genome = args.new_genome
    User_ID = int(args.User_ID)
    Interest_ID_new_genome = int(args.Interest_ID)
    Attributes = args.Attributes
    privacy = args.privacy
    if privacy == "True":
        privacy = 1
    else:
        privacy = 0
    eastern = timezone("EST")
    currenttime = eastern.localize(datetime.now())
    fmt_time_display = '%Y-%m-%d %H:%M:%S %Z%z'
    standardtime = currenttime.strftime(fmt_time_display)
    db = Connect('localhost','root')
    c = db.cursor()
    logging.info("Connecting to the database")
    c.execute('USE LINdb_Psy')
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("#####################################")
    logging.info("New task from User ID {0}".format(User_ID))

    original_folder  = '/home/linproject/Workspace/Psy_166/'
    workspace_dir = '/home/linproject/Workspace/New/workspace/'
    # InfoFile = "/home/linproject/Workspace/Zika/Attribute_full.csv"
    # Update the submission table
    c.execute("INSERT INTO Submission (User_ID, Time) VALUES ({0},'{1}')".format(User_ID, standardtime))
    db.commit()
    # And get the new Submittion ID
    c.execute("SELECT Submission_ID FROM Submission where User_ID={0} and Time='{1}'".format(User_ID,standardtime))
    tmp = c.fetchone()
    Submission_ID = int(tmp[0])
    # And we have the file name of the genome
    # Fetched from the front end
    new_GenomeName = new_genome.split('.')[0]
    contig_number = get_contig_number(original_folder+new_genome)
    # As well as its Interest_ID
    c.execute('INSERT INTO Genome (Interest_ID, Submission_ID, FilePath, GenomeName) values ({0}, {1}, "{2}", "{3}")'
              .format(Interest_ID_new_genome ,Submission_ID, original_folder+new_genome, new_GenomeName))
    db.commit()
    c.execute('select Genome_ID from Genome where GenomeName="{0}"'.format(new_GenomeName))
    tmp = c.fetchone()
    new_Genome_ID = str(tmp[0])
    # Load Attribute values
    # First we need to know which field is which
    c.execute("SELECT Attribute_IDs FROM Interest WHERE Interest_ID={0}".format(Interest_ID_new_genome))
    tmp = c.fetchone()
    Attribute_ID_list = tmp[0].split(",")
    Attribute_ID_list = [int(id) for id in Attribute_ID_list]
    Attributes = Attributes.split("^^")
    for i in range(len(Attribute_ID_list)):
        c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) "
                  "VALUES ({0}, {1}, {2}, '{3}', {4}, {5})".format(Attribute_ID_list[i], new_Genome_ID,
                                                                   Interest_ID_new_genome, Attributes[i],
                                                                   User_ID, privacy))
        db.commit()

    # And Also, for reading these attributes from front end, a table called Genome_to_Attribute is created with
    # columns of Genome_IDs and each Attributes, This maybe tricky because if we have a new attributes, it's hard
    # to extend the table, but it's always possible to extract the existing data and put it in a new table.
    # Attribute_Names = []
    # for i in Attribute_ID_list:
    #     c.execute("SELECT AttributeName FROM Attribute WHERE Attribute_ID={0}".format(i))
    #     tmp = c.fetchone()
    #     this_Attribute_Name = "_".join(tmp[0].split(" "))
    #     Attribute_Names.append(this_Attribute_Name)
    c.execute("SELECT Attribute_IDs FROM Interest WHERE Interest_ID={0}".format(Interest_ID_new_genome))
    tmp = c.fetchone()
    Attribute_ID_list_string = tmp[0]
    c.execute("Select AttributeName from Attribute WHERE Attribute_ID in ({0})".format(Attribute_ID_list_string))
    tmp = c.fetchall()
    # Attribute_Names = ['_'.join(i[0].split(' ')) for i in tmp]
    Attribute_Names = [i[0].replace(" ","_").replace("/","__").replace("-","_") for i in tmp]
    Attributes_in_Genome_to_Attribute = "','".join(Attributes)
    AttributeNames_in_Genome_to_Attribute = ",".join(Attribute_Names)
    Attributes_in_Genome_to_Attribute = "'" + Attributes_in_Genome_to_Attribute + "'"
    # print "INSERT INTO Genome_to_Attribute ({0}, Genome_ID) VALUES ({1}, {2})".format(AttributeNames_in_Genome_to_Attribute, Attributes_in_Genome_to_Attribute, new_Genome_ID)
    c.execute("INSERT INTO Genome_to_Attribute ({0}, Genome_ID) VALUES ({1}, {2})".
              format(AttributeNames_in_Genome_to_Attribute, Attributes_in_Genome_to_Attribute, new_Genome_ID))
    db.commit()

    ## For Zika virus case only, comment out when initialization is done.
    # LoadInfo(InfoFile,c,new_GenomeName,Interest_ID_new_genome)
    # db.commit()
    # # Fetch the file paths of all the genomes from the database that have the same interest ID
    # c.execute('SELECT FilePath FROM Genome WHERE Interest_ID = {0}'.format(Interest_ID_new_genome))
    # FilePaths = c.fetchall()
    # FilePaths = [i[0] for i in FilePaths]
    # We need first to read their k-mer frequencies, which, are calculated beforehand and store in the server.
    # This reminds me of adding one more table for the location of those frequency files.
    ## And supposedly, the frequency of new genome shuold also be added to this file once it's calculated.
    similarity = k_mer.generate_distance(queryfilepath=original_folder+new_genome,Genome_ID=new_Genome_ID,User_ID=User_ID)
    # Check the biggest value, if it is bigger than the bottomline of the cutoff being used
    # if similarity[new_GenomeName].max() < 0.6:
    #     print "No similar genome found, run ANIb calculation sequentially to all genome is recommended."
    #     sys.exit()
    # else:
    logging.info("Looking for similars genome from our database.")
    if len(similarity['Genome'])<=10:
        n_top = len(similarity['Genome'])
    else: # As for selecting numbers of clusters, 3 might not be very good since it is possible when we have a lot of
          # genomes to select while there might be way more than 10 genomes in the top cluster. So we are setting the
          # the number of clusters based on the total number of genomes to select from. Since 10 is a preferred number
          # of genomes to perform pairwise blasting, I guess we can do (M/10)+1
        n_clusters = 3
        km = KMeans(n_clusters=n_clusters)
        km.fit(similarity[new_Genome_ID].reshape(-1,1))
        centroid_idx = list(km.cluster_centers_).index(max(km.cluster_centers_))
        top_cluster_idx = [i for i,x in enumerate(km.labels_) if x==centroid_idx]
        n_top = len(top_cluster_idx)
        logging.info("We are comparing your genome with {0} genomes in our database.".format(n_top))
    top10 = similarity.head(n_top)['Genome'].values # top10 is a list of Genome_IDs in the database

    logging.info("Writing k-mer result.")
    IntermediateResult.write_kmer_result(top10=top10,db_cursor=c,User_ID=User_ID)
    logging.info("Sending k-mer result to the user, User_ID {0}.".format(User_ID))

    #################################################################################
    # IntermediateResult.send_email(file_source="kmer",User_ID=User_ID,db_cursor=c)
    #################################################################################

    # top10_LINs = [ExtractInfo.get_top10_LIN(i,c) for i in top10] # This can be used to send preliminary results
    # print top10_LINs
    # Get their file paths and copy them to the workspace
    similarities = pd.DataFrame()

    logging.info("Iteratively calculating ANIs.")
    for i in top10:
        c.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(i))
        target_filepath = c.fetchone()[0]
        target_filename_rename = target_filepath.split("/")[-1]
        target_filename_rename = "{0}.fasta".format(i)
        cmd = "cp {0} {1}".format(target_filepath, workspace_dir+target_filename_rename)
        os.system(cmd)
        c.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(int(new_Genome_ID)))
        query_filepath = c.fetchone()[0]
        query_filename_rename = query_filepath.split("/")[-1]
        query_filename_rename = "{0}.fasta".format(new_Genome_ID)
        os.system('cp {0} {1}'.format(query_filepath,workspace_dir+query_filename_rename))
        # Now we have all of them in the workspace
        ANIb_result = ANI_Wrapper_2.unified_anib(workspace_dir,User_ID)[new_Genome_ID]
        os.system('rm -rf {0}*'.format(workspace_dir))
        similarity = ANIb_result.loc[str(i)]
        similarities[str(i)]=[similarity]
    top1_genome = similarities.idxmax(axis=1)[0]
    top1_Genome_ID = int(top1_genome)
    top1_similarity = similarities.max(axis=1)[0]
    new_LIN_object = LIN_Assign.getLIN(Genome_ID=top1_Genome_ID, Scheme_ID=3, similarity=top1_similarity,c=c)
    logging.info("This most similar record is Genome ID {0}, whose LIN is {1}.".format(top1_Genome_ID,
                                                                             ','.join(new_LIN_object.LIN)))
    logging.info("The similarity to it is " + str(top1_similarity*100) + "%.")
    new_LIN = LIN_Assign.Assign_LIN(new_LIN_object,c=c).new_LIN
    logging.info("The LIN assigned to your genome is " + new_LIN)
    c.execute("INSERT INTO LIN (Genome_ID, Scheme_ID, LIN, SubjectGenome, ANI) values ({0}, 3, '{1}', '{2}', {3})"
              .format(new_Genome_ID, new_LIN, top1_Genome_ID, top1_similarity))
    db.commit()

    url = IntermediateResult.write_result_page(new_Genome_ID=new_Genome_ID,new_LIN_object=new_LIN_object,new_LIN=new_LIN,db_cursor=c, User_ID=User_ID,Interest_ID=Interest_ID_new_genome)
    IntermediateResult.write_ANI_result(new_Genome_ID=new_Genome_ID,new_LIN_object=new_LIN_object,new_LIN=new_LIN,db_cursor=c,User_ID=User_ID,url=url,Interest_ID=Interest_ID_new_genome)
    # IntermediateResult.send_email(file_source="ANI",User_ID=User_ID,db_cursor=c)

    c.close()
    db.close()
    logging.info("Task completed.")
    logging.info("#####################################")




if __name__ == '__main__':
    main()

