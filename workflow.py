#!/usr/bin/python
"""
This script implements functions of fast k-mer calculation to narrow search subjects, ANIb calculation for accurate
genomic percentage identity and LIN assignment to the new genome.
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
import filecmp
import uuid
import sendEmail
import shutil

# import ExtractInfo

sourmash_dir = "/home/linproject/Workspace/Sourmash/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash/rep_bac/"

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
    c.execute('USE LINdb_RefSeq')
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("#####################################")
    logging.info("New task from User ID {0}".format(User_ID))

    original_folder  = '/home/linproject/Workspace/Psy_166/'
    workspace_dir = '/home/linproject/Workspace/New/workspace/'
    # InfoFile = "/home/linproject/Workspace/Zika/Attribute_full.csv"
    # Update the submission table
    # c.execute("INSERT INTO Submission (User_ID, Time) VALUES ({0},'{1}')".format(User_ID, standardtime))
    # db.commit()
    # # And get the new Submittion ID
    # c.execute("SELECT Submission_ID FROM Submission where User_ID={0} and Time='{1}'".format(User_ID,standardtime))
    # tmp = c.fetchone()
    # Submission_ID = int(tmp[0])
    # # And we have the file name of the genome
    # # Fetched from the front end
    # new_GenomeName = ".".join(new_genome.split('.')[:-1])
    # contig_number = get_contig_number(original_folder+new_genome)
    # # As well as its Interest_ID
    # c.execute('INSERT INTO Genome (Interest_ID, Submission_ID, FilePath, GenomeName) values ({0}, {1}, "{2}", "{3}")'
    #           .format(Interest_ID_new_genome ,Submission_ID, original_folder+new_genome, new_GenomeName))
    # db.commit()
    # c.execute('select Genome_ID from Genome where GenomeName="{0}"'.format(new_GenomeName))
    # tmp = c.fetchone()
    # new_Genome_ID = str(tmp[0])
    # # Load Attribute values
    # # First we need to know which field is which
    # c.execute("SELECT Attribute_IDs FROM Interest WHERE Interest_ID={0}".format(Interest_ID_new_genome))
    # tmp = c.fetchone()
    # Attribute_ID_list = tmp[0].split(",")
    # Attribute_ID_list = [int(id) for id in Attribute_ID_list]
    # Attributes = Attributes.split("^^")
    # for i in range(len(Attribute_ID_list)):
    #     c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) "
    #               "VALUES ({0}, {1}, {2}, '{3}', {4}, {5})".format(Attribute_ID_list[i], new_Genome_ID,
    #                                                                Interest_ID_new_genome, Attributes[i],
    #                                                                User_ID, privacy))
    # First check if we ever have this file, roughly
    c.execute("SELECT Genome_ID, FilePath from Genome")
    tmp = c.fetchall()
    # FilePath_table = pd.DataFrame()
    # Genome_IDs = [int(i[0]) for i in tmp]
    # FilePath_table["FilePath"] = [i[1] for i in tmp]
    # FilePath_table.index=Genome_IDs
    # Identical_File = False
    # for each_recorded_genome in Genome_IDs:
    #     if filecmp.cmp(original_folder+new_genome, FilePath_table.get_value(each_recorded_genome,"FilePath")):
    #         Identical_File = each_recorded_genome
    #     else:
    #         continue
    # if Identical_File:
    #     identical_genome_id = int(Identical_File)
    # else:
    #     # Update the submission table
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

    new_LIN, SubjectGenome, ANIb_result,new_SigPath,conserved_LIN = mash_indexing.mash_indexing(cursor=c, new_Genome_ID=new_Genome_ID,User_ID=User_ID,conn=db)
    logging.info("The LIN assigned to your genome is " + new_LIN)
    # c.execute("INSERT INTO LIN (Genome_ID, Scheme_ID, LIN, SubjectGenome, ANI) values ({0}, 4, '{1}', '{2}', {3})"
    #           .format(new_Genome_ID, new_LIN, top1_Genome_ID, top1_similarity))
    # db.commit()
    new_LINgroup = ",".join(new_LIN.split(",")[:6])
    if not isdir(sourmash_dir + new_LINgroup):
        os.mkdir(sourmash_dir + new_LINgroup)
        shutil.copy(new_SigPath, sourmash_dir + "rep_bac/")
    shutil.copy(new_SigPath, sourmash_dir + new_LINgroup + "/")
    c.execute("select LIN_ID from LIN where Scheme_ID=4 and LIN='{0}'".format(new_LIN))
    LIN_ID = int(c.fetchone()[0])
    Job_uuid = str(uuid.uuid4())
    c.execute("insert into Job (LIN_ID,User_ID,Job_uuid,Conserved_LIN) values ({0},{1},'{2}','{3}')".format(LIN_ID,User_ID,Job_uuid,conserved_LIN))
    db.commit()

    c.execute("select Email from User where User_ID={0}".format(User_ID))
    user_email = c.fetchone()[0]
    email_cmd = "python /home/linproject/Projects/LIN_proto/sendEmail.py {0} Submission_result {1}".format(user_email,Job_uuid)
    os.system(email_cmd)
    c.close()
    db.close()
    logging.info("Task completed.")
    logging.info("#####################################")




if __name__ == '__main__':
    main()

