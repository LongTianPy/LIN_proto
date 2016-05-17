#!/usr/bin/python
"""
This script sends an email to the user with the result by k-mer calculation.
The result will include the genomes and corresponding LINs
"""

import email
import smtplib
from email.mime.text import MIMEText
import os

def write_kmer_result(top10,db_cursor):
    c = db_cursor
    f = open("/home/linproject/Workspace/email_content/kmer.txt","w")
    f.write("After analyzing the k-mer profile of your uploaded genome by fragmenting it into 12-mer, the following "
            "records are decided to be further analyzed by calculating the Average Nucleotide Identity (ANI) "
            "with your uploaded genome, which may take a while.\n\n")
    for i in top10:
        c.execute("SELECT LIN.LIN from LIN, Genome where LIN.Genome_ID=Genome.Genome_ID and Genome.GenomeName='{0}'".format(i))
        tmp = c.fetchone()
        tmp = tmp[0] # Now we get its LIN
        f.write(i+"\t\t"+tmp+"\n")
    f.close()
    os.system("chmod 777 /home/linproject/Workspace/email_content/kmer.txt")

def write_ANI_result(new_genomeID, new_LIN_object, new_LIN, db_cursor):
    new_GenomeName = new_genomeID
    db_cursor.execute("SELECT LIN.SubjectGenome, LIN.ANI FROM LIN,Genome WHERE LIN.Genome_ID=Genome.Genome_ID "
                      "and Genome.GenomeName='{0}'".format(new_GenomeName))
    tmp = db_cursor.fetchone()
    best_hit = tmp[0]
    ANI_best_hit = str(float(tmp[1])*100)+'%'
    LIN_best_hit = new_LIN_object.LIN
    f = open("/home/linproject/Workspace/email_content/ANI.txt","w")
    f.write("The final result of your recent submission is here, by calculateing the Average Nucleotide Identity (ANI) "
            "between your submission and those best hit candidates chosen according to k-mer profile.\n\n")
    f.write("The result of your submission:\n")
    f.write("Genome: {0}\t\tLIN: {1}\n".format(new_GenomeName, new_LIN))
    f.write("It is most similar with {0}, whose LIN is {1}, with the ANI of {2}.\n\n".format(best_hit, LIN_best_hit, ANI_best_hit))
    f.write("A result web-page is generating for your submission. We will notify you via E-mail once it is done.\n")
    f.close()
    os.system("chmod 777 /home/linproject/Workspace/email_content/ANI.txt")


def send_email(file_source, db_cursor, User_ID=1):
    """
    Send the user an email about the intermediate result
    :param User_ID: Should be able to read from front end, default = 2, which is me.
    :return:
    """
    file_switch = {"kmer":"/home/linproject/Workspace/email_content/kmer.txt","ANI":"/home/linproject/Workspace/email_content/ANI.txt"}
    filepath = file_switch[file_source]
    fp = open(filepath,"r")
    msg = MIMEText(fp.read())
    fp.close()
    me = "longtian@vt.edu"
    db_cursor.execute("SELECT LastName, Email from User where User_ID={0}".format(User_ID))
    tmp = db_cursor.fetchone()
    LastName = tmp[0]
    you = tmp[1]
    if file_source == "kmer":
        msg['Subject'] = "Mr. {0}, here's the preliminary result of your recently submission".format(LastName)
    else: # file_source == "ANI":
        msg['Subject'] = "Mr. {0}, here's the final result of your recently submission".format(LastName)
    s = smtplib.SMTP('smtp.gmail.com',587)
    s.ehlo_or_helo_if_needed()
    s.starttls()
    s.login(me, '!Lt19881120')
    s.sendmail(me, you, msg.as_string())
    s.close()
