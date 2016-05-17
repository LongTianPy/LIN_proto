#!/usr/bin/python
"""
This script sends an email to the user with the result by k-mer calculation.
The result will include the genomes and corresponding LINs
"""

import email
import smtplib
from email.mime.text import MIMEText

def write_kmer_result(top10,db_cursor):
    c = db_cursor
    f = open("/home/linproject/Workspace/email_content/tmp.txt","w")
    f.write("After analyzing the k-mer profile of your uploaded genome by fragmenting it into 12-mer, the following "
            "records are decided to be further analyzed by calculating the Average Nucleotide Identity (ANI) "
            "with your uploaded genome, which may take a while.\n")
    for i in top10:
        c.execute("SELECT LIN.LIN from LIN, Genome where LIN.Genome_ID=Genome.Genome_ID and Genome.GenomeName='{0}'".format(i))
        tmp = c.fetchone()
        tmp = tmp[0] # Now we get its LIN
        f.write(i+"\t\t"+tmp+"\n")
    f.close()

def send_email(User_ID=2, db_cursor):
    """
    Send the user an email about the intermediate result
    :param User_ID: Should be able to read from front end, default = 2, which is me.
    :return:
    """
    fp = open("/home/linproject/Workspace/email_content/tmp.txt","r")
    msg = MIMEText(fp.read())
    fp.close()

    me = "longtian@vt.edu"
    c = db_cursor
    c.execute("SELECT LastName, Email from User where User_ID={0}".format(User_ID))
    tmp = c.fetone()
    LastName = tmp[0]
    you = tmp[1]
    msg['Subject'] = "Mr. {0}, here's the preliminary result of your recently submission".format(LastName)

    s = smtplib.SMTP('smtp.gmail.com',587)
    s.ehlo_or_helo_if_needed()
    s.starttls()
    s.login(me, '!Lt19881120')
    s.sendmail(me, you, msg.as_string())
    s.close()
