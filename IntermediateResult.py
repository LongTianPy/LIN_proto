#!/usr/bin/python
"""
This script sends an email to the user with the result by k-mer calculation.
The result will include the genomes and corresponding LINs
"""

import email
import smtplib
from email.mime.text import MIMEText
import os
import logging

def write_kmer_result(top10,db_cursor,User_ID):
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    c = db_cursor
    f = open("/home/linproject/Workspace/email_content/kmer.txt","w")
    f.write("<html><body>\n")
    f.write("<p>After analyzing the k-mer profile of your uploaded genome by fragmenting it into 12-mer, the following "
            "records are decided to be further analyzed by calculating the Average Nucleotide Identity (ANI) "
            "with your uploaded genome, which may take a while.<p>\n\n")
    f.write("<table style='width:100%'>\n")
    f.write("<tr><th>Genus</th><th>Species</th><th>Strain</th>"
            "<th>A</th><th>B</th><th>C</th><th>D</th><th>E</th>"
            "<th>F</th><th>G</th><th>H</th><th>I</th><th>J</th>"
            "<th>K</th><th>L</th><th>M</th><th>N</th><th>O</th>"
            "<th>P</th><th>Q</th><th>R</th><th>S</th><th>T</th>"
            "</tr>\n")
    for i in top10:
        logging.info("SELECT AttributeValue.AttributeValue, LIN.LIN from LIN, AttributeValue where "
                  "LIN.Genome_ID=AttributeValue.Genome_ID and "
                  "LIN.Genome_ID={0} and AttributeValue.Attribute_ID in (1,4,5)".format(int(i)))
        c.execute("SELECT AttributeValue.AttributeValue, LIN.LIN from LIN, AttributeValue where "
                  "LIN.Genome_ID=AttributeValue.Genome_ID and "
                  "LIN.Genome_ID={0} and AttributeValue.Attribute_ID in (1,4,5)".format(int(i)))
        tmp = c.fetchall() # By which, we will get a list of 3 elements where for each element, 0 is an attribute, 1 is LIN
        try:
            Genus = tmp[1][0]
        except:
            Genus = "N/A"
        try:
            Species = tmp[2][0]
        except:
            Species = "N/A"
        try:
            Strain = tmp[0][0]
        except:
            Strain = "N/A"
        LIN = tmp[0][1] # Could also be tmp[1][1] or tmp[2][1]
        LIN = LIN.split(",")
        f.write("<tr><td>{0}</td><td>{1}</td><td>{2}</td>".format(Genus,Species,Strain))
        for lin in LIN:
            f.write("<td>{0}</td>".format(lin))
        f.write("</tr>\n")
    f.write("</table></body></html>")
    f.close()

def write_ANI_result(new_Genome_ID, new_LIN_object, new_LIN, db_cursor,User_ID):
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    db_cursor.execute("SELECT LIN.SubjectGenome, LIN.ANI FROM LIN,Genome WHERE LIN.Genome_ID=Genome.Genome_ID "
                      "and Genome.Genome_ID='{0}'".format(new_Genome_ID))
    tmp = db_cursor.fetchone()
    best_hit = tmp[0]
    ANI_best_hit = str(float(tmp[1])*100)+'%'
    LIN_best_hit = new_LIN_object.LIN

    # Get info of the new submission
    db_cursor.execute("SELECT AtttributeValue.AttributeValue from AttributeValue WHERE Genome_ID={0} AND Attribute_ID "
                      "IN (1,4,5)".format(new_Genome_ID))
    tmp=db_cursor.fetchall()
    try:
        Genus_new_Genome = tmp[1][0]
    except:
        Genus_new_Genome = "N/A"
    try:
        Species_new_Genome = tmp[2][0]
    except:
        Species_new_Genome = "N/A"
    try:
        Strain_new_Genome = tmp[0][0]
    except:
        Strain_new_Genome = "N/A"
    LIN_new_Genome = new_LIN

    # Get info of the best match
    Genome_ID_best_hit = new_LIN_object.Genome_ID
    db_cursor.execute("SELECT AtttributeValue.AttributeValue from AttributeValue WHERE Genome_ID={0} AND Attribute_ID "
                      "IN (1,4,5)".format(Genome_ID_best_hit))
    try:
        Genus_best_hit = tmp[1][0]
    except:
        Genus_best_hit = "N/A"
    try:
        Species_best_hit = tmp[2][0]
    except:
        Species_best_hit = "N/A"
    try:
        Strain_best_hit = tmp[0][0]
    except:
        Strain_best_hit = "N/A"
    LIN_best_hit = new_LIN_object.LIN  # This is a list already

    # Get the Genome_IDs of all those sharing the same conserved part of LINs
    db_cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE LIN LIKE '{0}%' AND Genome_ID <> {1} and Genome_ID <> {2}}".
                      format(new_LIN_object.conserved_LIN,new_Genome_ID,Genome_ID_best_hit))
    tmp = db_cursor.fetchall()
    Genome_IDs_related_hits = [int(i[0]) for i in tmp]
    LINs_related_hits = [i[1] for i in tmp]

    f = open("/home/linproject/Workspace/email_content/ANI.txt","w")
    f.write("<html><body>\n")
    f.write("<p>The final result of your recent submission is here, by calculateing the Average Nucleotide Identity (ANI) "
            "between your submission and those best hit candidates chosen according to k-mer profile.\nThe ANI between"
            "your submission and the best match is {0}.</p>\n\n".format(ANI_best_hit))
    f.write("<h2>The result of your submission:</h2>\n")
    f.write("<table style='width:100%>\n")
    f.write("<tr><th>Category</th><th>Genus</th><th>Species</th><th>Strain</th>"
            "<th>A</th><th>B</th><th>C</th><th>D</th><th>E</th>"
            "<th>F</th><th>G</th><th>H</th><th>I</th><th>J</th>"
            "<th>K</th><th>L</th><th>M</th><th>N</th><th>O</th>"
            "<th>P</th><th>Q</th><th>R</th><th>S</th><th>T</th>"
            "</tr>\n")
    f.write("<tr><td>New Submission</td><td>{0}</td><td>{1}</td><td>{2}</td>".format(Genus_new_Genome,
                                                                                     Species_new_Genome,
                                                                                     Strain_new_Genome))
    for lin in LIN_new_Genome.split(","):
        f.write("<td>{0}</td>".format(lin))
    f.write("</tr>\n")
    f.write("<tr><td>Best match</td><td>{0}</td><td>{1}</td><td>{2}</td>".format(Genus_best_hit,
                                                                                 Species_best_hit,
                                                                                 Strain_best_hit))
    for lin in LIN_best_hit:
        f.write("<td>{0}</td>".format(lin))
    f.write("</tr>\n")
    # Write the related hits, who share the conserved LINs
    for i in range(len(Genome_IDs_related_hits)):
        db_cursor.execute("SELECT AttributeValue FROM AttributeValue WHERE Genome_ID={0} AND Attribute_ID in (1,4,5)"
                          .format(Genome_IDs_related_hits[i]))
        tmp = db_cursor.fetchall()
        try:
            Genus_related_hit = tmp[1][0]
        except:
            Genus_related_hit = "N/A"
        try:
            Species_related_hit = tmp[2][0]
        except:
            Species_related_hit = "N/A"
        try:
            Strain_related_hit = tmp[0][0]
        except:
            Strain_related_hit = "N/A"
        f.write("<tr><td>Similar record</td><td>{0}</td><td>{1}</td><td>{2}</td>".format(Genus_related_hit,
                                                                                         Species_related_hit,
                                                                                         Strain_related_hit))
        for lin in LINs_related_hits[i].split(","):
            f.write("<td>{0}</td>".format(lin))
        f.write("</tr>\n")
    f.write("</table></body></html>")
    f.close()

def send_email(file_source, db_cursor, User_ID):
    """
    Send the user an email about the intermediate result
    :param User_ID: Should be able to read from front end, default = 2, which is me.
    :return:
    """
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    file_switch = {"kmer":"/home/linproject/Workspace/email_content/kmer.txt","ANI":"/home/linproject/Workspace/email_content/ANI.txt"}
    filepath = file_switch[file_source]
    fp = open(filepath,"r")
    msg = MIMEText(fp.read())
    fp.close()
    me = "kingdom586@hotmail.com"
    db_cursor.execute("SELECT LastName, Email from User where User_ID={0}".format(User_ID))
    tmp = db_cursor.fetchone()
    LastName = tmp[0]
    you = tmp[1]
    if file_source == "kmer":
        msg['Subject'] = "Mr. {0}, here's the preliminary result of your recent submission".format(LastName)
    else: # file_source == "ANI":
        msg['Subject'] = "Mr. {0}, here's the final result of your recent submission".format(LastName)
    s = smtplib.SMTP('smtp.live.com',587)
    s.ehlo_or_helo_if_needed()
    s.starttls()
    s.login(me, '1988112019881120')
    s.sendmail(me, you, msg.as_string())
    s.quit()
