#!/usr/bin/python
"""
This script sends an email to the user with the result by k-mer calculation.
The result will include the genomes and corresponding LINs
"""

import email
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import os
import logging

def write_kmer_result(top10,db_cursor,User_ID):
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    c = db_cursor
    f = open("/home/linproject/Workspace/email_content/kmer.txt","w")
    f.write("<html><head>\n")
    f.write("<style>table,td,th{border:1px solid black;}"
            "table{border-collapse:collapse;width:100%}"
            "th{height:20px;}"
            "</style></head>")
    f.write("<body><p>After analyzing the k-mer profile of your uploaded genome by fragmenting it into 12-mer, the following "
            "records are decided to be further analyzed by calculating the Average Nucleotide Identity (ANI) "
            "with your uploaded genome, which may take a while.</p>\n\n")
    f.write("<table>\n")
    f.write("<tr><th align='left>Genus</th><th align='left>Species</th><th align='left>Strain</th>"
            "<th align='left>A</th><th align='left>B</th><th align='left>C</th><th align='left>D</th><th align='left>E</th>"
            "<th align='left>F</th><th align='left>G</th><th align='left>H</th><th align='left>I</th><th align='left>J</th>"
            "<th align='left>K</th><th align='left>L</th><th align='left>M</th><th align='left>N</th><th align='left>O</th>"
            "<th align='left>P</th><th align='left>Q</th><th align='left>R</th><th align='left>S</th><th align='left>T</th>"
            "</tr>\n")
    for i in top10:
        db_cursor.execute("SELECT AttributeValue FROM AttributeValue WHERE Genome_ID={0} AND Attribute_ID IN (1,4,5)".format(int(i)))
        tmp = c.fetchall() # By which, we will get a list of 3 elements where for each element, 0 is an attribute, 1 is LIN
        if len(tmp) == 0:
            Genus = Species = Strain = "N/A"
        else:
            Genus = tmp[1][0]
            Species = tmp[2][0]
            Strain = tmp[0][0]
        db_cursor.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(int(i)))
        tmp = c.fetchone()
        LIN = tmp[0] # Could also be tmp[1][1] or tmp[2][1]
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
    db_cursor.execute("SELECT AttributeValue from AttributeValue WHERE Genome_ID={0} AND Attribute_ID "
                      "IN (1,4,5)".format(new_Genome_ID))
    tmp=db_cursor.fetchall()
    if len(tmp) == 0:
        Genus_new_Genome = Species_new_Genome = Strain_new_Genome = "N/A"
    else:
        Genus_new_Genome = tmp[1][0]
        Species_new_Genome = tmp[2][0]
        Strain_new_Genome = tmp[0][0]
    LIN_new_Genome = new_LIN

    # Get info of the best match
    Genome_ID_best_hit = new_LIN_object.Genome_ID
    db_cursor.execute("SELECT AttributeValue from AttributeValue WHERE Genome_ID={0} AND Attribute_ID "
                      "IN (1,4,5)".format(Genome_ID_best_hit))
    if len(tmp) == 0:
        Genus_best_hit = Species_best_hit = Strain_best_hit = "N/A"
    else:
        Genus_best_hit = tmp[1][0]
        Species_best_hit = tmp[2][0]
        Strain_best_hit = tmp[0][0]
    LIN_best_hit = new_LIN_object.LIN  # This is a list already

    # Get the Genome_IDs of all those sharing the same conserved part of LINs
    db_cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE LIN LIKE '{0}%' AND Genome_ID <> {1} and Genome_ID <> {2}".
                      format(",".join(new_LIN_object.conserved_LIN),new_Genome_ID,Genome_ID_best_hit))
    tmp = db_cursor.fetchall()
    Genome_IDs_related_hits = [int(i[0]) for i in tmp]
    LINs_related_hits = [i[1] for i in tmp]

    f = open("/home/linproject/Workspace/email_content/ANI.txt","w")
    f.write("<html><head>\n")
    f.write("<style>table,td,th{border:1px solid black;}"
            "table{border-collapse:collapse;width:100%}"
            "th{height:20px;}"
            "</style></head>")
    f.write("<body><p>The final result of your recent submission is here, by calculateing the Average Nucleotide Identity (ANI) "
            "between your submission and those best hit candidates chosen according to k-mer profile.\nThe ANI between"
            "your submission and the best match is {0}.</p>\n\n".format(ANI_best_hit))
    f.write("<h2>The result of your submission:</h2>\n")
    f.write("<table>\n")
    f.write("<tr><th align='left>Category</th><th align='left>Genus</th><th align='left>Species</th><th align='left>Strain</th>"
            "<th align='left>A</th><th align='left>B</th><th align='left>C</th><th align='left>D</th><th align='left>E</th>"
            "<th align='left>F</th><th align='left>G</th><th align='left>H</th><th align='left>I</th><th align='left>J</th>"
            "<th align='left>K</th><th align='left>L</th><th align='left>M</th><th align='left>N</th><th align='left>O</th>"
            "<th align='left>P</th><th align='left>Q</th><th align='left>R</th><th align='left>S</th><th align='left>T</th>"
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
        if len(tmp) == 0:
            Genus_related_hit = Species_related_hit = Strain_related_hit = "N/A"
        else:
            Genus_related_hit = tmp[1][0]
            Species_related_hit = tmp[2][0]
            Strain_related_hit = tmp[0][0]
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
    msg = MIMEMultipart('alternative')
    fp = open(filepath,"r")
    content = MIMEText(fp.read(),"html")
    fp.close()
    msg.attach(content)
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
