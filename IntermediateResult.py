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
import uuid
import pandas as pd

def write_kmer_result(top10,db_cursor,User_ID):
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    c = db_cursor
    f = open("/home/linproject/Workspace/email_content/kmer.txt","w")
    f.write("<html><head>\n")
    f.write("<style>"
            "table{border-collapse:collapse;width:100%}"
            "th{width:1.5em;height:20px;}"
            "</style></head>")
    f.write("<body><p>After analyzing the k-mer profile of your uploaded genome by fragmenting it into 12-mer, the following "
            "records are decided to be further analyzed by calculating the Average Nucleotide Identity (ANI) "
            "with your uploaded genome, which may take a while.</p>\n\n")
    f.write("<table>\n")
    f.write("<tr><th align='left'>Genus</th><th align='left'>Species</th><th align='left'>Strain</th>"
            "<th align='left'>A</th><th align='left'>B</th><th align='left'>C</th><th align='left'>D</th><th align='left'>E</th>"
            "<th align='left'>F</th><th align='left'>G</th><th align='left'>H</th><th align='left'>I</th><th align='left'>J</th>"
            "<th align='left'>K</th><th align='left'>L</th><th align='left'>M</th><th align='left'>N</th><th align='left'>O</th>"
            "<th align='left'>P</th><th align='left'>Q</th><th align='left'>R</th><th align='left'>S</th><th align='left'>T</th>"
            "</tr>\n")
    for i in top10:
        c.execute("SELECT AttributeValue FROM AttributeValue WHERE Genome_ID={0} AND Attribute_ID IN (1,2,3)".format(int(i)))
        tmp = c.fetchall() # By which, we will get a list of 3 elements where for each element, 0 is an attribute, 1 is LIN
        if len(tmp) == 0:
            Genus = Species = Strain = "N/A"
        else:
            Genus = tmp[0][0]
            Species = tmp[1][0]
            Strain = tmp[2][0]
        c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(int(i)))
        tmp = c.fetchone()
        LIN = tmp[0] # Could also be tmp[1][1] or tmp[2][1]
        LIN = LIN.split(",")
        f.write("<tr><td>{0}</td><td>{1}</td><td>{2}</td>".format(Genus,Species,Strain))
        for lin in LIN:
            f.write("<td>{0}</td>".format(lin))
        f.write("</tr>\n")
    f.write("</table></body></html>")
    f.close()

def write_ANI_result(new_Genome_ID, new_LIN_object, new_LIN, db_cursor,User_ID,url,Interest_ID):
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    # Get Attribute_IDs
    db_cursor.execute("SELECT Attribute_IDs FROM Interest WHERE Interest_ID={0}".format(Interest_ID))
    Attribute_IDs = db_cursor.fetchone()[0]
    # Get AttributeName
    db_cursor.execute("SELECT AttributeName FROM Attribute WHERE Attribute_ID in ({0})".format(Attribute_IDs))
    tmp = db_cursor.fetchall()
    AttributeName = [i[0] for i in tmp]
    AttributeName_in_table = ['_'.join(i[0].split(' ')) for i in tmp]
    AttributeName_in_table_w_o_gss = [i for i in AttributeName_in_table if
                                      i != 'Genus' and i != 'Strain' and i != 'Species']
    df = {}
    for i in AttributeName_in_table:
        df[i] = []
    keys = df.keys()  # In case the order changes
    AttributeName_string = ','.join(keys)

    # Let's first get the info of the new submission
    db_cursor.execute(
        "SELECT {0} FROM Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string, new_Genome_ID))
    new_genome_meta = db_cursor.fetchall()[0]
    for i in range(len(keys)):
        df[keys[i]].append(new_genome_meta[i])

    # The second row is the best match
    # But we need the Genome_ID of it
    db_cursor.execute(
        "SELECT SubjectGenome,ANI,LIN FROM LIN where Genome_ID={0}".format(new_Genome_ID)
    )
    tmp = db_cursor.fetchall()[0]
    Genome_ID_best_match = tmp[0]
    ANI = '{0:.5f}%'.format(float(tmp[1]) * 100)
    LIN_new_Genome = tmp[2]  # A string
    db_cursor.execute(
        "SELECT LIN FROM LIN WHERE Genome_ID={0}".format(Genome_ID_best_match)
    )
    LIN_best_match = db_cursor.fetchone()[0]  # A string
    db_cursor.execute(
        "SELECT {0} FROM Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string, Genome_ID_best_match)
    )
    best_match_meta = db_cursor.fetchall()[0]
    for i in range(len(keys)):
        df[keys[i]].append(best_match_meta[i])

    f = open("/home/linproject/Workspace/email_content/ANI.txt","w")
    f.write("<html><head>\n")
    f.write("<style>"
            "table{border-collapse:collapse;width:100%}"
            "th{width:1.5em;height:20px;}"
            "</style></head>")
    f.write("<body><p>The final result of your recent submission is here, by calculateing the Average Nucleotide Identity (ANI) "
            "between your submission and those best hit candidates chosen according to k-mer profile.\nThe ANI between"
            "your submission and the best match is {0}.</p>\n\n".format(ANI))
    f.write("<h2>The result of your submission:</h2>\n")
    f.write("<table>\n")
    f.write("<tr><th align='left'>Category</th><th align='left'>Genus</th><th align='left'>Species</th><th align='left'>Strain</th>"
            "<th align='left'>A</th><th align='left'>B</th><th align='left'>C</th><th align='left'>D</th><th align='left'>E</th>"
            "<th align='left'>F</th><th align='left'>G</th><th align='left'>H</th><th align='left'>I</th><th align='left'>J</th>"
            "<th align='left'>K</th><th align='left'>L</th><th align='left'>M</th><th align='left'>N</th><th align='left'>O</th>"
            "<th align='left'>P</th><th align='left'>Q</th><th align='left'>R</th><th align='left'>S</th><th align='left'>T</th>"
            "</tr>\n")
    f.write("<tr><td>New Submission</td><td>{0}</td><td>{1}</td><td>{2}</td>".format(df["Genus"][0],df["Species"][0],df["Strain"][0]))
    for each_single_LIN in LIN_new_Genome.split(','):
        f.write("<td>{0}</td>\n".format(each_single_LIN))
    f.write("</tr>\n")

    f.write("<tr><td>Best match</td><td>{0}</td><td>{1}</td><td>{2}</td>".format(df["Genus"][1],df["Species"][1],df["Strain"][1]))
    for each_single_LIN in LIN_best_match.split(','):
        f.write("<td>{0}</td>\n".format(each_single_LIN))
    f.write("</tr>\n")

    f.write("</table>")
    f.write("<p>you can visit the following page to check more details and add descriptions. {0}</p>".format(url))
    f.write("</body></html>")
    f.close()

def write_result_page(new_Genome_ID, new_LIN_object, new_LIN, db_cursor,User_ID,Interest_ID):
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    # Get Attribute_IDs
    db_cursor.execute("SELECT Attribute_IDs FROM Interest WHERE Interest_ID={0}".format(Interest_ID))
    Attribute_IDs = db_cursor.fetchone()[0]
    # Get AttributeName
    db_cursor.execute("SELECT AttributeName FROM Attribute WHERE Attribute_ID in ({0})".format(Attribute_IDs))
    tmp = db_cursor.fetchall()
    AttributeName = [i[0] for i in tmp]
    AttributeName_in_table = ['_'.join(i[0].split(' ')) for i in tmp]
    AttributeName_in_table_w_o_gss = [i for i in AttributeName_in_table if i != 'Genus' and i != 'Strain' and i != 'Species']
    df = {}
    for i in AttributeName_in_table:
        df[i] = []
    keys = df.keys() # In case the order changes
    AttributeName_string = ','.join(keys)

    # Let's first get the info of the new submission
    db_cursor.execute(
        "SELECT {0} FROM Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string, new_Genome_ID))
    new_genome_meta = db_cursor.fetchall()[0]
    for i in range(len(keys)):
        df[keys[i]].append(new_genome_meta[i])

    # The second row is the best match
    # But we need the Genome_ID of it
    db_cursor.execute(
        "SELECT SubjectGenome,ANI,LIN FROM LIN where Genome_ID={0}".format(new_Genome_ID)
    )
    tmp = db_cursor.fetchall()[0]
    Genome_ID_best_match = tmp[0]
    ANI = '{0:.5f}%'.format(float(tmp[1])*100)
    LIN_new_Genome = tmp[2] # A string
    db_cursor.execute(
        "SELECT LIN FROM LIN WHERE Genome_ID={0}".format(Genome_ID_best_match)
    )
    LIN_best_match = db_cursor.fetchone()[0] # A string
    db_cursor.execute(
        "SELECT {0} FROM Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string, Genome_ID_best_match)
    )
    best_match_meta = db_cursor.fetchall()[0]
    for i in range(len(keys)):
        df[keys[i]].append(best_match_meta[i])

    # And the related hits sharing the same part of LIN
    conserved_LIN = ','.join(new_LIN_object.conserved_LIN)
    db_cursor.execute(
        "SELECT Genome_ID, LIN from LIN where LIN like '{0}%' and Genome_ID NOT in ({1},{2}) ORDER BY LIN.LIN".format(conserved_LIN,new_Genome_ID,Genome_ID_best_match)
    )
    tmp = db_cursor.fetchall()
    if len(tmp) != 0:
        related_hits = True
        Genome_IDs_related_hits = [int(i[0]) for i in tmp]
        LINs_related_hits = [i[1] for i in tmp]
        for i in Genome_IDs_related_hits:
            db_cursor.execute(
                "SELECT {0} FROM Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string,i)
            )
            related_hit_meta = db_cursor.fetchall()[0]
            for j in range(len(keys)):
                df[keys[j]].append(related_hit_meta[j])
    else:
        related_hits = False


    unique_filename = str(uuid.uuid4())
    filename = '/var/www/html/CodeIgniter/resultPages/' + unique_filename + '.php'
    url = 'http://128.173.74.68/CodeIgniter/resultPages/' + unique_filename + '.php'
    # Read the static page header
    f = open('/var/www/html/CodeIgniter/resultPages/Before','r')
    static_header = f.read()
    f.close()
    # Read the static page footer
    f = open('/var/www/html/CodeIgniter/resultPages/After','r')
    static_footer = f.read()
    f.close()
    f = open(filename,'a')
    f.write(static_header)

    # Popup info
    for i in range(len(df['Genus'])):
        popup_div_id = 'detail_'+str(i)
        f.write("<div id='detail_{0}' style='display:none; position:absolute;background-color: #CAE1FF;' class='popup'>".format(str(i)))
        f.write("<p><b>Genus: </b>{0}</p>".format(df["Genus"][i]))
        f.write("<p><b>Species: </b>{0}</p>".format(df["Species"][i]))
        f.write("<p><b>Strain: </b>{0}</p>".format(df["Strain"][i]))
        for each_item in AttributeName_in_table_w_o_gss:
            item = ' '.join(each_item.split("_"))
            f.write("<p><b>{0}: </b>{1}</p>".format(item, df[each_item][i]))
        f.write("</div>")
    # New Submission
    f.write("<tr class='info'><td colspan='1'></td><td rowspan='1' colspan='23' style='text-align: left; font-weight: bold;'>New Submission</td></tr>")
    f.write("<tr id='0'><td style='padding-left: 10px'><input type='checkbox' id='row_checkbox[{0}]' class='row_checkbox' name='row_checkbox[{0}]'"
            "value='{0}'></td>"
            "<td style='padding-left: 1px; width: 100px; overflow-x: auto;'>{1}</td>"
            "<td style='padding-left: 1px; width: 80px; overflow-x: auto;'>{2}</td>"
            "<td style='padding-left: 1px; width: 80px; overflow-x: auto;'>{3}</td>".format(0,df['Genus'][0], df['Species'][0], df['Strain'][0]))
    for i in range(len(LIN_new_Genome.split(','))):
        f.write("<td class='0_LIN_{1}'>{0}</td>".format(LIN_new_Genome.split(',')[i],i))
    f.write("</tr>")
    # Best hit
    f.write("<tr class='info'><td colspan='1'></td><td rowspan='1' colspan='23' style='text-align: left; font-weight: bold;'>Best match, ANI: {0}</td></tr>".format(ANI))
    f.write("<tr id='1'><td style='padding-left: 10px'><input type='checkbox' id='row_checkbox[{0}]' class='row_checkbox' name='row_checkbox[{0}]'"
            "value='{0}'></td>"
            "<td style='padding-left: 1px; width: 100px; overflow-x: auto;'>{1}</td>"
            "<td style='padding-left: 1px; width: 80px; overflow-x: auto;'>{2}</td>"
            "<td style='padding-left: 1px; width: 80px; overflow-x: auto;'>{3}</td>".format(1, df['Genus'][1], df['Species'][1],
                                                                         df['Strain'][1]))
    for i in range(len(LIN_best_match.split(','))):
        f.write("<td class='1_LIN_{1}'>{0}</td>".format(LIN_best_match.split(',')[i],i))
    f.write("</tr>")
    # Related hits, if any
    if related_hits:
        f.write("<tr class='info'><td colspan='1'></td><td rowspan='1' colspan='23' style='text-align: left; font-weight: bold;'>Related record(s)</td></tr>")
        for i in range(len(df["Genus"])-2):
            real_i = i+2
            f.write(
                "<tr id='{0}'><td style='padding-left: 10px'><input type='checkbox' id='row_checkbox[{0}]' class='row_checkbox' name='row_checkbox["
                "{0}]'"
                "value='{0}'></td>"
                "<td style='padding-left: 1px; width: 100px; overflow-x: auto;'>{1}</td>"
                "<td style='padding-left: 1px; width: 80px; overflow-x: auto;'>{2}</td>"
                "<td style='padding-left: 1px; width: 80px; overflow-x: auto;'>{3}</td>".format(real_i, df['Genus'][real_i], df['Species'][real_i],
                                                                             df['Strain'][real_i]))
            for j in range(len(LINs_related_hits[i].split(','))):
                f.write("<td class='{1}_LIN_{2}'>{0}</td>".format(LINs_related_hits.split(',')[j], real_i, j))
            f.write("</tr>")
    f.write(static_footer)
    f.close()
    return url

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
    s = smtplib.SMTP('smtp-mail.outlook.com',587)
    s.ehlo_or_helo_if_needed()
    s.starttls()
    s.login(me, '!Tl19881120')
    s.sendmail(me, you, msg.as_string())
    s.quit()
