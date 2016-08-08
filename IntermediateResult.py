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
        db_cursor.execute("SELECT AttributeValue FROM AttributeValue WHERE Genome_ID={0} AND Attribute_ID IN (1,2,3)".format(int(i)))
        tmp = c.fetchall() # By which, we will get a list of 3 elements where for each element, 0 is an attribute, 1 is LIN
        if len(tmp) == 0:
            Genus = Species = Strain = "N/A"
        else:
            Genus = tmp[0][0]
            Species = tmp[1][0]
            Strain = tmp[2][0]
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

def write_ANI_result(new_Genome_ID, new_LIN_object, new_LIN, db_cursor,User_ID,url):
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
                      "IN (1,2,3)".format(Genome_ID_best_hit))
    if len(tmp) == 0:
        Genus_best_hit = Species_best_hit = Strain_best_hit = "N/A"
    else:
        Genus_best_hit = tmp[0][0]
        Species_best_hit = tmp[1][0]
        Strain_best_hit = tmp[2][0]
    LIN_best_hit = new_LIN_object.LIN  # This is a list already

    # Get the Genome_IDs of all those sharing the same conserved part of LINs
    db_cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE LIN LIKE '{0}%' AND Genome_ID <> {1} and Genome_ID <> {2}".
                      format(",".join(new_LIN_object.conserved_LIN),new_Genome_ID,Genome_ID_best_hit))
    tmp = db_cursor.fetchall()
    Genome_IDs_related_hits = [int(i[0]) for i in tmp]
    LINs_related_hits = [i[1] for i in tmp]

    f = open("/home/linproject/Workspace/email_content/ANI.txt","w")
    f.write("<html><head>\n")
    f.write("<style>"
            "table{border-collapse:collapse;width:100%}"
            "th{width:1.5em;height:20px;}"
            "</style></head>")
    f.write("<body><p>The final result of your recent submission is here, by calculateing the Average Nucleotide Identity (ANI) "
            "between your submission and those best hit candidates chosen according to k-mer profile.\nThe ANI between"
            "your submission and the best match is {0}.</p>\n\n".format(ANI_best_hit))
    f.write("<h2>The result of your submission:</h2>\n")
    f.write("<table>\n")
    f.write("<tr><th align='left'>Category</th><th align='left'>Genus</th><th align='left'>Species</th><th align='left'>Strain</th>"
            "<th align='left'>A</th><th align='left'>B</th><th align='left'>C</th><th align='left'>D</th><th align='left'>E</th>"
            "<th align='left'>F</th><th align='left'>G</th><th align='left'>H</th><th align='left'>I</th><th align='left'>J</th>"
            "<th align='left'>K</th><th align='left'>L</th><th align='left'>M</th><th align='left'>N</th><th align='left'>O</th>"
            "<th align='left'>P</th><th align='left'>Q</th><th align='left'>R</th><th align='left'>S</th><th align='left'>T</th>"
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
    f.write("</table>")
    f.write("<p>You can visit the following page to check more details and add descriptions.{0}</p>".format(url))
    f.write("</body></html>")
    f.close()

def write_result_page(new_Genome_ID, new_LIN_object, new_LIN, db_cursor,User_ID,Interest_ID):
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    db_cursor.execute("SELECT Attribute_IDs from Interest WHERE Interest_ID={0}".format(Interest_ID))
    Attribute_IDs = db_cursor.fetchone()[0]
    db_cursor.execute("select AttributeName from Attribute where Attribute_ID in ({0})".format(Attribute_IDs))
    tmp = db_cursor.fetchall()
    AttributeName = ['_'.join(i[0].split(' ')) for i in tmp]
    # AttributeName_string = ','.join(AttributeName)
    df = {}
    for i in AttributeName:
        df[i] = []
    keys = df.keys()
    AttributeName_string = ','.join(df.keys())
    db_cursor.execute("SELECT {0} FROM Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string,new_Genome_ID))
    tmp = db_cursor.fetchone()
    for i in range(len(df.keys())):
        df[df.keys()[i]].append(tmp[i])

    db_cursor.execute("SELECT LIN.SubjectGenome, LIN.ANI FROM LIN,Genome WHERE LIN.Genome_ID=Genome.Genome_ID "
                      "and Genome.Genome_ID='{0}'".format(new_Genome_ID))
    tmp = db_cursor.fetchone()
    best_hit = tmp[0]
    ANI_best_hit = str(float(tmp[1]) * 100) + '%'
    LIN_best_hit = new_LIN_object.LIN

    # Get info of the new submission
    db_cursor.execute("SELECT AttributeValue from AttributeValue WHERE Genome_ID={0} AND Attribute_ID "
                      "IN (1,2,3)".format(new_Genome_ID))
    tmp = db_cursor.fetchall()
    if len(tmp) == 0:
        Genus_new_Genome = Species_new_Genome = Strain_new_Genome = "N/A"
    else:
        Genus_new_Genome = tmp[0][0]
        Species_new_Genome = tmp[1][0]
        Strain_new_Genome = tmp[2][0]
    LIN_new_Genome = new_LIN # A string
    # Get info of the best match
    Genome_ID_best_hit = new_LIN_object.Genome_ID
    db_cursor.execute("SELECT AttributeValue from AttributeValue WHERE Genome_ID={0} AND Attribute_ID "
                      "IN (1,2,3)".format(Genome_ID_best_hit))
    if len(tmp) == 0:
        Genus_best_hit = Species_best_hit = Strain_best_hit = "N/A"
    else:
        Genus_best_hit = tmp[0][0]
        Species_best_hit = tmp[1][0]
        Strain_best_hit = tmp[2][0]
    db_cursor.execute("SELECT {0} from Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string,Genome_ID_best_hit))
    tmp = db_cursor.fetchone()
    for i in range(len(df.keys())):
        df[df.keys()[i]].append(tmp[i])
    # Get the Genome_IDs of all those sharing the same conserved part of LINs
    db_cursor.execute("SELECT Genome_ID, LIN FROM LIN WHERE LIN LIKE '{0}%' AND Genome_ID <> {1} and Genome_ID <> {2}".
                      format(",".join(new_LIN_object.conserved_LIN), new_Genome_ID, Genome_ID_best_hit))
    tmp = db_cursor.fetchall()
    Genome_IDs_related_hits = [int(i[0]) for i in tmp]
    LINs_related_hits = [i[1] for i in tmp]
    num_of_all_lines = len(Genome_IDs_related_hits)


    for each_related in Genome_IDs_related_hits:
        db_cursor.execute("SELECT {0} FROM Genome_to_Attribute WHERE Genome_ID={1}".format(AttributeName_string,each_related))
        tmp = db_cursor.fetchone()
        for i in range(len(df.keys())):
            df[df.keys()[i]].append(tmp[i])

    # A dataframe
    table = pd.DataFrame.from_dict(df)
    # table_withoutname = table.drop(["Genus","Species","Strain"],axis=1)
    keys_withoutname = [i for i in keys if i not in ["Genus","Species","Strain"]]
    LIN_all = [','.join(LIN_new_Genome)] + [','.join(LIN_best_hit)] + LINs_related_hits
    # Generate a random and unique result page file name
    unique_filename = str(uuid.uuid4())
    filename = '/var/www/html/linSite/ResultPage/' + unique_filename + '.php'
    url = 'http://128.173.74.68/linSite/ResultPage/' + unique_filename + '.php'
    # Read the static page header
    f = open('/var/www/html/linSite/ResultPage/static_header','r')
    static_header = f.read()
    f.close()
    # Read the static page footer
    f = open('/var/www/html/linSite/ResultPage/static_footer','r')
    static_footer = f.read()
    f.close()
    f = open(filename,'a')
    f.write(static_header)
    f.write("<div id='main-content'>")
    f.write("<form action='add_description.php' method='post' class='pure-form pure-form-aligned'>")
    f.write("<fieldset>")
    f.write("<legend>Results:</legend>")
    f.write('<table class="pure-table pure-table-horizontal" id="Table_main" role="grid" style="width: 100%">')
    f.write("<thead>")
    f.write("<tr id='header' role='row'>")
    f.write("<th rowspan='1' colspan='1' style='width:2%;'></th>")
    f.write("<th rowspan='1' colspan='1' style='width:6%;'>Category</th>")
    f.write("<th rowspan='1' colspan='1' style='width:4%;'>Genus</th>")
    f.write("<th rowspan='1' colspan='1' style='width:4%;'>Species</th>")
    f.write("<th rowspan='1' colspan='1' style='width:4%;'>Strain</th>")
    positions = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    for i in range(len(positions)):
        f.write("<th class='LIN' rowspan='1' colspan='1' style='width:2%;'>"
                "<input class='col_checkbox' type='checkbox' id='{0}' name='col_checkbox[{1}]' value='{2}'>"
                "{3}"
                "</th>".format(i,i,i,positions[i]))
    f.write("</tr>")
    f.write("</thead>")
    f.write("<tbody>")
    f.write("<tr id='0' role='row' class='record'>"
            "<td><input type='checkbox' id='row_checkbox[0]' name='row_checkbox[0]' value='0'></td>"
            "<td>New submission</td>"
            "<td>{0}</td><td>{1}</td><td>{2}</td>".format(Genus_new_Genome,Species_new_Genome,Strain_new_Genome))
    for i in LIN_new_Genome:
        f.write("<td class='LIN'>{0}</td>".format(i))
    f.write("<div id='detail_0' style='display:none; position:abosolute; background-color: #CAE1FF;' class='popup'>")
    f.write("<p><b>Genus: </b>{0}</p>".format(table['Genus'][0]))
    f.write("<p><b>Species: </b>{0}</p>".format(table['Species'][0]))
    f.write("<p><b>Strain: </b>{0}</p>".format(table['Strain'][0]))
    for i in keys_withoutname:
        f.write("<p><b>{0}: </b>{1}</p>".format(i, table[i][0]))
    f.write("</div>")
    f.write("</tr>")

    f.write("<tr id='0' role='row' class='record pure-table-odd'>"
            "<td><input type='checkbox' id='row_checkbox[0]' name='row_checkbox[0]' value='0'></td>"
            "<td>Best hit</td>"
            "<td>{0}</td><td>{1}</td><td>{2}</td>".format(Genus_best_hit, Species_best_hit, Strain_best_hit))
    for i in LIN_best_hit:
        f.write("<td class='LIN'>{0}</td>".format(i))
    f.write("<div id='detail_1' style='display:none; position:abosolute; background-color: #CAE1FF;' class='popup'>")
    f.write("<p><b>Genus: </b>{0}</p>".format(table['Genus'][1]))
    f.write("<p><b>Species: </b>{0}</p>".format(table['Species'][1]))
    f.write("<p><b>Strain: </b>{0}</p>".format(table['Strain'][1]))
    for i in keys_withoutname:
        f.write("<p><b>{0}: </b>{1}</p>".format(i, table[i][1]))
    f.write("</div>")
    f.write("</tr>")

    for i in range(num_of_all_lines):
        idx = i+2
        if idx % 2 == 0:
            f.write("<tr id='{0}' role='row' class='record'>".format(idx))
        else:
            f.write("<tr id='{0}' role='row' class='record pure-table-odd'>".format(idx))
        f.write("<td><input type='checkbox' id='row_checkbox[{0}]' name='row_checkbox[{1}]' value='{2}'></td>".format(idx,idx,idx))
        f.write("<td>Related hit</td>")
        f.write("<td>{0}</td><td>{1}</td><td>{2}</td>".format(table['Genus'][idx],table['Species'][idx],table['Strain'][idx]))
        lin = LINs_related_hits[i].split(',')
        for each_lin_position in lin:
            f.write("<td class='LIN'>{0}</td>".format(each_lin_position))
        f.write("<div id='detail_{0}' style='display:none; position:abosolute; background-color: #CAE1FF;' class='popup'>".format(idx))
        f.write("<p><b>Genus: </b>{0}</p>".format(table['Genus'][idx]))
        f.write("<p><b>Species: </b>{0}</p>".format(table['Species'][idx]))
        f.write("<p><b>Strain: </b>{0}</p>".format(table['Strain'][idx]))
        for i in keys_withoutname:
            f.write("<p><b>{0}: </b>{1}</p>".format(i, table[i][idx]))
        f.write("</div>")
        f.write("</tr>")
    f.write("</tbody>")
    f.write("</table>")
    f.write("</fieldset>")
    f.write("<br>")
    f.write("<fieldset>")
    f.write("<legend>Add description(s) to selected LINs</legend>")
    db_cursor.execute("select * from Description_Items")
    tmp = db_cursor.fetchall()
    Description = [i[1] for i in tmp]
    Description_Item_ID = [int(i[0]) for i in tmp]
    num_of_descriptions = len(Description)
    f.write("<a style='display: none;' id='attribute_counter'>{0}</a>".format(num_of_descriptions))
    f.write("<div id='counter' style='display: none;'>")
    f.write("<a id='count'>1</a>")
    f.write("</div>")
    f.write("<div id='default_option' class='pure-control'>")
    f.write("<select id='description_input[0]' name='description_input[0]'>")
    for i in range(num_of_descriptions):
        f.write("<option value='{0}'>{1}</option>".format(Description_Item_ID[i],Description[i]))
    f.write("</select>")
    f.write("<input id='description_input_text[0]' name='description_input_text[0]' type='text' value='' /></div>")
    for i in range(1,num_of_descriptions):
        idx = i+1
        f.write("<div id='MoreOption_{0}' style='display:none' class='pure-control-group'>".format(idx))
        f.write("<select id='description_input[{0}]' name='description_input[{0}]'>".format(i,i))
        f.write("<option value='0'>--Select--</option>")
        for j in range(num_of_descriptions):
            f.write("<option value='{0}'>{1}</option>".format(Description_Item_ID[j],Description[j]))
        f.write("</select>")
        f.write("<input id='description_input_text[{0}]' name='description_input_text[{1}]' type='text' value='' /></div>")
    f.write("<div class='pure-control-group'><button class='pure-button button-secondary' "
            "type='button' id='MoreOption_original' onclick='javascript:ExtendOption();'>Add an option</button></div>")
    f.write("<?php")
    f.write("session_start();")
    f.write("$_SESSION['LIN'] = {0};".format(LIN_all))
    f.write("$_SESSION['Description_Item_Name'] = {0};".format(Description))
    f.write("?>")
    f.write("</fieldset>")
    f.write("</form>")
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
    s = smtplib.SMTP('smtp.live.com',587)
    s.ehlo_or_helo_if_needed()
    s.starttls()
    s.login(me, '1988112019881120')
    s.sendmail(me, you, msg.as_string())
    s.quit()
