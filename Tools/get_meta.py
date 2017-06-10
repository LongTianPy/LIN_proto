#!/usr/bin/python
"""
This script extracts metadata from database and generate an excel file
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
from export_db import connect_to_db
import sys

# FUNCTIONS
def extract_meta(value,output_file):
    c = connect_to_db("LINdb")
    c.execute("SELECT Attribute_IDs FROM Interest WHERE InterestName='Plant pathogens'")
    tmp = c.fetchall()
    Attribute_IDs=tmp[0][0]
    Attribute_IDs_list = Attribute_IDs.split(",")
    c.execute("SELECT AttributeName FROM Attribute WHERE Attribute_ID IN ({0})".format(Attribute_IDs))
    tmp=c.fetchall()
    AttributeNames = [i[0] for i in tmp]
    Attribute_dict = {str(Attribute_IDs_list[i]):AttributeNames[i] for i in range(len(Attribute_IDs_list))}
    # Let's say we want to pull result of all Pseudomonas
    c.execute("SELECT Genome_ID FROM AttributeValue WHERE AttributeValue='{0}'".format(value))
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    df = pd.DataFrame(index=Genome_ID)
    for attributename in AttributeNames:
        df[attributename] = ['']*len(Genome_ID)
    df["LIN"] = ['']*len(Genome_ID)
    for each_genome in Genome_ID:
        for each_attribute in Attribute_IDs_list:
            c.execute("SELECT AttributeValue FROM AttributeValue where Attribute_ID={0} AND Genome_ID={1}".format(int(each_attribute),each_genome))
            tmp = c.fetchone()
            attributevalue = tmp[0]
            df.set_value(each_genome,Attribute_dict[str(each_attribute)],attributevalue)
        c.execute("select LIN from LIN where Genome_ID={0} and Scheme_ID=4".format(each_genome))
        tmp = c.fetchone()
        lin = tmp[0]
        df.set_value(each_genome,"LIN",lin)
    writer = pd.ExcelWriter(output_file)
    df.to_excel(excel_writer=writer)


# MAIN
if __name__ == '__main__':
    value = sys.argv[1]
    output_file = sys.argv[2]
    extract_meta(value,output_file)