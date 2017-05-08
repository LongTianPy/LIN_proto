#!/usr/bin/python
"""This script functions to assign new LINs in a different scheme based on LIN.SubjectGenome and LIN.ANI
    Remember to keep the LIN table of the new database only first row, or not. You don't have to.
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
import sys
sys.path.append("/home/linproject/Projects/LIN_proto/")
import LIN_Assign

# FUNCTIONS
def connect_to_db(db):
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use {0}".format(db))
    return conn, c

def fetch_data(old):
    conn, c = connect_to_db(old)
    c.execute("select Genome_ID, SubjectGenome, ANI from LIN where Scheme_ID=4")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    SubjectGenome = [int(i[1]) for i in tmp]
    ANI = [i[2] for i in tmp]
    df = pd.DataFrame({"SubjectGenome":SubjectGenome,"ANI":ANI},index=Genome_ID)
    print(len(df.index))
    return df

def reassign(new,df):
    conn, c = connect_to_db(new)
    Genome_ID = df.index
    sql = "insert into LIN (Genome_ID, Scheme_ID, SubjectGenome, ANI, LIN) value ({0},4,{1},{2},'{3}')"
    for each_genome in Genome_ID:
        if each_genome == 1:
            c.execute(sql.format(1,1,1,'0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0'))
            conn.commit()
        else:
            subjectgenome = df.get_value(each_genome,"SubjectGenome")
            ani = df.get_value(each_genome,"ANI")
            new_getLIN_object = LIN_Assign.getLIN(Genome_ID=subjectgenome,Scheme_ID=4,similarity=ani,c=c)
            new_LIN =  LIN_Assign.Assign_LIN(getLIN_object=new_getLIN_object,c=c,current_genome=each_genome).new_LIN
            c.execute(sql.format(each_genome,subjectgenome,ani,new_LIN))
            conn.commit()

def reassign2(new,df):
    conn,c=connect_to_db(new)
    Genome_ID=df.index
    sql = "update LIN set LIN='{0}' where Scheme_ID=4 and Genome_ID={1}"
    for query in Genome_ID[1:]:
        subject = df.get_value(query,"SubjectGenome")
        ani = df.get_value(query,"ANI")
        new_getLIN_object=LIN_Assign.getLIN(Genome_ID=subject,Scheme_ID=4,similarity=ani,c=c)
        new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_getLIN_object,c=c,current_genome=query).new_LIN
        c.execute(sql.format(new_LIN,query))
        conn.commit()


# MAIN
if __name__ == '__main__':
    df = fetch_data("LINdb")
    reassign2("LINdb",df)
