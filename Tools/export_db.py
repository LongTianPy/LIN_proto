#!/usr/bin/python
"""This script export bacteria info including names and LINs and LIN-associated
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
import sys

# FUNCTIONS
def connect_to_db(db):
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use {0}".format(db))
    return c

def export_table(c):
    c.execute("select Genome_ID,SubjectGenome,ANI,LIN from LIN WHERE Scheme_ID=4")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    SubjectGenome = [int(i[1]) for i in tmp]
    ANI = [i[2] for i in tmp]
    LIN = [i[3] for i in tmp]
    c.execute("select AttributeValue from AttributeValue where Attribute_ID=1 AND Genome_ID in ({0})".format(",".join([str(i) for i in Genome_ID])))
    tmp = c.fetchall()
    Genus = [i[0] for i in tmp]
    c.execute("select AttributeValue from AttributeValue where Attribute_ID=2 AND Genome_ID in ({0})".format(",".join([str(i) for i in Genome_ID])))
    tmp = c.fetchall()
    Species = [i[0] for i in tmp]
    # c.execute("select AttributeValue from AttributeValue where Attribute_ID=3")
    # tmp = c.fetchall()
    # Subspecies = [i[0] for i in tmp]
    c.execute("select AttributeValue from AttributeValue where Attribute_ID=4 AND Genome_ID in ({0})".format(",".join([str(i) for i in Genome_ID])))
    tmp = c.fetchall()
    Strain = [i[0] for i in tmp]
    df = pd.DataFrame({"SubjectGenome":SubjectGenome,"ANI":ANI,"LIN":LIN},index=Genome_ID)
    df["Genus"] = Genus
    df["Species"] = Species
    # df["Subspecies"] = Subspecies
    df["Strain"] = Strain
    return df

def main(db):
    c = connect_to_db(db)
    df = export_table(c)
    df.to_csv("{0}.csv".format(db))

# MAIN
if __name__ == "__main__":
    db = sys.argv[1]
    main(db)