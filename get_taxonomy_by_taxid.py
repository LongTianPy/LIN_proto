#!/usr/bin/python
"""
"""

# IMPORT
from Bio import Entrez
from MySQLdb import Connect
import pandas as pd

Entrez.email = "aaa@bbb.ccc"

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq")
    return conn, c

def extract_taxonomy(tax_id):
    tax_id = str(tax_id)
    handler = Entrez.efetch(db="taxonomy", id=tax_id, rettype="xml")
    record = Entrez.read(handler)[0]
    lineage = record["LineageEx"]  # Returns a list of dictionaries
    full_lineage = [i['ScientificName'] for i in lineage]
    full_tax_ids = [i['TaxId'] for i in lineage]
    full_rank = [i['Rank'] for i in lineage]
    return full_lineage, full_tax_ids, full_rank




# MAIN
if __name__ == '__main__':
    conn,c = connect_to_db()
    c.execute("select Genome_ID,AttributeValue from AttributeValue WHERE Attribute_ID=15 and AttributeValue<>'N/A'")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    Tax_ID = [str(i[1]) for i in tmp]
    df = pd.DataFrame()
    df["Taxonomy_ID"] = Tax_ID
    df.index = Genome_ID
    full_lineages = []
    full_tax_ids = []
    full_rank = []
    for i in Tax_ID:

