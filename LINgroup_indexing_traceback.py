#!/usr/bin/python
"""
This method traces back to the LINgroup indexing process to see how it reduced
the times of calculation
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd


# FUNCTIONS
def LINgroup_indexing_traceback():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use LINdb_Psy")
    c.execute("SELECT Genome_ID, SubjectGenome, LIN FROM LIN")
    tmp = c.fetchall()
    Genome_ID = [i[0] for i in tmp]
    SubjectGenome = [i[1] for i in tmp]
    LIN = [i[2] for i in tmp]
    df = pd.DataFrame()
    
# MAIN