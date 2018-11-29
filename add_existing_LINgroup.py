#!/usr/bin/python
"""
"""

# IMPORT
from MySQLdb import Connect

# VARIABLES


# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost", "LINbase","Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def add_LINgroup(c,conn):
    c.execute("SELECT Genome_ID,LIN FROM LIN WHERE Scheme_ID=4")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    LIN = [i[1] for i in tmp]
    c.execute("SELECT LINgroup_ID, LINgroup FROM LINgroup ORDER BY LINgroup ASC")
    tmp = c.fetchall()
    LINgroup_ID = [int(i[0]) for i in tmp]
    LINgroup = [i[1] for i in tmp]
    for i in range(len(Genome_ID)):
        for j in range(len(LINgroup_ID)):
            if LIN[i].startswith(LINgroup[j]):
                sql = "UPDATE Genome SET LINgroup=CONCAT(LINgroup,',{0}') WHERE Genome_ID={1}"
                c.execute(sql.format(str(LINgroup_ID[j]), Genome_ID[i]))
                conn.commit()

# MAIN
