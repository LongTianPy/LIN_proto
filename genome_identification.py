#!/usr/bin/python
""" Workflow for genome identification of LINbase
    No metadata needed, and only FastANI is calculated.
"""

# IMPORT
import os
import sys
from os.path import join
from MySQLdb import Connect

# VARIABLES
working_dir = "/home/linproject/Workspace/Genome_identification"
genome_dir = "/home/linproject/Workspace/Genome_identification/uploaded_genome"
reference_list_file = join(working_dir,"reference_list.txt")

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def make_reference_list(c):
    c.execute("select Genome.FilePath from Genome,LIN where Genome.Genome_ID=LIN.Genome_ID and LIN.Scheme_ID=4")
    tmp = c.fetchall()
    with open(reference_list_file,"w") as f:
        for i in tmp:
            f.write("{0}\n".format(i[0]))

# MAIN
if __name__ == '__main__':
    conn, c = connect_to_db()
    make_reference_list(c)
    input_genome = join(genome_dir,sys.argv[1])
    output_file = sys.argv[2]
    fastANI_cmd = "fastANI -q {0} -rl {1} -o {2}/{3}".format(input_genome,reference_list_file,working_dir,output_file)
    os.system(fastANI_cmd)

