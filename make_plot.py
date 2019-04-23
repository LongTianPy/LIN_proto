#!/usr/bin/python
"""This script generates plot based on pairwise Jaccard similarity of a given LINgroup
"""

# IMPORT
import os
import sys
from MySQLdb import Connect
from os.path import isdir,join
import shutil

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2.0/all_sketches/"
plot_dir = "/home/linproject/Workspace/Sourmash2.0/plot/"

# FUNCTIONS
def connect_to_db():
    conn = Connect("127.0.0.1", "LINbase","Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def get_genomes(LINgroup,c):
    if len(LINgroup.split(",")) == 20:
        LINgroup = LINgroup
    else:
        LINgroup = LINgroup + ","
    c.execute("SELECT Genome_ID FROM LIN WHERE LIN LIKE '{0}%'".format(LINgroup))
    tmp = c.fetchall()
    genome_ids = [str(i[0]) for i in tmp]
    return genome_ids

def compare(genome_ids,job_id):
    working_dir = join(plot_dir + job_id)
    if not isdir(working_dir):
        os.mkdir(working_dir)
    for i in genome_ids:
        shutil.copy(join(sourmash_dir,'{0}.sig'.format(i)), working_dir)
    os.system('sourmash compare -k 21 {0} -o {1} -q'.format(join(working_dir,'*.sig'), join(working_dir,job_id)))
    matrix = join(working_dir,job_id)
    labels = join(working_dir,job_id+'.labels.txt')
    return matrix,labels
def plot(matrix,labels,c,job_id):
    with open(labels,'r') as f:
        lines = [i.strip() for i in f.readlines()]
    with open(labels,'w') as f:
        for line in lines:
            c.execute('SELECT Genome.Genome_ID,LIN.LIN FROM Genome LEFT JOIN LIN ON Genome.Genome_ID=LIN.Genome_ID WHERE Genome.FilePath="{0}"'.format(line))
            [id, lin] = c.fetchone()
            c.execute('SELECT NCBI_Tax_ID.Taxon nt, Taxonomy.Taxon tt FROM Taxonomy LEFT JOIN NCBI_Tax_ID ON Taxonomy.NCBI_Tax_ID=NCBI_Tax_ID.NCBI_Tax_ID WHERE Taxonomy.Genome_ID={0} AND Taxonomy.Rank_ID=20'.format(id))
            [tt,nt] = c.fetchone()
            if tt is None:
                name = nt
            else:
                name = tt
            f.write('{:<10}    {:<45}\n'.format(name,lin))
    os.system("sourmash plot {0} --labels 2> /dev/null".format(matrix))
    return join(plot_dir, job_id, '{0}.matrix.png'.format(job_id))

def make_plot(LINgroup,job_id):
    conn, c = connect_to_db()
    genome_ids = get_genomes(LINgroup, c)
    matrix, labels = compare(genome_ids, job_id)
    figure = plot(matrix, labels, c, job_id)
    return figure

# MAIN
if __name__ == '__main__':
    LINgroup = sys.argv[1]
    job_id = sys.argv[2]
    conn, c = connect_to_db()
    genome_ids = get_genomes(LINgroup,c)
    matrix,labels = compare(genome_ids,job_id)
    figure = plot(matrix,labels,c,job_id)