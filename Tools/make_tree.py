#!/usr/bin/python
"""This script generate a phylogenetic tree from LINs, given a list of genome IDs
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
import string

# OBJECT
class tree_topo(object):
    def __init__(self,num,route):
        self.num = num
        self.route = route
        self.find_children()
    def find_chilren(self,route):


# FUNCTIONS
def connect_to_db(dbname):
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use {0}".format(dbname))
    return c

def make_dictionary(Genome_ID_list,c):
    dictionary = {}
    for each_genome in Genome_ID_list:
        c.execute("select LIN from LIN where Genome_ID={0} and Scheme_ID=4")
        dictionary[str(each_genome)] = c.fetchone()[0]
    return dictionary

def sort_genomes(dictionary):
    cutoff = [70,75,80,85,90,95,96,97,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999]
    cutoff = [i/100.00 for i in cutoff]
    genomes = dictionary.keys()
    rev_dict = {dictionary[key]:key for key in dictionary.keys()}
    lins = rev_dict.keys()
    tree = {}
    leaf_num = 0
    lin_length = len(lins[0])
    positions = list(string.ascii_uppercase[:lin_length])
    df = pd.DataFrame(columns=genomes,index=genomes)
    def compare_2_genomes(g1,g2,df,dictionary):
        lin1 = dictionary[g1]
        lin2 = dictionary[g2]
        i = 0
        while lin1[i] is lin2[i] and i<19:
            i += 1
        if i is 19:
            if lin1[i] is lin2[i]:
                similarity = 1
            else:
                similarity = cutoff[i]
        else:
            similarity = cutoff[i]
        df.loc[g1,g2] = similarity
        return df
    for i in genomes:
        for j in genomes:
            df = compare_2_genomes(i,j,df,dictionary)
    df.to_csv("similarity.csv")



# MAIN