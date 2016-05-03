#!/usr/bin/python

# This script consists of nearly all the queries we need to extract information from the database

# IMPORT


# FUNCTIONS
def get_top10_LIN(genome,cursor):
    c = cursor
    LIN_query = "SELECT LIN.LIN FROM Genome, LIN WHERE Genome.Genome_ID = LIN.Genome_ID and Genome.GenomeName = '{0}'".format(genome)
    c.execute(LIN_query)
    LIN = c.fetchone()[0]
    return LIN