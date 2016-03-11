#!/usr/bin/python

"""
This script aims to assign LIN to the new genome by reading and parsing the LINs of its top similar
genome. See the new LIN assignment strategy for more details.
"""

# IMPORT
import MySQLdb
from MySQLdb import Connect


# OBJECT
class getLIN(object):
    """
    Read the LIN of the top similar genome.
    """
    def __init__(self, genome, Scheme_ID, similarity):
        self.genome = genome
        self.Scheme_ID = Scheme_ID
        db = Connect('localhost', 'root')
        c = db.cursor()
        c.execute('use LINdb_test_2')
        c.execute("SELECT LabelNum from Scheme WHERE Scheme_ID=3")
        self.label_num = int(c.fetchone()[0])
        self.similarity = float(similarity)*100
        self.parse()
    def parse(self, genome = None, Scheme_ID = None, similarity = None):
        if not genome:
            genome = self.genome
        if not Scheme_ID:
            Scheme_ID = self.Scheme_ID
        if not similarity:
            similarity = self.similarity
        # Read the LIN of this genome
        db = Connect('localhost', 'root')
        c = db.cursor()
        c.execute('use LINdb_test_2')
        c.execute('SELECT LIN.LIN from LIN, Genome where LIN.Genome_ID=Genome.Genome_ID and Genome.FilePath LIKE "%{0}%" and LIN.Scheme_ID=3'.format(genome))
        lin = c.fetchone()[0].split(',')
        self.LIN = lin
        # Read the cutoff of this scheme
        c.execute('SELECT Cutoff from Scheme where Scheme_ID={0}'.format(Scheme_ID))
        cutoff = c.fetchone()[0].split(',')
        cutoff = [float(i) for i in cutoff]
        idx_to_change = 0
        if cutoff[0] > similarity:
            idx_to_change = 0
            self.idx_to_change = idx_to_change
            self.conserved_LIN = ''
        else:
            i = 0
            while i < len(cutoff)-1:
                if cutoff[i] < similarity and cutoff[i+1] >= similarity:
                    idx_to_change = i
                    break
                else:
                    idx_to_change = len(cutoff) - 1
                    i += 1
            self.idx_to_change = idx_to_change
            self.conserved_LIN = lin[:idx_to_change]

class Assign_LIN(object):
    """ Get the biggest number assigned to the idx_to_change with the same conserved part of LIN
    """
    def __init__(self, getLIN_object):
        self.idx_to_change = getLIN_object.idx_to_change
        self.conserved_LIN = ','.join(getLIN_object.conserved_LIN)
        self.label_num = getLIN_object.label_num
        self.assign()
    def assign(self, idx_to_change=None, conserved_LIN=None, label_num=None):
        if not idx_to_change:
            idx_to_change = self.idx_to_change
        if not conserved_LIN:
            conserved_LIN = self.conserved_LIN
        if not label_num:
            label_num = self.label_num
        db = Connect('localhost', 'root')
        c = db.cursor()
        c.execute('use LINdb_test_2')
        if conserved_LIN == '':
            c.execute("SELECT LIN.LIN FROM LIN")
            tmp = c.fetchall()
        else:
            c.execute('SELECT LIN.LIN from LIN WHERE LIN.LIN LIKE "{0}%"'.format(conserved_LIN))
            tmp = c.fetchall()
        LINs = [int(i[0].split(',')[idx_to_change]) for i in tmp]
        num_to_assign = str(max(LINs)+1)
        if idx_to_change != label_num - 1:
            tail = ['0'] * (label_num - 1 - idx_to_change)
            tail = ','.join(tail)
            new_LIN = conserved_LIN + ',%s,'%num_to_assign + tail
        else:
            new_LIN = conserved_LIN + ',%s'%num_to_assign
        if new_LIN.startswith(','):
            new_LIN = new_LIN[1:]
        else:
            new_LIN = new_LIN
        self.new_LIN = new_LIN









