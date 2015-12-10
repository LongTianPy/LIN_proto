#!/usr/bin/python

__author__ = 'longtian'

# IMPORT
from itertools import product
import random
# CLASS

# FUNCTIONS
def KMerLib(k): # Prepare K-mer library based on k. Time cost increases exponentially with k.
    """We set k start from 10, aka 10-mer, so we get 4^10 kmers in the library.
        """
    bases = ['A', 'T', 'C', 'G']
    kmers = [''.join(i) for i in product(bases,repeat=k)]
    # Sample 500 of them randomly
    kmers_500 = random.sample(kmers,500)
    return kmers_500



