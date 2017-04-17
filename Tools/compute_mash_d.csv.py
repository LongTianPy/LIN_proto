#!/usr/bin/python
"""For now sourmash only returns a distance matrix in binary file, not sure how to look into it...
This script computes the mash distance/similarities and return it as a data frame.
"""

# IMPORT
import pandas as pd
import os
from os import listdir
from os.path import isfile
import multiprocessing as mp
from functools import partial


# FUNCTIONS
def parse_result(current_sig):
    f = open("results_{0}.txt".format(current_sig),"r")
    lines = [i.strip().split(" \t ") for i in f.readlines()[3:]]
    f.close()
    mash_dict = {}
    for i in lines:
        sig = lines[2].split(".")[0]
        mash_dict[sig] = float(lines[1])
    return mash_dict

def run_sourmash(current_sig,total_number_of_sigs,df):
    cmd = "sourmash search {0} *.sig -n {1} > results_{0}.txt".format(current_sig,total_number_of_sigs)
    os.system(cmd)
    dict = parse_result(current_sig)
    current_col = current_sig.split(".")[0]
    for each_key in dict.keys():
        df.set_value(current_col,each_key,dict[each_key])


# MAIN
if __name__ == "__main__":
    sigs = [sig for sig in listdir("./") if isfile(sig) and sig.endswith("sig")]
    number_of_sigs = len(sigs)
    names = [sig.split(".")[0] for sig in sigs]
    df = pd.DataFrame(index=names,columns=names)
    pool = mp.Pool(8)
    partial_run_sourmash = partial(run_sourmash,total_number_of_sigs=number_of_sigs,df=df)
    pool.map(partial_run_sourmash,sigs)
    df.to_csv("pairwise_mash.csv")
