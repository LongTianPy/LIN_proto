#!/usr/bin/python
"""
"""

# IMPORT
import os
from os import listdir
from os.path import isfile,isdir,join
import sys
import pandas as pd
# import matplotlib.pyplot as plt
# from pandas.stats.api import ols
# import matplotlib as mpl
# mpl.rc('text', usetex=True)
# mpl.rcParams['text.latex.preamble']=[r"\usepackage{amssymb}", r"\usepackage{amsmath}"]

# FUNCTIONS
#working_dir = "/home/linproject/Workspace_playaround/bbmap_sketches"
#genome_dir = "/home/linproject/Workspace/Psy_166"

def create_sketches(working_dir, genome_dir):
    cmd = "sh sketch.sh in={0} out={1}"
    genomes = [join(genome_dir,f) for f in listdir(genome_dir) if isfile(join(genome_dir,f)) and f.endswith(".fasta")]
    genomes_nopath = [f for f in listdir(genome_dir) if isfile(join(genome_dir,f))]
    for i in range(len(genomes)):
        os.system(cmd.format(genomes[i], join(working_dir,genomes_nopath[i]+".sketch")))

def run_comparesketch(working_dir,result_dir):
    files = [file for file in listdir(working_dir) if isfile(join(working_dir, file)) and file.endswith(".fasta.sketch")]
    cmd = "sh comparesketch.sh in={0} *.sketch minid=0 records=4000 > {1}"
    for i in files:
        os.system(cmd.format(i,result_dir+"/"+i.split(".")[0]+".txt"))

def parse_bbmap_result(file):
    f = open(file,"r")
    lines = [i.strip().split("\t") for i in f.readlines()[3:-1]]
    f.close()
    wkid = []
    kid = []
    est_ani = []
    subject = []
    for i in lines:
        subject.append(i[-1].split(" ")[0].split(".")[0])
        wkid.append(float(i[0][:-1])/100)
        kid.append(float(i[1][:-1])/100)
        est_ani.append(float(i[2][:-1])/100)
        df = pd.DataFrame({"wkid":wkid,"kid":kid,"est_ani":est_ani},index=[int(i) for i in subject])
    return df

def fill_df(result_dir):
    os.chdir(result_dir)
    files = [file for file in listdir(result_dir) if
             isfile(join(result_dir, file)) and file.endswith(".txt")]
    prefix = [i.split(".")[0] for i in files]
    wkid_df = pd.DataFrame(0,index=[int(i) for i in prefix],columns=prefix)
    kid_df = pd.DataFrame(0,index=[int(i) for i in prefix],columns=prefix)
    est_ANI_df = pd.DataFrame(0,index=[int(i) for i in prefix],columns=prefix)
    for file in files:
        this = int(file.split(".")[0])
        df = parse_bbmap_result(file=file)
        ## print(df)
        idx = df.index
        for i in idx:
            wkid_df.loc[this,str(i)] = df.get_value(i,"wkid")
            kid_df.loc[this,str(i)] =  df.get_value(i,"kid")
            est_ANI_df.loc[this,str(i)] = df.get_value(i,"est_ani")
    wkid_df.to_csv("../wkid.csv")
    kid_df.to_csv("../kid.csv")
    est_ANI_df.to_csv("../est_ANI.csv")

def concat_dfs(df_dir):
    os.chdir(df_dir)
    est_ANI = pd.read_csv("est_ANI.csv", header=0, index_col=0)
    wkid = pd.read_csv("wkid.csv", header=0, index_col=0)
    kid = pd.read_csv("kid.csv", header=0, index_col=0)
    idxs = est_ANI.index
    cols = est_ANI.columns
    total = len(idxs) * len(cols)
    df = pd.DataFrame(columns=['wkid', 'kid', 'est_ANI'],index=[i+1 for i in range(total)])
    job_pairs = []
    for idx in idxs:
        for col in cols:
            job_pairs.append(str(idx) + "_" + str(col))
    idx_count = 1



    for idx in idxs:
        for col in cols:
            df.loc[idx_count, "wkid"] = wkid.get_value(idx, col)
            df.loc[idx_count, "kid"] = kid.get_value(idx, col)
            df.loc[idx_count, "est_ANI"] = est_ANI.get_value(idx, col)
            idx_count += 1
    df.to_csv("concat_df.csv")
    # df_ols = ols(x=df["est_ANI"], y=df["wkid"])
    # beta1 = df_ols.beta.x
    # beta0 = df_ols.beta.intercept
    # plt.plot(df["est_ANI"], df["wkid"], ".b")
    # plt.xlabel("Estimated ANI")
    # plt.ylabel("WKID(Weighted k-mer identity)")
    # plt.plot(df["est_ANI"], df["est_ANI"] * beta1 + beta0, "-r")
    # plt.title("r'Correlation between estimated ANI and weighted k-mer identity'")
    # plt.savefig("ANI_vs_wkid.pdf")

# MAIN
if __name__ == '__main__':
    working_dir = sys.argv[1]
    genome_dir = sys.argv[2]
    result_dir = sys.argv[3]
    df_dir = sys.argv[4]
    # create_sketches(working_dir=working_dir, genome_dir=genome_dir)
    # os.chdir(working_dir)
    # run_comparesketch(working_dir=working_dir,result_dir=result_dir)
    # fill_df(result_dir=result_dir)
    concat_dfs(df_dir=df_dir)
