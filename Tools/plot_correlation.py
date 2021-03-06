#!/usr/bin/python
"""This scripts takes the pairwise ANI and mash similarity and the outputs are scatter plots
"""

# IMPORT
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pandas.stats.api import ols
import sys
import matplotlib as mpl

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amssymb}",
                                     r"\usepackage{amsmath}"]

# FUNCTIONS
def concat_dfs(df_ani,df_mash):
    idxs = df_ani.index
    cols = df_ani.columns
    idx_count = 1
    df = pd.DataFrame(columns=["ANI","Mash similarity","Query Genome","Subject Genome"])
    for idx in idxs:
        for col in cols:
            df.set_value(idx_count,"ANI",df_ani.get_value(idx,col))
            df.set_value(idx_count,"Query Genome",idx)
            df.set_value(idx_count,"Subject Genome",col)
            if pd.isnull(df_mash.get_value(idx,col)):
                df.set_value(idx_count,"Mash similarity",0)
            else:
                df.set_value(idx_count,"Mash similarity",df_mash.get_value(idx,col))
            idx_count += 1
    df = df[df["ANI"] != 1]
    df_95up = df[df["ANI"]>0.95]
    df_95up.index = [i for i in range(len(df_95up.index))]
    df_zero = df_95up[df_95up["Mash similarity"]==0]
    df_zero.to_csv("mash_zero.csv")
    return df, df_95up

def make_single_plot(df):
    df_ols = ols(x=df["ANI"],y=df["Mash similarity"])
    beta1 = df_ols.beta.x
    beta0 = df_ols.beta.intercept
    plt.plot(df["ANI"], df["Mash similarity"], ".b")
    plt.xlabel("ANI")
    plt.ylabel("Mash similarity")
    plt.plot(df["ANI"], df["ANI"] * beta1 + beta0, "-r")

def make_plots(df, df_95up):
    scheme=[60,70,75,80,85,90,95,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999,99.9999]
    scheme_percentage = [str(i) for i in scheme if i>=95]
    scheme = [float(i)/100 for i in scheme if i>=95]
    plt.figure(figsize=(20,40))
    plt.suptitle(r"Correlation between ANI and Mash similarity, with Mash params k=31, sketch size=500",
                 fontsize="36", fontweight="bold")
    # plt.figure()
    plt.subplot(7,2,1)
    plt.plot(df["ANI"],df["Mash similarity"],".b")
    plt.xlabel("ANI")
    plt.ylabel("Mash similarity")
    title = r"Correlation between ANI and Mash similarity"
    plt.title(title)
    subplot_pos = 2
    plt.subplot(7,2,subplot_pos)
    df_95up_ols = ols(x=df_95up["ANI"],y=df_95up["Mash similarity"])
    beta1 = df_95up_ols.beta.x
    beta0 = df_95up_ols.beta.intercept
    r2 = "{0:.3f}".format(df_95up_ols.r2)
    plt.plot(df_95up["ANI"],df_95up["Mash similarity"],".b")
    plt.xlabel("ANI")
    plt.ylabel("Mash similarity")
    plt.plot(df_95up["ANI"],df_95up["ANI"]*beta1+beta0,"-r")
    title = r"Correlation between ANI ($\geqslant95\%$) and Mash similarity, $R^2$={0}".format(r2)
    plt.title(title)
    # annotation = r"y={0}$\times$x+{1}".format(beta1,beta0)
    # plt.annotate(annotation,xy=(4,4))
    # plt.savefig("95up.pdf")
    for i in range(1,len(scheme)):
        lower = scheme[i-1]
        upper = scheme[i]
        lower_percentage = scheme_percentage[i-1]
        upper_percentage = scheme_percentage[i]
        df_sub = df_95up[df_95up["ANI"]>=lower][df_95up["ANI"]<=upper]
        if len(df_sub.index)>1:
            subplot_pos += 1
            plt.subplot(7, 2, subplot_pos)
            # plt.figure()
            df_sub_ols = ols(x=df_sub["ANI"],y=df_sub["Mash similarity"])
            beta1 = df_sub_ols.beta.x
            beta0 = df_sub_ols.beta.intercept
            r2 = df_sub_ols.r2
            plt.plot(df_sub["ANI"],df_sub["Mash similarity"],".b")
            plt.xlabel("ANI")
            plt.ylabel("Mash similarity")
            plt.plot(df_sub["ANI"],df_sub["ANI"]*beta1+beta0,"-r")
            title = r"Correlation between ANI ({0}$\%\leqslant$ANI$\leqslant${1}$\%$) and Mash similarity, $R^2$={2:.3f}".format(lower_percentage,
                                                                                                         upper_percentage,
                                                                                                        r2)
            plt.title(title)
            # annotation = r"y={0}$\times$x+{1}".format(beta1, beta0)
            # plt.annotate(annotation,xy=(4,4))
            # plt.savefig("{0}%_{1}%.pdf".format(lower_percentage,upper_percentage))
    plt.savefig("all_n_scheme.pdf")

# MAIN
if __name__ == "__main__":
    ani = sys.argv[1]
    mash = sys.argv[2]
    df_ani = pd.read_table(ani,sep="\t",header=0,index_col=0)
    df_mash = pd.read_table(mash,sep=",",header=0,index_col=0)
    df, df_95up = concat_dfs(df_ani=df_ani,df_mash=df_mash)
    make_plots(df, df_95up)
