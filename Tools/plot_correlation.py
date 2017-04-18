#!/usr/bin/python
"""This scripts takes the pairwise ANI and mash similarity and the outputs are scatter plots
"""

# IMPORT
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pandas.stats.api import ols
import sys

# FUNCTIONS
def concat_dfs(df_ani,df_mash):
    idxs = df_ani.index
    cols = df_ani.columns
    idx_count = 1
    df = pd.DataFrame(columns=["ANI","Mash similarity"])
    for idx in idxs:
        for col in cols:
            df.set_value(idx_count,"ANI",df_ani.get_value(idx,col))
            if pd.isnull(df_mash.get_value(idx,col)):
                df.set_value(idx_count,"Mash similarity",0)
            else:
                df.set_value(idx_count,"Mash similarity",df_mash.get_value(idx,col))
            idx_count += 1
    df = df[df["ANI"] != 1]
    df_95up = df[df["ANI"]>0.95]
    df_95up.index = [i for i in range(len(df_95up.index))]
    return df_95up

def make_single_plot(df):
    df_ols = ols(x=df["ANI"],y=df["Mash similarity"])
    beta1 = df_ols.beta.x
    beta0 = df_ols.beta.intercept
    plt.plot(df["ANI"], df["Mash similarity"], ".b")
    plt.xlabel("ANI")
    plt.ylabel("Mash similarity")
    plt.plot(df["ANI"], df["ANI"] * beta1 + beta0, "-r")

def make_plots(df_95up):
    scheme=[60,70,75,80,85,90,95,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999,99.9999]
    scheme_percentage = [str(i) for i in scheme if i>=95]
    scheme = [float(i)/100 for i in scheme if i>=95]
    plt.figure(figsize=(8.27,11.69))
    subplot_pos = 1
    plt.subplot(2,7,subplot_pos)
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
    for i in range(1,len(scheme)):
        subplot_pos += 1
        lower = scheme[i-1]
        upper = scheme[i]
        lower_percentage = scheme_percentage[i-1]
        upper_percentage = scheme_percentage[i]
        df_sub = df_95up[df_95up["ANI"]>=lower][df_95up["ANI"]<=upper]
        plt.subplot(2,7,subplot_pos)
        df_sub_ols = ols(x=df_sub["ANI"],y=df_sub["Mash similarity"])
        beta1 = df_sub_ols.beta.x
        beta0 = df_sub_ols.beta.intercept
        plt.plot(df_sub["ANI"],df_sub["Mash similarity"],".b")
        plt.xlabel("ANI")
        plt.ylabel("Mash similarity")
        plt.plot(df_sub["ANI"],df_sub["ANI"]*beta1+beta0,"-r")
        title = r"Correlation between ANI ({0}$\%\leqslant$ANI$\leqslant${1}$\%$) and Mash similarity, $R^2$={2}".format(lower_percentage,
                                                                                                     upper_percentage,
                                                                                                    r2)
        plt.title(title)
    plt.savefig("overall_and_scheme.pdf")

# MAIN
if __name__ == "__main__":
    ani = sys.argv[1]
    mash = sys.argv[2]
    df_ani = pd.read_table(ani,sep="\t",header=0,index_col=0)
    df_mash = pd.read_table(mash,sep=",",header=0,index_col=0)
    df_95up = concat_dfs(df_ani=df_ani,df_mash=df_mash)
    make_plots(df_95up)
