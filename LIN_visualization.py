#!/usr/bin/python
"""
Thinking about how to visualize LIN, since tree structure usually is more straightforward, let's first try tree.
Our sample file is Psy_211.csv, where there are columns of Genus, Species, Subspecies, LIN (of course) and Genome_ID.
"""

# IMPORT
import pandas as pd
import sys

# FUNCTIONS
def update_LIN_table(current_position, updated_index,df,f,positions,LINgroup):
    current_number_undup = list(set(df[positions[current_position]]))
    if len(current_number_undup) == 1:
        return current_position-1, updated_index,LINgroup
    else:
        index_to_list = list(updated_index)
        sub_df = df.loc[updated_index,positions[:current_position]]
        sub_df["sub_LIN"] = [".".join(sub_df.loc[[each_index], ].values[0]) for each_index in updated_index]
        for each_index in sub_df.index:
            if df.get_value(each_index, positions[current_position]) != "0":
                this_LIN_leading = ".".join(df.get_value(each_index, "LIN").split(",")[:current_position])
                if this_LIN_leading not in LINgroup:
                    LINgroup.append(this_LIN_leading)
                    for each_step in range(1,len(this_LIN_leading.split("."))):
                        LINgroup.append(".".join(this_LIN_leading.split(".")[:each_step]))
                else:
                    continue
                for each_of_each_index in sub_df.index:
                    if this_LIN_leading == sub_df.get_value(each_of_each_index,"sub_LIN"):
                        f.write(this_LIN_leading+"."+df.get_value(each_of_each_index,"Genus")+"_"+
                                df.get_value(each_of_each_index,"Species")+"_"+
                                df.get_value(each_of_each_index,"Pathovar")+"_"+
                                df.get_value(each_of_each_index,"Strain")+"\n")
                        if each_of_each_index in index_to_list:
                            index_to_list.remove(each_of_each_index)
                        else:
                            continue
                    else:
                        continue
            else:
                continue
        return current_position-1, index_to_list, LINgroup

def write_structured_csv(file):
    df = pd.read_csv(file)
    df.index = df["Genome_ID"]
    df_index = df.index
    positions = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T"]
    for i in range(len(positions)):
        df[positions[i]] = [j.split(",")[i] for j in df["LIN"]]
    f = open("first_tryout.csv","w")
    f.write("id\n")
    current_position = 19
    updated_index = df_index
    LINgroup = []
    while current_position >= 0:
        current_position, updated_index,LINgroup = update_LIN_table(current_position=current_position, updated_index=updated_index,
                                                        df=df,f=f,positions=positions,LINgroup=LINgroup)
    undup_LINgroup = list(set(LINgroup))
    for each_LINgroup in undup_LINgroup:
        f.write(each_LINgroup+"\n")
    f.close()

# MAIN
if __name__ == "__main__":
    csv_file = sys.argv[1]
    write_structured_csv(csv_file)
