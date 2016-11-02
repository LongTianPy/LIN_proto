#!/usr/bin/python

"""
Prototype of LINgroup indexing, code might (very, very likely to) be ugly since it's the first edition.
"""

# IMPORT
import pandas as pd
import sys


# FUNCTIONS
# def process_current_LIN_table(current_LIN_table, current_level):


# MAIN
if __name__ == "__main__":
    ani_table = sys.argv[1]
    LIN_table = sys.argv[2]
    # Read ANI table and sort it
    ani = pd.read_table(ani_table, header=0, index_col=0)
    ani_sorted = ani.sort().sort(axis=1)
    # Read LIN table, sorted already
    LIN = pd.read_table(LIN_table, header=0, index_col=0)
    LIN["listed_LIN"] = [i.split(".") for i in LIN["LIN"]]
    genomes = ani_sorted.columns
    reverse_LIN_dict = {LIN.get_value(each_genome, "LIN"): each_genome for each_genome in genomes}
    for i in range(1, len(genomes)):
        print "#####################################"
        print "Current query: " + genomes[i] + "\t" + LIN.get_value(genomes[i], "LIN")
        # print "The whole set of genomes compared to this current query: "
        LIN_table_piece = LIN.loc[genomes[:i], ["LIN", "listed_LIN"]]
        # print LIN_table_piece["LIN"]
        ani_table_piece = ani_sorted.loc[genomes[:i], genomes[i]]
        # print "The ANI table of the current query:"
        # print ani_table_piece
        print "The best hit is " + ani_table_piece.idxmax() + ", with ANI of " + str(ani_table_piece.max())
        print "The best hit's LIN is " + LIN.get_value(ani_table_piece.idxmax(),"LIN")
        for j in range(19):
            print "Current position: " + str(j + 1)
            LIN_dictionary = {}
            previous_genomes = ani_table_piece.index
            for each_previous_genome in previous_genomes:
                each_LIN = LIN_table_piece.get_value(each_previous_genome, "listed_LIN")
                each_leading_part = ".".join(each_LIN[:j + 1])  # Leading part of LIN includes the current position
                if each_leading_part not in LIN_dictionary:
                    LIN_dictionary[each_leading_part] = {each_LIN[j + 1]: [each_previous_genome]}
                    # The sub-dictionary of the
                    # dictionary consists of the keys corresponding to the next number assigned and the genomes
                    # having the same
                    # leading_part + next number
                else:
                    if each_LIN[j + 1] not in LIN_dictionary[each_leading_part]:
                        LIN_dictionary[each_leading_part][each_LIN[j + 1]] = [each_previous_genome]
                    else:
                        LIN_dictionary[each_leading_part][each_LIN[j + 1]].append(each_previous_genome)
            print "Detail of the current level LIN table:"
            for each_LIN_dictionary_key in LIN_dictionary.keys():
                print each_LIN_dictionary_key
                for each_next_number in LIN_dictionary[each_LIN_dictionary_key].keys():
                    for each_genome_in_dictionary in LIN_dictionary[each_LIN_dictionary_key][each_next_number]:
                        print ".".join("+" * len(each_LIN_dictionary_key.split(
                            "."))) + "." + each_next_number + "\t" + each_genome_in_dictionary
            # Now we need to pre-set some conditions to reduce the time used
            print "Start iterating the current LIN table"
            LIN_ANI_storage = {}
            LIN_ANI_max_storage = {}
            for each_LIN_dictionary_key in LIN_dictionary.keys():
                LIN_ANI_storage[each_LIN_dictionary_key] = []

                for each_next_number in LIN_dictionary[each_LIN_dictionary_key].keys():
                    # Assume the trailing part of each target genome is all 0, which is the safest choice
                    subject_LIN = each_LIN_dictionary_key + "." + each_next_number + \
                                  "".join([".0"]*(20 - 1 - j - 1))
                    subject_genome = reverse_LIN_dict[subject_LIN]
                    subject_ANI = ani_table_piece[subject_genome]
                    LIN_ANI_storage[each_LIN_dictionary_key].append(subject_ANI)
                LIN_ANI_max_storage[each_LIN_dictionary_key] = max(LIN_ANI_storage[each_LIN_dictionary_key])
            leading_part_w_max_ANI = max(LIN_ANI_max_storage,key=LIN_ANI_max_storage.get)
            leading_part_real_LIN = ".".join(LIN.get_value(ani_table_piece.idxmax(), "listed_LIN")[:j+1])
            leading_part_assigned_LIN = ".".join(LIN.get_value(genomes[i],"listed_LIN")[:j+1])
            print "For now, the selected leading part of LIN as the LINgroup index is    " + leading_part_w_max_ANI
            print "The actual current leading part of the referece LIN is                " + leading_part_real_LIN
            print "The leading part of the later assigned LIN is                         " + leading_part_assigned_LIN
            if leading_part_w_max_ANI == leading_part_real_LIN:
                print "True \n\n"
            else:
                if leading_part_w_max_ANI[:-2] == leading_part_real_LIN[:-2]:
                    print "Done \n\n"
                    break
                else:
                    print LIN_ANI_max_storage
                    print "False \n\n"
            """
            According to the new evaluation report, route shifting happens when the real best route consists of multiple members at
            the next step hence the average value got dragged down sometimes and less than the predicted route's average. Solution to
            this issue would be
            1.
                a. keep track of such fluctuation, if a new leading part, different from the previous one, appears, number of members
                    of the leading part inherited from the previous best route is calculated
                b. if this number is larger than 1, the tracked fluctuation might possibily caused by the reason aforementioned
                c.
            """
