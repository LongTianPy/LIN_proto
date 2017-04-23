#!/usr/bin/python
"""A bunch of functions related to mash or sourmash
"""

# IMPORT
import os
from os import listdir
from Bio import SeqIO
from os.path import join, isfile
from MySQLdb import Connect
from LINgroup_indexing_comparison import fetch_current, fetch_genomes
from os.path import isdir
import shutil
import pandas as pd
import LIN_Assign
import multiprocessing as mp
from functools import partial

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use LINdb_Psy")
    return conn, c

def write_both_strand(Genome_ID,cursor,sourmash_dir):
    cursor.execute("select FilePath from Genome where Genome_ID={0}".format(Genome_ID))
    tmp = cursor.fetchone()
    filepath = tmp[0]
    f = open(filepath,"r")
    records = list(SeqIO.parse(f,"fasta"))
    f.close()
    f = open("{0}{1}.fasta".format(sourmash_dir,Genome_ID),"w")
    for record in records:
        f.write(">" + str(record.id) + "\n")
        f.write(str(record.seq) + "\n")
        f.write(">" + str(record.id) + "_rc\n")
        f.write(str(record.seq.reverse_complement()) + "\n")
    f.close()

def create_signature(Genome_ID,sourmash_dir,cursor,conn):
    sig_path = sourmash_dir+str(Genome_ID)+".sig"
    sourmash_cmd = "sourmash compute -o {0} {1}.fasta -k 10,31 -n 1000".format(sig_path,sourmash_dir+str(Genome_ID))
    os.system(sourmash_cmd)
    # cursor.execute("insert into Signature (Genome_ID,SignaturePath) values ({0}, '{1}')".format(Genome_ID, sig_path))
    # conn.commit()
    return sig_path

# def sourmash_indexing(sourmash_dir, LINgroup):
#     target_folder = sourmash_dir+LINgroup+"/"
#     cmd = "sourmash sbt_index {0}index {0}*.sig".format(target_folder)
#     os.system(cmd)

def sourmash_searching(sourmash_dir,LINgroup,current_sig_path,current_genome):
    target_folder = sourmash_dir + LINgroup + "/"
    # shutil.copy(current_sig_path,target_folder)
    # copied_sig = target_folder + str(current_genome) + ".sig"
    files = [i for i in listdir(target_folder) if isfile(join(target_folder,i))]
    size = len(files)
    print size
    cmd = "sourmash search {0} {1}*.sig -n {2} > {1}{3}result.txt".format(current_sig_path,target_folder,size,current_genome)
    # print cmd
    os.system(cmd)
    f = open("{0}result.txt".format(target_folder,current_genome),"r")
    lines = [i.strip().split(" \t ") for i in f.readlines()[3:]]
    f.close()
    # candidates = []
    # candidates.append(lines[0])
    # if len(lines) > 1:
    #     for i in range(1,len(lines)):
    #         if float(lines[i][1]) - float(lines[i-1][1]) < -0.1:
    #             break
    #         else:
    #             candidates.append(lines[i])
    df = pd.DataFrame()
    if len(lines) == 0:
        return df
    else:
        mash_d = [float(i[1]) for i in lines]
        df["mash_d"] = mash_d
        df.index = [i[0].split("/")[-1].split(".")[0] for i in candidates]
        # df = df[df["mash_d"] > (df.get_value(df.index[0], "mash_d") - 0.1)]
        return df

def write_LIN_to_db(current_genome,subject_genome,ani,new_LIN,conn,c):
    sql = "INSERT INTO LIN (Genome_ID, Scheme_ID, LIN, SubjectGenome, ANI) values ({0}, 3, '{1}', '{2}', {3})"
    c.execute(sql.format(current_genome,3,subject_genome,ani,new_LIN))
    conn.commit()
# def check_result():

def test_mash():
    sourmash_dir = "/home/linproject/Workspace/Sourmash/"
    output_handler = open(sourmash_dir+"test_result.txt","w")
    output_handler.write("Genome_ID\tAssigned_LIN\tNew_G_LINgroup\tG_LINgroup\tNew_LIN\tANI_calcs\n")
    line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
    conn, c = connect_to_db()
    full_df, All_genomes = fetch_genomes(c)
    full_df["G_LINgroup"] = [",".join(each.split(",")[:7]) for each in full_df["LIN"]]
    write_both_strand(All_genomes[0], c, sourmash_dir)
    sig_path0 = create_signature(All_genomes[0],sourmash_dir,c,conn)
    if not isdir(sourmash_dir+"0,0,0,0,0,0,0/"):
        os.mkdir(sourmash_dir+"0,0,0,0,0,0,0/")
    shutil.copy(sig_path0,sourmash_dir+"0,0,0,0,0,0,0/")
    if not isdir(sourmash_dir+"rep_bac/"):
        os.mkdir(sourmash_dir+"rep_bac/")
    shutil.copy(sig_path0, sourmash_dir + "rep_bac/")
    output_handler.write(line.format(All_genomes[0],",".join(['0']*20),"Y","0,0,0,0,0,0,0",",".join(['0']*20),0))
    # write_LIN_to_db(1,1,1,",".join(['0']*20),conn,c)
    for idx in range(1,len(All_genomes)):
        current_genome, current_df = fetch_current(full_df,All_genomes,idx)
        write_both_strand(current_genome,c,sourmash_dir)
        sig_path_current = create_signature(current_genome,sourmash_dir,c,conn)
        current_LIN = full_df.get_value(All_genomes[idx],"LIN") # A string
        current_G_LINgroup = ",".join(current_LIN.split(",")[:7])
        G_LINgroup = [",".join(each.split(",")[:7]) for each in current_df["LIN"]]
        # subject_genome = full_df.get_value(All_genomes[idx],"SubjectGenome")
        # ani = full_df.get_value(All_genomes[idx],"ANI")
        # Trick is, mash d may or may not respond to 90% ani.
        # We are using 95% ani.
        # First thing is to check if it correlates with representative genomes.
        # If yes, there is change that the ani is between 90% and 95%. Need to check its LIN to determine if it's
        # in any LINgroup or creates a new one.
        if current_G_LINgroup not in G_LINgroup:
            if not isdir(sourmash_dir+current_G_LINgroup+"/"):
                os.mkdir(sourmash_dir+current_G_LINgroup+"/")
            shutil.copy(sig_path_current,sourmash_dir+current_G_LINgroup+"/")
            shutil.copy(sig_path_current,sourmash_dir+"rep_bac/")
            # sourmash_indexing(sourmash_dir,"rep_bac")
            output_handler.write(
            line.format(current_genome, current_LIN, "Y", current_G_LINgroup, current_LIN, 1))
            # db = Connect("localhost","root")
            # result_rep = sourmash_searching(sourmash_dir,"rep_bac",sig_path_current,current_genome)
            # if len(result_rep) == 2:
            #     output_handler.write(line.format(current_genome,current_LIN,"Y",current_G_LINgroup,"rep_bac","NA","NA",
            #                                      subject_genome,ani))
            # else:
            #     for i in range(1,len(result_rep)):
            #         if result_rep[i][1].split("/")[-1].split(".")[0] != str(current_genome):
            #             top_rep_mash = result_rep[i][1].split("/")[-1].split(".")[0]
            #             top_rep_mash_d = result_rep[i][0]
            #             break
            #     output_handler.write(line.format(current_genome,current_LIN,"Y",current_G_LINgroup,"rep_bac",top_rep_mash,
            #                                      top_rep_mash_d,
            #                                      subject_genome,ani))

        else:
            result_rep = sourmash_searching(sourmash_dir,"rep_bac",sig_path_current,current_genome)
            top_rep_mash = result_rep.index[0]
            top_rep_mash_d = result_rep.get_value(result_rep.index[0],"mash_d")
            top_rep_mash_G_LINgroup = full_df.get_value(int(top_rep_mash),"G_LINgroup")
            result = sourmash_searching(sourmash_dir,top_rep_mash_G_LINgroup,sig_path_current,current_genome)
            shutil.copy(sig_path_current, sourmash_dir + top_rep_mash_G_LINgroup + "/")
            mash_candidates = result.index
            mash_d = result["mash_d"]
            similarity = {}
            for each_candidate in mash_candidates:
                ani = assign_LIN_based_on_mash(current_genome=current_genome,subject_genome=each_candidate,c=c)
                similarity[each_candidate] = ani
            top_hit = max(similarity,key=similarity.get)
            top_ani = similarity[top_hit]
            new_LIN_obj = LIN_Assign.getLIN(Genome_ID=top_hit,Scheme_ID=3,similarity=top_ani,c=c)
            new_LIN = LIN_Assign.Assign_LIN(new_LIN_obj,c=c,current_genome=current_genome).new_LIN
            output_handler.write(line.format(current_genome,current_LIN,"N",current_G_LINgroup,
                                             new_LIN,len(similarity.keys())))
            print(str(current_genome)+"\n"+current_LIN+"\n"+new_LIN+"\n")
    output_handler.close()

def assign_LIN_based_on_mash(current_genome,subject_genome,c):
    ANI_calc_dir = "/home/linproject/Workspace/ANI/"
    if not isdir(ANI_calc_dir):
        os.mkdir(ANI_calc_dir)
    c.execute("select FilePath from Genome where Genome_ID={0}".format(int(current_genome)))
    tmp = c.fetchone()
    filepath_current_genome = tmp[0]
    c.execute("select FilePath from Genome where Genome_ID={0}".format(int(subject_genome)))
    tmp = c.fetchone()
    filepath_subject_genome = tmp[0]
    shutil.copyfile(filepath_current_genome,ANI_calc_dir+"{0}.fasta".format(str(current_genome)))
    shutil.copyfile(filepath_subject_genome,ANI_calc_dir+"{0}.fasta".format(str(subject_genome)))
    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py -i {0} " \
                "-o {0}output -m ANIblastall --nocompress".format(ANI_calc_dir)
    os.system(pyani_cmd)
    ANIb_result = pd.read_table(ANI_calc_dir + "output/ANIblastall_percentage_identity.tab",
                                header=0, index_col=0).get_value(int(current_genome), str(subject_genome))
    # print ANIb_result
    os.system("rm -rf {0}*".format(ANI_calc_dir))
    # new_LIN_object = LIN_Assign.getLIN(Genome_ID=subject_genome, Scheme_ID=3, similarity=ANIb_result, c=c)
    # new_LIN = LIN_Assign.Assign_LIN(new_LIN_object, c=c,current_genome=current_genome).new_LIN
    return ANIb_result

def whatsgoingon():
    # This is happening in Sourmash dir
    sourmash_dir = "/home/linproject/Workspace/Sourmash/"
    conn,c = connect_to_db()
    f = open("/home/linproject/Workspace/Sourmash/unmatch.txt","r")
    Genome_ID = [i.strip() for i in f.readlines()]
    f.close()
    df = pd.read_table("/home/linproject/Workspace/Sourmash/test_result.txt",sep="\t",header=0,index_col=0)
    # sig_pool = {}
    for genome in Genome_ID:
        c.execute("SELECT SignaturePath FROM Signature WHERE Genome_ID={0}".format(genome))
        current_sig_path = c.fetchone()[0]
            # write_both_strand(genome,c,sourmash_dir)
            # current_sig_path = create_signature(genome,sourmash_dir,c,conn)
            # sig_pool[str(genome)] = current_sig_path
        LINgroup = str(df.get_value(int(genome),"G_LINgroup"))
        c.execute("SELECT Genome_ID FROM LIN"
                  " WHERE LIN LIKE '{0}%' "
                  "AND Genome_ID<{1}".format(LINgroup,genome))
        tmp = c.fetchall()
        subject_genomes = [i[0] for i in tmp]
        subject_sig_path = []
        for each_subject in subject_genomes:
            c.execute("SELECT SignaturePath FROM Signature where Genome_ID={0}".format(each_subject))
            sig_path = c.fetchone()[0]
            subject_sig_path.append(sig_path)
        if not isdir(sourmash_dir+"tmp"):
            os.mkdir(sourmash_dir+"tmp")
        for each in subject_sig_path:
            shutil.copy(each,sourmash_dir+"tmp")
        result = sourmash_searching(sourmash_dir,"tmp",current_sig_path,genome)
        os.system("rm -rf tmp/*")
        candidates = result.index
        similarity = []
        similarity_dict = {}
        for each_candidate in candidates:
            ani = assign_LIN_based_on_mash(current_genome=int(genome), subject_genome=each_candidate, c=c)
            similarity.append(ani)
        result["ANI"] = similarity
        result.to_csv("{0}_unmatch_result.csv".format(genome))


# MAIN
if __name__ == "__main__":
    whatsgoingon()
    # test_mash()
    # conn, c = connect_to_db()
    # df = pd.read_table("/home/linproject/Workspace/Sourmash/test_result.txt",sep="\t",header=0,index_col=0)
    # height = len(df.index)
    # mash_based_LIN = []
    # for each_genome in df.index:
    #     if df.get_value(each_genome,"Top_Mash") != df.get_value(each_genome,"SubjectGenome") \
    #             and df.get_value(each_genome,"New_G_LINgroup") == "N":
    #         print int(each_genome), int(df.get_value(each_genome,"Top_Mash"))
    #         new_LIN = assign_LIN_based_on_mash(each_genome,int(df.get_value(each_genome,"Top_Mash")),c,conn)
    #     else:
    #         new_LIN = df.get_value(each_genome,"Assigned_LIN")
    #     mash_based_LIN.append(new_LIN)
    # df["mash_based_LIN"] = mash_based_LIN
    # df.to_csv("/home/linproject/Workspace/Sourmash/mash_LIN.csv")
