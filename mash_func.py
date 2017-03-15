#!/usr/bin/python
"""A bunch of functions related to mash or sourmash
"""

# IMPORT
import os
from Bio import SeqIO
from os.path import join
from MySQLdb import Connect
from LINgroup_indexing_comparison import fetch_current, fetch_genomes
from os.path import isdir
import shutil

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
    sourmash_cmd = "sourmash compute -o {0} {1}.fasta".format(sig_path,sourmash_dir+str(Genome_ID))
    os.system(sourmash_cmd)
    cursor.execute("insert into Signature (Genome_ID,SignaturePath) values ({0}, '{1}')".format(Genome_ID, sig_path))
    conn.commit()
    return sig_path

def sourmash_indexing(sourmash_dir, LINgroup):
    target_folder = sourmash_dir+LINgroup+"/"
    cmd = "sourmash sbt_index {0}index {0}*.sig".format(target_folder)
    os.system()

def sourmash_searching(sourmash_dir,LINgroup,current_sig_path):
    target_folder = sourmash_dir + LINgroup + "/"
    cmd = "sourmash sbt_search {0} {1} > {0}result.txt".format(target_folder,current_sig_path)
    os.system(cmd)
    f = open("{0}result.txt".format(target_folder),"r")
    lines = [i.strip().split(" ") for i in f.readlines()]
    f.close()
    return lines

# def check_result():

def test_mash():
    sourmash_dir = "/home/linproject/Workspace/Sourmash/"
    output_handler = open(sourmash_dir+"test_result.txt","w")
    output_handler.write("Genome_ID\tAssigned_LIN\tNew_G_LINgroup\tG_LINgroup\tSearching_G_LINgroup\tTop_Mash\tTop_Mash_D\tSubjectGenome\tANI\n")
    line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n"
    conn, c = connect_to_db()
    full_df, All_genomes = fetch_genomes(c)
    write_both_strand(All_genomes[0], c, sourmash_dir)
    sig_path0 = create_signature(All_genomes[0],sourmash_dir,c,conn)
    if not isdir(sourmash_dir+"0,0,0,0,0,0,0/"):
        os.mkdir(sourmash_dir+"0,0,0,0,0,0,0/")
    shutil.copyfile(sig_path0,sourmash_dir+"0,0,0,0,0,0,0/")
    sourmash_indexing(sourmash_dir=sourmash_dir,LINgroup="0,0,0,0,0,0,0")
    if not isdir(sourmash_dir+"rep_bac/"):
        os.mkdir(sourmash_dir+"rep_bac/")
    shutil.copyfile(sig_path0, sourmash_dir + "rep_bac/")
    sourmash_indexing(sourmash_dir=sourmash_dir, LINgroup="rep_bac")
    output_handler.write(line.format(All_genomes[0],",".join([0]*20),"Y","0,0,0,0,0,0,0",1,1,1,1))
    for idx in range(1,len(All_genomes)):
        current_genome, current_df = fetch_current(full_df,All_genomes[idx],idx)
        write_both_strand(current_genome,c,sourmash_dir)
        sig_path_current = create_signature(current_genome,sourmash_dir,c,conn)
        current_LIN = full_df.get_values(current_genome,"LIN") # A string
        current_G_LINgroup = ",".join(current_LIN.split(",")[:7])
        current_df["G_LINgroup"] = [",".join(each.split(",")[:7]) for each in current_df["LIN"]]
        subject_genome = full_df.get_values(current_genome,"SubjectGenome")
        ani = full_df.get_values(current_genome,"ANI")
        # Trick is, mash d may or may not respond to 90% ani.
        # We are using 95% ani.
        # First thing is to check if it correlates with representative genomes.
        # If yes, there is change that the ani is between 90% and 95%. Need to check its LIN to determine if it's
        # in any LINgroup or creates a new one.
        if current_G_LINgroup not in current_df["G_LINgroup"]:
            os.mkdir(sourmash_dir+current_G_LINgroup+"/")
            shutil.copyfile(sig_path_current,sourmash_dir+current_G_LINgroup+"/")
            sourmash_indexing(sourmash_dir,current_G_LINgroup)
            shutil.copyfile(sig_path_current,sourmash_dir+"rep_bac/")
            sourmash_indexing(sourmash_dir,"rep_bac")
            result_rep = sourmash_searching(sourmash_dir,"rep_bac",sig_path_current)
            if len(result_rep) == 2:
                output_handler.write(line.format(current_genome,current_LIN,"Y",current_G_LINgroup,"rep_bac","NA","NA",
                                                 subject_genome,ani))
            else:
                top_rep_mash = result_rep[2][1].split("/")[-1].split(".")[0]
                top_rep_mash_d = result_rep[2][0]
                output_handler.write(line.format(current_genome,current_LIN,"Y",current_G_LINgroup,"rep_bac",top_rep_mash,
                                                 top_rep_mash_d,
                                                 subject_genome,ani))
        else:
            result_rep = sourmash_searching(sourmash_dir,"rep_bac",sig_path_current)
            if len(result_rep) == 1: # Means you are screwed at 95% level??? Hopefully this will not happen.
                output_handler.write(line.format(current_genome,current_LIN,"N",current_G_LINgroup,"NA","NA",
                                                 subject_genome,ani))
            else:
                top_rep_mash = result_rep[1][1].split("/")[-1].split(".")[0]
                top_rep_mash_G_LINgroup = current_df.get_values(top_rep_mash,"G_LINgroup")
                result = sourmash_searching(sourmash_dir,top_rep_mash_G_LINgroup,sig_path_current)
                top_mash = result[1][1].split("/")[-1].split(".")[0]
                top_mash_d = result[1][0]
                output_handler.write(line.format(current_genome,current_LIN,"N",current_G_LINgroup,
                                                 top_rep_mash_G_LINgroup,top_mash,top_mash_d,
                                                 subject_genome,ani))
            shutil.copyfile(sig_path_current,sourmash_dir+top_rep_mash_G_LINgroup+"/")
            sourmash_indexing(sourmash_dir,top_rep_mash_G_LINgroup)
    output_handler.close()

# MAIN
if __name__ == "__main__":
    test_mash()