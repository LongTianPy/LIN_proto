#!/usr/bin/python
"""New workflow wrapper, make the code and thoughts cleaner and straightforward.

"""

# IMPORT
import LIN_Assign
import sys
from MySQLdb import Connect
import pandas as pd
from pytz import timezone
from Bio import SeqIO, Entrez
import os
from os.path import isdir, isfile, join
import argparse
import filecmp
import uuid
from datetime import datetime
import shutil
import uuid

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2.0/all_sketches/"
rep_bac_dir = "/home/linproject/Workspace/Sourmash2.0/rep_bac/"
sourmash_tmp = "/home/linproject/Workspace/Sourmash2.0/tmp_2/"
sourmash_result = "/home/linproject/Workspace/Sourmash2.0/result/"
bbmap_bin = "/home/linproject/Projects/bbmap/"
bbmap_dir = "/home/linproject/Workspace/bbmap/"
bbmap_rep_bac_dir = "/home/linproject/Workspace/bbmap/rep_bac/"
bbmap_all_sketches = "/home/linproject/Workspace/bbmap/all_sketches/"
bbmap_tmp = "/home/linproject/Workspace/bbmap/tmp_sketches/"
bbmap_results = "/home/linproject/Workspace/bbmap/results/"
original_folder  = '/home/linproject/Workspace/LINdb/'
tmp_folder = '/var/www/html/LINbase/tmp_uploads/'
workspace_dir = '/home/linproject/Workspace/New/workspace/'
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']
# ranks_dict = {'superkingdom':1, 'phylum':2, 'class':3, 'order':4, 'family':5, 'genus':6, 'species':7, 'strain':20}
Entrez.email = "aaa@bb.cc"

# OBJECTS
class DecisionTree(object):
    def __init__(self,ANI,cov,wkid):
        self.decide(ANI=ANI,cov=cov,wkid=wkid)
    def decide(self,ANI,cov,wkid):
        if wkid == "N/A":
            if ANI > 0.7:
                if cov > 0.145:
                    self.same_family = True
                else:
                    self.same_family = False
            else:
                self.same_family = False
        elif wkid == 0:
            if cov > 0.046:
                self.same_family = True
            else:
                self.same_family = False
        elif wkid > 0:
            if ANI > 0.7:
                if cov > 0.229:
                    self.same_family = True
                else:
                    self.same_family = False
            else:
                self.same_family = False


# FUNCTIONS
### Parse arguments
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="LIN platform backend"
    )
    parser.add_argument("-i", dest="new_genome", help="xxxxxx.fasta")
    parser.add_argument("-u", dest="User_ID", help="An interger")
    parser.add_argument("-s", dest="Interest_ID", help="Interest ID")
    parser.add_argument("-t", dest="Taxonomy", help="Taxonomy")
    parser.add_argument("-a", dest="Attributes",help="Attributes")
    # parser.add_argument("-p", dest="privacy", help="Is it private information")
    args = parser.parse_args()
    return args

### Connect to database
def connect_to_db():
    conn = Connect("127.0.0.1", "LINbase","Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

### Get metadata from database
def extract_metadata(c):
    metadata = pd.DataFrame()
    c.execute("SELECT Genome.Genome_ID, Genome.FilePath,LIN.LIN FROM Genome,LIN WHERE Genome.Genome_ID=LIN.Genome_ID AND LIN.Scheme_ID=4")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    FilePath = [i[1] for i in tmp]
    LIN = [i[2] for i in tmp]
    metadata["FilePath"] = FilePath
    metadata["LIN"] = LIN
    metadata.index = Genome_ID
    return metadata

def extract_ranks(c):
    ranks_dict = pd.DataFrame()
    c.execute("select Rank_ID,Rank,Rank_order from Taxonomic_ranks")
    tmp = c.fetchall()
    Rank_ID = [int(i[0]) for i in tmp]
    Rank = [i[1] for i in tmp]
    Rank_order = [i[2] for i in tmp]
    ranks_dict["Rank_ID"] = Rank_ID
    ranks_dict["Rank_order"] = Rank_order
    ranks_dict.index=Rank
    return ranks_dict

def extract_attributes(c):
    c.execute("select * from Attribute")
    tmp = c.fetchall()
    attributes_dict = {i[1]:int(i[0]) for i in tmp}
    return attributes_dict

def extract_taxonomy_by_taxid(tax_id):
    name_list = {rank:[] for rank in ranks}
    handler = Entrez.efetch(db='taxonomy',id=str(tax_id),retmode='xml')
    record = Entrez.read(handler)[0]
    lineage_list = record["LineageEx"]
    for taxon in lineage_list:
        if taxon["Rank"] in ranks:
            name_list[taxon["Rank"]] = [taxon["ScientificName"],taxon["TaxId"]]
    species_name_full = name_list["species"][0]
    genus_name = name_list["genus"][0]
    species_name_simple = species_name_full[len(genus_name)+1:]
    strain_name_full = record['ScientificName']
    strain_name_simple = strain_name_full[len(species_name_full)+1:]
    name_list["species"][0] = species_name_full
    name_list["strain"] = [strain_name_simple,tax_id]
    for i in name_list.keys():
        if name_list[i] == []:
            name_list[i] = ["N/A","N/A"]
    current_strain_name = name_list["strain"][0].replace("=","")
    current_strain_name_list = current_strain_name.split(" ")
    remove_duplicate = []
    for i in current_strain_name_list:
        if i not in remove_duplicate:
            remove_duplicate.append(i)
    name_list["strain"][0] = " ".join(remove_duplicate)
    return name_list

def get_Tax_ID_by_entry(entry):
    handler = Entrez.esearch(db='taxonomy',term=entry)
    record = Entrez.read(handler)
    if record["Count"] != '0':
        return record["IdList"][0]
    else:
        return 'N/A'

def check_and_load(entry,c,conn,Rank_ID,Genome_ID):
    if entry != "":
        c.execute("select exists(select NCBI_Tax_ID from NCBI_Tax_ID where Taxon='{0}' and Rank_ID={1})".format(entry,Rank_ID))
        tmp = c.fetchone()[0]
        if tmp == 1:
            c.execute("select NCBI_Tax_ID from NCBI_Tax_ID where Taxon='{0}' and Rank_ID={1}".format(entry,Rank_ID))
            tax_id = c.fetchone()[0]
        else:
            tax_id = 'N/A'
        if tax_id != 'N/A':
            c.execute('insert into Taxonomy (Genome_ID,Rank_ID,NCBI_Tax_ID,Taxon) values ({0},{1},{2},"")'.format(Genome_ID,Rank_ID,int(tax_id)))
            conn.commit()
        else:
            c.execute("insert into Taxonomy (Genome_ID,Rank_ID,NCBI_Tax_ID,Taxon) values ({0},{1},0,'{2}')".format(Genome_ID,Rank_ID,entry))
            conn.commit()

def check_and_load_w_taxid(tax_list,c,conn,Rank_ID,Genome_ID):
    [taxon,taxid] = tax_list
    c.execute("select exists(select NCBI_Tax_ID from NCBI_Tax_ID where NCBI_Tax_ID={0} and Rank_ID={1})".format(taxid,Rank_ID))
    tmp = c.fetchone()[0]
    if tmp == 1:
        c.execute('insert into Taxonomy (Genome_ID,Rank_ID,NCBI_Tax_ID,Taxon) values ({0},{1},{2},"")'.format(Genome_ID, Rank_ID,
                                                                                                     int(taxid)))
        conn.commit()
    else:
        c.execute("insert into NCBI_Tax_ID (NCBI_Tax_ID,Taxon,Rank_ID) values ({0},'{1}',{2})".format(taxid,taxon,Rank_ID))
        conn.commit()
        c.execute("insert into Taxonomy (Genome_ID,Rank_ID,NCBI_Tax_ID,Taxon) values ({0},{1},{2},'')".format(Genome_ID,Rank_ID,taxid))
        conn.commit()

def fill_all_taxonomy_lineage(c,conn,Genome_ID):
    c.execute("SELECT Rank_ID,NCBI_Tax_ID FROM Taxonomy WHERE Genome_ID={0} AND NCBI_Tax_ID<>0 ORDER BY Rank_ID ASC".format(Genome_ID))
    tmp = c.fetchone()
    if tmp is not None:
        [rank_id, tax_id] = tmp
        # print(tmp)
        while int(rank_id)>1:
            c.execute("select * from NCBI_Tax_ID where NCBI_Tax_ID={0}".format(tax_id,rank_id))
            tmp = c.fetchone()
            c.execute("insert into Taxonomy (Genome_ID, Rank_ID,NCBI_Tax_ID,Taxon) values ({0}, {1}, {2},'')".format(Genome_ID,tmp[3]-1,tmp[1]))
            conn.commit()
            [rank_id,tax_id] = [tmp[3]-1, tmp[1]]


def load_attributes(c,db,Attributes,new_genome,Interest_ID,User_ID,attributes_dict):
    base_sql = "INSERT INTO AttributeValue (Genome_ID,Interest_ID,Attribute_ID,AttributeValue,User_ID) VALUES ({0},{1},{2},'{3}',{4})"
    for i in Attributes:
        sql = base_sql.format(new_genome,Interest_ID,attributes_dict[i],Attributes[i],User_ID)
        c.execute(sql)
        db.commit()

def load_new_metadata_newversion(c,db,Interest_ID,new_genome,Taxonomy,Attributes,User_ID,ranks_dict,standardtime):
    c.execute("INSERT INTO Submission (User_ID, Time) VALUES ({0},'{1}')".format(User_ID, standardtime))
    db.commit()
    c.execute("SELECT Submission_ID FROM Submission where User_ID={0} and Time='{1}'".format(User_ID, standardtime))
    Submission_ID = int(c.fetchone()[0])
    shutil.copy(tmp_folder + new_genome, original_folder + new_genome)
    c.execute("INSERT INTO Genome (Interest_ID, Submission_ID, FilePath,LINgroup) VALUES ({0}, {1}, '{2}','')"
              .format(Interest_ID, Submission_ID, original_folder + new_genome))
    db.commit()
    c.execute("SELECT Genome_ID FROM Genome WHERE Submission_ID={0}".format(Submission_ID))
    new_Genome_ID = int(c.fetchone()[0])
    Tax_ID = Attributes["NCBI Taxonomy ID"]
    try:
        Tax_ID = int(Tax_ID)
        lineage = extract_taxonomy_by_taxid(tax_id=Tax_ID)
        lineage['strain'] = [strain,Tax_ID]
        for rank in ranks:
            if lineage[rank] != ['N/A','N/A']:
                rank_id = ranks_dict.loc[rank,"Rank_ID"]
                # print(lineage[rank])
                check_and_load_w_taxid(lineage[rank],c,db,rank_id,new_Genome_ID)
    except:
        for i in Taxonomy:
            if i == "species":
                full_species_name = Taxonomy["genus"] + " " + Taxonomy["species"]
                check_and_load(full_species_name, c, db, ranks_dict.loc[i, "Rank_ID"], new_Genome_ID)
            else:
                check_and_load(Taxonomy[i],c,db,ranks_dict.loc[i,"Rank_ID"],new_Genome_ID)
        fill_all_taxonomy_lineage(c,db,new_Genome_ID)
    attributes_dict = extract_attributes(c)
    load_attributes(c,db,Attributes,new_Genome_ID,Interest_ID,User_ID,attributes_dict)
    return new_Genome_ID

# def create_sketch(filepath):
#     dest = sourmash_tmp+"tmp.sig"
#     cmd = "sourmash compute -o {0} {1} -k 31 -n 1000 > /dev/null 2>&1".format(dest,filepath)
#     os.system(cmd)

def create_sketch2(filepath,dest):
    cmd = "sourmash compute -o {0} {1} -k 21,31,51 -n 2000 -q".format(dest, filepath)
    os.system(cmd)
    time.sleep(3)
    return dest

def compare_sketch2(query,LINgroup,k):
    if LINgroup == "rep_bac":
        dest = rep_bac_dir
    else:
        dest = sourmash_dir + LINgroup + "/"
    folder_size = len([file for file in os.listdir(dest) if isfile(join(dest,file))])
    cmd = "sourmash search {0} {1}*.sig -n {2} -k 21 -q --threshold 0.0001 -o {3}"
    cmd = cmd.format(query, dest, folder_size, sourmash_result+"tmp_result.txt")
    os.system(cmd)


def parse_result2():
    df = pd.read_csv(sourmash_result + "tmp_result.txt",sep=",",header=0)
    if df.empty:
        return df
    else:
        ids = []
        for each in df['filename']:
            id = int(each.split('/')[-1].split('.')[0])
            ids.append(id)
        df.index = ids
        return df

def update_LINgroup(Genome_ID, c, new_LIN, conn):
    c.execute("select LINgroup_ID,LINgroup from LINgroup")
    tmp = c.fetchall()
    LINgroup_ID = [int(i[0]) for i in tmp]
    LINgroup = [i[1] for i in tmp]
    belongs_to = []
    for i in range(len(LINgroup_ID)):
        if new_LIN.startswith(LINgroup[i]):
            belongs_to.append(str(LINgroup_ID[i]))
    if belongs_to != []:
        c.execute("UPDATE Genome SET LINgroup='{0}' WHERE Genome_ID={1}".format(",".join(belongs_to), Genome_ID))
        conn.commit()

def check_belonged_LINgroups(conserved_LIN,c):
    c.execute("select LINgroup_ID,LINgroup from LINgroup")
    tmp = c.fetchall()
    LINgroup_ID = [int(i[0]) for i in tmp]
    LINgroup = [i[1] for i in tmp]
    belongs_to = []
    for i in range(len(LINgroup_ID)):
        if conserved_LIN.startswith(LINgroup[i]):
            belongs_to.append(LINgroup_ID[i])
    return belongs_to

def Genome_Submission(new_genome,Username,InterestName,Taxonomy,Attributes):
    eastern = timezone("EST")
    currenttime = eastern.localize(datetime.now())
    fmt_time_display = '%Y-%m-%d %H:%M:%S %Z%z'
    standardtime = currenttime.strftime(fmt_time_display)
    db, c = connect_to_db()
    c.execute("SELECT User_ID FROM User WHERE Username='{0}'".format(Username))
    User_ID = int(c.fetchone()[0])
    c.execute("SELECT Interest_ID FROM Interest WHERE InterestName='{0}'".format(InterestName))
    Interest_ID_new_genome = int(c.fetchone()[0])
    metadata = extract_metadata(c)
    ranks_dict = extract_ranks(c)
    new_genome_filepath = tmp_folder + new_genome
    file_duplication = 0
    for i in metadata.index:
        if filecmp.cmp(tmp_folder + new_genome,metadata.get_value(i,"FilePath")):
            file_duplication = 1
            SubjectGenome = int(i)
            break
    if file_duplication == 0:
        tmp_sig_filename = str(uuid.uuid4())+".sig"
        tmp_newgenome_sig = create_sketch2(tmp_folder + new_genome,sourmash_tmp+tmp_sig_filename)
        if metadata.empty:
            # if this is the first genome ever in the database
            new_LIN = ",".join(["0"] * 20)
            top1_Genome_ID = 1
            top1_similarity = 1
            top1_coverage = 1
            conserved_LIN = ""
            SubjectGenome = top1_Genome_ID
            ANIb_result = top1_similarity
            cov_result = top1_coverage
        else:
            compare_sketch2(tmp_newgenome_sig,"rep_bac",'21')
            df = parse_result2()
            if df.empty:
                # print("###########################################################")
                # print("System message:")
                # print("No Jaccard similarity detected, will use LINgroup indexing.")
                # print("###########################################################")
                ## LINgroup indexing
                # new_LIN, SubjectGenome, ANIb_result, cov_result, conserved_LIN = LINgroup_indexing(cursor=c,metadata=metadata,new_genome_filepath=new_genome_filepath)
                new_LIN_object = LIN_Assign.getLIN(Genome_ID=1, Scheme_ID=4,
                                                   similarity=0.6,
                                                   c=c)
                new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
                ANIb_result = 0.6
                cov_result = 0
                conserved_LIN = ""
                SubjectGenome = 1
            else:
                rep_bac_Genome_ID = int(df.index[0])
                rep_bac_LIN = metadata.get_value(rep_bac_Genome_ID,"LIN")
                rep_bac_LINgroup = ",".join(rep_bac_LIN.split(",")[:6])
                jaccard_similarity = df.loc[rep_bac_Genome_ID,'similarity']
                if jaccard_similarity > 0.2475:
                    compare_sketch2(tmp_newgenome_sig, rep_bac_LINgroup,'51')
                    time.sleep(2)
                    df = parse_result2()
                    # print("###########################################################")
                    # print("System message:")
                    # print("100% Jaccard similarity detected, checking duplication.")
                    # print("LIN will be assigned if new genome.")
                    # print("###########################################################")
                    # Same genome found
                    ANIb_result = 0
                    cov_result = 0
                    SubjectGenome = 0
                    # There is a table about same genome, better record it
                    #[new_LIN, ANIb_result,cov_result,conserved_LIN] = [None]*4
                    for each_subject_genome_ID in df.index[:3]:
                        subject_genome_filepath = metadata.get_value(int(each_subject_genome_ID),"FilePath")
                        sub_working_dir = workspace_dir + str(uuid.uuid4()) + "/"
                        if not isdir(sub_working_dir):
                            os.mkdir(sub_working_dir)
                        shutil.copyfile(new_genome_filepath, join(sub_working_dir, "tmp.fasta"))
                        shutil.copyfile(subject_genome_filepath, join(sub_working_dir,"{0}.fasta".format(each_subject_genome_ID)))
                        pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                                    "-i {0} -o {1} -m ANIb --nocompress -f".format(sub_working_dir,join(sub_working_dir,'output'))
                        os.system(pyani_cmd)
                        time.sleep(5)
                        this_ANIb_result = pd.read_table(join(sub_working_dir, "output","ANIb_percentage_identity.tab"), sep="\t",
                                                    header=0,
                                                    index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        this_cov_result = pd.read_table(join(sub_working_dir, "output", "ANIb_alignment_coverage.tab"), sep="\t",
                                                   header=0,
                                                   index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        # os.system("rm -rf {0}".format(sub_working_dir))
                        shutil.rmtree(sub_working_dir)
                        if this_ANIb_result > 0.99999:
                            ANIb_result = this_ANIb_result
                            cov_result = this_cov_result
                            SubjectGenome = each_subject_genome_ID
                            break
                        else:
                            if this_ANIb_result > ANIb_result:
                                ANIb_result = this_ANIb_result
                                cov_result = this_cov_result
                                SubjectGenome = each_subject_genome_ID
                            else:
                                continue
                    new_LIN_object = LIN_Assign.getLIN(Genome_ID=SubjectGenome, Scheme_ID=4,
                                                       similarity=ANIb_result, c=c)
                    new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
                    conserved_LIN = ",".join(new_LIN_object.conserved_LIN)
                elif jaccard_similarity <= 0.2475 and jaccard_similarity > 0.0025:
                    # print("###########################################################")
                    # print("System message:")
                    # print("Jaccard similarity detected, calculating ANIs.")
                    # print("LIN will be assigned.")
                    # print("###########################################################")
                    compare_sketch2(tmp_newgenome_sig, rep_bac_LINgroup, '21')
                    time.sleep(2)
                    df = parse_result2()
                    ANIb_result = 0
                    cov_result = 0
                    SubjectGenome = 0
                    for each_subject_genome_ID in df.index[:3]:
                        sub_working_dir = workspace_dir + str(uuid.uuid4()) + "/"
                        if not isdir(sub_working_dir):
                            os.mkdir(sub_working_dir)
                        subject_genome_filepath = metadata.get_value(int(each_subject_genome_ID), "FilePath")
                        shutil.copyfile(new_genome_filepath, join(sub_working_dir, "tmp.fasta"))
                        shutil.copyfile(subject_genome_filepath, join(sub_working_dir,"{0}.fasta".format(each_subject_genome_ID)))
                        pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                                    "-i {0} -o {1} -m ANIb --nocompress -f".format(sub_working_dir,join(sub_working_dir,'output'))
                        os.system(pyani_cmd)
                        time.sleep(3)
                        this_ANIb_result = pd.read_table(join(sub_working_dir, "output","ANIb_percentage_identity.tab"), sep="\t",
                                                    header=0,
                                                    index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        this_cov_result = pd.read_table(join(sub_working_dir, "output", "ANIb_alignment_coverage.tab"), sep="\t",
                                                   header=0,
                                                   index_col=0).get_value('tmp', str(each_subject_genome_ID))
                        # # os.system("rm -rf {0}".format(sub_working_dir))
                        predict = DecisionTree(ANI=ANIb_result, cov=cov_result, wkid=df.get_value(each_subject_genome_ID,"similarity"))
                        if predict.same_family:
                            break
                        else:
                            continue
                        if this_ANIb_result > ANIb_result:
                            ANIb_result = this_ANIb_result
                            cov_result = this_cov_result
                            SubjectGenome = each_subject_genome_ID
                    new_LIN_object = LIN_Assign.getLIN(Genome_ID=each_subject_genome_ID,Scheme_ID=4,similarity=ANIb_result,c=c)
                    new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
                    conserved_LIN = ",".join(new_LIN_object.conserved_LIN)
                else:
                    each_subject_genome_ID = int(df.index[0])
                    sub_working_dir = workspace_dir + str(uuid.uuid4()) + "/"
                    if not isdir(sub_working_dir):
                        os.mkdir(sub_working_dir)
                    subject_genome_filepath = metadata.get_value(int(each_subject_genome_ID), "FilePath")
                    shutil.copyfile(new_genome_filepath, join(sub_working_dir, "tmp.fasta"))
                    shutil.copyfile(subject_genome_filepath,
                                    join(sub_working_dir, "{0}.fasta".format(each_subject_genome_ID)))
                    pyani_cmd = "python3 /home/linproject/Projects/pyani/average_nucleotide_identity.py " \
                                "-i {0} -o {1} -m ANIb --nocompress -f".format(sub_working_dir,
                                                                               join(sub_working_dir, 'output'))
                    os.system(pyani_cmd)
                    time.sleep(3)
                    ANIb_result = pd.read_table(
                        join(sub_working_dir, "output", "ANIb_percentage_identity.tab"), sep="\t",
                        header=0,
                        index_col=0).get_value('tmp', str(each_subject_genome_ID))
                    cov_result = pd.read_table(join(sub_working_dir, "output", "ANIb_alignment_coverage.tab"),
                                                    sep="\t",
                                                    header=0,
                                                    index_col=0).get_value('tmp', str(each_subject_genome_ID))
                    new_LIN_object = LIN_Assign.getLIN(Genome_ID=each_subject_genome_ID, Scheme_ID=4, similarity=ANIb_result,
                                                       c=c)
                    new_LIN = LIN_Assign.Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
                    conserved_LIN = ""
                    SubjectGenome = each_subject_genome_ID
        c.execute("SELECT EXISTS(SELECT LIN FROM LIN WHERE LIN='{0}')".format(new_LIN))
        duplication = c.fetchone()[0]  # 0 = no, 1 = yes
        if duplication == 0:
            # print("###########################################################")
            # print("System message:")
            # print("New genome uploaded.")
            # print("LIN will be assigned.")
            # print("###########################################################")
            new_genome_ID = load_new_metadata_newversion(c=c,db=db,Interest_ID=Interest_ID_new_genome,
                                                         new_genome=new_genome,Taxonomy=Taxonomy,
                                                         Attributes=Attributes,User_ID=User_ID, ranks_dict=ranks_dict,
                                                         standardtime=standardtime,
                                                        )
            this_95_LINgroup = ",".join(new_LIN.split(",")[:6])
            this_95_LINgroup_path = sourmash_dir + this_95_LINgroup + "/"
            c.execute("INSERT INTO LIN (Genome_ID, Scheme_ID,SubjectGenome,ANI,Coverage,LIN) values "
                      "({0},4,{1},{2},{3},'{4}')".format(new_genome_ID, SubjectGenome, ANIb_result, cov_result,
                                                         new_LIN))
            db.commit()
            c.execute('SELECT FilePath FROM Genome WHERE Genome_ID={0}'.format(new_genome_ID))
            real_new_genome_filepath = c.fetchone()[0]
            new_genome_sig = create_sketch2(real_new_genome_filepath, sourmash_dir + str(new_genome_ID) + ".sig")
            if not isdir(this_95_LINgroup_path):  # It's a new rep_bac
                shutil.copyfile(new_genome_sig, rep_bac_dir + str(new_genome_ID) + ".sig")
                os.mkdir(this_95_LINgroup_path)
                shutil.copyfile(new_genome_sig, this_95_LINgroup_path + str(new_genome_ID) + ".sig")
            else:
                shutil.copyfile(new_genome_sig, this_95_LINgroup_path + str(new_genome_ID) + ".sig")
            update_LINgroup(Genome_ID=new_genome_ID, c=c, new_LIN=new_LIN, conn=db)
            c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
            best_LIN = c.fetchone()[0]
            belongs_to = check_belonged_LINgroups(conserved_LIN,c)
            result = {"new LIN":new_LIN, "best LIN":best_LIN,"ANI":ANIb_result,"LINgroup":conserved_LIN,"LINgroup_IDs":belongs_to}
        else:
            # print("###########################################################")
            # print("System message:")
            # print("Duplicate submission found, recording.")
            # print("###########################################################")
            c.execute("INSERT INTO Duplicated_upload (Reference_Genome_ID,Who_uploads_too) VALUES ({0},{1})".format(
                SubjectGenome, User_ID))
            db.commit()
            c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
            best_LIN = c.fetchone()[0]
            result = {"best LIN": best_LIN}
        if isfile(new_genome_sig):
            os.system("rm {0}".format(tmp_newgenome_sig))
    else:
        # duplicate genome file detected
        # print("###########################################################")
        # print("System message:")
        # print("Duplicate submission found, recording.")
        # print("###########################################################")
        c.execute(
                "INSERT INTO Duplicated_upload (Reference_Genome_ID,Who_uploads_too) VALUES ({0},{1})".format(
                    SubjectGenome,
                    User_ID))
        db.commit()
        c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0}".format(SubjectGenome))
        best_LIN = c.fetchone()[0]
        result = {"best LIN": best_LIN}
    db.commit()
    c.close()
    db.close()
    os.system('rm -rf {0}/*'.format(workspace_dir))
    return result
