#!/usr/bin/python
"""
I pulled bacteria genomes from NCBI refseq and got 6000+ genomes.
"""

# IMPORT
import filecmp
from os import listdir
from os.path import isfile, join
import shutil
from MySQLdb import Connect
import sys

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost","root")
    c = conn.cursor()
    c.execute("use LINdb_RefSeq")
    return conn, c

def check_identical(new_file,c):
    c.execute("select Genome_ID, FilePath from Genome")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    file_paths = [i[1] for i in tmp]
    for i in range(len(file_paths)):
        if not filecmp.cmp(new_file,file_paths[i]):
            duplicated_genome = i
            break
        else:
            continue
    if duplicated_genome:
        return duplicated_genome
    else:
        return None


def load_refseq():
    refseq_path = "/home/linproject/NCBI/"
    meta_data = "/home/linproject/NCBI/assembly_summary_type.txt"
    genome_path = "/home/linproject/NCBI/Genomes/"
    working_dir = "/home/linproject/Workspace/Psy_166/"
    conn, c = connect_to_db()
    refseq_genomes = [join(genome_path,f) for f in listdir(genome_path) if isfile(join(genome_path,f)) and f.endswith("fna")]
    f = open(meta_data,"r")
    meta = [i.strip().split("\t") for i in f.readlines()[2:]]
    f.close()
    accession = [i[0] for i in meta]
    type = []
    filename = []
    files = []
    for i in meta:
        if i[-2].startswith("ftp"):
            url = i[-2].split("/")[-1] + "_genomic.fna"
        else:
            url = i[-3].split("/")[-1] + "_genomic.fna"
        filename.append(join(genome_path,url))
        files.append(url)
    organism_name = [i[7] for i in meta]
    infra = [i[8].split("=")[-1] for i in meta]
    strain = []
    resource = []
    for i in infra:
        if "(" in infra:
            strain.appen(i.split("(")[0])
            resource.append(i.split("(")[1][:-1])
        else:
            strain.append(i)
            resource.append("N/A")
    asm_name = [i[15] for i in meta]
    # filename = [join(genome_path,i+"_genomic.fna") for i in ftp]

    for i in range(len(filename)):
        filepath = filename[i]
        duplicated_genome = check_identical(filepath,c)
        if not duplicated_genome:
            shutil.copy(filepath,working_dir)
            # if organism_name[i][0].isupper():
            #     if "(" in organism_name[i]:
            #         organism = organism_name[i].split("(")[0]
            #     else:
            #         organism = organism_name[i]
            #     if organism.startswith("Candidatus"):
            #         Genus = "_".join(organism.split(" ")[:2])
            #         Species = organism.split(" ")[2]
            #         Subspecies = "N/A"
            #         if strain[i] != "":
            #             Strain = strain[i].replace(" ","_")
            #         else:
            #             Strain = "_".join(organism.split(" ")[3:])
            #     else:
            #         Genus = organism.split(" ")[0]
            #         Species = organism.split(" ")[1]
            #         Subspecies = "N/A"
            #         if strain[i] != "":
            #             Strain = strain[i].replace(" ","_")
            #         else:
            #             Strain = "_".join(organism.split(" ")[2:])
            # else:
            #     Genus = "No_name"
            #     Species = "No_name"
            #     Subspecies = "N/A"
            #     if strain[i] != "":
            #         Strain = strain[i].replace(" ","_")
            #     else:
            #         Strain = organism_name[i].replace(" ","_")
            # Type_strain = "Yes"
            # NCBI_Accession = accession[i]
            # attribute = "^^".join([Genus,Species,Subspecies,Strain,Type_strain,NCBI_Accession])
            # cmd = "python /home/linproject/Projects/LIN_proto/workflow -i {0} -u 2 -s 6 -t {1}".format(files[i],attribute)
            # print cmd

def run_cmd():
    working_dir = "/home/linproject/Workspace/Psy_166/"
    files = [f for f in listdir(working_dir)]
    f = open("/home/linproject/NCBI/cmd.txt","r")
    lines = [i.strip() for i in f.readlines()]
    f.close()
    for i in lines:
        os.system(i)



if __name__ == '__main__':
    run_cmd()




            #     if "(" in organism_name[i]:
            #         organism = organism_name[i].split("(")[0]
            #     else:
            #         organism = organism_name[i]
            #     if organism.startswith("Candidatus"):
            #         Genus = "_".join(organism.split(" ")[:2])
            #         Species = organism.split(" ")[2]
            #         if "subsp." in organism:
            #             Subspecies = organism.split(" ")[4]
            #             if len(organism.split(" ")) > 5:
            #                 Strain = "_".join(organism.split(" ")[5:])
            #             else:
            #                 Strain = "_".join(strain[i])
            #         elif "pv."  in organism:
            #             Subspecies = "_".join(organism.split(" ")[3:5])
            #             if len(organism.split(" ")) > 5:
            #                 Strain = "_".join(organism.split(" ")[5:])
            #             else:
            #                 Strain = "_".join(strain[i])
            #         else:
            #             Subspecies = "N/A"
            #             if len(organism.split(" ")) > 3:
            #                 Strain = "_".join(organism.split(" ")[3:])
            #             else:
            #                 Strain = "_".join(strain[i])
            #     else:
            #         Genus = organism.split(" ")[0]
            #         Species = organism.split(" ")[1]
            #         if "subsp." in organism:
            #             Subspecies = organism.split(" ")[3]
            #             if len(organism.split(" ")) > 4:
            #                 Strain = "_".join(organism.split(" ")[4:])
            #             else:
            #                 Strain = "_".join(strain[i])
            #         elif "pv." in organism:
            #             Subspecies = "_".join(organism.split(" ")[2:4])
            #             if len(organism.split(" ")) > 4:
            #                 Strain = "_".join(organism.split(" ")[4:])
            #             else:
            #                 Strain = "_".join(strain[i])
            #         else:
            #             Subspecies = "N/A"
            #             if len(organism.split(" ")) > 2:
            #                 Strain = "_".join(organism.split(" ")[2:])
            #             else:
            #                 Strain = "_".join(strain[i])
            # else:









# MAIN