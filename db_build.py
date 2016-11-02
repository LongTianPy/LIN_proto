#!/usr/bin/python

import MySQLdb
from MySQLdb import Connect


def build_db():
    conn = Connect('localhost', 'root')
    c = conn.cursor()
    c.execute('CREATE DATABASE IF NOT EXISTS LINdb_Psy')
    c.execute('USE LINdb_Psy')
    c.execute('CREATE TABLE Genome (Genome_ID int NOT NULL AUTO_INCREMENT,'
              'Interest_ID int NOT NULL,'
              'Submission_ID int NOT NULL,'
              'FilePath text NOT NULL,'
              'GenomeName text NOT NULL,'
              'PRIMARY KEY (Genome_ID))')
    c.execute('CREATE TABLE Interest (Interest_ID int NOT NULL AUTO_INCREMENT,'
              'InterestName varchar(255) NOT NULL,'
              'Attribute_IDs text NOT NULL,'
              'PRIMARY KEY (Interest_ID))')
    c.execute('CREATE TABLE Attribute (Attribute_ID int NOT NULL AUTO_INCREMENT,'
              'AttributeName varchar(255) NOT NULL,'
              'PRIMARY KEY (Attribute_ID))')
    c.execute('CREATE TABLE LIN (LIN_ID int NOT NULL AUTO_INCREMENT,'
              'Genome_ID int NOT NULL,'
              'Scheme_ID int NOT NULL,'
              'SubjectGenome int NOT NULL,'
              'ANI DOUBLE NOT NULL,'
              'LIN text NOT NULL,'
              'PRIMARY KEY (LIN_ID))')
    c.execute('CREATE TABLE Scheme (Scheme_ID int NOT NULL AUTO_INCREMENT,'
              'Cutoff text(255) NOT NULL,'
              'LabelNum int NOT NULL,'
              'PRIMARY KEY (Scheme_ID))')
    c.execute('CREATE TABLE Submission (Submission_ID int NOT NULL AUTO_INCREMENT,'
              'User_ID int NOT NULL,'
              'Time text NOT NULL,'
              'PRIMARY KEY (Submission_ID))')
    c.execute(
        'CREATE TABLE User (User_ID int NOT NULL AUTO_INCREMENT,FirstName varchar(255) NOT NULL,LastName varchar(255) '
		'NOT NULL,Institute varchar(255) NOT NULL,RegistrationDate text NOT NULL,Username varchar(255) NOT NULL,'
		'Password text NOT NULL,Email text NOT NULL,PRIMARY KEY (User_ID))')
    c.execute(
        "CREATE TABLE AttributeValue (AttributeValue_ID int NOT NULL AUTO_INCREMENT,Genome_ID int NOT NULL,"
		"Interest_ID int NOT NULL, Attribute_ID int not null,AttributeValue text not null,User_ID int not null,"
		"Private boolean not null,PRIMARY KEY (AttributeValue_ID))")
    c.execute("CREATE TABLE Genome_to_Attribute (Genome_to_Attribute_ID int NOT NULL AUTO_INCREMENT, "
              "Genome_ID INT NOT NULL, Genus text, Species text,Subspecies__Pathovar text, Strain text, Type_strain text,"
              " GPS_Coordinates text, Date_of_isolation text,"
              "NCBI_Accession text, Country text, Link_to_peer_reviewed_paper text,Host_of_isolation__Environmental_source text,"
              "Disease text, Infectious_disease_name text,"
              "Host text, PRIMARY KEY (Genome_to_Attribute_ID))")
    c.execute(
        "CREATE TABLE Description (Description_ID int not null auto_increment, LIN_group text NOT NULL ,"
		"Description_Item_ID int not null, DescriptionValue text not null, User_ID INT NOT NULL , PRIMARY KEY (Description_ID))")
    c.execute(
        "CREATE TABLE Description_Items (Description_Item_ID int not null auto_increment, Description_Item_Name text "
		"not null, PRIMARY KEY (Description_Item_ID))")
    c.execute(
        "CREATE TABLE LIN_to_Description (LIN_to_Description_ID int NOT NULL AUTO_INCREMENT, LIN_group text NOT NULL, "
		"Genus text, Species text, Nickname text, Comment text, URL text, PRIMARY KEY(LIN_to_Description_ID))")
    ### INITIALIZATION
    # c.execute("insert into User (LastName, FirstName, Institute, RegistrationDate, Username, Password, Email) values
    #  ('Tian', 'Long', 'Virginia Tech', '2015-5-17', 'longtian', '123456', 'tianlongapp@gmail.com')")
    # Attributes
    # General entries
    c.execute("INSERT INTO Attribute (AttributeName) values ('Genus')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Species')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Subspecies/Pathovar')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Strain')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Type strain')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('GPS Coordinates')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Date of isolation')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('NCBI Accession')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Country')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Link to peer-reviewed paper')")
    # For Plant pathogen
    c.execute("INSERT INTO Attribute (AttributeName) values ('Host of isolation/Environmental source')")
    c.execute("INSERT INTO Attribute (AttributeName) VALUES ('Disease')")
    # For Virus
    c.execute("INSERT INTO Attribute (AttributeName) values ('Infectious disease name')")
    c.execute("INSERT INTO Attribute (AttributeName) values ('Host')")
    # Interest
    c.execute(
        "INSERT INTO Interest (InterestName, Attribute_IDs) values ('Plant pathogen','1,2,3,4,5,6,7,8,9,10,11,12')")
    c.execute(
        "INSERT INTO Interest (InterestName, Attribute_IDs) values ('Human & Animal viruses','1,2,4,5,6,7,9,10,13,14')")
    # Scheme
    c.execute("INSERT INTO Scheme (Cutoff, LabelNum) values ('70,75,80,85,90,95,96,97,98,99.99999', 10)")
    c.execute(
        "INSERT INTO Scheme (Cutoff, LabelNum) values ('60,70,80,85,90,95,98,99,99.5,99.6,99.7,99.8,99.9,99.91,99.92,"
		"99.93,99.94,99.95,99.96,99.97,99.98,99.99,99.999,99.9999', 24)")
    c.execute(
        "INSERT INTO Scheme (Cutoff, LabelNum) values ('60,70,75,80,85,90,95,98,98.5,99,99.25,99.5,99.75,99.9,99.925,"
		"99.95,99.975,99.99,99.999,99.9999', 20)")
    # Description Items
    c.execute("INSERT INTO Description_Items (Description_Item_Name) values ('Genus')")
    c.execute("INSERT INTO Description_Items (Description_Item_Name) values ('Species')")
    c.execute("INSERT INTO Description_Items (Description_Item_Name) values ('Nickname')")
    c.execute("INSERT INTO Description_Items (Description_Item_Name) values ('Comment')")
    c.execute("INSERT INTO Description_Items (Description_Item_Name) values ('URL')")

    # # Enter one genome to start

    # Start with Genome of DC3000, which is GCA_000007805.1.fasta
    c.execute("insert into User (LastName, FirstName, Institute, RegistrationDate, Username, Password, Email) values "
              "('Tian', 'Long', 'Virginia Tech', '2016-10-21', 'kingdom586@hotmail.com', '123456', 'kingdom586@hotmail.com')")
    c.execute("INSERT INTO Submission (User_ID, Time) VALUES (1,'2016-11-02')")
    c.execute(
        "INSERT INTO Genome (Interest_ID, Submission_ID, FilePath, GenomeName) VALUES (1, 1, '/home/linproject/Workspace/Psy_166/init/GCA_000007805.fasta', 'GCA_000007805')")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 1, 'Pseudomonas', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 2, 'syringae', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 3, 'pv. tomato', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 4, 'DC3000', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 5, 'No', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 6, '49.5,-2.5', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 7, 'NA', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 8, 'AE016853', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 9, 'United Kingdom', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 10, '10.1073/pnas.1731982100', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 11, 'Solanum lycopersicum', 1, False)")
    c.execute(
        "INSERT INTO AttributeValue (Genome_ID, Interest_ID, Attribute_ID, AttributeValue, User_ID, Private) VALUES (1, 1, 12, 'Plant rot, Speck disease', 1, False)")
    # c.execute(
    #     "INSERT INTO LIN (Genome_ID, Scheme_ID, SubjectGenome, ANI, LIN) VALUES (1, 3, 1, 1, '0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')")
    conn.commit()


if __name__ == '__main__':
    build_db()

#B728a	Pseudomonas	syringae	Pseudomonas syringae pv. syringae B728a (g-proteobacteria)		GCA_000012245.1	6093698	1	ASM1224v1	43.7555*-89.055234	snap bean leaflet	Plant rot		USA	CP000075.1
#M301072	Pseudomonas	syringae	Pseudomonas syringae pv. japonica str. M301072 (g-proteobacteria)		GCA_000145785.1	6516319	4661	ASM14578v1			Plant rot	10.1371/journal.ppat.1002132		AEAH01000001.1