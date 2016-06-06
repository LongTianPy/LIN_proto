#!/usr/bin/python

import MySQLdb
from MySQLdb import Connect

def build_db():
	conn = Connect('localhost','root')
	c = conn.cursor()
	c.execute('CREATE DATABASE IF NOT EXISTS LINdb_zika')
	c.execute('USE LINdb_zika')
	c.execute('CREATE TABLE Genome (Genome_ID int NOT NULL AUTO_INCREMENT,'
		'Interest_ID int NOT NULL,'
		'Submission_ID int NOT NULL,'
		'FilePath text NOT NULL,'
        	'GenomeName text NOT NULL,'
            'Contigs int NOT NULL,'
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
	c.execute('CREATE TABLE User (User_ID int NOT NULL AUTO_INCREMENT,LastName varchar(255) NOT NULL,FirstName varchar(255) NOT NULL,Institute varchar(255) NOT NULL,RegistrationDate text NOT NULL,Username varchar(255) NOT NULL,Password text NOT NULL,Email text NOT NULL,PRIMARY KEY (User_ID))')
	c.execute("CREATE TABLE AttributeValue (AttributeValue_ID int NOT NULL AUTO_INCREMENT,Genome_ID int NOT NULL,Interest_ID int NOT NULL, Attribute_ID int not null,AttributeValue text not null,User_ID int not null,Private boolean not null,PRIMARY KEY (AttributeValue_ID))")
    	### INITIALIZATION
	c.execute("insert into User (LastName, FirstName, Institute, RegistrationDate, Username, Password, Email) values ('Tian', 'Long', 'Virginia Tech', '2015-5-17', 'longtian', '123456', 'tianlongapp@gmail.com')")
    	# Attributes
    	# General entries
	c.execute("INSERT INTO Attribute (AttributeName) values ('Strain')")
	c.execute("INSERT INTO Attribute (AttributeName) values ('GPS Coordinates')")
	c.execute("INSERT INTO Attribute (AttributeName) values ('Date of isolation')")
    	c.execute("INSERT INTO Attribute (AttributeName) values ('Genus')")
    	c.execute("INSERT INTO Attribute (AttributeName) values ('Species')")
    	c.execute("INSERT INTO Attribute (AttributeName) values ('DOI')")
    	c.execute("INSERT INTO Attribute (AttributeName) values ('NCBI Accession')")
    	# For Plant pathogen
    	c.execute("INSERT INTO Attribute (AttributeName) values ('Host of isolation')")
	# For Virus
	c.execute("INSERT INTO Attribute (AttributeName) values ('Country discovered')")
	c.execute("INSERT INTO Attribute (AttributeName) values ('Infectious disease name')")
	c.execute("INSERT INTO Attribute (AttributeName) values ('Host')")
    	# Interest
    	c.execute("INSERT INTO Interest (InterestName, Attribute_IDs) values ('Plant pathogen','1,2,3,4,5,6,7,8')")
	c.execute("INSERT INTO Interest (InterestName, Attribute_IDs) values ('Human & Animal viruses','1,2,3,4,5,6,7,9,10,11')")
    	# Scheme
    	c.execute("INSERT INTO Scheme (Cutoff, LabelNum) values ('70,75,80,85,90,95,96,97,98,99.99999', 10)")
    	c.execute("INSERT INTO Scheme (Cutoff, LabelNum) values ('60,70,80,85,90,95,98,99,99.5,99.6,99.7,99.8,99.9,99.91,99.92,99.93,99.94,99.95,99.96,99.97,99.98,99.99,99.999,99.9999', 24)")
    	c.execute("INSERT INTO Scheme (Cutoff, LabelNum) values ('60,70,75,80,85,90,95,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999,99.9999', 20)")
    	# Enter one zika genome to start
    	c.execute("INSERT INTO Genome (Interest_ID, Submission_ID, FilePath, GenomeName) values (2,1,'/home/linproject/Workspace/Zika/init/AY632535.fasta','AY632535')")
    	c.execute("INSERT INTO LIN (Genome_ID, Scheme_ID, SubjectGenome, ANI, LIN) values (1,3,'AY632535', 1, '0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')")
    	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (1,1,2,'',1,True)") # Strain
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (2,1,2,'',1,True)") # GPS
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (3,1,2,'1947',1,True)") # Date
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (4,1,2,'Flavivirus',1,True)") # Genus
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (5,1,2,'Zika',1,True)") # Species
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (6,1,2,'',1,True)") # DOI
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (7,1,2,'AY632535',1,True)") # Accession
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (9,1,2,'Uganda',1,True)") # Country
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (10,1,2,'Microcephaly',1,True)") # Disease
	c.execute("INSERT INTO AttributeValue (Attribute_ID, Genome_ID, Interest_ID, AttributeValue, User_ID, Private) values (11,1,2,'Simiiformes',1,True)") # Host






    	conn.commit()








if __name__ == '__main__':
	build_db()
