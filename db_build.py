#!/usr/bin/python

import MySQLdb
from MySQLdb import Connect

def build_db():
	conn = Connect('localhost','root')
	c = conn.cursor()
	c.execute('CREATE DATABASE IF NOT EXISTS LINdb_test_Psy_3')
	c.execute('USE LINdb_test_Psy_3')
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
		'SubjectGenome TEXT NOT NULL,'
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
	c.execute('CREATE TABLE AttributeValue (AttributeValue_ID int NOT NULL AUTO_INCREMENT,Genome_ID int NOT NULL,Interest_ID int NOT NULL, Attribute_ID int not null,AttributeValue text not null,User_ID int not null,Private boolean not null,PRIMARY KEY (AttributeValue_ID))')
    	### INITIALIZATION
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
    	# Interest
    	c.execute("INSERT INTO Interest (InterestName, Attribute_IDs) values ('Plant pathogen','1,2,3,4,5,6,7,8,9')")
    	# Scheme
    	c.execute("INSERT INTO Scheme (Cutoff, LabelNum) values ('70,75,80,85,90,95,96,97,98,99.99999', 10)")
    	c.execute("INSERT INTO Scheme (Cutoff, LabelNum) values ('60,70,80,85,90,95,98,99,99.5,99.6,99.7,99.8,99.9,99.91,99.92,99.93,99.94,99.95,99.96,99.97,99.98,99.99,99.999,99.9999', 24)")
    	c.execute("INSERT INTO Scheme (Cutoff, LabelNum) values ('60,70,75,80,85,90,95,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999,99.9999', 20)")
    	conn.commit()








if __name__ == '__main__':
	build_db()
