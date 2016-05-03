#!/usr/bin/python

def LoadInfo(file, db_cursor,genome_name, interest_id, user=1, private=True):
    db_cursor.execute("SELECT Genome_ID from Genome where Genome.GenomeName = '{0}'".format(genome_name))
    Genome_ID = int(db_cursor.fetchone()[0])
    f = open(file,"r")
    lines = [line.split(",") for line in f.readlines()[1:]]
    f.close()
    accession = [line[6] for line in lines]
    idx = accession.index(genome_name)
    db_cursor.execute("SELECT Attribute_IDs FROM Interest WHERE Interest_ID={0}".format(interest_id))
    Attribute_IDs = db_cursor.fetchall()
    Attribute_IDs = [i[0] for i in Attribute_IDs]
    info_to_be_loaded = lines[idx]
    assert len(info_to_be_loaded) == len(Attribute_IDs),"The number of attributes to be added does not match with database record."
    for i in range(len(info_to_be_loaded)):
        for j in range(len(info_to_be_loaded[i])):
            query = "INSERT INTO AttributeValue (Attribute_ID, Genome_ID, AttributeValue, User_ID, Private) values" \
                    "{0}, {1}, {2}, {3}, {4}".format(Attribute_IDs[j],Genome_ID, info_to_be_loaded[i][j], user, private)
            db_cursor.execute(query)

