#!/usr/bin/python

"""
This script is a wrapper for the calculation reduction combining k-mer, k-means and ANI
"""

import k_mer
import ANI_Wrapper_2
import IntermediateResult
from sklearn.cluster import KMeans
import MySQLdb
import pandas as pd


def k_methods(queryfilepath, Genome_ID, User_ID, db_cursor, ):
    similarity = k_mer.generate_distance(queryfilepath=queryfilepath, Genome_ID=Genome_ID, User_ID=User_ID)
    if len(similarity['Genome']) <= 10:
        n_top = len(similarity['Genome'])
    else:
        n_clusters = 3
        km = KMeans(n_clusters=n_clusters)
        km.fit(similarity[Genome_ID].reshape(-1,1))
        centroid_idx = list(km.cluster_centers_).index(max(km.cluster_centers_))
        top_cluster_idx = [i for i, x in enumerate(km.labels_) if x == centroid_idx]
        n_top = len(top_cluster_idx)
    top10 = similarity.head(n_top)['Genome'].values
    IntermediateResult.write_kmer_result(top10=top10, db_cursor=db_cursor, User_ID=User_ID)
    IntermediateResult.send_email(file_source="kmer", User_ID=User_ID, db_cursor=db_cursor)
    similarities = pd.DataFrame()
    for i in top10:
        db_cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(i))
        target_filepath = c.fetchone()[0]
        target_filename_rename = target_filepath.split("/")[-1]
        target_filename_rename = "{0}.fasta".format(i)
