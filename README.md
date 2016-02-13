# LIN_proto

This repository mainly (at least for now) works for fast (I mean, really fast) generate similarity matrix of the newly uploaded genomes with original genomes in the database.

kPAL https://github.com/LUMC/kPAL is the software package we use to create k-mer profile for each genome.

At first, I tried to use the distance counting method from kPAL to create the similarity matrix, however, it works not fast as we expected for two different k-mer profiles -- we need to concatenate them together first, and this process would take a long time.

Then I thought about to calculate the distances without concatenate, i.e. using for loop and multithreading to calculate the pairwise distance one-by-one, and it would take at best around 60 min to do a 93*70-time calculation, which I think, there is still some methods to speed up.

So I came up with the idea that we read the k-mer counting profile, which is stored in HDF5-binary format, into memory and make everything calculatable in python. And this is what I'm working on.
