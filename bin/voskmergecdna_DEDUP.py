#! /usr/bin/env python
# -*- coding: UTF8 -*-

from optparse import OptionParser
import math
import sys
import datetime
from sys import stdout
import sqlite3
import os
from collections import *
from Bio import SeqIO
import string
import hashlib
# import numpy as np
# import matplotlib.mlab as mlab
# import matplotlib.pyplot as plts_SsD'
listfiles = os.popen("""find ./tomerge/ -maxdepth 1 -type f """).read().split("\n")[0:-1]
print listfiles

conn = sqlite3.connect(':memory:')
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
c.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar)''')
conn.commit()

os.system("rm -rf ./tomerge/merged.db")
conn3 = sqlite3.connect('./tomerge/merged.db')
c3 = conn3.cursor()
c3.execute('''PRAGMA synchronous = OFF''')
c3.execute('''PRAGMA journal_mode = OFF''')
c3.execute('''PRAGMA cache_size = 4096''')
c3.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar)''')
conn3.commit()


#print listfiles
dicohash = OrderedDict()
for i in listfiles:
    conn2 = sqlite3.connect(str(i))
    c2 = conn2.cursor()
    c2.execute('''PRAGMA synchronous = OFF''')
    c2.execute('''PRAGMA journal_mode = OFF''')
    c2.execute('''PRAGMA cache_size = 4096''')
    conn2.commit()

    c2.execute("""SELECT * FROM result ORDER BY gene_name,transcript_size""")
    conn2.commit()
    match = c2.fetchone()

    listdedup = []
    while match is not None :
        transid = str(match[0])
        geneid = str(match[1])
        protid = str(match[2])
        genename = str(match[3])
        namedesc = str(match[4])
        transtype = str(match[5])
        tssize = str(match[6])
        posdebts = str(match[7])
        posfints = str(match[8])
        posdebgn = str(match[9])
        posfingn = str(match[10])
        chromosome = str(match[11])
        specie = str(match[12])
        dbsource = str(match[13])
        sequence = str(match[14])
        hashid = hashlib.md5(sequence).hexdigest()
        #hashid = len(sequence)
        genespeciesid = geneid + specie
        try:
            test = dicohash[genespeciesid]
        except:
            dicohash[genespeciesid] = {}

        try:
            test = dicohash[genespeciesid][hashid]
        except:
            dicohash[genespeciesid][hashid] = {}

        try:
            dicohash[genespeciesid][hashid]["hits"] += 1
        except:
            dicohash[genespeciesid][hashid]["hits"] = 1
        #dicohash[geneid][hashid][transid] = sequence

        c.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(chromosome) + '''","''' + str(specie) + '''","''' + str(dbsource) + '''","''' + str(sequence) + '''")''')
        #print match
        match = c2.fetchone()
    conn.commit()
    print i


c.execute("""SELECT * FROM result ORDER BY species,gene_name,transcript_size""")
conn.commit()
match = c.fetchone()
countcontinue = 0
while match is not None:
    transid = str(match[0])
    protid = str(match[2])
    genename = str(match[3])
    namedesc = str(match[4])
    transtype = str(match[5])
    tssize = str(match[6])
    posdebts = str(match[7])
    posfints = str(match[8])
    posdebgn = str(match[9])
    posfingn = str(match[10])
    chromosome = str(match[11])
    specie = str(match[12])
    dbsource = str(match[13])
    geneid = str(match[1])
    sequence = str(match[14])
    genespeciesid = geneid + specie
    hashid = hashlib.md5(sequence).hexdigest()
    #hashid = len(sequence)
    duplivalue = dicohash[genespeciesid][hashid]["hits"]
    listdedup.append(duplivalue)
    if duplivalue == "continue":
        match = c.fetchone()
        #print "match = c.fetchone()"
        countcontinue += 1
        continue
    elif duplivalue > 1:
        #print dicohash[geneid][hashid]
        dbsource = "DEDUP"
        dicohash[genespeciesid][hashid]["hits"] = "continue"





    c3.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(chromosome) + '''","''' + str(specie) + '''","''' + str(dbsource) + '''","''' + str(sequence) + '''")''')
    match = c.fetchone()
conn3.commit()
print listdedup.count(1)
print len(listdedup)
print countcontinue