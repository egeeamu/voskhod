#! /usr/bin/env python
# -*- coding: UTF8 -*-

from collections import *
import sqlite3
import os
from random import shuffle
import random

dicoensdarg = OrderedDict()


listfiles = os.popen("""find ./tomerge/ -maxdepth 1 -type f """).read().split("\n")[0:-1]
print listfiles

ratenoise = 0.0

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



def nucchange(nuc):
    randomizer = random.randint(1, 3)
    if nuc == "A" or nuc == "a":
        if randomizer == 1:
            nucout = "T"
        if randomizer == 2:
            nucout = "G"
        if randomizer == 3:
            nucout = "C"

    elif nuc == "T" or nuc == "t":
        if randomizer == 1:
            nucout = "A"
        if randomizer == 2:
            nucout = "G"
        if randomizer == 3:
            nucout = "C"

    elif nuc == "G" or nuc == "g":
        if randomizer == 1:
            nucout = "T"
        if randomizer == 2:
            nucout = "A"
        if randomizer == 3:
            nucout = "C"

    elif nuc == "C" or nuc == "c":
        if randomizer == 1:
            nucout = "T"
        if randomizer == 2:
            nucout = "G"
        if randomizer == 3:
            nucout = "A"
    else:
        nucout = nuc

    return nucout


def brutos(sequence, rate):
    listseq = list(sequence)
    lenseq = len(sequence)
    numbernuc = int(round(lenseq * rate, 0))
    listposi = list(xrange(lenseq))
    shuffle(listposi)
    listprelev = listposi[0:numbernuc]
    for i in listprelev:
        listseq[i] = nucchange(listseq[i])
    return "".join(listseq)





#print listfiles

for i in listfiles:
    conn2 = sqlite3.connect(str(i))
    c2 = conn2.cursor()
    c2.execute('''PRAGMA synchronous = OFF''')
    c2.execute('''PRAGMA journal_mode = OFF''')
    c2.execute('''PRAGMA cache_size = 4096''')
    conn2.commit()

    c2.execute("""SELECT * FROM result ORDER BY gene_name,transcript_size DESC,transcript_id""")
    conn2.commit()
    match = c2.fetchone()
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
        sequence = brutos(str(match[14]), ratenoise)
        geneidspecies = geneid + "_" + specie
        try:
            dicoensdarg[geneidspecies]["geneidspecies"] = geneidspecies


        except:
            dicoensdarg[geneidspecies] = {"geneidspecies": geneidspecies}
            c.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(chromosome) + '''","''' + str(specie) + '''","''' + str(dbsource) + '''","''' + str(sequence) + '''")''')
        #print match
        match = c2.fetchone()
    conn.commit()
    print i

c.execute("""SELECT * FROM result ORDER BY species,gene_name,transcript_size""")
conn.commit()
match = c.fetchone()
while match is not None:
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
    #dbsource = str("Zv10biggest30pc")
    sequence = str(match[14])
    c3.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(chromosome) + '''","''' + str(specie) + '''","''' + str(dbsource) + '''","''' + str(sequence) + '''")''')
    match = c.fetchone()
conn3.commit()
