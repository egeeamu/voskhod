#! /usr/bin/python
# -*- coding: UTF8 -*-



from math import *
import datetime
from random import shuffle
import sqlite3
import os
from collections import *
import string
import random

# simulated reads come from cdna_infos.db, and must be in ./data_input
size_fragments = 100
coverage_value = 1
overfix = 50
# noise = "N" # "Y" or "N"
noise_rate = 0.00
species = "biggestZv10_30pcdb_cdna.db"
conn = sqlite3.connect("./data_input/cdna_infos.db")

mylist = []
today = datetime.date.today()
mylist.append(today)
today = str(mylist[0])

fastqR1 = str("./data_input/Simu_" + species + "_" + today + "_Size_" + str(size_fragments) + "_cov_" + str(
    int(coverage_value)) + "_overlap_" + str(int(overfix)) + "_Noise_" + str(noise_rate) + "_R1.fq")
os.system("rm -rfv " + fastqR1)
fastqR1w = open(fastqR1, "a+b")
fastqR2 = str("./data_input/Simu_" + species + "_" + today + "_Size_" + str(size_fragments) + "_cov_" + str(
    int(coverage_value)) + "_overlap_" + str(int(overfix)) + "_Noise_" + str(noise_rate) + "_R2.fq")
os.system("rm -rfv " + fastqR2)
fastqR2w = open(fastqR2, "a+b")


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


print """
______ _           _
| ___ \ |         | |
| |_/ / |__   ___ | |__   ___  ___
|  __/| '_ \ / _ \| '_ \ / _ \/ __|
| |   | | | | (_) | |_) | (_) \__ \\
\_|   |_| |_|\___/|_.__/ \___/|___/

¤ Illumina sequencing simulator
All sequences version

Version 20170721
Voskhod Pipeline version V1.2
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(GPL-3.0) 
Arnaud Ungaro contact@arnaud-ungaro.fr

"""
print "Depth     : " + str(coverage_value)
print "Overlap   : " + str(overfix)
print "FragSize  : " + str(size_fragments)
print "Noise rate: " + str(noise_rate)
print "Species   : " + str(species)

os.system("sleep 5")

#conn = sqlite3.connect("./data_input/cdna_infos.db")
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
conn.commit()

#dico_gene_id = OrderedDict()

tr = string.maketrans("acgtACGT", "tgcaTGCA")
#tr2 = string.maketrans("_", "-")

countgeneral = 0
ccc = 0
c.execute("""SELECT distinct(transcript_id) FROM result""")
conn.commit()
match = c.fetchone()
while match is not None:
    # gene_id = str(match[0])
    # dico_gene_id[gene_id] = {"gene_id": gene_id}
    # dico_gene_id[gene_id]["transcript_size"] = 0
    ccc += 1
    match = c.fetchone()

c.execute(
    """SELECT transcript_id,gene_id,gene_name,species,transcript_sequence, transcript_size FROM result ORDER BY RANDOM()""")
conn.commit()
match = c.fetchone()
while match is not None:
    # print dico_gene_id[gene_id]
    size = float(size_fragments)
    cov = coverage_value
    countfrag = 0
    countgeneral += 1
    print str(countgeneral) + "/" + str(ccc)
    Ensdart = str(match[0]).replace("#", "_")
    Ensdarg = str(match[1]).replace("#", "_")
    Gene_name = str(match[2]).replace("#", "_")
    Specie = str(match[3]).replace("#", "_")
    Sequence_ts = str(match[4])
    if noise_rate > 0:
        Seq = brutos(Sequence_ts, noise_rate)
    else:
        Seq = Sequence_ts
    lents = int(match[5])
    sizefragtocut = size + size - overfix

    nbrfrag = int(ceil((cov * lents) / (size * 1)))
    if sizefragtocut > lents:
        sizefragtocut = lents
        if size > sizefragtocut:
            size = sizefragtocut
    listseq = list(Seq)
    permutseq = Seq[::-1]
    permutseq = permutseq.translate(tr)
    permutlistseq = list(permutseq)

    for i in range(nbrfrag):
        countfrag += 1

        if countfrag == 1 or countfrag == 3:
            permut = 0
        elif countfrag == 2 or countfrag == 4:
            permut = 1
        else:
            permut = random.randint(0, 1)

        if permut == 1:
            listrun = list(permutlistseq)
        else:
            listrun = list(listseq)

        distmax = int(lents - sizefragtocut)
        posicut = int(random.randint(0, distmax))
        if countfrag == 1 or countfrag == 2:
            posicut = 0
        if countfrag == 3 or countfrag == 4:
            posicut = distmax
        sizefragtocut = int(sizefragtocut)
        size = int(size)
        frag = listrun[posicut:(posicut + sizefragtocut)]
        listR1 = frag[0:size]
        listR2 = frag[(sizefragtocut - size):sizefragtocut]
        seqR1 = "".join(listR1)
        seqR2 = "".join(listR2)
        seqR2r = seqR2[::-1]
        seqR2rc = seqR2r.translate(tr)

        fastqR1w.write(
            "@" + str(Ensdart) + "#" + str(Ensdarg) + "#" + str(Gene_name) + "#" + str(countfrag) + "#" + str(
                countgeneral) + "#" + str(Specie) + "#R1_" + str(permut) + "\n" + seqR1 + "\n" + "+" + "\n" + "H" * len(
                seqR1) + "\n")
        fastqR2w.write(
            "@" + str(Ensdart) + "#" + str(Ensdarg) + "#" + str(Gene_name) + "#" + str(countfrag) + "#" + str(
                countgeneral) + "#" + str(Specie) + "#R2_" + str(
                permut) + "\n" + seqR2rc + "\n" + "+" + "\n" + "H" * len(seqR2rc) + "\n")

    match = c.fetchone()
