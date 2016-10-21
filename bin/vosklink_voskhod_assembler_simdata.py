#! /usr/bin/env python
# -*- coding: UTF8 -*-

import datetime
import math
import sys
from sys import stdout
import sqlite3
import os
from collections import *
from Bio import SeqIO
import string

from optparse import OptionParser
os.system("rm -rf ./bin/SPAdes-3.6.2-Linux")
os.system("tar xf ./bin/SPAdes-3.6.2-Linux.tar.gz -C ./bin/")
os.system("rm -rf ./bin/CAP3")
os.system("tar xf ./bin/cap3.linux.x86_64.tar.gz -C ./bin/")


coveragewhished = 75


parser = OptionParser()
parser.add_option("-i", "--input", dest="inp", help="")
parser.add_option("-r", "--refspecies", dest="refsp", help="")
parser.add_option("-n", "--namespecies", dest="namespe", help="")
(options, args) = parser.parse_args()
inp = options.inp
refsp = options.refsp
namespe = options.namespe

#species = "Test_VV1"
#spe_model = "Danio_rerio"
species = namespe
spe_model = refsp


class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
os.system("mkdir -p ./logs")
logfile = str("./logs/Vosklink_" + species + "_" + str(datetime.datetime.now().strftime("%d-%m-%y")) + "_log.txt")
sys.stdout = Logger(logfile)



print """
       .
      /|\               .
     / | \            ./|\,
  ,-' \|/ `-.        <-=O=->
<'--==<O>==--`>       '\|/`
  `-. /|\ ,-'           '
     \ | /
      \|/
       '
 _    __           __   ___       __
| |  / /___  _____/ /__/ (_)___  / /__
| | / / __ \/ ___/ //_/ / / __ \/ //_/
| |/ / /_/ (__  ) ,< / / / / / / ,<
|___/\____/____/_/|_/_/_/_/ /_/_/|_|

Version 20160920A
Spades & CAP3 conductor
Part of the Voskhod project

(C) Arnaud Ungaro
contact@arnaud-ungaro.fr"
"""


conn = sqlite3.connect(inp)
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
conn.commit()

conn2 = sqlite3.connect("./data_input/cdna_infos.db")
c2 = conn2.cursor()
c2.execute('''PRAGMA synchronous = OFF''')
c2.execute('''PRAGMA journal_mode = OFF''')
c2.execute('''PRAGMA cache_size = 4096''')
conn2.commit()

namenew = str(inp.split("/")[-1]) + "_" + "cdna_infos.db"
os.system("rm -rf ./data_input/" + namenew)
conn3 = sqlite3.connect('./data_input/' + namenew)
c3 = conn3.cursor()
c3.execute('''PRAGMA synchronous = OFF''')
c3.execute('''PRAGMA journal_mode = OFF''')
c3.execute('''PRAGMA cache_size = 4096''')
c3.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id varchar , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar, reads int, coverage_avg_tsize float, cov_mintsize float, cov_maxtsize float)''')
c3.execute('''CREATE TABLE stats (geneid varchar, genename varchar, contigs varchar, ratiotriger varchar, ratiotrigermin , ratiotrigermax varchar, countcap3 varchar ,countspades varchar ,countread varchar)''')

conn3.commit()

c2.execute("""SELECT * FROM result WHERE species = '""" + str(spe_model) + """' ORDER BY gene_name,transcript_size""")
conn2.commit()
match = c2.fetchone()
dicoensdarg = OrderedDict()
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
    dicoensdarg[geneid] = {"ensdarg": geneid}
    dicoensdarg[geneid]["gene_name"] = genename
    dicoensdarg[geneid]["chromosome"] = chromosome
    dicoensdarg[geneid]["posdeb"] = posdebgn
    dicoensdarg[geneid]["posfin"] = posfingn
    dicoensdarg[geneid]["protid"] = protid
    dicoensdarg[geneid]["namedesc"] = namedesc
    match = c2.fetchone()

c.execute("""SELECT COUNT(gene_id) FROM Result WHERE Total_Hits != 0""")
conn.commit()
match = c.fetchone()
nb_genes = float(match[0])

c.execute("""SELECT gene_id, Total_hits, AvgRead, transcript_biotype, gene_size, transcript_size_Min, transcript_size_Max, AvgCov FROM result WHERE Total_Hits != 0 ORDER BY Gene_name, gene_id""")
conn.commit()
match = c.fetchone()
count = 0
coveragewhished = float(coveragewhished)
while match is not None:
    count += 1
    countass = 0
    countspades = 0
    countcap3 = 0
    countread = 0
    sumts = 0
    listsizes = []
    geneid = str(match[0])
    totalhits = int(match[1])
    avgread = float(match[2])
    tsbiotyp = str(match[3])
    gnsize = float(match[4])
    tsminsize = float(match[5])
    tsmaxsize = float(match[6])
    limitfq = int((gnsize * coveragewhished) / avgread)
    limitfa = int((gnsize * coveragewhished) / avgread)
    covdispo = (avgread * totalhits ) / gnsize
    posdeb = dicoensdarg[geneid]["posdeb"]
    posfin = dicoensdarg[geneid]["posfin"]
    chromosome = dicoensdarg[geneid]["chromosome"]
    gname = dicoensdarg[geneid]["gene_name"]
    namedesc = dicoensdarg[geneid]["namedesc"]
    protid = dicoensdarg[geneid]["protid"]
    namefile = str("./fastq/" + str(gname) + "_" + str(geneid) + ".fastq")
    os.system("rm -rf ./runass")
    os.system("mkdir ./runass")
    os.system("cp " + namefile + " ./runass/")
    print "########################################################"
    print str(count) + "/" + str(int(nb_genes)) + " " + str(gname) + " " + str(geneid)
    print "Available coverage : " + str(int(covdispo))
    print "Available reads    : " + str(totalhits)
    print "Wished coverage    : " + str(int(coveragewhished))
    print "Wished reads       : " + str(int(limitfa))
    if totalhits > limitfa:
        print "Resampling available reads.."
    if limitfa > totalhits:
        limitfa = totalhits
    os.system("python ./limitfq.py -i ./runass/" + str(gname) + "_" + str(geneid) + ".fastq -n " + str(int(limitfq)) + " -o ./runass/runfq.fq" )
    if totalhits > 1:
        try:
            print "SPADES Assembling starting.."
            os.system("./bin/SPAdes-3.6.2-Linux/bin/spades.py -t 8 --careful --disable-gzip-output --only-assembler --cov-cutoff auto -s ./runass/runfq.fq -o ./runass >/dev/null 2>/dev/null")
            #os.system("./bin/SPAdes-3.6.2-Linux/bin/spades.py -t 8 --careful --disable-gzip-output --cov-cutoff auto -s ./runass/runfq.fq -o ./runass >/dev/null 2>/dev/null")
            #os.system("./bin/SPAdes-3.6.1-Linux/bin/spades.py -t 8 --cov-cutoff auto --untrusted-contigs ./runass/runfa.fa.cap.contigs.fasta -s ./runass/runfq.fq -o ./runass >/dev/null 2>/dev/null")
            print "Done"
        except:
            print "ERROR WITH SPADES"
    else:
        print "Only 1 read !"
    os.system("touch ./runass/contigs.fasta")

    try:
        for record in SeqIO.parse("./runass/contigs.fasta", "fasta"):
            countass += 1
            countspades += 1
            sqname = "TS_" + str(countass) + "_" + geneid + "_" + species
            sq = str(record.seq)
            sumts += len(sq)
            listsizes.append(len(sq))
            avgcovcalc = round(((len(sq)) / gnsize), 3)
            avgcovTsmincalc = round(((len(sq)) / tsminsize), 3)
            avgcovTsmaxcalc = round(((len(sq)) / tsmaxsize), 3)
            c3.execute('''INSERT INTO result VALUES ("''' + str(sqname) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(gname) + '''","''' + str(namedesc) + '''","''' + str(tsbiotyp) + '''","''' + str(len(sq)) + '''","''' + str(posdeb) + '''","''' + str(posfin) + '''","''' + str(posdeb) + '''","''' + str(posfin) + '''","''' + str(chromosome) + '''","''' + str(species) + '''","''' + str("SPADES") + '''","''' + str(sq) + '''","''' + str(limitfa) + '''","''' + str(avgcovcalc) + '''","''' + str(avgcovTsmincalc) + '''","''' + str(avgcovTsmaxcalc) + '''")''')
    except:
        print "ERROR PARSING ./runass/contigs.fasta"
    try:
        if countspades == 0 or covdispo <= 20:
            print "Spades assembly failed or cov < 20, trying CAP3.."
            os.system("python ./limitfa.py -i ./runass/" + str(gname) + "_" + str(geneid) + ".fastq -n " + str(int(limitfa)) + " -o ./runass/runfa.fa" )
            if totalhits > 1:
                try:
                    print "CAP3 Assembling starting.."
                    os.system("./bin/CAP3/cap3 ./runass/runfa.fa >/dev/null 2>/dev/null")
                    #os.system("./bin/CAP3/cap3 ./runass/runfa.fa")
                    print "Done"
                except:
                    print "ERROR WITH CAP3"
            os.system("touch ./runass/runfa.fa.cap.contigs")
            os.system("cp ./runass/runfa.fa.cap.contigs ./runass/runfa.fa.cap.contigs.fasta")

            for record in SeqIO.parse("./runass/runfa.fa.cap.contigs.fasta", "fasta"):
                countass += 1
                countcap3 += 1
                sqname = "TS_" + str(countass) + "_" + geneid + "_" + species
                sq = str(record.seq)
                sumts += len(sq)
                listsizes.append(len(sq))
                avgcovcalc = round(((len(sq)) / gnsize), 3)
                avgcovTsmincalc = round(((len(sq)) / tsminsize), 3)
                avgcovTsmaxcalc = round(((len(sq)) / tsmaxsize), 3)
                c3.execute('''INSERT INTO result VALUES ("''' + str(sqname) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(gname) + '''","''' + str(namedesc) + '''","''' + str(tsbiotyp) + '''","''' + str(len(sq)) + '''","''' + str(posdeb) + '''","''' + str(posfin) + '''","''' + str(posdeb) + '''","''' + str(posfin) + '''","''' + str(chromosome) + '''","''' + str(species) + '''","''' + str("CAP3") + '''","''' + str(sq) + '''","''' + str(limitfa) + '''","''' + str(avgcovcalc) + '''","''' + str(avgcovTsmincalc) + '''","''' + str(avgcovTsmaxcalc) + '''")''')
        else:
            print "Spades assembly worked, skiping CAP3"
    except:
        print "ERROR PARSING ./runass/runfa.fa.cap.contigs.fasta"
    try:
        avgtssize = float(sum(listsizes))/len(listsizes)
        maxts = float(max(listsizes))
        mints = float(min(listsizes))
        ratiotriger = round(maxts/gnsize, 3)
        ratiotrigermin = round(maxts/tsminsize, 3)
        ratiotrigermax = round(maxts/tsmaxsize, 3)
    except:
        avgtssize = 0
        maxts = 0
        mints =0
        ratiotriger = 0
        ratiotrigermin = 0
        ratiotrigermax = 0
    #print ratiotriger
    # shunted
    # try:
    #     #if countass == 0 or ratiotrigermin <= 0.5:
    #     if countass == 0 or ratiotrigermin <= 0.5:
    #         countassread = countass
    #         for record in SeqIO.parse("./runass/runfa.fa", "fasta"):
    #             countread += 1
    #             countassread += 1
    #             sqname = "TS_" + str(countassread) + "_" + geneid + "_" + species
    #             sq = str(record.seq)
    #             avgcovcalc = round(((len(sq)) / gnsize), 3)
    #             avgcovTsmincalc = round(((len(sq)) / tsminsize), 3)
    #             avgcovTsmaxcalc = round(((len(sq)) / tsmaxsize), 3)
    #             c3.execute('''INSERT INTO result VALUES ("''' + str(sqname) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(gname) + '''","''' + str(namedesc) + '''","''' + str(tsbiotyp) + '''","''' + str(len(sq)) + '''","''' + str(posdeb) + '''","''' + str(posfin) + '''","''' + str(posdeb) + '''","''' + str(posfin) + '''","''' + str(chromosome) + '''","''' + str(species) + '''","''' + str("READ") + '''","''' + str(sq) + '''","''' + str(totalhits) + '''","''' + str(avgcovcalc) + '''","''' + str(avgcovTsmincalc) + '''","''' + str(avgcovTsmaxcalc) + '''")''')
    # except:
    #     print "ERROR PARSING ./runass/runfa.fa"
    print "Contigs Assembled by CAP3      : " + str(countcap3)
    print "Contigs Assembled by SPADES    : " + str(countspades)
    print "Avg contigs size               : " + str(avgtssize)
    print "Min contigs size               : " + str(mints)
    print "Max contigs size               : " + str(maxts)
    print "READs copied as contigs        : " + str(countread)
    print "Assembling rate MAX_TS/AVG_GN  : " + str(ratiotriger)
    print "Assembling rate MAX_TS/MIN_GN  : " + str(ratiotrigermin)
    print "Assembling rate MAX_TS/MAX_GN  : " + str(ratiotrigermax)
    print "DONE !"
    c3.execute('''INSERT INTO stats VALUES ("''' + str(geneid) + '''","''' + str(gname) + '''","''' + str(countass) + '''","''' + str(ratiotriger) + '''","''' + str(ratiotrigermin) + '''","''' + str(ratiotrigermax) + '''","''' + str(countcap3) + '''","''' + str(countspades) + '''","''' + str(countread) + '''")''')
    conn3.commit()
    match = c.fetchone()
#os.system("rm -rf ./runass")
os.system("rm -rf ./bin/SPAdes-3.6.2-Linux")
os.system("rm -rf ./bin/CAP3")
