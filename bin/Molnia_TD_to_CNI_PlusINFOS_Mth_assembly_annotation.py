#! /usr/bin/python
# -*- coding: UTF8 -*-

from collections import *
import sqlite3
import os
import sys
from sys import stdout
import hashlib
from Bio.Seq import Seq
from collections import OrderedDict
import time
from decimal import Decimal
import datetime


qq = str(sys.argv[1])
spe = str(sys.argv[2])

conn = sqlite3.connect(qq)
specie = spe
dbsource = "Voskhod"
speciesshort = specie

datetoday = str(time.strftime("%Y_%m_%d"))


logsall = ""
class Logger(object):

    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        global logsall
        logsall += message


#os.system("mkdir -p ./logs")
logfile = str("./" + str("logs_" + specie + "_" + dbsource + "_" + datetoday + ".db") + '_' + str(datetime.datetime.now().strftime("%y-%m-%d_%Hh%Mm%Ss")) + "_log.txt")
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
    __  ___      __      _
   /  |/  /___  / /___  (_)___ _
  / /|_/ / __ \/ / __ \/ / __ `/
 / /  / / /_/ / / / / / / /_/ /
/_/  /_/\____/_/_/ /_/_/\__,_/

Version 20160911A_MTH_PKTX
Convert Voskhod output to usable "cdnainfo.db" database
Part of the Voskhod project

Need reference cdnainfos.db & identified sequences (via Voskhod)

(C) Arnaud Ungaro
contact@arnaud-ungaro.fr"
"""

c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
conn.commit()

os.system("rm -rf ./cdna_infos_" + specie + "_" + dbsource + "_" + datetoday + ".db")
conn2 = sqlite3.connect("./cdna_infos_" + specie + "_" + dbsource + "_" + datetoday + ".db")
c2 = conn2.cursor()
c2.execute('''PRAGMA synchronous = OFF''')
c2.execute('''PRAGMA journal_mode = OFF''')
c2.execute('''PRAGMA cache_size = 4096''')
c2.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar)''')
conn2.commit()

conn9 = sqlite3.connect('./data_input/cdna_infos.db')
c9 = conn9.cursor()
c9.execute('''PRAGMA synchronous = OFF''')
c9.execute('''PRAGMA journal_mode = OFF''')
c9.execute('''PRAGMA cache_size = 4096''')
conn9.commit()


dict_genes = OrderedDict()



c9.execute("""SELECT gene_id, gene_description, transcript_biotype, gene_start_position, gene_end_position, chromosome FROM result""")
conn9.commit()
match = c9.fetchone()
while match is not None:
    gene_id = str(match[0])
    #print gene_id
    dict_genes[gene_id] = {"gene_id": gene_id}
    dict_genes[gene_id]["transcript_biotype"] = []
    match = c9.fetchone()

c9.execute("""SELECT gene_id, gene_description, transcript_biotype, gene_start_position, gene_end_position, chromosome, gene_name FROM result""")
conn9.commit()
match = c9.fetchone()
while match is not None:
    gene_id = str(match[0])
    gene_description = str(match[1])
    transcript_biotype = str(match[2])
    gene_start_position = str(match[3])
    gene_end_position = str(match[4])
    chromosome = str(match[5])
    gene_name = str(match[6])
    dict_genes[gene_id]["gene_description"] = gene_description
    if transcript_biotype not in dict_genes[gene_id]["transcript_biotype"]:
        dict_genes[gene_id]["transcript_biotype"].append(transcript_biotype)
    dict_genes[gene_id]["gene_start_position"] = gene_start_position
    dict_genes[gene_id]["gene_end_position"] = gene_end_position
    dict_genes[gene_id]["chromosome"] = chromosome
    dict_genes[gene_id]["gene_name"] = gene_name
    match = c9.fetchone()

for i in dict_genes.keys():
    dict_genes[i]["transcript_biotype"] = "#".join(dict_genes[i]["transcript_biotype"])

# for i in dict_genes.keys():
#     print dict_genes[i]["transcript_biotype"]



c.execute("""SELECT gene_name, gene_id, transcript_id, read_name, hsp, read_name, strand, read_sequence, multi_hits FROM result WHERE gene_name != '' AND (hsp_size * 1.0) / (read_size * 1.0) >= .7 AND (hsp_size * 1.0) / (read_size * 1.0) <= 1.3 AND identity >= 0.70 ORDER BY gene_name, read_size """)
conn.commit()
match = c.fetchone()
count = 0
megabase = 0
listtailleshsp = []
listtaillesraw = []
listgenid = []
while match is not None :
    tmplistmth = []
    mthlist = str(match[8]).split("|")
    for bob in mthlist:
        splited = bob.split("+")
        splitedgeneid = splited[2]
        if splitedgeneid not in tmplistmth:
            tmplistmth.append(splitedgeneid)
    for wywyw in tmplistmth:
        if len(tmplistmth) >= 1:


            gene_id = str(wywyw)
            geneid = str(wywyw)

            count += 1
            sens = int(match[6])

            sequence = str(match[4])
            rrrrr = str(match[7])
            listtaillesraw.append(len(rrrrr))
            listtailleshsp.append(len(sequence))
            seqseq = Seq(sequence)
            if sens == -1:
                sequence = str(seqseq.reverse_complement())

            transidname = str(match[3]) + "_" + str(count)

            listgenid.append(gene_id)
            protid = ""
            genename = dict_genes[gene_id]["gene_name"]
            namedesc = dict_genes[gene_id]["gene_description"]
            transtype = dict_genes[gene_id]["transcript_biotype"]
            tssize = len(str(match[4]))
            posdebts = dict_genes[gene_id]["gene_start_position"]
            posfints = dict_genes[gene_id]["gene_end_position"]
            posdebgn = dict_genes[gene_id]["gene_start_position"]
            posfingn = dict_genes[gene_id]["gene_end_position"]
            chromosome = dict_genes[gene_id]["chromosome"]

            #hashing = str(hashlib.md5(sequence).hexdigest())
            hashing = str(int(str(Decimal(time.time())).replace(".","")))[:16]
            #transid = str(dbsource + str(count) + "XX" + hashing).upper()
            transid = str(speciesshort + hashing + "SIH").upper()
            if len(tmplistmth) > 1:
                transid = str(speciesshort + hashing + "MTH").upper()
            c2.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(chromosome) + '''","''' + str(specie) + '''","''' + str(dbsource) + '''","''' + str(sequence) + '''")''')
            megabase += len(sequence)

    match = c.fetchone()
conn2.commit()
print "Species: " + specie + " Source: " + dbsource
print "Genes : " + str(len(list(OrderedDict.fromkeys(listgenid))))
print "Contigs: " + str(count) + "     contigs/gene: " + str(round(float(count) / len(list(OrderedDict.fromkeys(listgenid))), 2))
print "HSP: " + str(round(sum(listtailleshsp) / 1000000., 2)) + " megabases      avg size: " + str(sum(listtailleshsp) / len(listtailleshsp))
print "RAW: " + str(round(sum(listtaillesraw) / 1000000., 2)) + " megabases      avg size: " + str(sum(listtaillesraw) / len(listtaillesraw))
print
