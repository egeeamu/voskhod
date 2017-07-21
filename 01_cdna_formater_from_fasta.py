#! /usr/bin/env python
# -*- coding: UTF8 -*-

from Bio import SeqIO
import os
import sqlite3
from optparse import OptionParser
import sys


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

¤ Failsafe cdna formater from embl cdna

Caution : Annotations are not real but "fixed", fasta file do not contains all values ..

Version 20170721
Voskhod Pipeline version V1.2
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(GPL-3.0) 
Arnaud Ungaro contact@arnaud-ungaro.fr

This script allow you to work with non-standard dataset (not from ensembl) or problematic biomart dataset.

If you are working with ensembl (may work with plan,metazoa,etc..) dataset :

Get coding rna :
xx.xx.xx.cdna.all.fa.gz from :
ftp://ftp.ensembl.org/pub/release-XX/fasta/species_sp/cdna/

Get non-coding rna :
and  xx.xx.xx.ncrna.fa.gz from :
ftp://ftp.ensembl.org/pub/release-XX/fasta/species_sp/ncrna/

unpack (gzip -d) and put this in "failsafe_input"


If you are working with custom data-set, make sure your fasta is formated like this :

>transcriptname_1
sequence
>transcriptname_2
sequence

without space or special character in names, then put it in "failsafe_input".


Launch the sript with :

python 01b_failsafe_cdna_formater_from_fasta.py -s SpeciesName

"""


parser = OptionParser()
parser.add_option("-s", "--species-name", dest="spename", help="")
(options, args) = parser.parse_args()
species = filter(str.isalnum, str(options.spename))

conn = sqlite3.connect(':memory:')
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
c.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar)''')
conn.commit()

os.system("rm -rfv ./reference_ts/cdna_infos.db_" + species + ".db")
conn2 = sqlite3.connect('./reference_ts/cdna_infos.db_' + species + ".db")
c2 = conn2.cursor()
c2.execute('''PRAGMA synchronous = OFF''')
c2.execute('''PRAGMA journal_mode = OFF''')
c2.execute('''PRAGMA cache_size = 4096''')
c2.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar)''')
conn2.commit()

dbsource = "FailSafe"
listfiles = os.popen("""ls ./failsafe_input/*.f*a""").read().split("\n")[0:-1]
print listfiles
listtransid = []
count = 0
try:
    for i in listfiles:
        for record in SeqIO.parse(i, "fasta"):
            count += 1
            recsplit = record.description.split(" ")
            transidtmp = str(recsplit[0]).replace("-", "").replace(" ", "").replace("_", "").replace("+", "").replace(".", "").replace("/", "").replace("\\", "")

            if transidtmp not in listtransid:
                listtransid.append(transidtmp)
                transid = transidtmp
            else:
                transid = transidtmp + str(count)

            geneid = transid
            genename = geneid
            tbiotyp = ""

            for y in recsplit:
                findgene = y.find("gene:")
                if findgene != -1:
                    geneid = y
                    geneid = str(geneid.split("gene:")[-1]).replace("-", "").replace(" ", "").replace("_", "").replace("+", "").replace(".", "").replace("/", "_").replace("\\", "_")

            for y in recsplit:
                findgene = y.find("gene_symbol:")
                if findgene != -1:
                    genename = y
                    genename = str(genename.split("gene_symbol:")[-1]).replace("-", "-").replace(" ", "-").replace("_", "-").replace("+", "-").replace(".", "-").replace(":", "-").replace("/", "_").replace("\\", "_")

            for y in recsplit:
                findgene = y.find("transcript_biotype:")
                if findgene != -1:
                    tbiotyp = y
                    tbiotyp = str(tbiotyp.split("transcript_biotype:")[-1]).replace("-", "-").replace(" ", "-").replace("_", "-").replace("+", "-").replace(".", "-").replace(":", "-").replace("/", "_").replace("\\", "_")



            sequence = str(record.seq)


            # print line
            sequence = filter(str.isalnum, sequence)
            genename = filter(str.isalnum, genename)
            transid = filter(str.isalnum, transid)
            protid = ""
            namedesc = ""
            posdebgn = 100
            posfingn = len(sequence) + posdebgn
            posdebts = 100
            posfints = len(sequence) + posdebts
            transtype = tbiotyp
            geneid = geneid
            chromosome = "1"
            tssize = len(sequence)
            if genename == "":
                genename = geneid
            c.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(
                protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(
                transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(
                posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(
                chromosome) + '''","''' + str(species) + '''","''' + str(dbsource) + '''","''' + str(
                sequence) + '''")''')
        conn.commit()

except:
    print "ERROR"

conn.commit()

c.execute("""SELECT * FROM result ORDER BY gene_name,transcript_size""")
conn.commit()
match = c.fetchone()
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
    c2.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(chromosome) + '''","''' + str(specie) + '''","''' + str(dbsource) + '''","''' + str(sequence) + '''")''')
    match = c.fetchone()

conn2.commit()
