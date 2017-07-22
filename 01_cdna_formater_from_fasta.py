#! /usr/bin/env python
# -*- coding: UTF8 -*-

from Bio import SeqIO
import os
import sqlite3
from optparse import OptionParser
import mygene
import re
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

¤ Cdna formater from fasta file(s) cdna

Caution : Annotations may be not real but "fixed", like positions, chromosomes , etc.. if not found, fasta file do not contains all values ..


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
ftp://ftp.ensembl.org/pub/current_fasta/species_sp/cdna/

Get non-coding rna :
and  xx.xx.xx.ncrna.fa.gz from :
ftp://ftp.ensembl.org/pub/current_fasta/species_sp/ncrna/

unpack (gzip -d) and put this in "fasta_ref_input"


If you are working with custom data-set, make sure your fasta is formated like this :

>transcriptname_1
sequence
>transcriptname_2
sequence

without space or special character in names, then put it in "fasta_ref_input".


Launch the sript with :

01_cdna_formater_from_fasta.py -s SpeciesName

"""
# see http://mygene.info/doc/query_service.html#available-fields for choose scope
scopemg = "ensembl.transcript"


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

dbsource = "V1.2_20170721"
listfiles = os.popen("""ls ./fasta_ref_input/*.f*a""").read().split("\n")[0:-1]
print listfiles
listtransid = []
count = 0
listtoqmyge = []
dico_mginfos = {}


for i in listfiles:
    for record in SeqIO.parse(i, "fasta"):
        strtoadd = record.description.split(" ")[0].split(".")[0]
        strtoadd = re.sub('[^0-9a-zA-Z]+', '_', strtoadd)
        if strtoadd not in listtoqmyge:
            listtoqmyge.append(strtoadd)
        else:
            print ""
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print ""
            print "DUPLICATED TRANSCRIPT ID, ABORT"
            print i
            print str(strtoadd)
            print ""
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

            sys.exit()

mg = mygene.MyGeneInfo()
ginfo = mg.querymany(listtoqmyge, scopes=scopemg, species='all', fields='all')

for g in ginfo:
    try:
        gquery = re.sub('[^0-9a-zA-Z]+', '_', str(g["query"]))
    except:
        gquery = "NA"
    try:
        gene_id = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["ensemblgene"]))

    except:
        gene_id = gquery

    try:
        chromosome = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["chr"]))
    except:
        chromosome = "NA"

    try:
        posstart = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["start"]))
    except:
        posstart = "0"

    try:
        posend = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["end"]))
    except:
        posend = "0"

    try:
        transtype = re.sub('[^0-9a-zA-Z]+', '_', str(g["type_of_gene"]))
    except:
        transtype = "NA"

    try:
        gsymbol = re.sub('[^0-9a-zA-Z]+', '_', str(g["symbol"]))
    except:
        gsymbol = gquery

    try:
        gdesc = re.sub('[^0-9a-zA-Z]+', '_', str(g["name"]))
    except:
        gdesc = gsymbol

    dico_mginfos[gquery] = {"TransID": gquery, "GeneID": gene_id, "Gene_name": gsymbol, "Chromosome": chromosome,
                            "Posstart": posstart, "Posend": posend, "Transtype": transtype, "Gdesc": gdesc}

try:
    for i in listfiles:

        for record in SeqIO.parse(i, "fasta"):
            count += 1

            transid = record.description.split(" ")[0].split(".")[0]
            transid = re.sub('[^0-9a-zA-Z]+', '_', transid)

            gene_id = dico_mginfos[transid]["GeneID"]
            geneid = gene_id
            chromosome = dico_mginfos[transid]["Chromosome"]
            genename = dico_mginfos[transid]["Gene_name"]
            protid = ""
            namedesc = dico_mginfos[transid]["Gdesc"]
            posdebgn = dico_mginfos[transid]["Posstart"]
            posfingn = dico_mginfos[transid]["Posend"]
            posdebts = posdebgn
            posfints = posfingn
            transtype = dico_mginfos[transid]["Transtype"]
            sequence = str(record.seq)
            sequence = re.sub('[^0-9a-zA-Z]+', '_', sequence)
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

secondlistverifuniq = []

c.execute("""SELECT * FROM result ORDER BY gene_name,transcript_size""")
conn.commit()
match = c.fetchone()
while match is not None :
    transid = str(match[0])
    if transid not in secondlistverifuniq:
        secondlistverifuniq.append(transid)
    else:
        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print ""
        print "DUPLICATED TRANSCRIPT ID, ABORT"
        print ""
        print str(transid)
        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

        sys.exit()

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


