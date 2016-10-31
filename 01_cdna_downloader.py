#! /usr/bin/python
# -*- coding: UTF8 -*-


import sqlite3
import os
import string
from biomart import BiomartDataset
import math
from sys import stdout

specie = "drerio_gene_ensembl" # dyakuba_eg_gene , drerio_gene_ensembl, dmelanogaster_eg_gene, pformosa_gene_ensembl, gmorhua_gene_ensembl, trubripes_gene_ensembl
dbsource = "www_ensembl_org" # metazoa_ensembl_org or www_ensembl_org or www.biomart.org


if dbsource == "www_ensembl_org":
    dbinput = "http://www.ensembl.org/biomart"
# elif dbsource == "metazoa_ensembl_org":
#     dbinput = "http://metazoa.ensembl.org/biomart/"
else:
    dbinput = "http://www.biomart.org/biomart"

countts = 0
count = 0.
countprint = 0.
tr = string.maketrans("'()*\/:|, ", "---------_")
os.system("mkdir ./reference_ts >/dev/null 2>/dev/null")
#os.system("rm -rf ./reference_ts/cdna_infos.db.tmp")
conn = sqlite3.connect(':memory:')
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
c.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar)''')
conn.commit()

os.system("rm -rf ./reference_ts/cdna_infos_" + specie + ".db")
conn2 = sqlite3.connect("./reference_ts/cdna_infos_" + specie + ".db")
c2 = conn2.cursor()
c2.execute('''PRAGMA synchronous = OFF''')
c2.execute('''PRAGMA journal_mode = OFF''')
c2.execute('''PRAGMA cache_size = 4096''')
c2.execute('''CREATE TABLE result (transcript_id varchar, gene_id varchar, peptide_id , gene_name varchar,gene_description varchar,transcript_biotype varchar, transcript_size int,  transcript_start_position int, transcript_end_position int, gene_start_position int, gene_end_position int, chromosome varchar, species varchar, source varchar, transcript_sequence varchar)''')
conn2.commit()

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

Version 20160911
Automatic cdna_downloader & formater
Part of the Voskhod project

(C) Arnaud Ungaro
contact@arnaud-ungaro.fr"
"""

print "Connecting to " + dbinput + ".."

martquery = BiomartDataset( dbinput, name = specie )
#martquery = BiomartDataset( "http://www.biomart.org/biomart", name = specie )
print "fetching transcripts count.."

response = martquery.search({'filters': {},'attributes': [ 'ensembl_transcript_id']})

for line in response.iter_lines():
    countts += 1
    line.split("\t")
print "Done, " + str(countts) + " trancripts found for " + specie
print "Retrieving data and filling temporary database.."
if dbsource == "www_ensembl_org":
    response = martquery.search({'filters': {},'attributes': [ 'cdna', 'external_gene_name', 'ensembl_transcript_id', 'ensembl_peptide_id', 'description', 'start_position', 'end_position', 'transcript_start', 'transcript_end', 'transcript_biotype', 'ensembl_gene_id', 'chromosome_name']})
# if dbsource == "metazoa_ensembl_org":
else:
    response = martquery.search({'filters': {},'attributes': [ 'cdna', 'external_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'description', 'start_position', 'end_position', 'transcript_start', 'transcript_end', 'transcript_biotype', 'ensembl_gene_id', 'chromosome_name']})

for line in response.iter_lines():
    #print line
    count += 1
    countprintold = countprint
    countprint = int(math.floor(round((count / countts * 100.), 0)))
    #if countprint > countprintold:
    stdout.write("\r%s" % str(int(countprintold)) + "%  " + str(int(count)) + "/" + str(int(countts)) + " transcripts")
    stdout.flush()
    sequence = line.split("\t")[0]
    genename = str(line.split("\t")[1].split(" (")[0]).translate(tr)
    transid = line.split("\t")[2]
    protid = line.split("\t")[3]
    namedesc = str(line.split("\t")[4].split(" [")[0]).translate(tr)
    posdebgn = line.split("\t")[5]
    posfingn = line.split("\t")[6]
    posdebts = line.split("\t")[7]
    posfints = line.split("\t")[8]
    transtype = line.split("\t")[9]
    geneid = line.split("\t")[10]
    chromosome = str(line.split("\t")[11]).translate(tr)
    tssize = len(sequence)
    if genename == "":
        genename = geneid
    c.execute('''INSERT INTO result VALUES ("''' + str(transid) + '''","''' + str(geneid) + '''","''' + str(protid) + '''","''' + str(genename) + '''","''' + str(namedesc) + '''","''' + str(transtype) + '''","''' + str(tssize) + '''","''' + str(posdebts) + '''","''' + str(posfints) + '''","''' + str(posdebgn) + '''","''' + str(posfingn) + '''","''' + str(chromosome) + '''","''' + str(specie) + '''","''' + str(dbsource) + '''","''' + str(sequence) + '''")''')
conn.commit()
stdout.write("\r%s" % str(100) + "%  " + str(int(countts)) + "/" + str(int(countts)) + " transcripts")
stdout.flush()
print "\n"
print "Done.."
print "Sorting transcript by Gene_name and filling final database.."
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
os.system("rm -rf ./reference_ts/cdna_infos.db.tmp")
conn2.commit()
print "Done..!"
