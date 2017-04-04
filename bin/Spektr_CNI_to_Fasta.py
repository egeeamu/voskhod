#! /usr/bin/env python
# -*- coding: UTF8 -*-

import sqlite3
import sys

count = 0

qq = str(sys.argv[1])
conn3 = sqlite3.connect(qq)
c3 = conn3.cursor()
c3.execute('''PRAGMA synchronous = OFF''')
c3.execute('''PRAGMA journal_mode = OFF''')
c3.execute('''PRAGMA cache_size = 4096''')
conn3.commit()

c3.execute("""SELECT gene_name, gene_id, transcript_sequence, source FROM result WHERE source != 'READ' ORDER BY gene_name,transcript_size""")
conn3.commit()
match = c3.fetchone()

while match is not None :
    count += 1
    #print match
    name = ">" + str(match[1]) + "_" + str(match[0]) + "_" + str(count) + "_" + str(match[3]) + "\n" + str(match[2]) + "\n"
    #print name
    with open("export.fasta", "a") as myfile:
        myfile.write(name)
    match = c3.fetchone()

