#! /usr/bin/env python
# -*- coding: UTF8 -*-

from optparse import OptionParser
import math
import datetime
import sys
from sys import stdout
import sqlite3
import os
from collections import *
from Bio import SeqIO
from collections import defaultdict
import hashlib
from bashplotlib.histogram import plot_hist



parser = OptionParser()
parser.add_option("-i", "--input", dest="inp", help="")
parser.add_option("-r", "--refspecies", dest="refsp", help="")
(options, args) = parser.parse_args()
inp = options.inp
refsp = options.refsp

#inp = "./results/table_data_Simu_mtco1_2016-02-29_Size_70_cov_5_overlap_8_Noise_0.0_R1.fq.db"

indiv = inp.split("/")[-1]

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


os.system("mkdir -p ./logs")
logfile = str("./logs/Voskhits_" + str(indiv) + '_' + str(datetime.datetime.now().strftime("%y-%m-%d_%Hh%Mm%Ss")) + "_log.txt")
sys.stdout = Logger(logfile)

# inp = "./results/toto3.db"
#inp = "./results/table_data_Simu_mtco1_2016-02-29_Size_70_cov_5_overlap_8_Noise_0.0_R1.fq.db"

quality_check = False
cope_merged = False
assemblymode = True
len_head_deduplication = 50
hspminsize = 30.
evaluemin = 0.05
identmin = 0.70
ratio_mininimum = 0.70
ratio_mini_ts = 0.70
#spe_model = "Danio_rerio"
spe_model = refsp

if assemblymode == True:
    deduplucate = True # use True for assembly only
    output_fasta = True
    quality_check = False
else:
    deduplucate = False
    output_fasta = False

if output_fasta == True:
    os.system("rm -rf ./fastq")
    os.system("mkdir ./fastq")


serv = os.popen("""hostname""").read().split("\n")
serv.pop(-1)
dat = os.popen("""date""").read().split("\n")
dat.pop(-1)
distrib = os.popen("""cat /etc/issue""").read().split("\n")
distrib.pop(-1)
kernel = os.popen("""uname -a""").read().split("\n")
kernel.pop(-1)
user = os.popen("""whoami""").read().split("\n")
user.pop(-1)
pwd = os.popen("""pwd""").read().split("\n")
pwd.pop(-1)
cpu = os.popen('''cat /proc/cpuinfo | grep "model name" | cut -d ":" -f2''').read().split("\n")
cpu.pop(-1)


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
 _    __           __   __    _ __
| |  / /___  _____/ /__/ /_  (_) /______
| | / / __ \/ ___/ //_/ __ \/ / __/ ___/
| |/ / /_/ (__  ) ,< / / / / / /_(__  )
|___/\____/____/_/|_/_/ /_/_/\__/____/

¤ Hits counter & accuracy tester

Version 20170721
Voskhod Pipeline version V1.2
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(GPL-3.0) 
Arnaud Ungaro contact@arnaud-ungaro.fr

"""


print "\n"
print "File name:       " + str(inp)
print "Dir:             " + str(pwd[0])
print "Server:          " + str(serv[0])
print "Distro:          " + str(distrib[0])
print "Kernel:          " + str(kernel[0])
print "CPU:             " + str(cpu[0])
print "Date:            " + str(dat[0])
print "User:            " + str(user[0])
print "\n"
print "Quality_check:   " + str(quality_check)
print "Cope_merged:     " + str(cope_merged)
print "Output_fasta:    " + str(output_fasta)
print "Hspminsize:      " + str(hspminsize)
print "Evaluemin:       " + str(evaluemin)
print "Identmin:        " + str(identmin)
print "Ratio_mininimum: " + str(ratio_mininimum)
print "Ratio_mini_ts:   " + str(ratio_mini_ts)
print "Ref-Species:     " + str(spe_model)
print "Deduplicate:     " + str(deduplucate)
print "Len_head_dedup:  " + str(len_head_deduplication)
print "\n"



def countorder(xs, top=100000000000000):
    counts = defaultdict(int)
    for x in xs:
        counts[x] += 1
    return sorted(counts.items(), reverse=True, key=lambda tup: tup[1])[:top]

conn = sqlite3.connect("./data_input/cdna_infos.db")
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
conn.commit()

dicoensdarg = OrderedDict()
c.execute("""SELECT gene_id,gene_name, chromosome, species, transcript_size,transcript_start_position, transcript_end_position , transcript_size, gene_start_position, gene_end_position, transcript_biotype FROM result WHERE species = '""" + str(spe_model) + """' ORDER BY gene_name""")
conn.commit()
match = c.fetchone()
while match is not None:
    ensdarg = str(match[0])
    gene_name = str(match[1])
    chromosome = str(match[2])
    sizets = int(match[4])
    posdeb = int(match[5])
    posfin = int(match[6])
    posdebgn = int(match[8])
    posfingn = int(match[9])
    biotype = str(match[10])
    try:
        dicoensdarg[ensdarg]["posdeb"] = posdebgn
        dicoensdarg[ensdarg]["posfin"] = posfingn
        sizetsmin = dicoensdarg[ensdarg]["sizetsmin"]
        sizetsmax = dicoensdarg[ensdarg]["sizetsmax"]
        dicoensdarg[ensdarg]["tssizelist"].append(sizets)
        if str(biotype) not in dicoensdarg[ensdarg]["transtype"]:
            dicoensdarg[ensdarg]["transtype"].append(str(biotype))
        if sizets > sizetsmax:
            dicoensdarg[ensdarg]["sizetsmax"] = sizets
        if sizets < sizetsmin:
            dicoensdarg[ensdarg]["sizetsmin"] = sizets
    except:
        dicoensdarg[ensdarg] = {"ensdarg": ensdarg}
        dicoensdarg[ensdarg]["gene_name"] = gene_name
        dicoensdarg[ensdarg]["chromosome"] = chromosome
        dicoensdarg[ensdarg]["posdeb"] = posdebgn
        dicoensdarg[ensdarg]["posfin"] = posfingn
        dicoensdarg[ensdarg]["sizetsmax"] = sizets
        dicoensdarg[ensdarg]["sizetsmin"] = sizets
        dicoensdarg[ensdarg]["transtype"] = [str(biotype)]
        dicoensdarg[ensdarg]["tssizelist"] = [sizets]
        dicoensdarg[ensdarg]["ParaNoise"] = 0
        dicoensdarg[ensdarg]["ParaNoiseList"] = []
    match = c.fetchone()

for iii in dicoensdarg.keys():
    try:
        dicoensdarg[iii]["sizegn"] = sum(dicoensdarg[iii]["tssizelist"]) / len(dicoensdarg[iii]["tssizelist"])
    except:
        pass

conn2 = sqlite3.connect(inp)
c2 = conn2.cursor()
c2.execute('''PRAGMA synchronous = OFF''')
c2.execute('''PRAGMA journal_mode = OFF''')
c2.execute('''PRAGMA cache_size = 4096''')
conn2.commit()

notfound = 0
resize_minsize = 0
count = 0.
countprint = 0.
dicospecies = OrderedDict()
c2.execute("""SELECT Count(*) FROM result""")
conn2.commit()
match = c2.fetchone()
nbrtotal = int(match[0])
conn2.commit()

c2.execute("""SELECT read_name,gene_id,species,evalue,hsp_size,identity,read_size,quality,read_sequence,multi_hits FROM result""")

conn2.commit()
match = c2.fetchone()
countok = 0
countbad = 0
countparanoise = 0
countduplicates = 0
countresizehsp = 0
file_handles = []
dicohash = OrderedDict()
dicohashheadtail = OrderedDict()
listident = []
listlenhsp = []
listlenread = []
listavgreadsize = []
listavghspsize = []
countignorednotdico = 0
countnotfound = 0
countshared = 0

while match is not None:
    triggerparanoise = False
    trigerhash = True
    count += 1
    countprintold = countprint
    countprint = int(math.floor(round((count / nbrtotal * 100.), 0)))
    if output_fasta == True:
        dict_qual = SeqIO.QualityIO._phred_to_sanger_quality_str
        quality = str(match[7])
        readseq = str(match[8])
        liste_quality = quality.split(" ")
        qual_fastaq = ""
        for i in liste_quality:
            qual_fastaq += dict_qual[int(i)]
        if len(readseq) != len(qual_fastaq):
            print "ERROR SIZE !="
    if countprint > countprintold:
        stdout.write("\r%s" % str(int(countprintold)) + "%  " + str(int(count)) + "/" + str(int(nbrtotal)) + " reads")
        stdout.flush()
    if quality_check == True:
        try:
            readname = str(match[0])
            # print readname.split("#")
            if cope_merged == False:
                endarggen = readname.split("#")[1]
                spegen = readname.split("#")[5]
            else:
                endarggen = readname.split("#")[2]
                spegen = readname.split("#")[6]
        except:
            endarggen = "notfound"
            pass
    ensdarg = str(match[1])
    hspif = hspminsize

    if ensdarg != "":
        try:
            sizeminits = dicoensdarg[ensdarg]["sizetsmin"]
            gname = dicoensdarg[ensdarg]["gene_name"]
        except:
            match = c2.fetchone()
            countignorednotdico += 1
            continue
        if hspif >= sizeminits:
            hspif = int((sizeminits * ratio_mini_ts))
            countresizehsp += 1
        gene_id = ensdarg
        resize_minsize += 1
        specie = str(match[2])
        evalue = float(match[3])
        hspsize = float(match[4])
        identity = float(match[5])
        readsize = float(match[6])
        ratiohsp = hspsize / readsize
        headtailvalue = (len_head_deduplication / 2)
        headseq = str(match[8])[0:len_head_deduplication]
        headseq2 = str(match[8])[0:headtailvalue]
        tailseq = str(match[8])[-headtailvalue:]
        hashid = hashlib.md5(headseq).hexdigest()
        headtail = headseq2 + tailseq
        hashidheadtail = hashlib.md5(headtail).hexdigest()

        # print headseq
        # print headseq2
        # print tailseq
        # print headtail
        # print hashid
        # print hashidheadtail
    if ensdarg == '':
        countnotfound += 1

    mastertriger = False
    if ensdarg != "" and evalue <= evaluemin and hspsize >= hspif and identity >= identmin and ratiohsp >= ratio_mininimum:
        mastertriger = True

    if deduplucate == True:
        if mastertriger == True:
            try:
                test = dicohash[gene_id][hashid]["hits"]
                trigerhash = False
                countduplicates += 1
            except:
                trigerhash = True
    else:
        trigerhash = True

    try:
        if mastertriger == True:

            multivalue = str(match[9])
            multivalue.split("|")
            listmultispe = []
            lismultigene = []
            for i in multivalue.split("|"):
                multielements = i.split("+")
                speciesmulti = multielements[0]
                transidmulti = multielements[1]
                geneidsmulti = multielements[2]

                test = dicoensdarg[geneidsmulti]["gene_name"]
                if geneidsmulti not in lismultigene:
                    lismultigene.append(geneidsmulti)

                if speciesmulti not in listmultispe:
                    listmultispe.append(speciesmulti)


            if len(listmultispe) > 1:
                specie = "MultiSpecies"

            if len(lismultigene) > 1:
                triggerparanoise = True
                # if assemblymode == True:
                #     triggerparanoise = False
                #     specie = "ParaNoise"
                #     #print "bob"
                #     countparanoise -= 1
                countparanoise += 1
                for iii in lismultigene:
                    #print iii
                    if assemblymode == True:
                        countshared += 1
                        genametmp = dicoensdarg[iii]["gene_name"]
                        filequal = str("./fastq/" + str(genametmp) + "_" + str(iii) + ".fastq")
                        qual = str("""@""" + str(genametmp) + "_" + str(iii) + "_" + str(int(count)) + "_" + "ParaNoiseShare" + "\n" + str(readseq) + "\n+\n" + str(qual_fastaq) + "\n")
                        filework = open(filequal, "a+b")
                        #file_handles.append(filework)
                        filework.write(qual)
                        filework.close()
                    try:
                        dicoensdarg[iii]["ParaNoise"] += 1
                        for zzz in lismultigene:
                            #if zzz not in dicoensdarg[iii]["ParaNoiseList"] and zzz != iii:
                            if zzz != iii:
                                dicoensdarg[iii]["ParaNoiseList"].append(zzz)
                    except:
                        dicoensdarg[iii]["ParaNoise"] = 1
    except:
        countignorednotdico += 1
        #print "PROBLEM MULTISPECIES PARSING"
        match = c2.fetchone()
        continue

    if mastertriger == True and trigerhash == True and triggerparanoise == False:

        if deduplucate == False:
            try :
                test = dicohash[gene_id][hashid]["hits"]
                countduplicates += 1
            except:
                pass

        listident.append(round(identity,2))
        listlenhsp.append(int(hspsize))
        listlenread.append((readsize))
        try:
            test = dicohash[gene_id]
        except:
            dicohash[gene_id] = {}

        try:
            test = dicohash[gene_id][hashid]
        except:
            dicohash[gene_id][hashid] = {}

        try:
            dicohash[gene_id][hashid]["hits"] += 1
        except:
            dicohash[gene_id][hashid]["hits"] = 1

        # try:
        #     dicohash[gene_id][hashid][specie] += 1
        # except:
        #     dicohash[gene_id][hashid][specie] = 1





        try:
            test = dicohashheadtail[gene_id]
        except:
            dicohashheadtail[gene_id] = {}


        try:
            test = dicohashheadtail[gene_id][hashidheadtail]
        except:
            dicohashheadtail[gene_id][hashidheadtail] = {}

        try:
            dicohashheadtail[gene_id][hashidheadtail]["hits"] += 1
        except:
            dicohashheadtail[gene_id][hashidheadtail]["hits"] = 1

        # try:
        #     dicohashheadtail[gene_id][hashidheadtail][specie] += 1
        # except:
        #     dicohashheadtail[gene_id][hashidheadtail][specie] = 1






        countok += 1

        # if countok > 1000000:
        #     break

        if output_fasta == True:
            filequal = str("./fastq/" + str(gname) + "_" + str(ensdarg) + ".fastq")
            qual = str("""@""" + str(gname) + "_" + str(ensdarg) + "_" + str(int(count)) + "_" + specie + "\n" + str(readseq) + "\n+\n" + str(qual_fastaq) + "\n")
            filework = open(filequal, "a+b")
            #file_handles.append(filework)
            filework.write(qual)
            filework.close()
        if quality_check == True:
            try:
                dicoensdarg[endarggen]["ReadsGen"] += 1
            except:
                dicoensdarg[endarggen]["ReadsGen"] = 1
            try:
                dicoensdarg[endarggen]["ReadsGen_spe_" + str(spegen)] += 1
            except:
                dicoensdarg[endarggen]["ReadsGen_spe_" + str(spegen)] = 1
            if ensdarg == endarggen and specie == spegen:
                try:
                    dicoensdarg[ensdarg]["HitsGood"] += 1
                except:
                    dicoensdarg[ensdarg]["HitsGood"] = 1
                try:
                    dicoensdarg[ensdarg]["HitsGood_spe_" + str(spegen)] += 1
                except:
                    dicoensdarg[ensdarg]["HitsGood_spe_" + str(spegen)] = 1
            else:
                if endarggen == ensdarg:
                    try:
                        dicoensdarg[ensdarg]["HitsGood"] += 1
                    except:
                        dicoensdarg[ensdarg]["HitsGood"] = 1
                else:
                    try:
                        dicoensdarg[ensdarg]["HitsBad"] += 1
                    except:
                        dicoensdarg[ensdarg]["HitsBad"] = 1
                        pass
                    try:
                        dicoensdarg[ensdarg]["HitsBad_spe_" + str(spegen)] += 1
                    except:
                        dicoensdarg[ensdarg]["HitsBad_spe_" + str(spegen)] = 1
        try:
            dicoensdarg[ensdarg]["Hits"] += 1
        except:
            dicoensdarg[ensdarg]["Hits"] = 1
        try:
            dicoensdarg[ensdarg]["SumReadSize"] += readsize
        except:
            dicoensdarg[ensdarg]["SumReadSize"] = readsize
        try:
            dicoensdarg[ensdarg]["SumHSPSize"] += hspsize
        except:
            dicoensdarg[ensdarg]["SumHSPSize"] = hspsize
        try:
            dicoensdarg[ensdarg]["Hits_" + str(specie)] += 1
        except:
            dicoensdarg[ensdarg]["Hits_" + str(specie)] = 1

    else:
        notfound += 1
        countbad += 1
        if quality_check == True:
            try:
                dicoensdarg[endarggen]["ReadsGen"] += 1
            except:
                dicoensdarg[endarggen]["ReadsGen"] = 1
            try:
                dicoensdarg[endarggen]["ReadsGen_spe_" + str(spegen)] += 1
            except:
                dicoensdarg[endarggen]["ReadsGen_spe_" + str(spegen)] = 1
        match = c2.fetchone()
        continue

    try:
        dicospecies[specie]["TotalHitsSpeRaw"] += 1
    except:
        dicospecies[specie] = {"specie": specie}
        dicospecies[specie]["TotalHitsSpeRaw"] = 1

    match = c2.fetchone()
stdout.write("\r%s" % str(100) + "%  " + str(int(nbrtotal)) + "/" + str(int(nbrtotal)) + " reads")
stdout.flush()
print "\n"

for fic in file_handles:
    fic.close()

listspecies = []
for i in dicospecies:
    listspecies.append(str(dicospecies[i]["specie"]))
listspecies = sorted(listspecies)
#print listspecies

strspe = ""
listfly = ["HitsGood", "HitsBad", "ReadsGen", "SumHSPSize", "SumReadSize", "Hits"]
for iii in listspecies:
    replaced = str(iii).replace(".", "_")
    strspe = strspe + ",Hits_" + replaced + " varchar"
    listfly.append(str("Hits_" + replaced))
#print strspe

if quality_check == True:
    strspeqc = ",HitsGood int,HitsBad int,ReadsGen int"
    for iii in listspecies:
        replaced = str(iii).replace(".", "_")
        strspeqc =  strspeqc + ",HitsGood_spe_" + replaced + " varchar"
        listfly.append(str("HitsGood_spe_" + replaced))
        strspeqc =  strspeqc + ",HitsBad_spe_" + replaced + " varchar"
        listfly.append(str("HitsBad_spe_" + replaced))
        strspeqc =  strspeqc + ",ReadsGen_spe_" + replaced + " varchar"
        listfly.append(str("ReadsGen_spe_" + replaced))
    #print strspeqc


for iii in dicoensdarg.keys():
    #print iii
    #print dicoensdarg[iii]
    for i in listfly:
        try:
            test = dicoensdarg[iii][i]
        except:
            dicoensdarg[iii][i] = 0
            test = dicoensdarg[iii][i]
        #print str(i) + " " + str(test)
    try:
        if int(dicoensdarg[iii]["Hits"]) >= 1:
            dicoensdarg[iii]["AvgRead"] = round(float(dicoensdarg[iii]["SumReadSize"]) / float(dicoensdarg[iii]["Hits"]),3)
        else:
            dicoensdarg[iii]["AvgRead"] = 0
    except:
        dicoensdarg[iii]["AvgRead"] = 0

    try:
        if int(dicoensdarg[iii]["Hits"]) >= 1:
            dicoensdarg[iii]["AvgHSP"] = round(float(dicoensdarg[iii]["SumHSPSize"]) / float(dicoensdarg[iii]["Hits"]),3)
        else:
            dicoensdarg[iii]["AvgHSP"] = 0
    except:
        dicoensdarg[iii]["AvgRead"] = 0


    try:
        if int(dicoensdarg[iii]["Hits"]) >= 1:
            dicoensdarg[iii]["AvgCov"] = round(float(dicoensdarg[iii]["SumReadSize"]) / float(dicoensdarg[iii]["sizegn"]),3)
        else:
            dicoensdarg[iii]["AvgCov"] = 0
    except:
        dicoensdarg[iii]["AvgCov"] = 0


    #print dicoensdarg[iii]




dbstr = "CREATE TABLE result (gene_id varchar,gene_name varchar, transcript_biotype varchar ,gene_size varchar, transcript_size_Min varchar, transcript_size_Max varchar, gene_start_position varchar, gene_end_position varchar, chromosome varchar,AvgHSP float, AvgRead float, AvgCov float, Total_Hits int, Non_uniq_Sequences_Head int, Duplicates_Head_Distrib varchar, Non_uniq_Sequences_Head_Tail int, Duplicates_Head_Tail_Distrib varchar, ParaNoise int, ParaNoiseList varchar"
dbstr = dbstr + strspe
if quality_check == True:
    dbstr = dbstr + strspeqc
dbstr = dbstr + ")"
#print dbstr
savedb = "./results/hits/" + str(inp.split("/")[-1]) + "_HITS.db"
os.system("mkdir -p ./results/hits/")
#savedb = "./results/test_hits.db"
os.system("rm -rfv " + savedb)
conn3 = sqlite3.connect(savedb)
c3 = conn3.cursor()
c3.execute('''PRAGMA synchronous = OFF''')
c3.execute('''PRAGMA journal_mode = OFF''')
c3.execute('''PRAGMA cache_size = 4096''')
print dbstr
c3.execute(dbstr)
conn3.commit()


# iii = "ENSDARG00000063924"
# print dicoensdarg[iii]

sommedupli = 0
sommedupliht = 0
listcov = []
listheaddup = []
listheadtaildup = []
listhitsgenes = []
listhitsparanoise = []

for iii in dicoensdarg.keys():
    listhash = []

    dupli = 0
    dupliht = 0
    try:
        #paranoiselistegenes = str("|".join(dicoensdarg[iii]["ParaNoiseList"]))
        paranoiselistegenes = str(countorder(dicoensdarg[iii]["ParaNoiseList"]))
    except:
        print "parabug1"
    try:
        for y in dicohash[iii]:
            listhash.append(int(dicohash[iii][y]["hits"]))
        distribduplicates = str(countorder(listhash))
        dupli = int(dicoensdarg[iii]["Hits"]) - listhash.count(1)
        sommedupli += dupli
    except:
        distribduplicates = str("")
    listhash = []

    try:
        for y in dicohashheadtail[iii]:
            listhash.append(int(dicohashheadtail[iii][y]["hits"]))
        distribduplicatesht = str(countorder(listhash))
        dupliht = int(dicoensdarg[iii]["Hits"]) - listhash.count(1)
        sommedupliht += dupliht
    except:
        distribduplicatesht = str("")


    try:
        paranoise = dicoensdarg[iii]["ParaNoise"]
    except:
        paranoise = 0
        #print "error para"

    try:
        g = '"'
        gv = '","'
        dbinsert = """INSERT INTO result VALUES (""" + g + iii + gv + str(dicoensdarg[iii]["gene_name"]) + gv + str("|".join(dicoensdarg[iii]["transtype"])) + gv + str(dicoensdarg[iii]["sizegn"]) + gv + str(dicoensdarg[iii]["sizetsmin"]) + gv + str(dicoensdarg[iii]["sizetsmax"])
        dbinsert = dbinsert + gv + str(dicoensdarg[iii]["posdeb"]) + gv + str(dicoensdarg[iii]["posfin"]) + gv + str(dicoensdarg[iii]["chromosome"]) + gv + str(dicoensdarg[iii]["AvgHSP"])
        dbinsert = dbinsert + gv + str(dicoensdarg[iii]["AvgRead"]) + gv + str(dicoensdarg[iii]["AvgCov"]) + gv + str(dicoensdarg[iii]["Hits"]) + gv + str(dupli) + gv + str(distribduplicates) + gv + str(dupliht) + gv + str(distribduplicatesht) + gv + str(paranoise) + gv + paranoiselistegenes

        if dicoensdarg[iii]["Hits"] != 0:
            listcov.append(round(dicoensdarg[iii]["AvgCov"],3))
            listheaddup.append(int(dupli))
            listheadtaildup.append(int(dupliht))
            listhitsgenes.append(int(dicoensdarg[iii]["Hits"]))
            listavghspsize.append(float(dicoensdarg[iii]["AvgHSP"]))
            listavgreadsize.append(float(dicoensdarg[iii]["AvgRead"]))
            listhitsparanoise.append(float(dicoensdarg[iii]["ParaNoise"]))

        stradd = ""
        for yyy in listspecies:
            replaced = str(yyy).replace(".", "_")
            stradd = stradd + gv + str(dicoensdarg[iii]["Hits_" + replaced])
        dbinsert = dbinsert + stradd

        if quality_check == True:
            dbinsert = dbinsert + gv + str(dicoensdarg[iii]["HitsGood"]) + gv + str(dicoensdarg[iii]["HitsBad"]) + gv + str(dicoensdarg[iii]["ReadsGen"])

        if quality_check == True:
            stradd2 = ""
            for yyy in listspecies:
                replaced = str(yyy).replace(".", "_")
                stradd2 = stradd2 + gv + str(dicoensdarg[iii]["HitsGood_spe_" + replaced])
                stradd2 = stradd2 + gv + str(dicoensdarg[iii]["HitsBad_spe_" + replaced])
                stradd2 = stradd2 + gv + str(dicoensdarg[iii]["ReadsGen_spe_" + replaced])
            dbinsert = dbinsert + stradd2
        dbinsert = dbinsert + g + ")"
        c3.execute(dbinsert)

    except:
        print "Error " + str(iii)
conn3.commit()
#print dicoensdarg["ENSDARG00000070003"]

if deduplucate == True:
    rejectedprint = (countbad - countnotfound - countparanoise - countduplicates)
else:
    rejectedprint = (countbad - countnotfound - countparanoise)

print "Total                            : " + str(nbrtotal)
print "Total OK                         : " + str(countok)
print "Not identified                   : " + str(countnotfound)
print "Rejected  (filters)              : " + str(rejectedprint)
print "Ignored (ParaNoise)              : " + str(countparanoise)
print "Count Shared Reads (Ass)         : " + str(countshared)
print "Ignored (Not in ref)             : " + str(countignorednotdico)
print "Deduplicates Head (Assembly mode): " + str(countduplicates)
print "Non uniq sequences (head)        : " + str(sommedupli) + " uniqs : " + str(countok - sommedupli)
print "Non uniq sequences(head+tail)    : " + str(sommedupliht) + " uniqs : " + str(countok - sommedupliht)
print "Hsp resize                       : " + str(countresizehsp)

print "Done.."

# for i in dicohash:
#     print i,dicohash[i]
# print "\n"
# for i in dicohash:
#     for y in dicohash[i]:
#         print i,y,dicohash[i][y]["hits"]

c3.execute('''CREATE TABLE logs (logs varchar)''')
conn3.commit()
#c3.execute('''INSERT INTO result VALUES ("''' + logsall + '''")''')
#print logsall.split("\n")

#print logsall








try:
    plot_hist(listlenhsp,showSummary=True,title="Reads_Size",xlab=True,bincount=100,pch="*")
    plot_hist(listlenread,showSummary=True,title="HSPs_Size",xlab=True,bincount=100,pch="*")
    plot_hist(listavgreadsize,showSummary=True,title="Reads_Size_By_Gene",xlab=True,bincount=100,pch="*")
    plot_hist(listavghspsize,showSummary=True,title="HSPs_Size_By_Gene",xlab=True,bincount=100,pch="*")
    plot_hist(listident,showSummary=True,title="Identity",xlab=True,bincount=100,pch="*")
    plot_hist(listhitsgenes,showSummary=True,title="Hits_By_Gene",xlab=True,bincount=100,pch="*")
    plot_hist(listcov,showSummary=True,title="Coverage_By_Gene",xlab=True,bincount=100,pch="*")
    plot_hist(listheaddup,showSummary=True,title="Head_Duplication_By_Gene",xlab=True,bincount=100,pch="*")
    plot_hist(listheadtaildup,showSummary=True,title="Head_&_Tail_Duplication_By_Gene",xlab=True,bincount=100,pch="*")
    plot_hist(listhitsparanoise,showSummary=True,title="ParaNoise_By_Gene",xlab=True,bincount=100,pch="*")
except:
    pass

try:
    for i in logsall.split("\n"):
        c3.execute('''INSERT INTO logs VALUES ("''' + str(i) + '''")''')

    conn3.commit()

except:
    pass
