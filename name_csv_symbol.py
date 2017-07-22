import csv
import mygene
import re
import sys

#qq = str(sys.argv[1])


ens = ['ENSDART00000093606', 'ENSDART00000093609', "ENSDART00000163904", "ENSDART00000164314"]


# outputtsv = list(csv.reader(open(qq, 'rb'), delimiter=';'))
#
# listgeneid = []
#
# for i in outputtsv:
#     try:
#         #print i[0]
#         gentmp = str(i[0])
#         if gentmp not in listgeneid:
#             listgeneid.append(gentmp)
#     except:
#         pass
#
# ens = listgeneid

ginfo = mg.querymany(ens, scopes='ensembl.transcript' ,species='all', fields='all')

for g in ginfo:
    try:
        gquery = re.sub('[^0-9a-zA-Z]+', '_', str(g["query"]))
    except:
        gquery = "NA"
    try:
        gene_id = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["ensemblgene"]))

    except:
        gene_id = "NA"

    try:
        chromosome = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["chr"]))
    except:
        chromosome = "NA"

    try:
        posstart = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["start"]))
    except:
        posstart = "NA"

    try:
        posend = re.sub('[^0-9a-zA-Z]+', '_', str(g["genomic_pos"]["end"]))
    except:
        posend = "NA"

    try :
        transtype = re.sub('[^0-9a-zA-Z]+', '_', str(g["type_of_gene"]))
    except:
        transtype = "NA"

    try:
        gdesc = re.sub('[^0-9a-zA-Z]+', '_', str(g["name"]))
    except:
        gdesc = gquery

    try:
        gsymbol = re.sub('[^0-9a-zA-Z]+', '_', str(g["symbol"]))
    except:
        gsymbol = gquery

    print gene_id
    print chromosome
    print posstart
    print posend
    print transtype
    print gdesc
    print gsymbol