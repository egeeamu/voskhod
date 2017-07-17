#!/bin/bash
set -e

cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ Voskhod CDNA merger

Version 20170331
Voskhod Pipeline version V1.1
Part of the Voskhod project
https://github.com/egeeamu/voskhod

GPL-3.0
Arnaud Ungaro contact@arnaud-ungaro.fr

Input must be in ./assembly/tomerge in sqlite file format.
Output will be in ./reference_ts in sqlite file format.

Type the whished name of the validated assembly (without space ex: drer_validated), followed by [ENTER]:
EndOfMessage

echo ""

read -n1 -r -p "Place all XX_cdna_infos.db you want merge in ./assembly/tomerge Press 'C' to continue..." key

echo ""

if [ "$key" = 'C' ]; then
	echo ""
	cp ./bin/cdna_merger.py ./assembly
	mkdir -p ./assembly/merged/
	echo "Type the name of the merged cdna_infos.db (without space), followed by [ENTER]:"
	read name
	cd ./assembly/
	python ./cdna_merger.py
	rm -rf ../reference_ts/$name.db
	mv ./tomerge/merged.db ../reference_ts/$name.db
	rm -rf ./cdna_merger.py
	

else
    echo "ABORT"
    exit
    # Anything else pressed, do whatever else.
    # echo [$key] not empty
fi
