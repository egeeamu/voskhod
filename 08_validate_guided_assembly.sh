#!/bin/bash
set -e

cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ Voskhod assembly validator

Version 20170331
Voskhod Pipeline version V1.1
Part of the Voskhod project
https://github.com/egeeamu/voskhod

GPL-3.0
Arnaud Ungaro contact@arnaud-ungaro.fr

Input must be in ./assembly/raw in fasta file format.
Output will be in ./assembly/validated in sqlite file format.

Type the whished name of the validated assembly (without space ex: drer_validated), followed by [ENTER]:
EndOfMessage
read name
echo ""


echo ""
echo ""

prompt="Please select Reference species:"
options=( $(find ./reference_ts -maxdepth 1 -type f -iregex '.*\.\(db\|db\)$' -print0 | xargs -0) )

PS3="$prompt "
select transref in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  "You picked $opt which is file $REPLY"
        break

    else
        echo "Invalid option. Try another one."
    fi
done    

echo ""
echo ""

prompt="Please select transcriptome to validate:"
options=( $(find ./assembly/raw/ -maxdepth 2 -type f -iregex '.*\.\(fq\|fastq\)$' -print0 | xargs -0) )

PS3="$prompt "
select transtovalid in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  "You picked $opt which is file $REPLY"
        break

    else
        echo "Invalid option. Try another one."
    fi
done    

echo ""
echo ""
rm -rf ./voskhod_part"$name"_validate
mkdir -p ./logs
mkdir -p ./voskhod_part"$name"_validate
mkdir -p ./voskhod_part"$name"_validate/validate_assembly
mkdir -p ./voskhod_part"$name"_validate/validate_assembly/data_input
mkdir -p ./voskhod_part"$name"_validate/validate_assembly/tomerge
mkdir -p ./assembly/validated
cp ./bin/voskhod_validate_assembly.py ./voskhod_part"$name"_validate/validate_assembly/
cp ./bin/voskload_preload_database_validate_contigs.py ./voskhod_part"$name"_validate/validate_assembly/
cp ./bin/Molnia_TD_to_CNI_PlusINFOS_Mth_assembly_annotation.py ./voskhod_part"$name"_validate/validate_assembly/
cp ./bin/voskmergecdna_DEDUP.py ./voskhod_part"$name"_validate/validate_assembly/


cd ./voskhod_part"$name"_validate/validate_assembly/data_input/
ln -s ../../."$transref" ./cdna_infos.db
ln -s ../../."$transtovalid"  ./"$name".fastq
cd ../../../

cd ./voskhod_part"$name"_validate/validate_assembly/
python ./voskload_preload_database_validate_contigs.py
python voskhod_validate_assembly.py ./data_input/"$name".fastq
mv  ./logs/* ../../logs

python Molnia_TD_to_CNI_PlusINFOS_Mth_assembly_annotation.py ./results/table_data_"$name".fastq.db "$name"
mv  *.txt ../../logs
mv cdna_infos_"$name"*.db ./tomerge
python voskmergecdna_DEDUP.py
rm ./tomerge/cdna_infos_"$name"*.db
rm -rf ../../assembly/validated/"$name"_validated.db
mv ./tomerge/merged.db ../../assembly/validated/"$name"_validated.db
exit
#mv  cdna_infos_"$name"*.db ../../assembly/validated


