#!/bin/bash
set -e



cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

Version 20160920
¤ Voskhod assembler
Part of the Voskhod project

(C) Arnaud Ungaro
contact@arnaud-ungaro.fr

Input must be in ./cleaned_input/assembly
Output will be in ./assembly/raw

Type the name assembled species (without space ex: drer_vosk), followed by [ENTER]:
EndOfMessage

read name
echo ""
echo ""

prompt="Please select cdna_infos.db (will be the blast database) (should be merged trinity(s) and MUST countain ref species!!)"
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

prompt="Please select input fastq file (single fastq, merged fastq, concatenated 'R1+R2') :"
options=( $(find ./cleaned_input/assembly -maxdepth 2 -type f -iregex '.*\.\(fq\|fastq\)$' -print0 | xargs -0) )

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

mkdir -p ./logs
rm -rfv ./voskhod_part_assembly_"$name"
mkdir -p ./assembly/raw/voshod/
mkdir -p ./voskhod_part_assembly_"$name"
mkdir -p ./voskhod_part_assembly_"$name"/ident_reads
mkdir -p ./voskhod_part_assembly_"$name"/ident_reads/data_input
mkdir -p ./voskhod_part_assembly_"$name"/ident_reads/bin
mkdir -p ./voskhod_part_assembly_"$name"/ident_reads/tomerge
mkdir -p ./assembly/validated

cp ./bin/voskhod_identreads.py ./voskhod_part_assembly_"$name"/ident_reads/
cp ./bin/cap3.linux.x86_64.tar.gz ./voskhod_part_assembly_"$name"/ident_reads/bin/
cp ./bin/SPAdes-3.6.2-Linux.tar.gz ./voskhod_part_assembly_"$name"/ident_reads/bin/
cp ./bin/limitfa.py ./voskhod_part_assembly_"$name"/ident_reads/
cp ./bin/limitfq.py ./voskhod_part_assembly_"$name"/ident_reads/
cp ./bin/vosklink_voskhod_assembler_realdata.py ./voskhod_part_assembly_"$name"/ident_reads/vosklink_voskhod_assembler.py
cp ./bin/voskmergecdna_DEDUP.py ./voskhod_part_assembly_"$name"/ident_reads/
cp ./bin/Spektr_CNI_to_Fasta.py ./voskhod_part_assembly_"$name"/ident_reads/
cp ./bin/convert_fasta_to_fastq.py ./voskhod_part_assembly_"$name"/ident_reads/


cp ./bin/voskhits_mth_share_assemb.py ./voskhod_part_assembly_"$name"/ident_reads/
cp ./bin/voskload_preload_database_validate_contigs.py ./voskhod_part_assembly_"$name"/ident_reads/
cp ./bin/Molnia_TD_to_CNI_PlusINFOS_Mth_assembly_annotation.py ./voskhod_part_assembly_"$name"/ident_reads/

#ln -rs "$transref" ./voskhod_part_assembly_"$name"/ident_reads/data_input/cdna_infos.db
cd ./voskhod_part_assembly_"$name"/ident_reads/data_input/
ln -s ../../."$transref" ./cdna_infos.db
ln -s ../../."$transtovalid"  ./"$name".fastq
cd ../../../
cd ./voskhod_part_assembly_"$name"/ident_reads/


mkdir list_species
sqlite3 ./data_input/cdna_infos.db "select species from result" | sort -u > listspe.txt
xargs -a listspe.txt -I name touch ./list_species/name

echo "#########################################################################"
echo "#########################################################################"

prompt="Please select the name of the reference species in this list :"
options=( $(find ./list_species/ -maxdepth 1 -type f -print0 | xargs -0) )

PS3="$prompt "
select spemodel in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  "You picked $opt which is file $REPLY"
        break

    else
        echo "Invalid option. Try another one."
    fi
done    

echo "#########################################################################"
echo "#########################################################################"

STR="$spemodel"  
IFS='/' read -ra NAMES <<< "$STR" 
refspecies="${NAMES[-1]}"


python ./voskload_preload_database_validate_contigs.py
python voskhod_identreads.py ./data_input/"$name".fastq
echo "$refspecies"
python voskhits_mth_share_assemb.py -i ./results/table_data_"$name".fastq.db -r "$refspecies"

python vosklink_voskhod_assembler.py -i ./results/hits/table_data_"$name".fastq.db_HITS.db -r "$refspecies" -n "$name"_Voskhod
mv ./data_input/table_data_"$name".fastq.db_HITS.db_cdna_infos.db ./tomerge
python voskmergecdna_DEDUP.py
rm ./tomerge/table_data_"$name".fastq.db_HITS.db_cdna_infos.db
python Spektr_CNI_to_Fasta.py ./tomerge/merged.db
python convert_fasta_to_fastq.py ./export.fasta export.fastq
mv ./export.fastq ../../assembly/raw/voshod/"$name"_Voskhod.fastq
mv ./export.fasta ../../assembly/raw/voshod/"$name"_Voskhod.fasta

#mv ./data_input/table_data_"$name".fastq.db_HITS.db_cdna_infos.db ../../assembly/raw/"$name"_Voskhod.db

mv ./logs/* ../../logs
#~ 
#~ python Molnia_TD_to_CNI_PlusINFOS_Mth_assembly_annotation.py ./results/table_data_"$name".fastq.db "$name"
#~ mv -v *.txt ../../logs
#~ mv -v cdna_infos_"$name"*.db ../../assembly/validated


