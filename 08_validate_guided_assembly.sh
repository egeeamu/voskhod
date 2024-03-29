#!/bin/bash
set -e

cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ Voskhod assembly validator

Version 20170809
Voskhod Pipeline version V1.2.1
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
cp ./bin/Spektr_CNI_to_Fasta.py ./voskhod_part"$name"_validate/validate_assembly/
cp ./bin/convert_fasta_to_fastq.py ./voskhod_part"$name"_validate/validate_assembly/
cp ./bin/Demeter_Biggest_Noise.py ./voskhod_part"$name"_validate/validate_assembly/




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

#deduplication
python voskmergecdna_DEDUP.py
rm ./tomerge/cdna_infos_"$name"*.db

#convert sqlite to fasta & fastq
python Spektr_CNI_to_Fasta.py ./tomerge/merged.db
python convert_fasta_to_fastq.py ./export.fasta export.fastq
mv ./export.fastq ../../assembly/validated/"$name"_validated_guided.fastq
mv ./export.fasta ../../assembly/validated/"$name"_validated_guided.fasta


rm -rf ../../assembly/validated/"$name"_validated.db
cp ./tomerge/merged.db ../../assembly/validated/"$name"_validated_guided.db
mv ./tomerge ./tomergededup
mkdir -p ./tomerge
mv ./tomergededup/merged.db ./tomerge/tofindbiggest.db

rm -rfv ./export.fastq
rm -rfv ./export.fasta

# biggest contig by genes
python Demeter_Biggest_Noise.py
cp ./tomerge/merged.db ../../assembly/validated/"$name"_validated_guided_biggest_contig_by_gene.db

python Spektr_CNI_to_Fasta.py ./tomerge/merged.db
python convert_fasta_to_fastq.py ./export.fasta export.fastq

mv -f ./export.fastq ../../assembly/validated/"$name"_validated_guided_biggest_contig_by_gene.fastq
mv -f ./export.fasta ../../assembly/validated/"$name"_validated_guided_biggest_contig_by_gene.fasta



exit
#mv  cdna_infos_"$name"*.db ../../assembly/validated
