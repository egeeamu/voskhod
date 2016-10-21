#!/bin/bash
set -e



cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

Version 20160920
¤ Voskhod Expression
Part of the Voskhod project

(C) Arnaud Ungaro
contact@arnaud-ungaro.fr

Input must be in ./cleaned_input/expression in fastq file format.
Output will be in ./expression_results/ in sqlite and csv format.

Caution : all files in ./expression_results/"STUDYNAME" will be erased automatically before process.

Type the name of the study (without space ex: drer_liver_exp1), followed by [ENTER]:
EndOfMessage
read name



echo ""
echo ""

prompt="Please select cdna database (all species(final) merged + ref):"
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

#~ prompt="Please select transcriptom to validate:"
#~ options=( $(find ./assembly/raw/ -maxdepth 2 -type f -iregex '.*\.\(fq\|fastq\)$' -print0 | xargs -0) )
#~ 
#~ PS3="$prompt "
#~ select transtovalid in "${options[@]}" "Quit" ; do 
    #~ if (( REPLY == 1 + ${#options[@]} )) ; then
        #~ exit
#~ 
    #~ elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        #~ echo  "You picked $opt which is file $REPLY"
        #~ break
#~ 
    #~ else
        #~ echo "Invalid option. Try another one."
    #~ fi
#~ done    

echo ""
echo ""
rm -rfv ./voskhod_part_"$name"_expression
mkdir -p ./expression_results/
rm -rfv ./expression_results/"$name"
mkdir -p ./expression_results/"$name"/db
mkdir -p ./expression_results/"$name"/csv
mkdir -p ./logs
mkdir -p ./voskhod_part_"$name"_expression
mkdir -p ./voskhod_part_"$name"_expression/espression_hits
mkdir -p ./voskhod_part_"$name"_expression/espression_hits/data_input
mkdir -p ./voskhod_part_"$name"_expression/espression_hits/tomerge
mkdir -p ./assembly/validated
cp ./bin/voskhod_expression.py ./voskhod_part_"$name"_expression/espression_hits/
cp ./bin/voskhits_mth_share_expression.py ./voskhod_part_"$name"_expression/espression_hits/
cp ./bin/voskload_preload_database_validate_contigs.py ./voskhod_part_"$name"_expression/espression_hits/




cd ./voskhod_part_"$name"_expression/espression_hits/data_input/
ln -s ../../."$transref" ./cdna_infos.db
cd ../../../




cd ./voskhod_part_"$name"_expression/espression_hits/


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
echo ""
echo "Work in progress, may take long time, do not exit.."
echo ""
time parallel --ungroup --gnu --progress --max-procs 1 'time python voskhod_expression.py {}' ::: ../../cleaned_input/expression/*.f*q
time parallel --ungroup --gnu --progress --max-procs 1 "time python voskhits_mth_share_expression.py -r '$refspecies' -i {}" ::: ./results/*.db
mv -v ./logs/* ../../logs
cd ./results/hits
parallel --ungroup --gnu --progress --max-procs 8 'sqlite3 -header -csv {} "select * from result;" > {}.csv' ::: ./*.db
mv *.db ../../../../expression_results/"$name"/db
mv *.csv ../../../../expression_results/"$name"/csv

echo ""
echo ""
echo ""

echo "DONE !!"




