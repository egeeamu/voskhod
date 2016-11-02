#!/bin/bash
# Bash Menu Script Example


set -e
mkdir -p logs
mkdir -p ./cleaned_input
mkdir -p ./cleaned_input/assembly
cp ./bin/vosklean_sequences_cleaner_assembly_keep_pw.py ./raw_input/assembly
cp ./bin/afterclean_remove_bad_pw.py ./raw_input/assembly
cd ./raw_input/assembly
cat << EndOfMessage
 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

Version 20160920
¤ Cleaner
Part of the Voskhod project


(CC-BY-NC-ND 4.0 International license) 
Arnaud Ungaro contact@arnaud-ungaro.fr


¤ Option 1 :
  You have a single end library with only one fastq file
  This step will filter this file based on the quality and length

¤ Option 2 :
  You have a paired end library with R1 & R2 fastq files
  This step will filter the files then sync R1 and R2 into two filtered and synchronyzed fastq files

The input file(s) must be in ./raw_input/assembly and in fastq format.
The filtered files will be in ./cleaned_input/assembly in fastq format.

EndOfMessage
PS3='Please enter your choice: '
options=("Option 1" "Option 2" "Quit")
select opt in "${options[@]}"
do
    case $opt in
        "Option 1")
            argu=1
            break
            ;;
        "Option 2")
            argu=2
            break
            ;;
        "Quit")
            exit
            ;;
        *) echo invalid option;;
    esac
done


if [ "$argu" = "1" ]; then


echo ""
echo ""

prompt="Please select your file:"
options=( $(find ./ -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

PS3="$prompt "
select R1 in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  "You picked $opt which is file $REPLY"
        break

    else
        echo "Invalid option. Try another one."
    fi
done    

rm -rf ./"$R1"_query.fq
rm -rf ./cleantrinityR1.fq


echo ""
echo "Filtering in progress ..."
echo ""

time parallel --gnu --progress --max-procs 2 'nice -n 19 python vosklean_sequences_cleaner_assembly_keep_pw.py {}' ::: "$R1"
#time parallel --gnu --progress --max-procs 8 'nice -n 19 python vosklean_sequences_cleaner_assembly_keep_pw.py {}' ::: ./*.fastq

echo ""
echo "Filtering done"
echo ""

rm -rf ../../cleaned_input/assembly/"$R1"_cleaned_SE_sync.fastq
mv "$R1"_query.fq ../../cleaned_input/assembly/"$R1"_cleaned_SE_sync.fastq

rm -rf python afterclean_remove_bad_pw.py
mv ./logs/* ../../logs/
rm -rf ./logs
rm -rf ./python vosklean_sequences_cleaner_assembly_keep_pw.py


fi

if [ "$argu" = "2" ]; then

echo ""
echo ""

prompt="Please select R1 file (caution, files may be in any order in the list):"
options=( $(find ./ -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

PS3="$prompt "
select R1 in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  "You picked $opt which is file $REPLY"
        break

    else
        echo "Invalid option. Try another one."
    fi
done    
prompt="Please select R2 file (caution, files may be in any order in the list):"
options=( $(find ./ -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

PS3="$prompt "
select R2 in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  "You picked $opt which is file $REPLY"
        break

    else
        echo "Invalid option. Try another one."
    fi
done    

rm -rf ./"$R1"_query.fq
rm -rf ./"$R2"_query.fq
rm -rf ./cleantrinityR1.fq
rm -rf ./cleantrinityR2.fq

echo ""
echo "Filtering in progress ..."
echo ""

time parallel --gnu --progress --max-procs 2 'nice -n 19 python vosklean_sequences_cleaner_assembly_keep_pw.py {}' ::: "$R1" "$R2"
#time parallel --gnu --progress --max-procs 8 'nice -n 19 python vosklean_sequences_cleaner_assembly_keep_pw.py {}' ::: ./*.fastq

echo ""
echo "Filtering done, syncing in progress ..."
echo ""

python afterclean_remove_bad_pw.py -f "$R1"_query.fq -r "$R2"_query.fq | tee ../../logs/logs_syncs_R1_R2.txt

echo ""
echo "Syncing done !"
echo ""


rm -rf ../../cleaned_input/assembly/"$R1"_cleaned_R1_sync.fastq
rm -rf ../../cleaned_input/assembly/"$R2"_cleaned_R2_sync.fastq
mv cleantrinityR1.fq ../../cleaned_input/assembly/"$R1"_cleaned_R1_sync.fastq
mv cleantrinityR2.fq ../../cleaned_input/assembly/"$R2"_cleaned_R2_sync.fastq

#mv *_query.fq ../../cleaned_input/assembly/

rm -rf "$R1"_query.fq
rm -rf "$R2"_query.fq
rm -rf python afterclean_remove_bad_pw.py
mv ./logs/* ../../logs/
rm -rf ./logs
rm -rf ./python vosklean_sequences_cleaner_assembly_keep_pw.py


fi
