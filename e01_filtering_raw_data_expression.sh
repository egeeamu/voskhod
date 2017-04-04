#!/bin/bash
set -e


cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ Cleaner & merger expression

Version 20170331
Voskhod Pipeline version V1.1
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(CC-BY-NC-ND 4.0 International license) 
Arnaud Ungaro contact@arnaud-ungaro.fr

¤ Option 1)

Input : Single end library with only one fastq file by experiment
____________________

Output : One filtered fastq file 
____________________


##################################################

¤ Option 2)

Input : Paired end library with R1 & R2 fastq files with overlap exptected by experiment
____________________
             ____________________

Output : One filtered & merged (with Pear) fastq file
_________________________________


##################################################

¤ Option 3)

Input : Paired end library with R1 & R2 fastq files without overlap exptected by experiment
____________________
                                 ____________________

Output : One filtered & concatenated fastq file
____________________
____________________


This script will process all fastq files ./raw_input/expression without questions.

The input file(s)'s names must countain "_R1_" & "_R2_" for steps 2 & 3 
The input file(s) must be in ./raw_input/expression and in fastq format.
The cleaned files will be in ./cleaned_input/expression in fastq format.

Caution : all files in ./cleaned_input/expression will be erased automatically before process.

EndOfMessage




PS3='Please enter your choice: '
options=("Option 1" "Option 2" "Option 3" "Quit")
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
        "Option 3")
            argu=3
            break
            ;;
        "Quit")
            exit
            ;;
        *) echo invalid option;;
    esac
done

rm -rfv ./cleaned_input/expression
if [ "$argu" = "1" ]; then

cp ./bin/vosklean_sequences_cleaner_expression.py ./raw_input/expression/
mkdir -p ./logs
mkdir -p ./cleaned_input
mkdir -p ./cleaned_input/expression

cd ./raw_input/expression
time parallel --gnu --progress --max-procs 8 -u 'python vosklean_sequences_cleaner_expression.py {}' ::: *.f*q
mv -v ./logs/* ../../logs
rename  "s/query/cleaned/" *_query.f*q
mv -v *cleaned.fq ../../cleaned_input/expression
rm -rfv *_query.fq *cleaned.fq logs *.py

fi

if [ "$argu" = "2" ]; then



cp ./bin/automated_pear.py ./raw_input/expression/.
cp ./bin/pear/pear ./raw_input/expression/.
#cp ./bin/vosklean_sequences_cleaner_assembly_keep_pw.py ./raw_input/assembly
cp ./bin/vosklean_sequences_cleaner_expression.py ./raw_input/expression/
mkdir -p ./logs
mkdir -p ./cleaned_input
mkdir -p ./cleaned_input/expression

cd ./raw_input/expression
python automated_pear.py
rm -rfv *_PEAR_Merged_*.unassembled.*.f*q *_PEAR_Merged_*.discarded.f*q automated_pear.py pear
mv -v *_logs.txt ../../logs
time parallel --gnu --progress --max-procs 8 -u 'python vosklean_sequences_cleaner_expression.py {}' ::: *_PEAR_Merged_*.f*q.assembled*
mv -v ./logs/* ../../logs
mv -v *PEAR_Merged_*.assembled.f*q_query.fq ../../cleaned_input/expression
rm -rfv *PEAR_Merged_*.assembled.f*q logs *.py

fi


if [ "$argu" = "3" ]; then


cp ./bin/automated_cat.py ./raw_input/expression/.
cp ./bin/pear/pear ./raw_input/expression/.
#cp ./bin/vosklean_sequences_cleaner_assembly_keep_pw.py ./raw_input/assembly
cp ./bin/vosklean_sequences_cleaner_expression.py ./raw_input/expression/
mkdir -p ./logs
mkdir -p ./cleaned_input
mkdir -p ./cleaned_input/expression

cd ./raw_input/expression
python automated_cat.py
rm -rfv *_concatenated_*.unassembled.*.f*q *_concatenated_*.discarded.f*q automated_pear.py pear
#mv -v *_logs.txt ../../logs
time parallel --gnu --progress --max-procs 8 -u 'python vosklean_sequences_cleaner_expression.py {}' ::: *_concatenated_*.f*q
mv -v ./logs/* ../../logs
rename  "s/query/cleaned/" *_query.fq
mv -v *cleaned*.f*q ../../cleaned_input/expression
rm -rfv *concatenated_*.f*q logs *.py



fi
