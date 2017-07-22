#!/bin/bash

set -e 


cat << EndOfMessage
 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ De novo assembly / Trinity wrapper

Version 20170721
Voskhod Pipeline version V1.2
Part of the Voskhod project
https://github.com/egeeamu/voskhod

GPL-3.0
Arnaud Ungaro contact@arnaud-ungaro.fr

If you have an unstranded library :

> Choose Option 1 if:
  You have a single-end unstranded library with only one fastq file
  This step will launch Trinity with relevant parameters for single input file and an unstranded library.

> Choose Option 2 if:
  You have a paired-end unstranded library with R1 & R2 fastq files
  This step will launch Trinity with relevant parameters for input files and an unstranded library.
  
If you have a stranded library :

> Choose Option 3 if:
  You have a single-end stranded library with only one fastq file
  This step will launch Trinity with relevant parameters for a single input file and a stranded library.

> Choose Option 4 if:
  You have a paired-end stranded library with R1 & R2 fastq files
  This step will launch Trinity with relevant parameters for input files and a stranded library.


There are four stranded library types:

> Paired:
¤ RF: first read (/1) of fragment pair is sequenced as anti-sense (reverse(R)), and second read (/2) is
in the sense strand (forward(F)); typical of the dUTP/UDG sequencing method.
¤ FR: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) is in the
antisense strand (reverse)

> single reads:
¤ F: the single read is in the sense (forward) orientation
¤ R: the single read is in the antisense (reverse) orientation 

This script assumes you use the dUTP/UDG sequencing method, and if you select stranded options (3 or 4), RF (paired) & R (single) parameter is given to Trinity.
If you use another sequencing method, modify this option in the script, or if unknown, use unstranded options (1 or 2). 


The input file(s) must be in ./cleaned_input/assembly and in fastq format.
The assembly files will be in ./assembly/raw/denovo

EndOfMessage

# --SS_lib_type RF

PS3='Please enter your choice: '
options=("Option 1" "Option 2" "Option 3" "Option 4" "Quit")
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
        "Option 4")
            argu=4
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
echo "Type the name of the assembled species (without space ex: drer_trinity), followed by [ENTER]:"
read name
echo ""
echo ""
echo "Type the number of cores to use, followed by [ENTER]:"
read cores
echo ""
echo ""
echo "Type the max memory tu use (in GB), followed by [ENTER]:"
read maxram
echo ""
echo ""

prompt="Please select the input file:"
options=( $(find ./cleaned_input/assembly -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

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
 
   
echo ""
echo ""

mkdir -p ./logs
rm -rfv trinity_part_"$name"
mkdir -p trinity_part_"$name"
mkdir -p assembly/raw/denovo
cd ./trinity_part_"$name"
../bin/trinityrnaseq/Trinity --seqType fq --single ."$R1" --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt
#../bin/trinityrnaseq/Trinity --seqType fq --single ."$R1" --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt

#--normalize_reads
#nice -n 19 ../bin/trinityrnaseq/Trinity --seqType fq --left ../cleaned_input/assembly/R1_cleaned_sync.fastq --right ../cleaned_input/assembly/R2_cleaned_sync.fastq --normalize_reads --CPU 8 --max_memory 12G | tee ../logs/logs_trinity_R1_R2.txt


cd trinity_out_dir
cp ../../bin/convert_fasta_to_fastq.py ./
mv Trinity.fasta "$name"_Trinity.fasta
python convert_fasta_to_fastq.py ./"$name"_Trinity.fasta "$name"_Trinity.fastq
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fasta
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fastq
mv -v "$name"_Trinity.* ../../assembly/raw/denovo/

fi

if [ "$argu" = "2" ]; then

echo ""
echo ""
echo "Type the name of the assembled species (without space ex: drer_trinity), followed by [ENTER]:"
read name
echo ""
echo ""
echo "Type the number of cores to use, followed by [ENTER]:"
read cores
echo ""
echo ""
echo "Type the max memory tu use (in GB), followed by [ENTER]:"
read maxram
echo ""
echo ""

prompt="Please select R1 file:"
options=( $(find ./cleaned_input/assembly -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

PS3="$prompt "
select R1 in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  ""
        break

    else
        echo "Invalid option. Try another one."
    fi
done   
 
echo ""
echo ""

prompt="Please select R2 file:"
options=( $(find ./cleaned_input/assembly -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

PS3="$prompt "
select R2 in "${options[@]}" "Quit" ; do 
    if (( REPLY == 1 + ${#options[@]} )) ; then
        exit

    elif (( REPLY > 0 && REPLY <= ${#options[@]} )) ; then
        echo  ""
        break

    else
        echo "Invalid option. Try another one."
    fi
done    

echo ""
echo ""

mkdir -p ./logs
rm -rfv trinity_part_"$name"
mkdir -p trinity_part_"$name"
mkdir -p assembly/raw/denovo
cd ./trinity_part_"$name"
../bin/trinityrnaseq/Trinity --seqType fq --left ."$R1" --right ."$R2" --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt
#nice -n "$nicevalue" ../bin/trinityrnaseq/Trinity --seqType fq --left ."$R1" --right ."$R2" --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt

#--normalize_reads
#nice -n 19 ../bin/trinityrnaseq/Trinity --seqType fq --left ../cleaned_input/assembly/R1_cleaned_sync.fastq --right ../cleaned_input/assembly/R2_cleaned_sync.fastq --normalize_reads --CPU 8 --max_memory 12G | tee ../logs/logs_trinity_R1_R2.txt


cd trinity_out_dir
cp ../../bin/convert_fasta_to_fastq.py ./
mv Trinity.fasta "$name"_Trinity.fasta
python convert_fasta_to_fastq.py ./"$name"_Trinity.fasta "$name"_Trinity.fastq
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fasta
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fastq
mv -v "$name"_Trinity.* ../../assembly/raw/denovo/

fi

if [ "$argu" = "3" ]; then

echo ""
echo ""
echo "Type the name of the assembled species (without space ex: drer_trinity), followed by [ENTER]:"
read name
echo ""
echo ""
echo "Type the number of cores to use, followed by [ENTER]:"
read cores
echo ""
echo ""
echo "Type the max memory tu use (in GB), followed by [ENTER]:"
read maxram
echo ""
echo ""

prompt="Please select the input file:"
options=( $(find ./cleaned_input/assembly -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

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
 
   
echo ""
echo ""

mkdir -p ./logs
rm -rfv trinity_part_"$name"
mkdir -p trinity_part_"$name"
mkdir -p assembly/raw/denovo
cd ./trinity_part_"$name"
../bin/trinityrnaseq/Trinity --seqType fq --single ."$R1" --SS_lib_type R --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt
#../bin/trinityrnaseq/Trinity --seqType fq --single ."$R1" --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt

#--normalize_reads
#nice -n 19 ../bin/trinityrnaseq/Trinity --seqType fq --left ../cleaned_input/assembly/R1_cleaned_sync.fastq --right ../cleaned_input/assembly/R2_cleaned_sync.fastq --normalize_reads --CPU 8 --max_memory 12G | tee ../logs/logs_trinity_R1_R2.txt


cd trinity_out_dir
cp ../../bin/convert_fasta_to_fastq.py ./
mv Trinity.fasta "$name"_Trinity.fasta
python convert_fasta_to_fastq.py ./"$name"_Trinity.fasta "$name"_Trinity.fastq
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fasta
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fastq
mv -v "$name"_Trinity.* ../../assembly/raw/denovo/

fi

if [ "$argu" = "4" ]; then

echo ""
echo ""
echo "Type the name of the assembled species (without space ex: drer_trinity), followed by [ENTER]:"
read name
echo ""
echo ""
echo "Type the number of cores to use, followed by [ENTER]:"
read cores
echo ""
echo ""
echo "Type the max memory tu use (in GB), followed by [ENTER]:"
read maxram
echo ""
echo ""

prompt="Please select R1 file:"
options=( $(find ./cleaned_input/assembly -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

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
 
echo ""
echo ""

prompt="Please select R2 file:"
options=( $(find ./cleaned_input/assembly -maxdepth 1 -type f -iregex '.*\.\(fastq\|fq\)$' -print0 | xargs -0) )

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

echo ""
echo ""

mkdir -p ./logs
rm -rfv trinity_part_"$name"
mkdir -p trinity_part_"$name"
mkdir -p assembly/raw/denovo
cd ./trinity_part_"$name"
../bin/trinityrnaseq/Trinity --seqType fq --left ."$R1" --right ."$R2" --SS_lib_type RF --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt
#nice -n "$nicevalue" ../bin/trinityrnaseq/Trinity --seqType fq --left ."$R1" --right ."$R2" --normalize_reads --CPU "$cores" --max_memory "$maxram"G | tee ../logs/logs_trinity_R1_R2.txt

#--normalize_reads
#nice -n 19 ../bin/trinityrnaseq/Trinity --seqType fq --left ../cleaned_input/assembly/R1_cleaned_sync.fastq --right ../cleaned_input/assembly/R2_cleaned_sync.fastq --normalize_reads --CPU 8 --max_memory 12G | tee ../logs/logs_trinity_R1_R2.txt


cd trinity_out_dir
cp ../../bin/convert_fasta_to_fastq.py ./
mv Trinity.fasta "$name"_Trinity.fasta
python convert_fasta_to_fastq.py ./"$name"_Trinity.fasta "$name"_Trinity.fastq
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fasta
rm -rfv ../../assembly/raw/denovo/"$name"_Trinity.fastq
mv -v "$name"_Trinity.* ../../assembly/raw/denovo/

fi
