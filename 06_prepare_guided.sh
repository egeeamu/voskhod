#!/bin/bash
set -e


cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ Cleaner & merger

Version 20170331
Voskhod Pipeline version V1.1
Part of the Voskhod project
https://github.com/egeeamu/voskhod

GPL-3.0
Arnaud Ungaro contact@arnaud-ungaro.fr


¤ Option 1)

Input : Single end library with only one fastq file
____________________

Output : One filtered fastq file 
____________________


##################################################

¤ Option 2)

Input : Paired end library with R1 & R2 fastq files with overlap exptected
____________________
             ____________________

Output : One filtered & merged (with Pear) fastq file
_________________________________


##################################################

¤ Option 3)

Input : Paired end library with R1 & R2 fastq files without overlap exptected
____________________
                                 ____________________

Output : One filtered & concatenated fastq file
____________________
____________________

The input file(s) must be in ./raw_input/assembly and in fastq format.
The cleaned files will be in ./cleaned_input/assembly in fastq format.


EndOfMessage






cp ./bin/pear/pear ./raw_input/assembly/
cp ./bin/vosklean_sequences_cleaner_assembly_keep_pw.py ./raw_input/assembly
mkdir -p ./cleaned_input/assembly
cd ./raw_input/assembly/




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

echo "Type a name for the filtered fastq (without space ex: drer_merged), followed by [ENTER]:"
read name

if [ "$argu" = "1" ]; then

prompt="Please select the fastq file:"
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

time parallel --gnu --progress --max-procs 1 'nice -n 19 python vosklean_sequences_cleaner_assembly_keep_pw.py {}' ::: "$R1"

rm -rfv ../../cleaned_input/assembly/"$name".assembled.fastq_query.fq
mv "$R1"_query.fq ../../cleaned_input/assembly/"$name".assembled.fastq_query.fq
rm -rfv "$R1"_query.fq 


mv -v ./logs/* ../../logs
rm -rfv logs "$name".assembled.fastq "$name".discarded.fastq "$name".unassembled.forward.fastq "$name".unassembled.reverse.fastq vosklean_sequences_cleaner_assembly_keep_pw.py pear afterclean_remove_bad_pw.py

fi


if [ "$argu" = "3" ]; then

prompt="Please select R1 file:"
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
prompt="Please select R2 file:"
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

echo ""
echo "Concatenating progress ... it may take long time.. "
echo ""



cat "$R2"_query.fq >> "$R1"_query.fq
rm -rfv ../../cleaned_input/assembly/"$name".assembled.fastq_query.fq
mv "$R1"_query.fq ../../cleaned_input/assembly/"$name".assembled.fastq_query.fq
rm -rfv "$R1"_query.fq "$R2"_query.fq


echo ""
echo "Concatenating done ! "
echo ""

mv -v ./logs/* ../../logs
rm -rfv logs "$name".assembled.fastq "$name".discarded.fastq "$name".unassembled.forward.fastq "$name".unassembled.reverse.fastq vosklean_sequences_cleaner_assembly_keep_pw.py pear afterclean_remove_bad_pw.py



fi

if [ "$argu" = "2" ]; then



prompt="Please select R1 file:"
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
 
echo ""
echo ""

prompt="Please select R2 file:"
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

echo ""
echo ""

./pear -j 8 -n 30 -v 8 -f "$R1" -r "$R2" -o "$name"
python ./vosklean_sequences_cleaner_assembly_keep_pw.py ./"$name".assembled.fastq
rm -rfv ../../cleaned_input/assembly/"$name".assembled.fastq_query.fq
mv -v "$name".assembled.fastq_query.fq ../../cleaned_input/assembly/"$name".assembled.fastq_query.fq
mv -v ./logs/* ../../logs
rm -rfv logs "$name".assembled.fastq "$name".discarded.fastq "$name".unassembled.forward.fastq "$name".unassembled.reverse.fastq vosklean_sequences_cleaner_assembly_keep_pw.py pear 

fi
