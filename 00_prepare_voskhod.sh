#!/bin/bash
set -e




cat << EndOfMessage

 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ Preloader and downloader

Version 20170331
Voskhod Pipeline version V1.1
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(CC-BY-NC-ND 4.0 International license) 
Arnaud Ungaro contact@arnaud-ungaro.fr

This script will download and unpack: Spades 3.10.1 ,Blast+ 2.6.0, Pear 0.9.10, Trinity v2.2.0 & compile Trinity
Unpack Cap3 v02/10/15 x64

Before continue, make sure you have on your system :

build-essential git python-pip  pbzip2 pigz zlib1g-dev libncurses5-dev bowtie parallel python-biopython sqlite3
java jre (at least 8) and install the following python packages : biomart (v0.5.0!!) & bashplotlib

If not :

aptitude install build-essential git python-pip pbzip2 pigz zlib1g-dev libncurses5-dev bowtie parallel python-biopython sqlite3 software-properties-common
add-apt-repository ppa:webupd8team/java
aptitude update
aptitude install oracle-java8-installer
aptitude install oracle-java8-set-default
pip install -Iv http://arnaud-ungaro.fr/voskhod/redist/lib/biomart-0.5.0.tar.gz
pip install bashplotlib

The pipeline is tested on a fresh Debian 7 x64 and fresh Ubuntu 16.04 LTS x64


EndOfMessage

read -n1 -r -p "If you read the message above, press C " key

echo ""

if [ "$key" = 'C' ]; then

mkdir -p ./logs
mkdir -p ./raw_input/
mkdir -p ./raw_input/expression
mkdir -p ./raw_input/assembly
mkdir -p ./assembly/
mkdir -p ./assembly/tomerge
mkdir -p ./failsafe_input

cd ./bin/
rm -rfv ./bin/trinityrnaseq-2.2.0

#rm -rfv ./SPAdes-3.6.2-Linux.tar.gz
#wget -O ./SPAdes-3.6.2-Linux.tar.gz "http://spades.bioinf.spbau.ru/release3.6.2/SPAdes-3.6.2-Linux.tar.gz"
#rm -rfv ./SPAdes-3.10.1-Linux.tar.gz
#wget -O ./SPAdes-3.10.1-Linux.tar.gz "http://spades.bioinf.spbau.ru/release3.10.1/SPAdes-3.10.1-Linux.tar.gz"



rm -rfv ./SPAdes-3.10.1-Linux.tar.gz
rm -rfv ./SPAdes-3.10.1-Linux
wget -O ./SPAdes-3.10.1-Linux.tar.gz "http://spades.bioinf.spbau.ru/release3.10.1/SPAdes-3.10.1-Linux.tar.gz"
tar xvf SPAdes-3.10.1-Linux.tar.gz
rm -rfv ./SPAdes-3.10.1-Linux.tar.gz
mv ./SPAdes-3.10.1-Linux SPAdes
tar cvf SPAdes.tar.gz SPAdes
rm -rfv ./SPAdes


rm -rfv cap3.linux.x86_64.tar
rm -rfv cap3.linux.x86_64.tar.gz
wget -O ./cap3.linux.x86_64.tar "http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar"
gzip cap3.linux.x86_64.tar


#rm -rfv ./pear
#rm -rfv ./pear-0.9.10-bin-64.tar.gz
#wget -O ./pear-0.9.10-bin-64.tar.gz "http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz"

tar xvf ./pear-0.9.10-bin-64.tar.gz
#rm -rfv ./pear-0.9.10-bin-64.tar.gz
mv ./pear-0.9.10-bin-64 ./pear
mv ./pear/pear-0.9.10-bin-64 ./pear/pear 
chmod +x ./pear/pear

rm -rfv ./trinity_v2.2.0.tar.gz
wget -O ./trinity_v2.2.0.tar.gz "https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.2.0.tar.gz"
tar xvf ./trinity_v2.2.0.tar.gz
mv trinityrnaseq-2.2.0 trinityrnaseq

#rm -rfv ./ncbi-blast-2.4.0+-x64-linux.tar.gz
#wget -O ./ncbi-blast-2.4.0+-x64-linux.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz"
#tar xvf ./ncbi-blast-2.4.0+-x64-linux.tar.gz 
rm -rfv ./ncbi-blast
rm -rfv ./ncbi-blast-2.6.0+-x64-linux.tar.gz
wget -O ./ncbi-blast-2.6.0+-x64-linux.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz"
tar xvf ./ncbi-blast-2.6.0+-x64-linux.tar.gz
rm -rfv ./ncbi-blast-2.6.0+-x64-linux.tar.gz
mv ./ncbi-blast-2.6.0+ ./ncbi-blast                                      

#create link to the latest_blast without breaking the pipeline
#ln -s ncbi-blast-2.6.0+ ncbi-blast


cd trinityrnaseq
make -j 8 | tee ../../logs/logs_compil_trinity.txt


else
echo "ABORT"
exit
# Anything else pressed, do whatever else.
# echo [$key] not empty
fi







