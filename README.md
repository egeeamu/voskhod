# voskhod
repository for the Transcriptome inference bioniformatic pipeline Voskhod

Warning : never use special characters or space in the names/values/answers
When asked "Please select Reference species" please select ensembl species !!

Warning : 
If test with sim data : use "sim_data.sh" scripts
If test with real data : use "real_data.sh" scripts


Install needed tools with root user :
aptitude install build-essential git python-pip  pbzip2 pigz zlib1g-dev libncurses5-dev bowtie parallel python-biopython sqlite3

aptitude install software-properties-common

add-apt-repository ppa:webupd8team/java

apt-get update

aptitude install oracle-java8-installer

pip install biomart

pip install bashplotlib

once  : 1. Prepare the pipeline directory:
(compil trinity & unpack blast)

> chmod +x *.sh
> ./00_unpack_blast_make_trinity.sh


2. Download the reference species.
go to http://www.ensembl.org/biomart/    and select Ensemble Genes 86 then your reference species in the dropdown list then select the show query as perl script and get the data file name in the setdataset variable (drerio_gene_ensembl in the case of Danio rerio)

Put in the script the correct value for the species
The reference cdna will be in ./reference_ts (with the correct sqlite format as a db file)
launch the script 
> python 00_lunokod_cdna_downloader.py

2. Clean and sync the R1 and R2 (fastq files) you want assembly with trinity.
If you have manys R1 and R2 files, concatenate them into one "big" R1 & one "big" R2.  (ex cat *R1* > bigR1.fastq)
The files must be in "./raw_input/assembly"

> ./01_01_cleaning_and_sync_R1_R2.sh

this will generate cleande files in ./cleaned_input/assembly

3. Lunch Trinity assembly (selecting the correct value while asking):
-Cores to use : the number of core on your computer (type 8 if you dont know)
-Ram to use : total ram on your computer minus 2 (ex 16Gb - 2 > 14GB)
-Nice : the priority of the process with 19 the lower and 0 the maximum

> ./01_trinity_02_assembling_R1_R2_de_novo.sh



4. Validate and anotate the assembly:

The assembly to validate must be in ./assembly/raw in fastq format.
The validated assembly will be in ./assembly/validated in sqlite format.

> ./99_voskhod_validate_assembly_.sh
sortie dans assembly/validated


Do steps 2 , 3 & 4 for all your species before continue.

5. Merge all de novo assemblys and reference cdna :
Copy your validated Trinity's assemblys (from ./assembly/validated)  and your reference cdna (from ./reference_ts) in ./assembly/tomerge
caution : put only the cdna(s) you want merge into ./assembly/tomerge

> ./XX_merge_cdna.sh

the result will be in ./reference_ts
	
Now we can identify reads using ref species and de novo assembly.


5. Merge (with Pear) and cleaning R1 & R2 before Voskhod assembly:
Using the same files as step 2.

> ./02_01_pear_and_clean_R1_R2.sh

6. Identify (blast) and make Voskhod assembly:
When asked, select your merged cdna (all species + ref), and the merged fastq file of your species.
attention donner le ficier correspondant à la sequece de reference, pas le mergé
> ./02_02_voskhod_matching_and_assembling.sh 

7. Validate and anotate the assembly (like step 4):
When asked, select the assembly in the list.

> ./99_voskhod_validate_assembly_.sh

Repeat steps 5 , 6 & 7 for all your species before continue.


8. Merge all de novo assemblys and reference cdna :
Copy your validated Voskhod's assemblys (from ./assembly/validated)  and your reference cdna (from ./reference_ts) in ./assembly/tomerge
caution : put only the cdna(s) you want merge into ./assembly/tomerge

> ./XX_merge_cdna.sh


9. Clean & merge all your R1/R2 for expression :
The inputs files must be in "./raw_input/expression", the result will be in "./cleaned_input/expression/" 
Caution, the names must countain "_R1_"/"_R2_"

> ./03_01_pear_and_clean_many_R1_R1_expression.sh

10. Expression step (blast) :
When asked, select the merged cdna (step 8)

> ./03_02_voskhod_expression.sh

The output will be in ./expression_results/xxx/  in csv and sqlite format.
(with xxx the name of the study)
