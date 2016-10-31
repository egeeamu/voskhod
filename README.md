# voskhod
This is the repository for the Transcriptome inference and RNAseq expression bioniformatics pipeline Voskhod.

Voskhod allows conducting the transcriptome inference and gene expression analyses as detailed in our paper Challenges and solutions for transcriptome assembly in non-model organisms with an application to hybrid specimens

Arnaud Ungaro, Nicolas Pech, Jean-Francois Martin, Scott RJ McCairns, Remi Chappaz, Andre Gilles

doi: http://dx.doi.org/10.1101/084145

http://www.biorxiv.org/content/early/2016/10/28/084145

#Installation

Voskhod has being extensively tested (in particular with paired end data) on Xubuntu 16.04 LTS and all instructions refer to this setup.
It surely works on other flavours of Linux OS but we will not be able to provide technical support on this kind of questions.

If using a fresh Ubuntu install, you will have to install dependencies and linux packages with administrator privileges :

> aptitude install build-essential git python-pip  pbzip2 pigz zlib1g-dev libncurses5-dev bowtie parallel python-biopython sqlite3 software-properties-common

> add-apt-repository ppa:webupd8team/java

> aptitude update

> aptitude insta ll oracle-java8-installer

> pip install biomart

> pip install bashplotlib

and finally get Voskhod!


> git clone https://github.com/egeeamu/voskhod

Make sure you have all permissions on the voskhod folders (else all commands would require using sudo permissions).
The first time you use the pipeline you have to download associated tools and compile them. The following commands achieve this task automatically. In the voskhod folder type:

> chmod +x *.sh

> ./00_prepare_voskhod.sh

Once this step is done you are ready to use Voskhod!

#Transcriptome assembly


1- Download and format the reference transcriptome from Ensembl.

From https://github.com/egeeamu/voskhod/blob/master/dataset_ensembl.txt select the Ensembl Genes 86 dataset name for your prefered  reference species (e.g. drerio_gene_ensembl in the case of Danio rerio) and edit the 01_cdna_downloader.py script with the correct dataset name for your reference transcriptome. Once you execute the python script (see below), the reference trancriptome will be downloaded and formated in ./reference_ts (with the correct sqlite format as a db file):

> python 01_cdna_downloader.py

2- Filter raw input data (uncompressed fastq files) before denovo assembly.
If you have multiple R1 and R2 files, concatenate them into one "large" R1 & one "large" R2.  (e.g. cat *R1* > largeR1.fastq) unless you need to keep track of the source for further analyses.

The input files must be in ./raw_input/assembly

> ./02_filtering_raw_data_denovo.sh

this will generate output files filtered for quality in ./cleaned_input/assembly/ (synchronized R1 and R2 files when Paired-End sequencing is selected)

3- Launch Denovo Trinity assembly (selecting the correct value when requested):
-Cores to use : the number of core on your computer available for the analysis
-Ram to use : total ram on your computer minus 2Gb (ex 16Gb - 2 > 14Gb)

> ./03_denovo_assembly.sh

this will generate output files ./assembly/raw/

4- Validate and annotate the assembly:

The assembly to validate must be in ./assembly/raw in fastq format. The validated and annotated assembly will be written in ./assembly/validated in sqlite format.

> ./04_validate_annotate_assembly.sh

If you have multiple sources for transcriptomes, repeat steps 2 , 3 & 4 for all your sources before continue.
If you only want denovo assemblies, you are good to go and can explore the results with any database browser as DB Browser for SQLite (http://sqlitebrowser.org/). 

5- Merging transcriptomes for further combined approach
If you want to use the combined Denovo + DGM approach as described in the paper, you have to merge the reconstructed denovo transcriptome(s) with the reference species before running the DGM step on the combined refrence.
To merge  denovo assembly(ies) and reference transcriptome :

Copy your validated denovo assembly(ies) (from ./assembly/validated) and your reference trascriptome (from ./reference_ts) in ./assembly/tomerge

caution : put only the transcriptome(s) you want merged into ./assembly/tomerge as the whole content of the folder will be used in merging step.

> ./05_merge_transcriptomes.sh

the combined output will be stored in ./reference_ts

Now we can identify reads using the combination of reference transcriptome species and the denovo assembly.


6- Prepare reads for Voskhod assembly:
Using the same raw input fastq files as step 2, filter reads for quality and optionally merge paired-end reads if relevant. The input fastq file(s) must be in ./raw_input/assembly/

> ./06_prepare_voskhod.sh

The filtered ans merged (when PE sequencing is used) reads will be found in ./cleaned_input/assembly/

7- Identify (blast) and make Voskhod assembly:
When requested, select your merged transcriptomes (all species + ref), and the filtered fastq file for your source to be inferred.

> ./07_voskhod_matching_and_assembling.sh

The script writes its output in ./assembly/raw/voskhod/. The next validation steps requires that you move (or copy) yourself the produced fastq file to ./assembly/raw before processing.


8- Validate and anotate the assembly (like step 4):
When asked, select the assembly from the corresponding fastq file found at ./assembly/raw

> ./08_voskhod_validate_assembly_.sh

Repeat steps 6, 7 & 8 for all your sources before continue if you need to merge transcriptomes.

if you are interrested in one single transcriptome you are good to go and can explore the results with any database browser as DB Browser for SQLite (http://sqlitebrowser.org/). 

9- Merge all  assemblies and reference cdna :
Copy your validated Voskhod's assemblys (from ./assembly/validated)  and your reference cdna (from ./reference_ts) in ./assembly/tomerge
caution : put only the transcriptome(s) you want merged into ./assembly/tomerge

> ./09_merge_transcriptomes.sh

The assembly is done and you can explore the results with any database browser as DB Browser for SQLite (http://sqlitebrowser.org/). 
It is usually a better idea to merge fastq files from multiples sources rather than merging the inferred transcriptomes as you get more coverage during the analyses. This although means you loose track of the source during the process. Merging transcriptomes in a single database allows identifying the transcript source.

#RNAseq Expression

For expression analysis, we postulate that your reference transcriptome(s) is/are obtained and correctly formatted in a sqlite database.

Clean & merge all your R1/R2 for expression :
The inputs files must be in "./raw_input/expression", the result will be in "./cleaned_input/expression/" 
Caution, the names must countain "_R1_"/"_R2_"

> ./03_01_pear_and_clean_many_R1_R1_expression.sh

10. Expression step (blast) :
When asked, select the merged cdna (step 8)

> ./03_02_voskhod_expression.sh

The output will be in ./expression_results/xxx/  in csv and sqlite format.
(with xxx the name of the study)
