Voskhod Pipeline version V1.2
Part of the Voskhod project
https://github.com/egeeamu/voskhod

GNU General Public License v3.0
Arnaud Ungaro contact@arnaud-ungaro.fr

# Voskhod
This is the repository for the Transcriptome inference and RNAseq expression bioniformatics pipeline Voskhod.

Voskhod allows conducting the transcriptome inference and gene expression analyses as detailed in our paper Challenges and solutions for transcriptome assembly in non-model organisms with an application to hybrid specimens

Arnaud Ungaro, Nicolas Pech, Jean-Francois Martin, Scott RJ McCairns, Remi Chappaz, Andre Gilles

doi: http://dx.doi.org/10.1101/084145

http://www.biorxiv.org/content/early/2016/10/28/084145

#Installation

Voskhod has being extensively tested (in particular with paired end data) on Xubuntu 16.04 LTS and all instructions refer to this setup.
It surely works on other flavours of Linux OS but we will not be able to provide technical support on this kind of questions.

If using a fresh Ubuntu install, you will have to install dependencies and linux packages with administrator privileges :

> apt-get install build-essential git python-pip  pbzip2 pigz zlib1g-dev libncurses5-dev bowtie parallel python-biopython sqlite3 software-properties-common

> add-apt-repository ppa:webupd8team/java (validate with pressing Enter)

> apt-get update

> apt-get install oracle-java8-installer (accept terms)

> pip install bashplotlib

and finally get Voskhod (as non-root user)!

> git clone https://github.com/egeeamu/voskhod

Make sure you have all permissions on the voskhod folders (else all commands would require using sudo permissions).
The first time you use the pipeline you have to download associated tools and compile them. The following commands achieve this task automatically. In the voskhod folder type:

> cd voskhod/

> chmod +x *.sh

> ./00_prepare_voskhod.sh

Once this step is done you are ready to use Voskhod!

#Transcriptome assembly

1- Download and format the reference transcriptome (one at a time).
The script uses uncompressed fasta files from voskhod/fail_safe/ folder.  
This script transforms a fasta file (preferably downloaded from Ensembl as it provides metadata)

If you are working with ensembl (may work with plan,metazoa,etc..) dataset : 
Get coding rna :
xx.xx.xx.cdna.all.fa.gz from for the species of your choice from :
ftp://ftp.ensembl.org/pub/current_fasta/
for example : ftp://ftp.ensembl.org/pub/current_fasta/danio_rerio/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz

if you like you can also get the corresponding non-coding rna :
and  xx.xx.xx.ncrna.fa.gz for the species of your choice from :
ftp://ftp.ensembl.org/pub/current_fasta/
for example ftp://ftp.ensembl.org/pub/current_fasta/danio_rerio/ncrna/Danio_rerio.GRCz10.ncrna.fa.gz

If compressed, unpack your reference fasta file(s) in "failsafe_input" within the voskhod folder. For example
> gzip -d Danio_rerio.GRCz10.cdna.all.fa.gz and put this in "failsafe_input"

Run the conversion script that converts the reference fasta into a SQlite database
> python 01_cdna_formater_from_fasta.py -s SpeciesName

If you are working with custom data-set, make sure your fasta is formated like this :

>transcriptname_1
sequence

>transcriptname_2
sequence

without space or special character in names, then put it in "failsafe_input".

Run the conversion script that converts the reference fasta into a SQlite database and precise the Species name
> python 01_cdna_formater_from_fasta.py -s SpeciesName

2- Filter raw input data (uncompressed fastq files) before denovo assembly.
If you have multiple R1 and R2 files, concatenate them into one "large" R1 & one "large" R2.  (e.g. cat *R1* > largeR1.fastq) unless you need to keep track of the source for further analyses.

The input files must be in ./raw_input/assembly

> ./02_filtering_raw_data_denovo.sh

this will generate output files filtered for quality in ./cleaned_input/assembly/ (synchronized R1 and R2 files when Paired-End sequencing is selected). If you already cleaned (and paired for PE-reads) your reads you can directly imort your reads to the cleaned_input/assembly folder of course. 

3- Launch Denovo assembly (selecting the correct value when requested):
-Cores to use : the number of cores on your computer available for the analysis
-Ram to use : total ram on your computer minus 2Gb (ex 16Gb - 2 > 14Gb)

> ./03_denovo_assembly.sh

this will generate output files ./assembly/raw/trinity/

4- Validate and annotate the assembly:

The assembly to validate must be in ./assembly/raw/trinity/ in fasta format. The validated and annotated assembly will be written in ./assembly/validated in sqlite format.

> ./04_validate_annotate_assembly.sh

If you have multiple sources for transcriptomes, repeat steps 2 , 3 & 4 for all your sources before continuing.
If you only want denovo assemblies, you are good to go and can explore the results with any database browser as DB Browser for SQLite (http://sqlitebrowser.org/). 

5- Merging transcriptomes for further combined approach
If you want to use the combined approach implemented in Voskohd as described in the paper, you have to merge the reconstructed denovo transcriptome(s) with the reference species before running the transcriptome-guided step on the combined refrence.
To merge  denovo assembly(ies) and reference transcriptome :

Copy your validated denovo assembly(ies) (from ./assembly/validated) and your reference trascriptome (from ./reference_ts) in ./assembly/tomerge

caution : put only the transcriptome(s) you want merged into ./assembly/tomerge as the whole content of the folder will be used in merging step.

> ./05_merge_transcriptomes.sh

the combined output will be stored in ./reference_ts

Now we can identify reads using the combination of reference transcriptome species and the denovo assembly.


6- Prepare reads for transcriptome-guided assembly:
Using the same raw input fastq files as step 2, filter reads for quality and optionally merge paired-end reads if relevant. The input fastq file(s) must be in ./raw_input/assembly/

> ./06_prepare_voskhod.sh

The filtered ans merged (when PE sequencing is used) reads will be found in ./cleaned_input/assembly/

7- Assign (blast) and transcriptome-guided assembly:
When requested, select your merged transcriptomes (all species + ref), and the filtered fastq file for your source to be inferred.

> ./07_assign_and_assemble.sh

The script writes its output in ./assembly/raw/voskhod/. 


8- Validate and anotate the assembly (like step 4):
When asked, select the assembly from the corresponding fastq file found at ./assembly/raw/voskhod/ The reference species required in this step is the reference transcriptome to avoid circularity with the denovo inference.


> ./08_validate_guided_assembly.sh

Repeat steps 6, 7 & 8 for all your sources before continuing if you need to merge transcriptomes.

if you are interrested in one single transcriptome you are good to go and can explore the results with any database browser as DB Browser for SQLite (http://sqlitebrowser.org/). 

9- To prepare for expression studies, merge all  assemblies and reference cdna :
Copy your validated Voskhod's assemblies (from ./assembly/validated)  and your reference cdna (from ./reference_ts) in ./assembly/tomerge
caution : put only the transcriptome(s) you want merged into ./assembly/tomerge

> ./09_merge_transcriptomes_for_expression.sh

It is usually a better idea to merge fastq files from multiples sources rather than merging the inferred transcriptomes as you get more coverage during the analyses. This although means you loose track of the source during the process. Merging transcriptomes in a single database allows identifying the transcript source.

#RNAseq Expression

For expression analysis, we postulate that your reference transcriptome(s) is/are obtained and correctly formatted in a sqlite database.

e1- filter on quality & merge paired end reads if relevant:
The input files must be in "./raw_input/expression", the result will be in "./cleaned_input/expression/" 
Caution, the file names must countain \_R1_ & \_R2_

> ./e01_filtering_raw_data_expression.sh

e2- Computing Expression (blast) :
When requested, select the relevant reference transcriptome (usually the transcriptome from step 8 or 9 in the assembly section).

> ./e02_expression_analysis.sh

The output will be in ./expression_results/STUDYNAME/  in csv and sqlite format.
(with STUDYNAME being the name of the study file)
