# voskhod
This is the repository for the Transcriptome inference bioniformatics pipeline Voskhod.

Voskhod allows conducting the transcriptome inference and gene expression analyses as detailed in our paper Challenges and solutions for transcriptome assembly in non-model organisms with an application to hybrid specimens

Arnaud Ungaro, Nicolas Pech, Jean-Francois Martin, Scott RJ McCairns, Remi Chappaz, Andre Gilles

doi: http://dx.doi.org/10.1101/084145

http://www.biorxiv.org/content/early/2016/10/28/084145

Voskhod has being extensively tested on Xubuntu 16.04 LTS and all instructions refer to this setup.

If using a fresh install, you will have to install dependencies and linux packages with administrator privileges :

> aptitude install build-essential git python-pip  pbzip2 pigz zlib1g-dev libncurses5-dev bowtie parallel python-biopython sqlite3 software-properties-common

> add-apt-repository ppa:webupd8team/java

> aptitude update

> aptitude insta ll oracle-java8-installer

> pip install biomart

> pip install bashplotlib

and finally get Voskhod!


> git clone https://github.com/egeeamu/voskhod

The first time you use the pipeline you have to download associated tools and compile them. The following commands achieve this task automatically. In the voskhod folder type:

> chmod +x *.sh

> ./00_prepare_voskhod.sh

Once this step is done you are ready to use Voskhod!

1- Download the reference species.

From https://github.com/egeeamu/voskhod/blob/master/dataset_ensembl.txt select the Ensembl Genes 86 dataset name for your prefered  reference species (e.g. drerio_gene_ensembl in the case of Danio rerio) and edit the 01_cdna_downloader.py script with the correct dataset name for your reference transcriptome. Once you execute the python script (see below), the reference trancriptome will be downloaded and formated in ./reference_ts (with the correct sqlite format as a db file):

> python 01_cdna_downloader.py

2- Filter raw input data (uncompressed fastq files) before denovo assembly.
If you have multiple R1 and R2 files, concatenate them into one "large" R1 & one "large" R2.  (e.g. cat *R1* > largeR1.fastq)
The input files must be in "./raw_input/assembly"

> ./02_filtering_raw_data_denovo.sh

this will generate input files filtered for quality in ./cleaned_input/assembly (synchronized R1 and R2 files when Paired-End sequencing is selected)

3- Launch Denovo Trinity assembly (selecting the correct value when requested):
-Cores to use : the number of core on your computer available for the analysis
-Ram to use : total ram on your computer minus 2Gb (ex 16Gb - 2 > 14Gb)

> ./03_denovo_assembly.sh


4- Validate and annotate the assembly:

The assembly to validate must be in ./assembly/raw in fastq format. This is precisely where the 03_denovo_assembly.sh script writes its output. The validated and annotated assembly will be written in ./assembly/validated in sqlite format.

> ./validate_annotate_assembly.sh

If you have multiple sources for transcriptomes, repeat steps 2 , 3 & 4 for all your sources before continue.
If you only want denovo assemblies, you are good to go and can explore the results with any database browser as DB Browser for SQLite (http://sqlitebrowser.org/). 
If you want to use the combined Denovo + DGM approach as described in the paper, you have to merge the reconstructed denovo transcriptome(s) with the reference transcriptome before running the DGM step on the combined refrence.
To merge  denovo assembly(ies) and reference transcriptome :

Copy your validated denovo assembly(ies) (from ./assembly/validated)  and your reference trascriptome (from ./reference_ts) in ./assembly/tomerge

caution : put only the transcriptome(s) you want merged into ./assembly/tomerge as the whole content of the folder will be used in merging step.

> ./merge_transcriptomes.sh

the output will be stored in ./reference_ts


Now we can identify reads using the comination of reference transcriptome species and the denovo assembly.


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
