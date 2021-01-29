## once you have conda installed, close the shell, reopen it and paste this 
## following line :
conda config --set auto_activate_base false

########################################################################
#STEP 1 : Create a new environment obitools

ENVYAML=./dada2_and_obitools/obitools_env_conda.yaml
conda env create -f $ENVYAML

########################################################################
#STEP 2 : Pair-ended merging

## unzip your data if you need :
unzip mullus_surmuletus_data.zip

## activate your environment :
conda activate obitools

## use the command "illuminapairedend" to make the pair-ended merging
## from the forward and reverse sequences you have in your data :
illuminapairedend --score-min=40 -r mullus_surmuletus_data/Aquarium_2_R1.fastq mullus_surmuletus_data/Aquarium_2_R2.fastq > Aquarium_2.fastq

## this command creates a new ".fastq" file which will contain the 
## sequences after the merging of forward and reverse strands

## alignments which have a quality score higher than 40
## (-- score-min=40) are merged and annotated "aligned", while
## alignemnts with a lower quality score are concatenated and
## annotated "joined"

## to only conserve the sequences which have been aligned, use "obigrep" :
obigrep -p 'mode!="joined"' Aquarium_2.fastq > Aquarium_2.ali.fastq

## "-p" requires a python expression

## python creates a new dataset (".ali.fastq") which only contains the
## sequences annotated "aligned"

########################################################################
#STEP 3 : Demultiplexing

## to compare the sequences next, you need to remove the tags and the
## primers, by using the "ngsfilter" command
ngsfilter -t mullus_surmuletus_data/Med_corr_tags.txt -u Aquarium_2.unidentified.fastq Aquarium_2.ali.fastq > Aquarium_2.ali.assigned.fastq

## new files are created :
## ".unidentified.fastq" file contains the sequences that were not 
## assigned whith a correct tag
## ".ali.assigned.fastq" file contains the sequences that were assigned 
## with a correct tag, so it contains only the barcode sequences

## separate your ".ali.assigned.fastq" files depending on their samples, 
## in placing them in a  dedicated folder (useful for next steps) :
mkdir samples

## create the folder

mv -t samples Aquarium_2.ali.assigned.fastq

## place the latests ".fastq" files in the folder

cd samples
obisplit -t sample --fastq Aquarium_2.ali.assigned.fastq

## separate the files depending on their sample

mv -t ./dada2_and_obitools Aquarium_2.ali.assigned.fastq

## remove the original files from the folder