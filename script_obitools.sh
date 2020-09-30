## once you have conda installed, close the shell, reopen it and paste this 
## following line :
conda config --set auto_activate_base false

########################################################################
#STEP 1 : Create a new environment obitools

ENVYAML=./dada2_and_obitools/obitools_env_conda.yaml
conda env create -f $ENVYAML

########################################################################
#STEP 2 : Pair-end sequencing

## unzip your data if you need :
unzip mullus_surmuletus_data.zip

## activate your environment :
conda activate obitools

## use the function "illuminapairedend" to make the pair-end sequencing from
## the forward and reverse sequences you have in your data :
illuminapairedend --score-min=40 -r mullus_surmuletus_data/Med_R1.fastq mullus_surmuletus_data/Med_R2.fastq > Med.fastq
illuminapairedend --score-min=40 -r mullus_surmuletus_data/Atl_R1.fastq mullus_surmuletus_data/Atl_R2.fastq > Atl.fastq

## this function creates a new ".fastq file" which will contain the sequences
## after the pair-end of forward and reverse sequences which have a quality 
## score higher than 40 ("-- score-min=40")

## to only conserve the sequences which have been aligned, use "obigrep" :
obigrep -p 'mode!="joined"' Med.fastq > Med.ali.fastq
obigrep -p 'mode!="joined"' Atl.fastq > Atl.ali.fastq

## "-p" requires a python expression

## the unaligned sequences are notified with mode="joined" by "illuminapairedend"
## whereas the aligned sequences are notified with mode="aligned"

## so here python creates new datasets (".ali.fastq") which only contain the
## sequences notified "aligned"

########################################################################
#STEP 3 : Demultiplexing

## to be able to compare the sequences next, you need to remove tags and primers,
## and to use the function ngsfilter :
ngsfilter -t mullus_surmuletus_data/Med_corr_tags.txt -u Med.unidentified.fastq Med.ali.fastq > Med.ali.assigned.fastq
ngsfilter -t mullus_surmuletus_data/Atl_corr_tags.txt -u Atl.unidentified.fastq Atl.ali.fastq > Atl.ali.assigned.fastq

## new files are created :
## ".unidentified.fastq" files contain the sequences that were not assigned
## whith a correct tag
## ".ali.assigned.fastq" files contain the sequences that were assigned with
## a correct tag, so it contains only the barcode sequences

## separate your ".ali.assigned.fastq" files depending on their samples, 
## in placing them in a  dedicated folder (useful for next steps) :
mkdir samples

## create the folder

mv -t samples Med.ali.assigned.fastq Atl.ali.assigned.fastq

## place the latests ".fastq" files in the folder

cd samples
obisplit -t sample --fastq Med.ali.assigned.fastq
obisplit -t sample --fastq Atl.ali.assigned.fastq

## separation of the files depending on their sample

mv -t ./dada2_and_obitools Med.ali.assigned.fastq Atl.ali.assigned.fastq

## removing the original files from the folder