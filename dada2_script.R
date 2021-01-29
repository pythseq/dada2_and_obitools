#STEP 1 : Be prepared

## load the package :
library("dada2")

## create a path to your ".fastq" files :
path <- "./dada2_and_obitools/samples"

## select the ".fastq" files you want to analyze :
fns <- sort(list.files(path, pattern = ".fastq", full.names = T))

## "sort" is a function that can be used to extract some files, from
## the list of the path files here

## the function only extracts files that end with the chosen pattern
## and they are extracted with their whole path

## then you can only keep the part of your files name you want :
sample.names <- sapply(strsplit(basename(fns), ".fastq"), '[', 1)

## "sapply" permits to apply a function to a vector, returning a vector

## here, the function is applied to the only vector "fns"

## the function "basename" removes all the path up to the file name

## the function "strsplit" removes the pattern written

########################################################################
#STEP 2 : Inspect the quality profiles of the reads

## display the plot of the reads quality profile for one file :
plotQualityProfile(fns[10])

## the plot gives different information :
## the grey-scale represents the quality score frequency at each base
## position on the sequences : darker is the plot, higher is the
## frequency
## the lines show summary statistics : mean in green, median in orange,
## and first and third quartiles in dashed orange
## the red line indicates the percentage of reads that extend to at least
## the position corresponding to the abscissa on the horizontal axis

########################################################################
#STEP 3 : Filtering & Trimming

## begin the creation of the new files and folder :
filts <- file.path(path, "filtered", paste0(sample.names, ".filt.fastq.gz"))

## "file.path" builds the path to the new folder, which will be located
## in the path already used and which name will be "filtered"

## the files will be named as described before with "sample.names", and
## the pattern ".filt.fastq.gz" will be added

## from the ".fastq files" of "fns", create the new ".fastq" files of
## "filts" after filtering and trimming :
out <- filterAndTrim(fns, filts,
                     truncLen = 235,
                     maxN = 0,
                     maxEE = 1,
                     compress = T,
                     verbose = T)

## "truncLen" value is chosen considering the marker length and define
## were the reads will be trimmed
## reads which are shorten than this value are filtered

## "maxN" is the number of N tolerated in the sequences after
## filtering

## "maxEE" defines the maximal number of expected errors tolerated in a
## read, based on the quality score (EE = sum(10^(-Q/10)))

## "compress = T" means that the files will be gzipped

## "verbose = T" means that information concerning the number of sequences after
## sequences after filtering will be given

########################################################################
#STEP 4 : Dereplication

## "derepFastq" function eliminates all the replications of each sequence in the files
derep <- derepFastq(filts)

## the function annotates each sequence with his abundance

########################################################################
#STEP 5 : Error rates learning

## generate the error model :
err <- learnErrors(derep, randomize = T)

## plot the estimated error rates :
plotErrors(err, nominalQ = T)

## eliminate the false sequences according to the model :
dadas <- dada(derep, err)

########################################################################
# STEP 6 : Generate your final files

## build the sequence table :
seqtab_with_bimeras <- makeSequenceTable(dadas)

## create the ".fasta" file :
uniqueSeqs <- getUniques(seqtab_with_bimeras)

## extract the vectors from the table

uniquesToFasta(uniqueSeqs, "./dada2_and_obitools/ASVs_with_bimeras.fasta")

########################################################################
# STEP 7 : Bimeras removal

## remove bimeras and repeat step 6 :
seqtab_without_bimeras <- removeBimeraDenovo(seqtab_with_bimeras, verbose = T)

uniqueSeqs1 <- getUniques(seqtab_without_bimeras)

uniquesToFasta(uniqueSeqs, "./dada2_and_obitools/ASVs_without_bimeras.fasta")