#STEP 1 : Be prepared

## load the package :
library("dada2")

## create a path to your .fastq files :
path1 <- "./pipeline_obitools_and_dada2/samples"

## list the files in your path to verify the .fastq files are in it :
list.files(path1)

## select the .fastq files you want to analyze :
fns <- sort(list.files(path1, pattern = ".fastq", full.names = T))

## sort is a function that can be used to extract some files, from the list of 
## the path files here

## the function will only extract files that end with the pattern chosen and
## they will be extracted their full names, that is to say with their path

## then you can only keep the part of your files name you want :
sample.names <- sapply(strsplit(basename(fns), ".fastq"), '[', 1)

## sapply permits to apply a function to a vector, returning a vector

## here, the function is applied to the only vector fns

## the function base name removes all the path up to the file name

## the function strsplit removes the pattern written

########################################################################
#STEP 2 : Inspect the quality profiles of the reads

## display the plot of the reads quality profile for one file :
plotQualityProfile(fns[10])

## the plot gives different information :
## the grey-scale represents the quality score frequency at base each position
## on the sequences : darker is the plot, higher is the frequency
## the lines show summary statistics : mean in green, median in orange, and
## first and third quartiles in dashed orange
## the red line indicates the percentage of reads that extend to at least the
## position corresponding to the abscissa on the horizontal axis

########################################################################
#STEP 3 : Filtering & Trimming

## begin the creation of the new files and folder :
filts1 <- file.path(path1, "filtered", paste0(sample.names, ".filt.fastq.gz"))

## file.path builds the path to the new folder, which will be located in the
## path already used and which name will be "filtered"

## the files will be named as described before with sample.names, and the
## pattern ".filt.fastq.gz" will be added


## from the .fastq files of fns, create the new .fastq files of filts after
## filtering and trimming :
out <- filterAndTrim(fns, filts1,
                     truncLen = 234,
                     maxN = 0,
                     maxEE = 1,
                     compress = T,
                     verbose = T)

## truncLen value is chosen considering the marker length and define were the
## reads will be trimmed
## reads which are shorten than this value are filtered

## maxN is the number of N tolerated in the sequences after filtering

## maxEE defines the maximal number of expected errors tolerated in a read,
## based on the quality score (EE = sum(10^(-Q/10)))

## compress = T means that the files will be gzipped

## verbose = T means that information concerning the number of sequences after
## filtering will be given

# create a new path to these files, because some files of filts1 may have been
# completely filtered :

path2 <- "./pipeline_obitools_and_dada2/samples/filtered"
filts2 <- sort(list.files(path2, pattern = ".fastq", full.names = T))

########################################################################
#STEP 4 : Dereplication

## derepFastq function eliminates the amplicons of each sequence in the files
derep <- derepFastq(filts2)

## the function annotates each sequence with his abundance

########################################################################
#STEP 5 : Error rates learning

## generate the error model :
err1 <- learnErrors(derep, randomize = T, errorEstimationFunction = loessErrfun)
err2 <- learnErrors(derep, randomize = T, errorEstimationFunction = noqualErrfun)

## errorEstimationFunction can take the value "loessErrfun" if you want to use
## the quality score of your sequence to build the model, or "noqualErrfun" if
## you don't want it

## plot the estimated error rates :
plotErrors(err1, nominalQ = T)
plotErrors(err2, nominalQ = T)

## eliminate the false sequences according to the model :
dadas1 <- dada(derep, err1)
dadas2 <- dada(derep, err2)

########################################################################
# STEP 6 : Observing the ASVs returned

## build the sequence tables :
seqtab_avec_chimeres1 <- makeSequenceTable(dadas1)
seqtab_avec_chimeres2 <- makeSequenceTable(dadas2)

## create the .csv files :
write.csv(seqtab_avec_chimeres1, "./pipeline_obitools_and_dada2/ASVs_avec_chimeres_loessErrfun.csv")
write.csv(seqtab_avec_chimeres2, "./pipeline_obitools_and_dada2/ASVs_avec_chimeres_noqualErrfun.csv")

## create the .fasta files :
uniqueSeqs1 <- getUniques(seqtab_avec_chimeres1)
uniqueSeqs2 <- getUniques(seqtab_avec_chimeres2)

## extract the vectors from the table

uniquesToFasta(uniqueSeqs1, ".pipeline_obitools_and_dada2/ASVs_avec_chimeres_loessErrfun.fasta")
uniquesToFasta(uniqueSeqs2, ".pipeline_obitools_and_dada2/ASVs_avec_chimeres_noqualErrfun.fasta")

## create the fasta file with these vectors

########################################################################
# STEP 7 : Removing chimeras

## remove bimeras ans repeat step 6 :
seqtab_sans_chimeres1 <- removeBimeraDenovo(seqtab_avec_chimeres1, verbose = T)
seqtab_sans_chimeres2 <- removeBimeraDenovo(seqtab_avec_chimeres2, verbose = T)

write.csv(seqtab_sans_chimeres1, ".pipeline_obitools_and_dada2/ASVs_sans_chimeres_loessErrfun.csv")
write.csv(seqtab_sans_chimeres2, ".pipeline_obitools_and_dada2/ASVs_sans_chimeres_noqualErrfun.csv")

uniqueSeqs1 <- getUniques(seqtab_sans_chimeres1)
uniqueSeqs2 <- getUniques(seqtab_sans_chimeres2)

uniquesToFasta(uniqueSeqs1, ".pipeline_obitools_and_dada2/ASVs_sans_chimeres_loessErrfun.fasta")
uniquesToFasta(uniqueSeqs2, ".pipeline_obitools_and_dada2/ASVs_sans_chimeres_noqualErrfun.fasta")

########################################################################
# STEP 8 : Summary of the filtering

## build the recapitulative tables :

getN <- function(x) sum(getUniques(x))

track1 <- cbind(out, sapply(dadas1, getN), rowSums(seqtab_sans_chimeres1))
colnames(track1) <- c("input", "filtered", "denoised", "nonchim")
rownames(track1) <- sample.names
track1

track2 <- cbind(out, sapply(dadas2, getN), rowSums(seqtab_sans_chimeres2))
colnames(track2) <- c("input", "filtered", "denoised", "nonchim")
rownames(track2) <- sample.names
track2