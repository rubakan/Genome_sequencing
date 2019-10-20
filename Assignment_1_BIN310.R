#!/usr/bin/env R

#*************************************************************************
#
#   Program:    Assignment_1.
#   File:       Assigntment_1.R
#   
#   Version:    R
#   Date:       19.10.2019
#   
#   
#   Copyright:  (c) NMBU / Rubakan Thurupan
#   Author:     Rubakan Thurupan
#   Address:    Universitetstunet 3, 1430 ??s
#   EMail:      rubakan.thurupan@nmbu.no

#*************************************************************************************************************************
#
#
#This is first assignment from course genome analysis(BIN310) at NMBU.
#
#********************************************************************************************************************
#
#   Description:
#   ============
#
#********************************************************************************************************************
#
#   Usage:
#   ======
#   
#
#********************************************************************************************************************
#
#   Revision History:
#   =================
#   V1.0   07.10.19  Original
#   V1.1   18.06.19  Final version
#                    
#
#
#********************************************************************************************************************


library(RLinuxModules)
moduleInit(modulesHome = "/local/genome/Modules/3.2.10")
module("load jellyfish")

### Settings that we may change
K <- 15
hash.size <- "16M"
threads <- 4
canonical <- "-C"   # use "" if not canonical K-mers
cnt.file <- "jelly/counts.jf"
hst.file <- "jelly/histo.txt"
in.files <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R1.fastq /mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R2.fastq" # edit this

### Building the command line for counting
cmd <- paste("jellyfish count",
             "-m", K,
             "-s", hash.size,
             "-t", threads,
             canonical,
             "-o", cnt.file,
             in.files)
system(cmd)

### Building the command line for histogram
cmd <- paste("jellyfish histo",
             "-o", hst.file,
             cnt.file)
system(cmd)

library(tidyverse)
Kmers.tb1 <- read_delim(hst.file, delim =" ", col_names = c("Kmer.count", "Frequency"))
ggplot(Kmers.tb1) +
  geom_col(aes(x = Kmer.count, y = Frequency), width = 1.0) +
  #scale_y_log10()
  xlim(0,200) + ylim(0, 150000)

############################################################################################
#3- assembley

###########################################################################################
#Should data be trimmed

library(RLinuxModules)
moduleInit(modulesHome = "/local/genome/Modules/3.2.10")
module("load trimmomatic")

### Settings
R1.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R1.fastq"
R2.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R2.fastq"
R1.paired   <- "/mnt/users/student10/BIN310/1_Assignment/trimmo/paired_R1.fastqc.zip"         # edit this
R2.paired   <- "/mnt/users/student10/BIN310/1_Assignment/trimmo/paired_R2.fastqc.zip"         # edit this
R1.unpaired <- "/mnt/users/student10/BIN310/1_Assignment/trimmo/unpaired_R1.fastqc.zip"       # edit this
R2.unpaired <- "/mnt/users/student10/BIN310/1_Assignment/trimmo/unpaired_R2.fastqc.zip."       # edit this
illuminaclip.adapters <- "/local/genome/packages/trimmomatic/0.36/adapters/TruSeq2-PE.fa"
illuminaclip.options <- ":2:30:10:3:TRUE"
maxinfo <- "MAXINFO:50:0.25"
leading <- "LEADING:10"
trailing <- "TRAILING:10"
slidingwindow <- "SLIDINGWINDOW:5:10"
minlen <- "MINLEN:32"
crop <- ""      # "CROP:25"
headcrop <- ""  # "HEADCROP:10"
threads <- 4

### Command line
cmd <- paste( "java -jar /local/genome/packages/trimmomatic/0.36/trimmomatic-0.36.jar PE",
              "-threads", threads,
              "-quiet",
              R1.file, R2.file,
              R1.paired, R1.unpaired, R2.paired, R2.unpaired,
              paste0("ILLUMINACLIP:", illuminaclip.adapters, illuminaclip.options),
              maxinfo,
              leading,
              trailing,
              slidingwindow,
              crop,
              headcrop,
              minlen)
system(cmd)



#########################################################################################

# How many contigs did you get? 57
#   What is the N50 value? 498161 
#   What is meant by NG50 and NA50? What are their values?
#The NG50 statistic is the same as N50 except that it is 50% of the known or estimated genome size that must be of the NG50 length or longer. This allows for meaningful comparisons between different assemblies. In the typical case that the assembly size is not more than the genome size, the NG50 statistic will not be more than the N50 statistic.

#   What is a misassembly? How many misassemblies are there, and how many bases does this cover?
#   How large fraction of the reference genome is re-constructed in the contigs?


library(RLinuxModules)
moduleInit(modulesHome = "/local/genome/Modules/3.2.10")
module("load spades")


### Settings
K.range <- "21,33,55,77,99,127"
threads <- 4
R1.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R1.fastq"      # edit this
R2.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R2.fastq"      # edit this
out.dir <- "/mnt/users/student10/BIN310/1_Assignment/spades/"   # edit this. Create folder


### Command line
cmd <- paste("spades.py",
             "-k", K.range,
             "--careful",
             "--threads", threads,
             "--pe1-1", R1.file,
             "--pe1-2", R2.file,
             "-o", out.dir)
system(cmd)

#Quast software

library(RLinuxModules)
moduleInit(modulesHome="/local/genome/Modules/3.2.10")
module("load quast")

### Settings - same for all cases
ref.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Bacillus_cereus_ATCC14579.fna"
contigs.file <- "/net/fs-1/home01/student10/BIN310/1_Assignment/spades/contigs.fasta"
out.dir <- "/net/fs-1/home01/student10/BIN310/1_Assignment/quast/"
threads <- 4


cmd <- paste("quast.py",
             "--threads", threads,
             "--pe1", R1.file,         # comment out if evaluating canu+nanopore
             "--pe2", R2.file,         # comment out if evaluating canu+nanopore
             "-r", ref.file,
             "-o", out.dir,
             contigs.file)
system(cmd)

################################### 4 - READ MAPPING ####################################################

#Bowtie

library(RLinuxModules)
moduleInit(modulesHome="/local/genome/Modules/3.2.10")
module("load bowtie2")
module("load samtools")

ref.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Bacillus_cereus_ATCC14579.fna"   #ference
db.name <- "/mnt/users/student10/BIN310/1_Assignment/bowtie/reference"            # bowtie2 will create this database
R1.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R1.fastq"
R2.file <- "/mnt/users/larssn/BIN310/coursedata/assignment1/Illumina_MiSeq_R2.fastq"
sam.file <- "/mnt/users/student10/BIN310/1_Assignment/bowtie/contigs_mapping.sam" # output in sam-format here
threads <- 4

cmd <- paste("bowtie2-build",   # builds database
             "-q",              # -q mean quiet (no output)
             ref.file,          # fasta-file with reference (contigs, genome etc)
             db.name,R1.file,R2.file)           # prefix of database
system(cmd)


cmd <- paste("bowtie2",         # mapping
             "-q",              # quiet
             "-p", threads,     # use this many threads
             "-x", db.name,     # name of database created above
             "--minins 0 --maxins 1500",  # rough estimates of insert size
             "-1", R1.file,     # the R1-reads
             "-2", R2.file,     # the R2-reads
             "-S", sam.file)    # the output sam-file
system(cmd)


coverage.file <- "~/BIN310/1_Assignment/bowtie/read_coverage.txt"
cmd <- paste("samtools view",             # convert sam-file to bam-file
             "--threads", threads,
             "-b", sam.file,
             "|",                         # ...and we pipe into...
             "samtools sort",             # sort by reference positions
             "--threads", threads,
             "|",                         # ...and we pipe into...
             "samtools depth",            # compute read coverage (=depth)
             "-a -a",
             "/dev/stdin",
             ">", coverage.file)               # results in a text file
system(cmd)


library(tidyverse)
tbl <- read_delim("~/BIN310/1_Assignment/bowtie/read_coverage.txt", delim = "\t", col_names = F)
p <- ggplot(tbl) +
  geom_histogram(aes(bindwidt = 1)(x = X3))  
print(p)



################5) Simulated data#################################
library(tidyverse)
library(microseq)
module("load art")

tbl <- readFasta("/mnt/users/larssn/BIN310/coursedata/assignment1/Bacillus_cereus_ATCC14579.fna")
tbl.new <- tbl[1,]
writeFasta(tbl.new, out.file = "reference_new")

library(RLinuxModules)
moduleInit(modulesHome = "/local/genome/Modules/3.2.10")
module("load art")

### Settings that we may change
reference.genome <- "/mnt/users/student10/BIN310/1_Assignment/reference_new"
sequencing.technology <- "MSv3"
insert.size <- 1500
insert.size.std <- 100
read.coverage <- 50
read.length <- 250
read.ID <- "my.ID"
out.prefix <- "feilfire_reads_R"   # existing folders

cmd <- paste("art_illumina",
             "--rndSeed", sample(1:1000000, 1),
             "--seqSys", sequencing.technology,
             "--paired",
             "--errfree",
             "--noALN",
             "--mflen", insert.size,
             "--sdev", insert.size.std,
             "--fcov", read.coverage,
             "--len", read.length,
             "--in", reference.genome,
             "--id", read.ID,
             "--out", out.prefix)
system(cmd)   ### running art

cmd <- paste("samtools fastq",
             "-1", str_c(out.prefix, "1.fq"),
             "-2", str_c(out.prefix, "2.fq"),
             str_c(out.prefix, ".sam"))
system(cmd)


#processes 1082350 reads


#Now are going to re-run spades and quast with the new reference
#We start with the spades

library(RLinuxModules)
moduleInit(modulesHome = "/local/genome/Modules/3.2.10")
module("load spades")



K.range <- "21,33,55,77,99,127"
threads <- 4
R1.file <- "/net/fs-1/home01/student10/BIN310/1_Assignment/feilfire_reads_R1.fq"      # This is the faultless reads R1.
R2.file <- "/net/fs-1/home01/student10/BIN310/1_Assignment/feilfire_reads_R2.fq"      #This is the faultless reads R2.                            #
out.dir <- "/net/fs-1/home01/student10/BIN310/1_Assignment/errorfree_rerun/spades_errorfree"   # output


### Command line
cmd <- paste("spades.py",
             "-k", K.range,
             "--careful",
             "--threads", threads,
             "--pe1-1", R1.file,
             "--pe1-2", R2.file,
             "-o", out.dir)
system(cmd)


# Quast

library(RLinuxModules)
moduleInit(modulesHome="/local/genome/Modules/3.2.10")
module("load quast")

### Settings - same for all cases
ref.file <- "/mnt/users/student10/BIN310/1_Assignment/reference_new"
contigs.file <- "/net/fs-1/home01/student10/BIN310/1_Assignment/errorfree_rerun/spades_errorfree/contigs.fasta"
out.dir <- "/net/fs-1/home01/student10/BIN310/1_Assignment/errorfree_rerun/quast_errorfree/"
threads <- 4


cmd <- paste("quast.py",
             "--threads", threads,
             "--pe1", R1.file,         # comment out if evaluating canu+nanopore
             "--pe2", R2.file,         # comment out if evaluating canu+nanopore
             "-r", ref.file,
             "-o", out.dir,
             contigs.file)
system(cmd)