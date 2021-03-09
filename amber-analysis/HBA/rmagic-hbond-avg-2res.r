## Run this with "Rscript rmagic-hbond-avg-2res.r"
## (Assuming you've already installed R...)

#-----------------------------------------------#
#--Specify the paths to the Files from cpptraj--#
#-----------------------------------------------#

## This script has been pre-built for a system with 3 replicates
## More or less than 3 reps (up to 5) can be achieved through
## Changing the `sets` number below.
## For more than 5 reps, add in an `else if` statement in the `Read in`
## block below.

## Paths to the Hbond-avg files
## Set A (system 1)
infile1A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-1/WT_protein_system_hbond_avg.dat")
infile2A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-2/WT_protein_system_hbond_avg.dat")
infile3A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-3/WT_protein_system_hbond_avg.dat")
#infile4A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-4/WT_protein_system_hbond_avg.dat")
#infile5A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-5/WT_protein_system_hbond_avg.dat")

#-----------------------------#
#--Define your outfile names--#
#-----------------------------#

## A is for infiles labeled A
## Each system gets an averaged file
## That is then used with the H bond analysis script

A_avg <- "/absolute/path/to/the/analysis/files/for/HBA/WT-System_total_hbond_avg_double.dat"

## Explicitly set the number of the SNP/changed residue to match
## fix_string should be the number you need to change (Format: _XXX)
## fixed_string should be the thing you see in the output
fix_stringA <- "_001"
fixed_stringA <-"Changed_ResA"

fix_stringB <- "_002"
fixed_stringB <-"Changed_ResB"

## How many data sets to add (use 4 after decimal)
sets <- 3.0000
#sets <- 5.0000

#----------------------------------------------------------------------#
#---------Behind the Curtain: No Need to Modify Past This Line---------#
#----------------------------------------------------------------------#

## Use the data tables package to read in data frames
## Remove comment to install locally
#install.packages("data.table")
library(data.table)

## Use the tidyverse package to perform string replacement
## Remove comment to install locally
#install.packages("tidyverse")
library(tidyverse)

## Turn off scientific notation
options(scipen = 999)

#-------------------------#
#--Read in Hbond Scripts--#
#-------------------------#

## Deal with files if X present
if (sets <= 2) { # Case if 1 or 2 files
  stop(sprintf("You gave me %s file(s), which is less than 3 input files.
       This is not the script for you because it does standard deviations.
       Stopping now...", sets))
} else if (sets == 3) { # Case of 3 files
  ## Reading each file as a data.table.
  ## Bonus - fread is much faster than read.csv
  read1A <- fread(infile1A, header=TRUE)
  read2A <- fread(infile2A, header=TRUE)
  read3A <- fread(infile3A, header=TRUE)
  ## Force name the columns for ease of access
  colnames(read1A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read2A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read3A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  ## Combine all the datasets into 1
  bound <- rbind(read1A, read2A, read3A)
} else if (sets == 4) { # Case of 4 files
  read1A <- fread(infile1A, header=TRUE)
  read2A <- fread(infile2A, header=TRUE)
  read3A <- fread(infile3A, header=TRUE)
  read4A <- fread(infile4A, header=TRUE)
  colnames(read1A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read2A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read3A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read4A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  bound <- rbind(read1A, read2A, read3A, read4A)
} else if (sets == 5) { # Case of 5 files
  read1A <- fread(infile1A, header=TRUE)
  read2A <- fread(infile2A, header=TRUE)
  read3A <- fread(infile3A, header=TRUE)
  read4A <- fread(infile4A, header=TRUE)
  read5A <- fread(infile5A, header=TRUE)
  colnames(read1A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read2A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read3A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read4A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  colnames(read5A) <- c("Acceptor", "DonorH", "Donor", "Frames", "Frac", "AvgDist", "AvgAng")
  bound <- rbind(read1A, read2A, read3A, read4A, read5A)
} else { # > 5 files
  stop(sprintf("Sorry, this is for 3-5 files. You gave me %s.
       You'll have to modify the if statements for the file reads yourself.
       Stopping now...", sets))
}

bound$Acceptor <- as.character(bound$Acceptor)
bound$DonorH <- as.character(bound$DonorH)
bound$Donor <- as.character(bound$Donor)
bound$Frac <- as.numeric(bound$Frac)
bound$Frames <- as.numeric(bound$Frames)
bound$AvgDist <- as.numeric(bound$AvgDist)
bound$AvgAng <- as.numeric(bound$AvgAng)

## Collapse repeat lines into themselves (i.e. add numbers together)
superbound <- aggregate(data=bound, cbind(Frames,Frac,AvgDist,AvgAng)~Acceptor+Donor, FUN=sum)

## If DonorH matters, then use:
#superbound2 <- aggregate(data=bound, cbind(Frames,Frac,AvgDist,AvgAng)~Acceptor+DonorH+Donor, FUN=sum)

## Get average based on number of sets combined [This if for 3]
superbound$AvgFrame <- format(superbound$Frames / sets, digits=4, format="f")
superbound$AvgFrac <- format(superbound$Frac / sets, digits=4, format="f")
superbound$AAvgDist <- format(superbound$AvgDist / sets, digits=4, format="f")
superbound$AAvgAng <- format(superbound$AvgAng / sets, digits=4, format="f")

save_cols_AH <- superbound[,c("Acceptor", "Donor", "AvgFrac")]

## Change the first SNP or Varied position name in Acceptor with something checkable
findmeAccA <- grep(fix_stringA, save_cols_AH$Acceptor)
save_cols_AH$Acceptor <- replace(save_cols_AH$Acceptor, findmeAccA, fixed_stringA)

## Change the first SNP or Varied position name in Donor with something checkable
findmeDonA <- grep(fix_stringA, save_cols_AH$Donor)
save_cols_AH$Donor <- replace(save_cols_AH$Donor, findmeDonA, fixed_stringA)

## Change the second SNP or Varied position name in Acceptor with something checkable
findmeAccB <- grep(fix_stringB, save_cols_AH$Acceptor)
save_cols_AH$Acceptor <- replace(save_cols_AH$Acceptor, findmeAccB, fixed_stringB)

## Change the second SNP or Varied position name in Donor with something checkable
findmeDonB <- grep(fix_stringB, save_cols_AH$Donor)
save_cols_AH$Donor <- replace(save_cols_AH$Donor, findmeDonB, fixed_stringB)

## Limit to 4 sig figs after decimal
save_cols_clean_AH <- format(save_cols_AH, digits=4)

colnames(save_cols_clean_AH) <- c("Acceptor", "Donor     ", "AvgFrac")

## Now write a tab-delimited outfile!
## Don't care about the index rownames because that's the residue number
#write.table(save_cols_clean_AH, file = A_avg, sep="\t", row.names=FALSE, quote=FALSE)

## Write a pseudo-fixed width outfile --> this keeps index rownames
capture.output( print(save_cols_clean_AH, print.gap=3), file = A_avg)
