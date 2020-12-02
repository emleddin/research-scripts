## Run this with "Rscript byres-processing.r"
## (Assuming you've already installed R...)

## cpptraj rmsd byres data comes as each residue with its own column
## with rows as the RMSD value per frame.
## This script gets an average RMSD by time for that residue, resulting in
## an output file with 2 columns: residue number (x) and RMSD (y)

#--------------------------------------------#
#--Specify the paths to the Files from RMSD--#
#--------------------------------------------#

## This script is for a single file

## Path to the rmsd_byres.dat file
## You're not *required* to use the absolute path (and too long a path will
## cause issues...)
infile1Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile2Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile3Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile4Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile5Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile6Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile7Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile8Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile9Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")
#infile10Rb <- Sys.glob("/absolute/path/to/the/analysis/files/for/rmsd_byres.dat")

## How many files are you doing this for? Assumes 1 by default
## Can do up to 10 as written
num_files = 1

#-----------------------------#
#--Define your outfile names--#
#-----------------------------#

## DO NOT use sys.glob here, it will NOT work!!!
## You're not *required* to use the absolute path (and too long a path will
## cause issues...)
RMSD_byres1 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres2 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres3 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres4 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres5 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres6 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres7 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres8 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres9 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"
#RMSD_byres10 <- "/absolute/path/to/the/analysis/files/for/rmsd_byres_timeavg.dat"

#----------------------------------------------------------------------#
#---------Behind the Curtain: No Need to Modify Past This Line---------#
#----------------------------------------------------------------------#

## Use the data tables package to read in data frames
## Remove comment to install locally
#install.packages("data.table")
library(data.table)

## Use the dplyr and tidyr packages (in tidyverse) to do col means
## Hence, you can either load tidyverse, or load dplyr and tidyr
## Remove comment to install locally
#install.packages("dplyr")
#install.packages("tidyr")
library(dplyr)
library(tidyr)

## Alternatively...
#install.packages("tidyverse")
##library(tidyverse)

## Turn off scientific notation
options(scipen = 999)

#-------------------------#
#--Read in the data file--#
#-------------------------#

## Reading each file as a data.table.
## Bonus - fread is much faster than read.csv
Read1Rb <- fread(infile1Rb, header=TRUE)

## read1Rb[[1]] refers to the Frame column

## Remove the Frame Column
NoFrame1Rb <- select(Read1Rb, -c("#Frame"))

## Get the mean of each residue (dplyr)
Avg1Rb <- summarise_all(NoFrame1Rb, mean)

## Reshape the data table (tidyr)
Shaped1Rb <- gather(Avg1Rb, key="Residue", value="Avg")

## Keep residue number, not name ["SER:1234" becomes just "1234"]
Shaped1Rb <- Shaped1Rb %>% mutate(Residue = gsub(".*:","",Residue))

## Write a whitespace-delimited outfile!
sink(RMSD_byres1, type=c("output"))
print(Shaped1Rb, row.names=FALSE)
sink()

## Deal with files 2-10 if present
if (num_files > 1) {
## Remove the previous data
rm(Read1Rb, NoFrame1Rb, Avg1Rb, Shaped1Rb)
Read2Rb <- fread(infile2Rb, header=TRUE)
NoFrame2Rb <- select(Read2Rb, -c("#Frame"))
Avg2Rb <- summarise_all(NoFrame2Rb, mean)
Shaped2Rb <- gather(Avg2Rb, key="Residue", value="Avg")
Shaped2Rb <- Shaped2Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres2, type=c("output"))
print(Shaped2Rb, row.names=FALSE)
sink()
}

if (num_files > 2) {
rm(Read2Rb, NoFrame2Rb, Avg2Rb, Shaped2Rb)
Read3Rb <- fread(infile3Rb, header=TRUE)
NoFrame3Rb <- select(Read3Rb, -c("#Frame"))
Avg3Rb <- summarise_all(NoFrame3Rb, mean)
Shaped3Rb <- gather(Avg3Rb, key="Residue", value="Avg")
Shaped3Rb <- Shaped3Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres3, type=c("output"))
print(Shaped3Rb, row.names=FALSE)
sink()
}

if (num_files > 3) {
rm(Read3Rb, NoFrame3Rb, Avg3Rb, Shaped3Rb)
Read4Rb <- fread(infile4Rb, header=TRUE)
NoFrame4Rb <- select(Read4Rb, -c("#Frame"))
Avg4Rb <- summarise_all(NoFrame4Rb, mean)
Shaped4Rb <- gather(Avg4Rb, key="Residue", value="Avg")
Shaped4Rb <- Shaped4Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres4, type=c("output"))
print(Shaped4Rb, row.names=FALSE)
sink()
}

if (num_files > 4) {
rm(Read4Rb, NoFrame4Rb, Avg4Rb, Shaped4Rb)
Read5Rb <- fread(infile5Rb, header=TRUE)
NoFrame5Rb <- select(Read5Rb, -c("#Frame"))
Avg5Rb <- summarise_all(NoFrame5Rb, mean)
Shaped5Rb <- gather(Avg5Rb, key="Residue", value="Avg")
Shaped5Rb <- Shaped5Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres5, type=c("output"))
print(Shaped5Rb, row.names=FALSE)
sink()
}

if (num_files > 5) {
rm(Read5Rb, NoFrame5Rb, Avg5Rb, Shaped5Rb)
Read6Rb <- fread(infile6Rb, header=TRUE)
NoFrame6Rb <- select(Read6Rb, -c("#Frame"))
Avg6Rb <- summarise_all(NoFrame6Rb, mean)
Shaped6Rb <- gather(Avg6Rb, key="Residue", value="Avg")
Shaped6Rb <- Shaped6Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres6, type=c("output"))
print(Shaped6Rb, row.names=FALSE)
sink()
}

if (num_files > 6) {
rm(Read6Rb, NoFrame6Rb, Avg6Rb, Shaped6Rb)
Read7Rb <- fread(infile7Rb, header=TRUE)
NoFrame7Rb <- select(Read7Rb, -c("#Frame"))
Avg7Rb <- summarise_all(NoFrame7Rb, mean)
Shaped7Rb <- gather(Avg7Rb, key="Residue", value="Avg")
Shaped7Rb <- Shaped7Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres7, type=c("output"))
print(Shaped7Rb, row.names=FALSE)
sink()
}

if (num_files > 7) {
rm(Read7Rb, NoFrame7Rb, Avg7Rb, Shaped7Rb)
Read8Rb <- fread(infile8Rb, header=TRUE)
NoFrame8Rb <- select(Read8Rb, -c("#Frame"))
Avg8Rb <- summarise_all(NoFrame8Rb, mean)
Shaped8Rb <- gather(Avg8Rb, key="Residue", value="Avg")
Shaped8Rb <- Shaped8Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres8, type=c("output"))
print(Shaped8Rb, row.names=FALSE)
sink()
}

if (num_files > 8) {
rm(Read8Rb, NoFrame8Rb, Avg8Rb, Shaped8Rb)
Read9Rb <- fread(infile9Rb, header=TRUE)
NoFrame9Rb <- select(Read9Rb, -c("#Frame"))
Avg9Rb <- summarise_all(NoFrame9Rb, mean)
Shaped9Rb <- gather(Avg9Rb, key="Residue", value="Avg")
Shaped9Rb <- Shaped9Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres9, type=c("output"))
print(Shaped9Rb, row.names=FALSE)
sink()
}

if (num_files > 9) {
rm(Read9Rb, NoFrame9Rb, Avg9Rb, Shaped9Rb)
Read10Rb <- fread(infile10Rb, header=TRUE)
NoFrame10Rb <- select(Read10Rb, -c("#Frame"))
Avg10Rb <- summarise_all(NoFrame10Rb, mean)
Shaped10Rb <- gather(Avg10Rb, key="Residue", value="Avg")
Shaped10Rb <- Shaped10Rb %>% mutate(Residue = gsub(".*:","",Residue))
sink(RMSD_byres10, type=c("output"))
print(Shaped10Rb, row.names=FALSE)
sink()
}
