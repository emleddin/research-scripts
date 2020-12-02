## Run this with "Rscript system-hbond-table.r"
## (Assuming you've already installed R...)

#----------------------------------------------------------#
#--Specify the paths to the Files from rmagic-hbond-avg.r--#
#----------------------------------------------------------#

## This script has been pre-built for a comparing 3 separate systems
## More or less than 3 reps (up to 5) can be achieved through
## Commenting or uncommenting
## 1a with 1b | 2a with 2b | 3a with 3b | 4a with 4b | 5a with 5b

## Paths to the Hbond-avg files
## Set A (system 1)
infile1a <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System/WT_protein_system_total_hbond_avg.dat")
infile1b <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A/MUT_A_system_total_hbond_avg.dat")
infile2a <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System/WT_protein_system_total_hbond_avg.dat")
infile2b <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-B/MUT_B_system_total_hbond_avg.dat")
infile3a <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System/WT_protein_system_total_hbond_avg.dat")
infile3b <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-C/MUT_C_system_total_hbond_avg.dat")
#infile4a <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System/MUT_E_system_total_hbond_avg.dat")
#infile4b <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-D/MUT_E_system_total_hbond_avg.dat")
#infile5a <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System/WT_protein_system_total_hbond_avg.dat")
#infile5b <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-E/MUT_E_system_total_hbond_avg.dat")

#-----------------------------#
#--Define your outfile names--#
#-----------------------------#

## A is for infiles labeled A
## Each system gets an averaged file
## That is then used with the H bond analysis script

A_diff <- "/absolute/path/to/the/output/files/WT-A-HBA-table.dat"
B_diff <- "/absolute/path/to/the/output/files/WT-B-HBA-table.dat"
C_diff <- "/absolute/path/to/the/output/files/WT-C-HBA-table.dat"
#D_diff <- "/absolute/path/to/the/output/files/WT-D-HBA-table.dat"
#E_diff <- "/absolute/path/to/the/output/files/WT-E-HBA-table.dat"
Tot_diff <- "/absolute/path/to/the/output/files/WT-ABCDE-HBA-table.dat"

## How many data pairs are there (a/b couples)? (use 1 after decimal)
## You can use 1.0, 2.0, 3.0, 4.0 or 5.0
sets <- 3.0
#sets <- 5.0

## For column names, limit A & B to 7 char
## Limit C to 10 char

## Syst 1 Column names:
## tag1a and tag1b (AbsFrac of each system)
## tag1c (AbsDiff of tag1a and tag1b)
tag1a = "AF_WT"
tag1b = "AF_MUTA"
tag1c = "AD_WT_MUTA"

## Syst 2 Column names:
## tag1a and tag1b (AbsFrac of each system)
## tag1c (AbsDiff of tag1a and tag1b)
tag2a = "AF_WT"
tag2b = "AF_MUTB"
tag2c = "AD_WT_MUTB"

## Syst 3 Column names:
## tag3a and tag3b (AbsFrac of each system)
## tag3c (AbsDiff of tag1a and tag1b)
tag3a = "AF_WT"
tag3b = "AF_MUTC"
tag3c = "AD_WT_MUTC"

## Syst 4 Column names:
## tag4a and tag4b (AbsFrac of each system)
## tag4c (AbsDiff of tag1a and tag1b)
#tag4a = "AF_WT"
#tag4b = "AF_MUTD"
#tag4c = "AD_WT_MUTD"

## Syst 5 Column names:
## tag5a and tag5b (AbsFrac of each system)
## tag5c (AbsDiff of tag1a and tag1b)
#tag5a = "AF_WT"
#tag5b = "AF_MUTE"
#tag5c = "AD_WT_MUTE"

## Difference cutoff (ex: 20% different)
cutoff = 20

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

## Reading each file as a data.table.
## Bonus - fread is much faster than read.csv
read1a <- fread(infile1a, header=TRUE)
read1b <- fread(infile1b, header=TRUE)
read2a <- fread(infile2a, header=TRUE)
read2b <- fread(infile2b, header=TRUE)
read3a <- fread(infile3a, header=TRUE)
read3b <- fread(infile3b, header=TRUE)
#read4a <- fread(infile4a, header=TRUE)
#read4b <- fread(infile4b, header=TRUE)
#read5a <- fread(infile5a, header=TRUE)
#read5b <- fread(infile5b, header=TRUE)

colnames(read1a) <- c("Index", "Acceptor", "Donor", "AvgFrac")
colnames(read1b) <- c("Index", "Acceptor", "Donor", "AvgFrac")
colnames(read2a) <- c("Index", "Acceptor", "Donor", "AvgFrac")
colnames(read2b) <- c("Index", "Acceptor", "Donor", "AvgFrac")
colnames(read3a) <- c("Index", "Acceptor", "Donor", "AvgFrac")
colnames(read3b) <- c("Index", "Acceptor", "Donor", "AvgFrac")
#colnames(read4a) <- c("Index", "Acceptor", "Donor", "AvgFrac")
#colnames(read4b) <- c("Index", "Acceptor", "Donor", "AvgFrac")
#colnames(read5a) <- c("Index", "Acceptor", "Donor", "AvgFrac")
#colnames(read5b) <- c("Index", "Acceptor", "Donor", "AvgFrac")

## Drop the Index column
read1a <- within(read1a, rm(Index))
read1b <- within(read1b, rm(Index))

read2a <- within(read2a, rm(Index))
read2b <- within(read2b, rm(Index))

read3a <- within(read3a, rm(Index))
read3b <- within(read3b, rm(Index))

#read4a <- within(read4a, rm(Index))
#read4b <- within(read4b, rm(Index))

#read5a <- within(read5a, rm(Index))
#read5b <- within(read5b, rm(Index))

## Rename the Average Frac Columns
read1a <- plyr:::rename(read1a, c("AvgFrac"=tag1a))
read1b <- plyr:::rename(read1b, c("AvgFrac"=tag1b))

read2a <- plyr:::rename(read2a, c("AvgFrac"=tag2a))
read2b <- plyr:::rename(read2b, c("AvgFrac"=tag2b))

read3a <- plyr:::rename(read3a, c("AvgFrac"=tag3a))
read3b <- plyr:::rename(read3b, c("AvgFrac"=tag3b))

#read4a <- plyr:::rename(read4a, c("AvgFrac"=tag4a))
#read4b <- plyr:::rename(read4b, c("AvgFrac"=tag4b))

#read5a <- plyr:::rename(read5a, c("AvgFrac"=tag5a))
#read5b <- plyr:::rename(read5b, c("AvgFrac"=tag5b))

## Merge the tables together
tab1ab <- merge(read1a, read1b, by=c("Acceptor","Donor"))
tab2ab <- merge(read2a, read2b, by=c("Acceptor","Donor"))
tab3ab <- merge(read3a, read3b, by=c("Acceptor","Donor"))
#tab4ab <- merge(read4a, read4b, by=c("Acceptor","Donor"))
#tab5ab <- merge(read5a, read4b, by=c("Acceptor","Donor"))


## Treat the AvgFrac columns numerically
tab1ab[[tag1a]] <- as.numeric(as.character(tab1ab[[tag1a]]))
tab1ab[[tag1b]] <- as.numeric(as.character(tab1ab[[tag1b]]))

tab2ab[[tag2a]] <- as.numeric(as.character(tab2ab[[tag2a]]))
tab2ab[[tag2b]] <- as.numeric(as.character(tab2ab[[tag2b]]))

tab3ab[[tag3a]] <- as.numeric(as.character(tab3ab[[tag3a]]))
tab3ab[[tag3b]] <- as.numeric(as.character(tab3ab[[tag3b]]))

#tab4ab[[tag4a]] <- as.numeric(as.character(tab4ab[[tag4a]]))
#tab4ab[[tag4b]] <- as.numeric(as.character(tab4ab[[tag4b]]))

#tab5ab[[tag5a]] <- as.numeric(as.character(tab5ab[[tag5a]]))
#tab5ab[[tag5b]] <- as.numeric(as.character(tab5ab[[tag5b]]))

## Explicitly set all NA values as "0"
tab1ab[is.na(tab1ab)] <- 0
tab2ab[is.na(tab2ab)] <- 0
tab3ab[is.na(tab3ab)] <- 0
#tab4ab[is.na(tab4ab)] <- 0
#tab5ab[is.na(tab5ab)] <- 0


## Get the absolute difference between columns
tab1ab[[tag1c]] <- abs(tab1ab[[tag1a]] - tab1ab[[tag1b]])
tab2ab[[tag2c]] <- abs(tab2ab[[tag2a]] - tab2ab[[tag2b]])
tab3ab[[tag3c]] <- abs(tab3ab[[tag3a]] - tab3ab[[tag3b]])
#tab4ab[[tag4c]] <- abs(tab4ab[[tag4a]] - tab4ab[[tag4b]])
#tab5ab[[tag5c]] <- abs(tab5ab[[tag5a]] - tab5ab[[tag5b]])


## Change to percentages
tab1ab[[tag1a]] <- tab1ab[[tag1a]]*100
tab1ab[[tag1b]] <- tab1ab[[tag1b]]*100
tab1ab[[tag1c]] <- tab1ab[[tag1c]]*100

tab2ab[[tag2a]] <- tab2ab[[tag2a]]*100
tab2ab[[tag2b]] <- tab2ab[[tag2b]]*100
tab2ab[[tag2c]] <- tab2ab[[tag2c]]*100

tab3ab[[tag3a]] <- tab3ab[[tag3a]]*100
tab3ab[[tag3b]] <- tab3ab[[tag3b]]*100
tab3ab[[tag3c]] <- tab3ab[[tag3c]]*100

#tab4ab[[tag4a]] <- tab4ab[[tag4a]]*100
#tab4ab[[tag4b]] <- tab4ab[[tag4b]]*100
#tab4ab[[tag4c]] <- tab4ab[[tag4c]]*100

#tab5ab[[tag5a]] <- tab5ab[[tag5a]]*100
#tab5ab[[tag5b]] <- tab5ab[[tag5b]]*100
#tab5ab[[tag5c]] <- tab5ab[[tag5c]]*100


## Get values > cutoff
## UQ removes the string quotes
tab1ab_cut <- filter(tab1ab, UQ(as.name(tag1c)) > cutoff)
tab2ab_cut <- filter(tab2ab, UQ(as.name(tag2c)) > cutoff)
tab3ab_cut <- filter(tab3ab, UQ(as.name(tag3c)) > cutoff)
#tab4ab_cut <- filter(tab4ab, UQ(as.name(tag4c)) > cutoff)
#tab5ab_cut <- filter(tab5ab, UQ(as.name(tag5c)) > cutoff)


## Limit to 4 sig figs after decimal
tab1ab_cleancut <- format(tab1ab_cut, digits=4)
tab2ab_cleancut <- format(tab2ab_cut, digits=4)
tab3ab_cleancut <- format(tab3ab_cut, digits=4)
#tab4ab_cleancut <- format(tab4ab_cut, digits=4)
#tab5ab_cleancut <- format(tab5ab_cut, digits=4)


## Now write a tab-delimited outfile!
## Don't care about the index rownames because that's the residue number
#write.table(save_cols_clean_AH, file = A_avg, sep="\t", row.names=FALSE, quote=FALSE)

## Write a pseudo-fixed width outfile
#capture.output( print(tab1ab_cleancut, print.gap=3, row.names=FALSE), file = A_diff)
#capture.output( print(tab2ab_cleancut, print.gap=3, row.names=FALSE), file = B_diff)
#capture.output( print(tab3ab_cleancut, print.gap=3, row.names=FALSE), file = C_diff)
#capture.output( print(tab4ab_cleancut, print.gap=3, row.names=FALSE), file = D_diff)
#capture.output( print(tab5ab_cleancut, print.gap=3, row.names=FALSE), file = E_diff)


## Merge all the datasets, keeping any blank columns
#### Coment all this out if you're only doing 1 set ####

## For 2 sets
#tabcombo <- merge(tab1ab_cleancut, tab2ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)

## For 3 sets
tabcomboA <- merge(tab1ab_cleancut, tab2ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)
tabcombo <- merge(tabcomboA, tab3ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)

## For 4 sets
#tabcomboA <- merge(tab1ab_cleancut, tab2ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)
#tabcomboB <- merge(tabcomboA, tab3ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)
#tabcombo <- merge(tabcomboB, tab4ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)

## For 5 sets
#tabcomboA <- merge(tab1ab_cleancut, tab2ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)
#tabcomboB <- merge(tabcomboA, tab3ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)
#tabcomboC <- merge(tabcomboB, tab4ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)
#tabcombo <- merge(tabcomboc, tab5ab_cleancut, by=c("Acceptor","Donor"), all=TRUE)

## Change the <NA> values to dashes
tabcombo[is.na(tabcombo)] <- "-"

#### End of block to comment out for skipping if 1 ####

## Write a pseudo-fixed width outfile
#capture.output( print(tabcombo, print.gap=3, row.names=FALSE), file = Tot_diff)
write.table(tabcombo, file = Tot_diff, sep="\t", row.names=FALSE, quote=FALSE)


## Determine How to format the columns with awk
## You've defined sets above, this selects the right terminology
## This section defines the spacing for printing columns
if (sets == 1.0){
thingA <- "%-9s %-9s %-12s"
} else if (sets == 2.0){
thingA <- "%-9s %-9s %-12s %-9s %-9s %-12s"
} else if (sets == 3.0){
thingA <- "%-9s %-9s %-12s %-9s %-9s %-12s %-9s %-9s %-12s"
} else if (sets == 4.0){
thingA <- "%-9s %-9s %-12s %-9s %-9s %-12s %-9s %-9s %-12s %-9s %-9s %-12s"
} else {
thingA <- "%-9s %-9s %-12s %-9s %-9s %-12s %-9s %-9s %-12s %-9s %-9s %-12s %-9s %-9s %-12s"
}

## This section determines the number of columns to be printed
if (sets == 1.0){
thingB <- ""
} else if (sets == 2.0){
thingB <- ", $6, $7, $8"
} else if (sets == 3.0){
thingB <- ", $6, $7, $8, $9, $10, $11"
} else if (sets == 4.0){
thingB <- ", $6, $7, $8, $9, $10, $11, $12, $13, $14"
} else {
thingB <- ", $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17"
}

## Format the file into nice columns with awk
## Note: you must escape double quotes (\")and escape the escape for newline (\\n)
## The appropriate number of %-8s is dependent on your data (this is for 2 sets)
## Acceptor, Donor, Read1a, Read1b, Tab1ab, Read2a, Read2b, Tab2ab
#system(paste("awk '{printf \"%-14s %-14s %-8s %-8s %-8s %-8s %-8s %-8s\\n\", $1, $2, $3, $4, $5, $6, $7, $8}'", Tot_diff, "> attempt.dat"), intern=TRUE)
system(paste("awk '{printf \"%-14s %-14s", as.name(thingA), "\\n\", $1, $2, $3, $4, $5", as.name(thingB), "}'", Tot_diff, "> attempt.dat"), intern=TRUE)
## Based on awk '{printf "%-14s %-14s %-8s %-8s %-8s %-8s %-8s %-8s\n", $1, $2, $3, $4, $5, $6, $7, $8}' thing.dat
