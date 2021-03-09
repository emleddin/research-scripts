## Run this with "Rscript rmagic-EDA-avg-diffs.r"
## (Assuming you've already installed R...)

#--------------------------------------------------------------#
#-----Specify the paths to the Files from rmagic-EDA-avg.r-----#
#--------------------------------------------------------------#

## This script has been pre-built for 2 systems that have gone through 
## `rmagic-EDA-avg.r` (meaning there were replicates originally)

## Paths to the -tot- files
## Set A (system 1)
infileACV <- Sys.glob("/absolute/path/to/the/avgeraging/output/WT_protein_system_EDA_resX_tot_avg.dat")

##Set B (system 2)
infileBCV <- Sys.glob("/absolute/path/to/the/avgeraging/output/MUT_A_system_EDA_resX_tot_avg.dat")

#-----------------------------#
#--Define your outfile names--#
#-----------------------------#

## A - B
TOTAB <- "/absolute/path/to/the/difference/output/WT-MUTA_total_interaction_resX_avg.dat"

## This is X in coul-X
## Y and Z are X+1 and X-1
## Other scripts call this the ROI
X_val <- "100"

#----------------------------------------------------------------------#
#---------Behind the Curtain: No Need to Modify Past This Line---------#
#----------------------------------------------------------------------#

## Use the data tables package to read in data frames
## Remove comment to install locally
#install.packages("data.table")
library(data.table)

## Use the abind package to combine data frames
## Remove comment to install locally
#install.packages("abind")
library(abind)

## Turn off scientific notation
options(scipen = 999)

#-------------------#
#--Begin with COUL--#
#-------------------#

## Reading each file as a data.table.
## Bonus - fread is much faster than read.csv
combineACV <- fread(infileACV, header=TRUE)
colnames(combineACV) <- c("R1", "R2", "TotAvg", "TotStdD")

combineBCV <- fread(infileBCV, header=TRUE)
colnames(combineBCV) <- c("R1", "R2", "TotAvg", "TotStdD")

## Redefine as a data frame
combineACV <- as.data.frame(combineACV)

combineBCV <- as.data.frame(combineBCV)

## They're not numbers, so make them numbers
combineACV$TotAvg <- as.numeric(as.character(combineACV$TotAvg))
combineACV$TotStdD <- as.numeric(as.character(combineACV$TotStdD))

combineBCV$TotAvg <- as.numeric(as.character(combineBCV$TotAvg))
combineBCV$TotStdD <- as.numeric(as.character(combineBCV$TotStdD))

## Combine A res numbers, tot average, tot average, tot stdev, tot stev
combineTotCV <- abind(combineACV[,1:3], combineBCV[,3], combineACV[,4], combineBCV[,4], along=2)

## Rename the columns
colnames(combineTotCV) <- c("R1", "R2", "ATotalE", "BTotalE", "AAvgStd", "BAvgStd")

## Redefine as a data frame
combineTotCV <- as.data.frame(combineTotCV)

## If the R1 column doesn't equal X_val, use R1. Else, use R2.
combineTotCV$Residue <- ifelse((combineTotCV$R1 != X_val), as.numeric(as.character(combineTotCV$R1)), as.numeric(as.character(combineTotCV$R2)))

## They're not numbers, so make them numbers
combineTotCV$ATotalE <- as.numeric(as.character(combineTotCV$ATotalE))
combineTotCV$BTotalE <- as.numeric(as.character(combineTotCV$BTotalE))
combineTotCV$AAvgStd <- as.numeric(as.character(combineTotCV$AAvgStd))
combineTotCV$BAvgStd <- as.numeric(as.character(combineTotCV$BAvgStd))

## Multiply B * -1
## THIS WILL DO A - B!!
combineTotCV$BTotalE <- (combineTotCV$BTotalE*(-1.0000000000))
combineTotCV$DiffE <- rowSums(combineTotCV[, c("ATotalE", "BTotalE")])

## Get the Avg Stdev
combineTotCV$AvgSTDEV <- rowMeans(combineTotCV[,5:6])

## Create a new variable with just Residue, DiffE, and AvgSTDEV
save_cols_total_CV <- combineTotCV[,c("Residue", "DiffE", "AvgSTDEV")]

## Limit to 8 sig figs after decimal
save_cols_clean_total_CV <- format(save_cols_total_CV, digits=8)

## Explicitly remove the two residues matched next to the residue of interest
## This is because it's more than interaction energy (stuff like bond E too)
## (Note: | is the or operator)
#save_cols_clean_total_CV <- save_cols_clean_total_CV[!(save_cols_clean_total_CV$Residue == as.numeric(X_val)+1 | #save_cols_clean_total_CV$Residue == as.numeric(X_val)-1),]

#---------------------------------------------------------------------#
#--------------------------TOT OUTFILES-------------------------------#
#---------------------------------------------------------------------#

## Now write a tab-delimited outfile!
## Don't care about the index rownames because that's the frame
#write.table(save_cols_clean_total_CV, file = TOTAB, sep="\t", row.names=FALSE, quote=FALSE)

## Write a whitespace-delimited outfile!
sink(TOTAB, type=c("output"))
print(save_cols_clean_total_CV, row.names=FALSE)
sink()
