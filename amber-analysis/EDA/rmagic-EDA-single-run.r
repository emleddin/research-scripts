## Run this with "Rscript rmagic-EDA-avg.r"
## (Assuming you've already installed R...)

#-------------------------------------------#
#--Specify the paths to the Files from EDA--#
#-------------------------------------------#

## This script has been pre-built for a system with 3 replicates
## More or less than 3 reps (up to 5) can be achieved through
## Commenting or uncommenting

## Paths to the fort.803 (Coul) files
## Set A (system 1)
infile1Ac <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-1/fort.803")

## Paths to the fort.806 (VdW) files
## Set A (system A)
infile1Av <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-1/fort.806")

#-----------------------------#
#--Define your outfile names--#
#-----------------------------#

## A is for infiles labeled A
## Each system gets an averaged file
## Have one for Coulomb, one for vdW, and one for Coul+vdW (total)

A_coul <- "/absolute/path/to/the/output/WT_protein_system_EDA_resX_coul_avg.dat"
A_vdw <- "/absolute/path/to/the/output/WT_protein_system_EDA_resX_vdw_avg.dat"
A_tot <- "/absolute/path/to/the/output/WT_protein_system_EDA_resX_tot_avg.dat"

## Residue of interest (A matched with B, which matches do you care about?) (use 4 after decimal)
## Ex. If mutant is residue 100, this would be 100
## This script will remove the matches directly surrounding ROI for you
## Which is good because by not being just Coul and vdW, they're too dominant
ROI <- 100

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

#----------------------------#
#--Read in Coul EDA Scripts--#
#----------------------------#

## First line of file is number of frames used for EDA
## This is skipped by R's fread by default to avoid
## Irregular header information

## Reading each file as a data.table.
## Bonus - fread is much faster than read.csv
read1Ac <- fread(infile1Ac, header=FALSE)

colnames(read1Ac) <- c("Index", "ResidueA", "ResidueB", "Coulomb", "StdErr")

## Combine all the datasets into 1
bound <- read1Ac

## Add in a blank row of the match for future plotting needs
extra <- data.frame(0, ROI, ROI, 0, 0)
bound <- rbind(bound, setNames(extra, names(read1Ac)))

#bound$Index <- as.numeric(bound$Index)
bound$Index <- as.numeric(bound$Index)
bound$ResidueA <- as.numeric(bound$ResidueA)
bound$ResidueB <- as.numeric(bound$ResidueB)
bound$Coulomb <- as.numeric(bound$Coulomb)
bound$StdErr <- as.numeric(bound$StdErr)

## Sort with that new zero row
bound <- bound[order(ResidueB),]

save_cols_Ac <- bound[,c("ResidueA", "ResidueB", "Coulomb")]

only_ROI_rows_Ac <- filter(save_cols_Ac, ResidueA == ROI | ResidueB == ROI)

## Change the NA from standard deviation to 0
only_ROI_rows_Ac[(ROI),4] <- 0

## Create a copy of the parsed data to format
clean_rows_Ac <- data.frame(only_ROI_rows_Ac)

## Pattern searching converted it to a character string, so back to numeric
clean_rows_Ac$Coulomb <- as.numeric(clean_rows_Ac$Coulomb)

## Limit to 4 sig figs after decimal
clean_rows_Ac$Coulomb <- formatC(clean_rows_Ac$Coulomb, digits=4, format="f")

## Set the two residues surrounding the ROI to zero
## This is because energy is overpowering due to other energy terms
## So if ROI=100, you remove matches between 99 & 100 as well as 100 & 101
if (ROI != 1) {clean_rows_Ac[(ROI-1),3] <- 0
}
clean_rows_Ac[(ROI+1),3] <- 0

#------------------------------------------------------------------------#
#-----------------------------COUL OUTFILE-------------------------------#
#------------------------------------------------------------------------#

## Now write a tab-delimited outfile!
## Don't care about the index rownames
#write.table(clean_rows_Ac, file = A_coul, sep="\t", row.names=FALSE, quote=FALSE)

## Write a whitespace-delimited outfile!
sink(A_coul, type=c("output"))
print(clean_rows_Ac, row.names=FALSE)
sink()

#---------------------------#
#--Read in VDW EDA Scripts--#
#---------------------------#

## First line of file is number of frames used for EDA
## This is skipped by R's fread by default to avoid
## Irregular header information

## Reading each file as a data.table.
## Bonus - fread is much faster than read.csv
read1Av <- fread(infile1Av, header=FALSE)

colnames(read1Av) <- c("Index", "ResidueA", "ResidueB", "VdW", "StdErr")

## Combine all the datasets into 1
bound <- read1Av

## Add in a blank row of the match for future plotting needs
extra <- data.frame(0, ROI, ROI, 0, 0)
bound <- rbind(bound, setNames(extra, names(read1Av)))

## Sort with that new zero row
bound <- bound[order(ResidueB),]

#bound$Index <- as.numeric(bound$Index)
bound$Index <- as.numeric(bound$Index)
bound$ResidueA <- as.numeric(bound$ResidueA)
bound$ResidueB <- as.numeric(bound$ResidueB)
bound$VdW <- as.numeric(bound$VdW)
bound$StdErr <- as.numeric(bound$StdErr)


save_cols_Av <- bound[,c("ResidueA", "ResidueB", "VdW")]

only_ROI_rows_Av <- filter(save_cols_Av, ResidueA == ROI | ResidueB == ROI)

## Change the NA from standard deviation to 0
only_ROI_rows_Av[(ROI),4] <- 0

## Create a copy of the parsed data to format
clean_rows_Av <- data.frame(only_ROI_rows_Av)

## Pattern searching converted it to a character string, so back to numeric
clean_rows_Av$VdW <- as.numeric(clean_rows_Av$VdW)

## Limit to 4 sig figs after decimal
clean_rows_Av$VdW <- formatC(clean_rows_Av$VdW, digits=4, format="f")

## Set the two residues surrounding the ROI to zero
## This is because energy is overpowering due to other energy terms
## So if ROI=100, you remove matches between 99 & 100 as well as 100 & 101
if (ROI != 1) {clean_rows_Av[(ROI-1),3] <- 0
}
clean_rows_Av[(ROI+1),3] <- 0

#------------------------------------------------------------------------#
#-----------------------------VDW OUTFILES-------------------------------#
#------------------------------------------------------------------------#

## Now write a tab-delimited outfile!
## Don't care about the index rownames
#write.table(clean_rows_Av, file = A_vdw, sep="\t", row.names=FALSE, quote=FALSE)

## Write a whitespace-delimited outfile!
sink(A_vdw, type=c("output"))
print(clean_rows_Av, row.names=FALSE)
sink()

#---------------------------------------#
#--Create the TOTAL (Coul + vdW) files--#
#---------------------------------------#

## Combine into one dataset
## Use all columns from _Ac and the VdW and VdWSD columns from _Av
## Note: this makes it a matrix
combine_Acv = cbind(clean_rows_Ac, clean_rows_Av[,3:4])

## Formatting the rows converted it to a character string, so back to numeric again!
combine_Acv$Coulomb <- as.numeric(combine_Acv$Coulomb)
combine_Acv$VdW <- as.numeric(combine_Acv$VdW)

## Your data are now ResidueA ResidueB AvgCoul AvgCoulSD AvgVdW AvgVdWSD
## Append a column called AvgIntTot thats the sum of AvgCoul and AvgVdW
combine_Acv$IntTot <- (combine_Acv$Coulomb + combine_Acv$VdW)

# ## Create a new variable that's just ResidueA ResidueB AvgIntTot AvgStdDev
# save_cols_tot <- combine_Acv[,c("ResidueA", "ResidueB", "IntTot")]

## If the ResidueA column doesn't equal ROI, use ResidueA. Else, use ResidueB.
combine_Acv$Residue <- ifelse((combine_Acv$ResidueA != ROI), as.numeric(as.character(combine_Acv$ResidueA)), as.numeric(as.character(combine_Acv$ResidueB)))

## Create a new variable that's just Residue IntTot
save_cols_tot <- combine_Acv[,c("Residue", "IntTot")]

## Sanity Check!
## Set the two residues surrounding the ROI to zero
## This is because energy is overpowering due to other energy terms
## So if ROI=100, you remove matches between 99 & 100 as well as 100 & 101
if (ROI != 1) {save_cols_tot[(ROI-1),2] <- 0
}
save_cols_tot[(ROI+1),2] <- 0

#------------------------------------------------------------------------#
#----------------------TOTAL INTERACTION OUTFILES------------------------#
#------------------------------------------------------------------------#

## Now write a tab-delimited outfile!
## Don't care about the index rownames
#write.table(save_cols_tot, file = A_tot, sep="\t", row.names=FALSE, quote=FALSE)

## Write a whitespace-delimited outfile!
sink(A_tot, type=c("output"))
print(save_cols_tot, row.names=FALSE)
sink()
