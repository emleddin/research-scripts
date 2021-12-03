## Run this with "Rscript rmagic.r"
## (Assuming you've already installed R...)

#-----------------------------------------------#
#--Specify the paths to the Files from cpptraj--#
#-----------------------------------------------#

## This script has been pre-built for 2 systems with 3 replicates
## More or less than 3 reps (up to 5) can be achieved through
## Commenting or uncommenting

## Paths to the RMSD files
## Set A (system 1)
infile1A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-1/WT_protein_system_total_bb_rms.dat")
infile2A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-2/WT_protein_system_total_bb_rms.dat")
infile3A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-3/WT_protein_system_total_bb_rms.dat")
#infile4A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-4/WT_protein_system_total_bb_rms.dat")
#infile5A <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-5/WT_protein_system_total_bb_rms.dat")

##Set B (system 2)
infile1B <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-1/MUT_A_system_total_bb_rms.dat")
infile2B <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-2/MUT_A_system_total_bb_rms.dat")
infile3B <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-3/MUT_A_system_total_bb_rms.dat")
#infile4B <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-4/MUT_A_system_total_bb_rms.dat")
#infile5B <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-5/MUT_A_system_total_bb_rms.dat")

## Paths to the RMSF files
## Set A (system 1)
infile1AF <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-1/WT_protein_system_rmsf_byres.dat")
infile2AF <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-2/WT_protein_system_rmsf_byres.dat")
infile3AF <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-3/WT_protein_system_rmsf_byres.dat")
#infile4AF <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-4/WT_protein_system_rmsf_byres.dat")
#infile5AF <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-5/WT_protein_system_rmsf_byres.dat")

infile1BF <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-1/MUT_A_system_rmsf_byres.dat")
infile2BF <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-2/MUT_A_system_rmsf_byres.dat")
infile3BF <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-3/MUT_A_system_rmsf_byres.dat")
#infile4BF <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-4/MUT_A_system_rmsf_byres.dat")
#infile5BF <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-5/MUT_A_system_rmsf_byres.dat")

## Paths to the H-Bond files
## Set A (system 1)
infile1AH <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-1/WT_protein_system_hbond.dat")
infile2AH <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-2/WT_protein_system_hbond.dat")
infile3AH <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-3/WT_protein_system_hbond.dat")
#infile4AH <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-4/WT_protein_system_hbond.dat")
#infile5AH <- Sys.glob("/absolute/path/to/the/analysis/files/for/WT-System-5/WT_protein_system_hbond.dat")

##Set B (system 2)
infile1BH <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-1/MUT_A_system_hbond.dat")
infile2BH <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-2/MUT_A_system_hbond.dat")
infile3BH <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-3/MUT_A_system_hbond.dat")
#infile4BH <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-4/MUT_A_system_hbond.dat")
#infile5BH <- Sys.glob("/absolute/path/to/the/analysis/files/for/MUT-A-System-5/MUT_A_system_hbond.dat")

## How many data sets to add (use 4 after decimal)
sets_rmsd  <- 3.0000
sets_rmsf  <- 3.0000
sets_hbond <- 3.0000
#sets_rmsd  <- 5.0000
#sets_rmsf  <- 5.0000
#sets_hbond <- 5.0000

#-----------------------------#
#--Define your outfile names--#
#-----------------------------#

## A is for infiles labeled A
## Each system gets an averaged file
## That is then used to plot with gnuplot
RMSDA <- "WT_protein_system_total_bb_rms_avgd.dat"
RMSDB <- "MUT_A_system_total_bb_rms_avgd.dat"

RMSFA <- "WT_protein_system_rmsf_byres_avgd.dat"
RMSFB <- "MUT_A_system_rmsf_byres_avgd.dat"

HBONDA <- "WT_protein_system_hbond_3trial_avgd.dat"
HBONDB <- "MUT_A_system_hbond_3trial_avgd.dat"

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

#-------------------#
#--Begin with RMSD--#
#-------------------#

## Deal with sets if X present
if (sets_rmsd <= 2) { # Case if 1 or 2 files
    stop(sprintf("You gave me %s file(s), which is less than 3 input files.
        This is not the script for you because it does standard deviations.
        Stopping now...", sets))
} else if (sets_rmsd == 3) { # Case of 3 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1A <- fread(infile1A, header=TRUE)
    read2A <- fread(infile2A, header=TRUE)
    read3A <- fread(infile3A, header=TRUE)
    ## Set 2
    read1B <- fread(infile1B, header=TRUE)
    read2B <- fread(infile2B, header=TRUE)
    read3B <- fread(infile3B, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1A and second column only for 2A and 3A
    ## Note: this makes it a matrix
    combineA = abind(read1A, read2A[,2], read3A[,2], along=2)
    combineB = abind(read1B, read2B[,2], read3B[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1A RMSD, infile2A RMSD, infile3A RMSD
    colnames(combineA) <- c("Frame", "RMSD1", "RMSD2", "RMSD3")
    colnames(combineB) <- c("Frame", "RMSD1", "RMSD2", "RMSD3")
    ## Redefine as a data frame
    combineA <- as.data.frame(combineA)
    combineB <- as.data.frame(combineB)
    ## Append a column that's the average of RMSD cols
    combineA$Average <- rowMeans(combineA[,2:4])
    combineB$Average <- rowMeans(combineB[,2:4])
    ## Append a column that's the STDEV of RMSD cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineA$STDEV <- apply(combineA[,2:4], 1, sd)
    combineB$STDEV <- apply(combineB[,2:4], 1, sd)
} else if (sets_rmsd == 4) { # Case of 4 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1A <- fread(infile1A, header=TRUE)
    read2A <- fread(infile2A, header=TRUE)
    read3A <- fread(infile3A, header=TRUE)
    read4A <- fread(infile4A, header=TRUE)
    ## Set 2
    read1B <- fread(infile1B, header=TRUE)
    read2B <- fread(infile2B, header=TRUE)
    read3B <- fread(infile3B, header=TRUE)
    read4B <- fread(infile4B, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1A and second column only for 2A and 3A
    ## Note: this makes it a matrix
    combineA = abind(read1A, read2A[,2], read3A[,2], read4A[,2], along=2)
    combineB = abind(read1B, read2B[,2], read3B[,2], read4B[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1A RMSD, infile2A RMSD, infile3A RMSD
    colnames(combineA) <- c("Frame", "RMSD1", "RMSD2", "RMSD3", "RMSD4")
    colnames(combineB) <- c("Frame", "RMSD1", "RMSD2", "RMSD3", "RMSD4")
    ## Redefine as a data frame
    combineA <- as.data.frame(combineA)
    combineB <- as.data.frame(combineB)
    ## Append a column that's the average of RMSD cols
    combineA$Average <- rowMeans(combineA[,2:5])
    combineB$Average <- rowMeans(combineB[,2:5])
    ## Append a column that's the STDEV of RMSD cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineA$STDEV <- apply(combineA[,2:5], 1, sd)
    combineB$STDEV <- apply(combineB[,2:5], 1, sd)
} else if (sets_rmsd == 5) { # Case of 5 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1A <- fread(infile1A, header=TRUE)
    read2A <- fread(infile2A, header=TRUE)
    read3A <- fread(infile3A, header=TRUE)
    read4A <- fread(infile4A, header=TRUE)
    read5A <- fread(infile5A, header=TRUE)
    ## Set 2
    read1B <- fread(infile1B, header=TRUE)
    read2B <- fread(infile2B, header=TRUE)
    read3B <- fread(infile3B, header=TRUE)
    read4B <- fread(infile4B, header=TRUE)
    read5B <- fread(infile5B, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1A and second column only for 2A and 3A
    ## Note: this makes it a matrix
    combineA = abind(read1A, read2A[,2], read3A[,2], read4A[,2], read5A[,2], along=2)
    combineB = abind(read1B, read2B[,2], read3B[,2], read4B[,2], read5B[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1A RMSD, infile2A RMSD, infile3A RMSD
    colnames(combineA) <- c("Frame", "RMSD1", "RMSD2", "RMSD3", "RMSD4", "RMSD5")
    colnames(combineB) <- c("Frame", "RMSD1", "RMSD2", "RMSD3", "RMSD4", "RMSD5")
    ## Redefine as a data frame
    combineA <- as.data.frame(combineA)
    combineB <- as.data.frame(combineB)
    ## Append a column that's the average of RMSD cols
    combineA$Average <- rowMeans(combineA[,2:6])
    combineB$Average <- rowMeans(combineB[,2:6])
    ## Append a column that's the STDEV of RMSD cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineA$STDEV <- apply(combineA[,2:6], 1, sd)
    combineB$STDEV <- apply(combineB[,2:6], 1, sd)
} else { # > 5 files
  stop(sprintf("Sorry, this is for 3-5 files. You gave me %s.
       You'll have to modify the if statements for the file reads yourself.
       Stopping now...", sets_rmsd))
}

## Create a new variable with just Frame, Average, and STDEV
save_cols_A <- combineA[,c("Frame", "Average", "STDEV")]

save_cols_B <- combineB[,c("Frame", "Average", "STDEV")]

## Limit to 4 sig figs after decimal
save_cols_clean_A <- format(save_cols_A, digits=4)

save_cols_clean_B <- format(save_cols_B, digits=4)

#----------------------------------------------------------------------#
#--------------------------RMSD OUTFILES-------------------------------#
#----------------------------------------------------------------------#

## Now write a tab-delimited outfile!
## Don't care about the index rownames because that's the frame
write.table(save_cols_clean_A, file = RMSDA, sep="\t", row.names=FALSE, quote=FALSE)

write.table(save_cols_clean_B, file = RMSDB, sep="\t", row.names=FALSE, quote=FALSE)

#----------------------#
#--CONTINUE WITH RMSF--#
#----------------------#

## Deal with sets if X present
if (sets_rmsf <= 2) { # Case if 1 or 2 files
    stop(sprintf("You gave me %s file(s), which is less than 3 input files.
        This is not the script for you because it does standard deviations.
        Stopping now...", sets))
} else if (sets_rmsf == 3) { # Case of 3 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1AF <- fread(infile1AF, header=TRUE)
    read2AF <- fread(infile2AF, header=TRUE)
    read3AF <- fread(infile3AF, header=TRUE)
    ## Set 2
    read1BF <- fread(infile1BF, header=TRUE)
    read2BF <- fread(infile2BF, header=TRUE)
    read3BF <- fread(infile3BF, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1AF and second column only for 2AF and 3AF
    ## Note: this makes it a matrix
    combineAF = abind(read1AF, read2AF[,2], read3AF[,2], along=2)
    combineBF = abind(read1BF, read2BF[,2], read3BF[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1AF RMSF, infile2AF RMSF, infile3AF RMSF
    colnames(combineAF) <- c("Residue", "RMSF1", "RMSF2", "RMSF3")
    colnames(combineBF) <- c("Residue", "RMSF1", "RMSF2", "RMSF3")
    ## Redefine as a data frame
    combineAF <- as.data.frame(combineAF)
    combineBF <- as.data.frame(combineBF)
    ## Append a column that's the average of RMSF cols
    combineAF$Average <- rowMeans(combineAF[,2:4])
    combineBF$Average <- rowMeans(combineBF[,2:4])
    ## Append a column that's the STDEV of RMSF cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineAF$STDEV <- apply(combineAF[,2:4], 1, sd)
    combineBF$STDEV <- apply(combineBF[,2:4], 1, sd)
} else if (sets_rmsf == 4) { # Case of 4 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1AF <- fread(infile1AF, header=TRUE)
    read2AF <- fread(infile2AF, header=TRUE)
    read3AF <- fread(infile3AF, header=TRUE)
    read4AF <- fread(infile4AF, header=TRUE)
    ## Set 2
    read1BF <- fread(infile1BF, header=TRUE)
    read2BF <- fread(infile2BF, header=TRUE)
    read3BF <- fread(infile3BF, header=TRUE)
    read4BF <- fread(infile4BF, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1AF and second column only for 2AF and 3AF
    ## Note: this makes it a matrix
    combineAF = abind(read1AF, read2AF[,2], read3AF[,2], read4AF[,2], along=2)
    combineBF = abind(read1BF, read2BF[,2], read3BF[,2], read4BF[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1AF RMSF, infile2AF RMSF, infile3AF RMSF
    colnames(combineAF) <- c("Residue", "RMSF1", "RMSF2", "RMSF3", "RMSF4")
    colnames(combineBF) <- c("Residue", "RMSF1", "RMSF2", "RMSF3", "RMSF4")
    ## Redefine as a data frame
    combineAF <- as.data.frame(combineAF)
    combineBF <- as.data.frame(combineBF)
    ## Append a column that's the average of RMSF cols
    combineAF$Average <- rowMeans(combineAF[,2:5])
    combineBF$Average <- rowMeans(combineBF[,2:5])
    ## Append a column that's the STDEV of RMSF cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineAF$STDEV <- apply(combineAF[,2:5], 1, sd)
    combineBF$STDEV <- apply(combineBF[,2:5], 1, sd)
} else if (sets_rmsf == 5) { # Case of 5 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1AF <- fread(infile1AF, header=TRUE)
    read2AF <- fread(infile2AF, header=TRUE)
    read3AF <- fread(infile3AF, header=TRUE)
    read4AF <- fread(infile4AF, header=TRUE)
    read5AF <- fread(infile5AF, header=TRUE)
    ## Set 2
    read1BF <- fread(infile1BF, header=TRUE)
    read2BF <- fread(infile2BF, header=TRUE)
    read3BF <- fread(infile3BF, header=TRUE)
    read4BF <- fread(infile4BF, header=TRUE)
    read5BF <- fread(infile5BF, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1AF and second column only for 2AF and 3AF
    ## Note: this makes it a matrix
    combineAF = abind(read1AF, read2AF[,2], read3AF[,2], read4AF[,2], read5AF[,2], along=2)
    combineBF = abind(read1BF, read2BF[,2], read3BF[,2], read4BF[,2], read5BF[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1AF RMSF, infile2AF RMSF, infile3AF RMSF
    colnames(combineAF) <- c("Residue", "RMSF1", "RMSF2", "RMSF3", "RMSF4", "RMSF5")
    colnames(combineBF) <- c("Residue", "RMSF1", "RMSF2", "RMSF3", "RMSF4", "RMSF5")
    ## Redefine as a data frame
    combineAF <- as.data.frame(combineAF)
    combineBF <- as.data.frame(combineBF)
    ## Append a column that's the average of RMSF cols
    combineAF$Average <- rowMeans(combineAF[,2:6])
    combineBF$Average <- rowMeans(combineBF[,2:6])
    ## Append a column that's the STDEV of RMSF cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineAF$STDEV <- apply(combineAF[,2:6], 1, sd)
    combineBF$STDEV <- apply(combineBF[,2:6], 1, sd)
} else { # > 5 files
  stop(sprintf("Sorry, this is for 3-5 files. You gave me %s.
       You'll have to modify the if statements for the file reads yourself.
       Stopping now...", sets_rmsf))
}

## Create a new variable with just Frame, Average, and STDEV
save_cols_AF <- combineAF[,c("Residue", "Average", "STDEV")]

save_cols_BF <- combineBF[,c("Residue", "Average", "STDEV")]

## Limit to 4 sig figs after decimal
save_cols_clean_AF <- format(save_cols_AF, digits=4)

save_cols_clean_BF <- format(save_cols_BF, digits=4)

#----------------------------------------------------------------------#
#--------------------------RMSF OUTFILES-------------------------------#
#----------------------------------------------------------------------#

## Now write a tab-delimited outfile!
## Don't care about the index rownames because that's the residue number
write.table(save_cols_clean_AF, file = RMSFA, sep="\t", row.names=FALSE, quote=FALSE)

write.table(save_cols_clean_BF, file = RMSFB, sep="\t", row.names=FALSE, quote=FALSE)


#------------------------#
#--CONTINUE WITH H-BOND--#
#------------------------#

## Deal with sets if X present
if (sets_hbond <= 2) { # Case if 1 or 2 files
    stop(sprintf("You gave me %s file(s), which is less than 3 input files.
        This is not the script for you because it does standard deviations.
        Stopping now...", sets))
} else if (sets_hbond == 3) { # Case of 3 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1AH <- fread(infile1AH, header=TRUE)
    read2AH <- fread(infile2AH, header=TRUE)
    read3AH <- fread(infile3AH, header=TRUE)
    ## Set 2
    read1BH <- fread(infile1BH, header=TRUE)
    read2BH <- fread(infile2BH, header=TRUE)
    read3BH <- fread(infile3BH, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1AH and second column only for 2AH and 3AH
    ## Note: this makes it a matrix
    combineAH = abind(read1AH, read2AH[,2], read3AH[,2], along=2)
    combineBH = abind(read1BH, read2BH[,2], read3BH[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1AH HB, infile2AH HB, infile3AH HB
    colnames(combineAH) <- c("Frame", "HB1", "HB2", "HB3")
    colnames(combineBH) <- c("Frame", "HB1", "HB2", "HB3")
    ## Redefine as a data frame
    combineAH <- as.data.frame(combineAH)
    combineBH <- as.data.frame(combineBH)
    ## Append a column that's the average of HB cols
    combineAH$Average <- rowMeans(combineAH[,2:4])
    combineBH$Average <- rowMeans(combineBH[,2:4])
    ## Append a column that's the STDEV of HB cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineAH$STDEV <- apply(combineAH[,2:4], 1, sd)
    combineBH$STDEV <- apply(combineBH[,2:4], 1, sd)
} else if (sets_hbond == 4) { # Case of 4 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1AH <- fread(infile1AH, header=TRUE)
    read2AH <- fread(infile2AH, header=TRUE)
    read3AH <- fread(infile3AH, header=TRUE)
    read4AH <- fread(infile4AH, header=TRUE)
    ## Set 2
    read1BH <- fread(infile1BH, header=TRUE)
    read2BH <- fread(infile2BH, header=TRUE)
    read3BH <- fread(infile3BH, header=TRUE)
    read4BH <- fread(infile4BH, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1AH and second column only for 2AH and 3AH
    ## Note: this makes it a matrix
    combineAH = abind(read1AH, read2AH[,2], read3AH[,2], read4AH[,2], along=2)
    combineBH = abind(read1BH, read2BH[,2], read3BH[,2], read4BH[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1AH HB, infile2AH HB, infile3AH HB
    colnames(combineAH) <- c("Frame", "HB1", "HB2", "HB3", "HB4")
    colnames(combineBH) <- c("Frame", "HB1", "HB2", "HB3", "HB4")
    ## Redefine as a data frame
    combineAH <- as.data.frame(combineAH)
    combineBH <- as.data.frame(combineBH)
    ## Append a column that's the average of HB cols
    combineAH$Average <- rowMeans(combineAH[,2:5])
    combineBH$Average <- rowMeans(combineBH[,2:5])
    ## Append a column that's the STDEV of HB cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineAH$STDEV <- apply(combineAH[,2:5], 1, sd)
    combineBH$STDEV <- apply(combineBH[,2:5], 1, sd)
} else if (sets_hbond == 5) { # Case of 5 files
    ## Reading each file as a data.table.
    ## Bonus - fread is much faster than read.csv
    read1AH <- fread(infile1AH, header=TRUE)
    read2AH <- fread(infile2AH, header=TRUE)
    read3AH <- fread(infile3AH, header=TRUE)
    read4AH <- fread(infile4AH, header=TRUE)
    read5AH <- fread(infile5AH, header=TRUE)
    ## Set 2
    read1BH <- fread(infile1BH, header=TRUE)
    read2BH <- fread(infile2BH, header=TRUE)
    read3BH <- fread(infile3BH, header=TRUE)
    read4BH <- fread(infile4BH, header=TRUE)
    read5BH <- fread(infile5BH, header=TRUE)
    ## Combine into one dataset
    ## Use both columns 1AH and second column only for 2AH and 3AH
    ## Note: this makes it a matrix
    combineAH = abind(read1AH, read2AH[,2], read3AH[,2], read4AH[,2], read5AH[,2], along=2)
    combineBH = abind(read1BH, read2BH[,2], read3BH[,2], read4BH[,2], read5BH[,2], along=2)
    ## Change the column names so future life makes sense
    ## Your data are now Frame, infile1AH HB, infile2AH HB, infile3AH HB
    colnames(combineAH) <- c("Frame", "HB1", "HB2", "HB3", "HB4", "HB5")
    colnames(combineBH) <- c("Frame", "HB1", "HB2", "HB3", "HB4", "HB5")
    ## Redefine as a data frame
    combineAH <- as.data.frame(combineAH)
    combineBH <- as.data.frame(combineBH)
    ## Append a column that's the average of HB cols
    combineAH$Average <- rowMeans(combineAH[,2:6])
    combineBH$Average <- rowMeans(combineBH[,2:6])
    ## Append a column that's the STDEV of HB cols
    ## The 1 means you're applying one function
    ## Which is the standard deviation function, sd
    combineAH$STDEV <- apply(combineAH[,2:6], 1, sd)
    combineBH$STDEV <- apply(combineBH[,2:6], 1, sd)
} else { # > 5 files
  stop(sprintf("Sorry, this is for 3-5 files. You gave me %s.
       You'll have to modify the if statements for the file reads yourself.
       Stopping now...", ssets_hbond))
}

## Create a new variable with just Frame, Average, and STDEV
save_cols_AH <- combineAH[,c("Frame", "Average", "STDEV")]

save_cols_BH <- combineBH[,c("Frame", "Average", "STDEV")]

## Limit to 4 sig figs after decimal
save_cols_clean_AH <- format(save_cols_AH, digits=4)

save_cols_clean_BH <- format(save_cols_BH, digits=4)

#------------------------------------------------------------------------#
#--------------------------H-BOND OUTFILES-------------------------------#
#------------------------------------------------------------------------#

## Now write a tab-delimited outfile!
## Don't care about the index rownames because that's the frame

write.table(save_cols_clean_AH, file = HBONDA, sep="\t", row.names=FALSE, quote=FALSE)

write.table(save_cols_clean_BH, file = HBONDB, sep="\t", row.names=FALSE, quote=FALSE)
