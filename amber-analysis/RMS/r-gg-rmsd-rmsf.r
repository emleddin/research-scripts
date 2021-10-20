## Run this with "Rscript r-gg-rmsd-rmsf.r"
## (Assuming you've already installed R...)
## Render RMSD and RMSF plots in R

## Load packages and install missing
## You can either load tidyverse, or load dplyr, tidyr, and ggplot2 separately
## If you want a specific font, you might want extrafont
packages <- c("data.table", "ggplot2", "dplyr", "tidyr", "ggpubr")
install.packages(setdiff(packages, rownames(installed.packages())))
invisible(lapply(packages, library, character.only = TRUE))

## Path to the total_bb_rm.dat file
## You're not *required* to use the absolute path (and too long a path will
## cause issues...)
infile_rmsd_wt <- Sys.glob("/absolute/path/to/wt/WT_total_bb_rms.dat")
infile_rmsf_wt <- Sys.glob("/absolute/path/to/wt/WT_total_bb_rms.dat")

infile_rmsd_snp1 <- Sys.glob("/absolute/path/to/snp1/SNP1_total_bb_rms.dat")
infile_rmsf_snp1 <- Sys.glob("/absolute/path/to/snp1/SNP1_total_bb_rms.dat")

infile_rmsd_snp2 <- Sys.glob("/absolute/path/to/snp2/SNP1_total_bb_rms.dat")
infile_rmsf_snp2 <- Sys.glob("/absolute/path/to/snp2/SNP1_total_bb_rms.dat")

infile_rmsd_snp3 <- Sys.glob("/absolute/path/to/snp3/SNP1_total_bb_rms.dat")
infile_rmsf_snp3 <- Sys.glob("/absolute/path/to/snp3/SNP1_total_bb_rms.dat")

## Variables
font_size = 12
## Background transparency in individual RMSD/RMSF plots you may make?
transparency=TRUE
## Frames/div_frame_by becomes nanoseconds (typically 500 or 1000)
div_frame_by = 500

## Colors for the data!
# https://www.schemecolor.com/high-contrast.php
# wt_color = "#A2117B" #pink
# snp1_color = "#0000FF" #blue
# snp2_color = "#8931EF" # purple
# snp3_color = "#3C8033" # green

## Define the systems
System = c("WT", "SNP1", "SNP2", "SNP3")
## Define the colors to use (same order as systems)
Color  = c("#A2117B", "#0000FF", "#8931EF", "#3C8033")
## Define the shapes to use (same order as systems)
Shape  = as.integer(c(15, 16, 17, 18))
## Create a DataFrame of the systems, colors, and shapes to use
## This is necessary so the legend on the combined image is correct!
prop_df    = data.frame(Color, Shape, System)
## Save the systems as rownames for accessing based on System
rownames(prop_df) <- System

## Plot between these for axes (Need 'c'!!!)
## You can comment them out if you want
## RMSD plotting limits
plot_x_ax_rmsd = c(0, 100*div_frame_by)
plot_y_ax_rmsd = c(0, 5)

## RMSF plotting limits
plot_x_ax_rmsf = c(1, 301)
plot_y_ax_rmsf = c(0, 7)

## Set up RMSF residue axis label Positions, Labels)
rmsf_pos = c(1, 51, 101, 151, 201, 251, 295)
rmsf_labs = c("16", "66", "GS", "166", "216", "266", "DNA")

##--------------------------##
##---- Set up Functions ----##
##--------------------------##
## These set up the common command types so that figures with similar data all
##  appear the same.
## From here, skip to "Creating the Plots" to read in the data/make plots

##--------------------------- Data Reading Function---------------------------##
read_file <- function(infile, head_row=FALSE, c1="Col1", c2="Col2", sys_lab="SYS") {
  ## Reading each file as a data.table.
  ## Bonus - fread is much faster than read.csv
  data <- fread(infile, header=head_row)
  ## Rename columns -- index starts at 1
  names(data)[1] <- c1
  names(data)[2] <- c2
  data$System <- sys_lab
  return(data)
}

##---------------------------- Set Up Plot Theme ----------------------------##
## Set up background transparency
if (transparency) {
  back_col = "transparent"
} else {
  back_col = "white"
}

## Create my own theme, based on the classic theme (kinda like theme_pubr)
## May want to modify theme_article https://rdrr.io/cran/egg/src/R/theme.r
theme_eml <- function (backgd_col=back_col, accent_col="black") {
  theme_classic(base_size=font_size, base_family="Arial") %+replace%
    theme(
      plot.title = element_text(face = "bold", color=accent_col, vjust=2),
      ## Top Right Bottom Left margins
      plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "in"),
      plot.background = element_rect(fill = backgd_col, color=backgd_col),
      panel.background = element_rect(fill = backgd_col,
                                      ## color is for the outline around plot
                                      color = backgd_col),
      panel.border = element_rect(fill = NA, color = accent_col, size = 1),
      panel.spacing =  unit(c(0.05, 0.05, 0.05, 0.05), "in"),

      ## Color x-ticks
      axis.ticks.x = element_line(color=accent_col),
      axis.ticks.y = element_line(color=accent_col),

      ## Remove axis line so x/y matches border
      axis.line = element_blank(),

      axis.title.x = element_text(color = accent_col, face="bold",
                                  margin = margin(5, 0, 0, 0)), #vjust=-0.5),
      axis.title.y = element_text(angle = 90, #hjust=0.5, vjust=-2,
                                  color = accent_col, face="bold",
                                  margin = margin(0, 10, 0, 0)),
      axis.text = element_text(color = accent_col, #face="bold",
                               size=(font_size-3)),

      # legend.position="bottom",
      legend.position="right",
      legend.text=element_text(color = accent_col, size=(font_size-3)),
      legend.margin=margin(0,0,0,0),
      # legend.border = element_blank(),
      # legend.box.margin=margin(-10,-10,-10,-10),
      # legend.key.size = unit(0.5, 'in'),
      legend.background = element_rect(fill = backgd_col, color=NA)
    )
}

##---------------------------- Set Up RMSD Plot ----------------------------##
## Uses Generalized Additive Model for smoothing!
plot_rmsd <- function(data, x_ax=NULL, y_ax=NULL,
  div_frame_by=1, pt=1, colors="black") {
  plot <- ggplot(data=data, aes(x=Frame, y=RMSD)) +
    ## Plot the data as points
    geom_point(size = 0.5, shape = pt, aes(color=System)) +
    ## Plot a line smoothing function over the points
    geom_smooth(color = "white", show.legend = FALSE) +
    ## Change the axis labelsgg_col
    labs(x = "Time (ns)", y = "RMSD (\u212b)") +
    ## Plot between given x-axes, but label using x / div_frame_by
    {if(exists("plot_x_ax_rmsd"))scale_x_continuous(limits = x_ax,
      labels=function(x)x/div_frame_by)} +
    ## Plot between given y-axes
    {if(exists("plot_y_ax_rmsd"))scale_y_continuous(limits = y_ax)} +
    ## Change color scheme
    scale_color_manual(values = colors, limits = force) +
    ## Increase size of points in legend
    guides(color = guide_legend(override.aes = list(size=10)),
           shape = pt) +
    ## Apply my theme
    theme_eml()
  return(plot)
}

##---------------------------- Set Up RMSF Plot ----------------------------##
plot_rmsf <- function(data, x_ax=NULL, y_ax=NULL, pt=1, colors="black") {
  plot <- ggplot(data=data, aes(x=Residue, y=RMSF)) +
    ## Plot the data as points
    geom_point(size = 1, shape = pt, aes(color=System)) +
    ## Change the axis labels
    labs(x = "Residue", y = "RMSF (\u212b)") +
    ## Plot between given x-axes
    {if(exists("plot_x_ax_rmsf"))scale_x_continuous(limits = plot_x_ax_rmsf,
      breaks = rmsf_pos, label = rmsf_labs)} +
    ## Distinct axis labels for biology residue labelling
    # {if(exists("rmsf_labs"))scale_x_discrete(labels=rmsf_labs)} +
    ## Plot between given y-axes
    {if(exists("plot_y_ax_rmsf"))scale_y_continuous(limits = plot_y_ax_rmsf)} +
    ## Change color scheme
    scale_color_manual(values = colors, limits = force) +
    ## Increase size of points in legend
    guides(color = guide_legend(override.aes = list(size=10))) +
    ## Apply my theme
    theme_eml()
  return(plot)
}

##---------------------------- Set Up for Legend ----------------------------##

plot_leg <- function(data, colors, shapes) {
  plot <- ggplot(data=data, aes(x=Shape, y=Color)) +
    ## Plot the data as points
    geom_point(size = 0.5, aes(color=System, shape=System)) +
    ## Change color scheme -- use unique to preserve order passed in DF
    ## Otherwise it'll use alphabetical and be wrong
    scale_color_manual(values = colors, breaks = unique(System)) +
    ## Assign points based on passed shapes, and use unique to match DF
    scale_shape_manual(values = shapes, breaks = unique(System)) +
    ## Increase size of points in legend
    guides(color = guide_legend(override.aes = list(size=10))) +
    ## Apply my theme & force legend to bottom
    theme_eml() + theme(legend.position="bottom")
  return(plot)
}

##----------------------------##
##---- Creating the Plots ----##
##----------------------------##
##----------------------------------------------------------------------------##
## ---------------- WT
## `sys_lab` label MUST match what was defined in `System`
read_rmsd_wt <- read_file(infile_rmsd_wt, head_row=TRUE,
  c1="Frame", c2="RMSD", sys_lab="WT")
read_rmsf_wt <- read_file(infile_rmsf_wt, head_row=TRUE,
  c1="Residue", c2="RMSF", sys_lab="WT")

## Plot RMSD
p_rmsd_wt <- plot_rmsd(read_rmsd_wt, x_ax=plot_x_ax_rmsd, y_ax=plot_y_ax_rmsd,
  div_frame_by, pt=prop_df["WT", "Shape"], colors=prop_df["WT", "Color"])

## Plot RMSF
p_rmsf_wt <- plot_rmsf(read_rmsf_wt, pt=prop_df["WT", "Shape"],
  colors=prop_df["WT", "Color"])

## Save an individual RMSD plot for WT
ggsave("rmsd.png", plot=p_rmsd_wt, width=5, height=3, units="in",
    bg = back_col, dpi=300)

## Save an individual RMSF plot for WT
ggsave("rmsf.png", plot=p_rmsf_wt, width=5, height=3, units="in",
    bg = back_col, dpi=300)

## ---------------- SNP1
read_rmsd_snp1 <- read_file(infile_rmsd_snp1, head_row=TRUE,
  c1="Frame", c2="RMSD", sys_lab="SNP1")
read_rmsf_snp1 <- read_file(infile_rmsf_snp1, head_row=TRUE,
  c1="Residue", c2="RMSF", sys_lab="SNP1")

## Plot RMSD
p_rmsd_snp1 <- plot_rmsd(read_rmsd_snp1,
  x_ax=plot_x_ax_rmsd, y_ax=plot_y_ax_rmsd, div_frame_by,
  pt=prop_df["SNP1", "Shape"], colors=prop_df["SNP1", "Color"])

## Plot RMSF
p_rmsf_snp1 <- plot_rmsf(read_rmsf_snp1,
  pt=prop_df["SNP1", "Shape"], colors=prop_df["SNP1", "Color"])

## ---------------- SNP2
read_rmsd_snp2 <- read_file(infile_rmsd_snp2, head_row=TRUE,
  c1="Frame", c2="RMSD", sys_lab="SNP2")
read_rmsf_snp2 <- read_file(infile_rmsf_snp2, head_row=TRUE,
  c1="Residue", c2="RMSF", sys_lab="SNP2")

## Plot RMSD
p_rmsd_snp2 <- plot_rmsd(read_rmsd_snp2,
  x_ax=plot_x_ax_rmsd, y_ax=plot_y_ax_rmsd, div_frame_by,
  pt=prop_df["SNP2", "Shape"], colors=prop_df["SNP2", "Color"])

## Plot RMSF
p_rmsf_snp2 <- plot_rmsf(read_rmsf_snp2,
  pt=prop_df["SNP1", "Shape"], colors=prop_df["SNP2", "Color"])

## ---------------- SNP3
read_rmsd_snp3 <- read_file(infile_rmsd_snp3, head_row=TRUE,
  c1="Frame", c2="RMSD", sys_lab="SNP3")
read_rmsf_snp3 <- read_file(infile_rmsf_snp3, head_row=TRUE,
  c1="Residue", c2="RMSF", sys_lab="SNP3")

## Plot RMSD
p_rmsd_snp3 <- plot_rmsd(read_rmsd_snp3,
  x_ax=plot_x_ax_rmsd, y_ax=plot_y_ax_rmsd, div_frame_by,
  pt=prop_df["SNP3", "Shape"], colors=prop_df["SNP3", "Color"])

## Plot RMSF
p_rmsf_snp3 <- plot_rmsf(read_rmsf_snp3,
  pt=prop_df["SNP3", "Shape"], colors=prop_df["SNP3", "Color"])

## ---------------- Legend
## Add a fake legend that grabs the shape/color combos
# my_leg <- plot_leg(prop_df, colors=unique(prop_df[, "Color"]),
#           shapes=unique(prop_df[, "Shape"]))
my_leg <- plot_leg(prop_df, colors=prop_df[, "Color"],
          shapes=prop_df[, "Shape"])

## Save the test legend
# ggsave("leg.png", plot=my_leg, width=5, height=3, units="in",
#     bg = "white", dpi=300)

##------------------ All together

## Plot all of the individual plots together in 1 massive figure with subplot
## labels (`labels`).
figure <- ggarrange(p_rmsd_wt, p_rmsf_wt, p_rmsd_snp1, p_rmsf_snp1,
      p_rmsd_snp2, p_rmsf_snp2, p_rmsd_snp3, p_rmsf_snp3,
      labels = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(G)", "(H)"),
      ncol = 2, nrow = 4, legend="bottom",
      ## Use that legend you made as a common legend
      legend.grob=get_legend(my_leg))

## Save the massive figure!
ggsave("plots.png", plot=figure, width=8.5, height=11, units="in",
    bg = "white", dpi=300)
