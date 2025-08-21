library(ParallelStructure)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run ParallelStructure                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Usage: nohup Rscript 6.1_run_structure.R &

options(scipen=999)

# Genotypes

gen <- fread("data/out/structure_sub/rhina_handmade_rad_drop_sort.stru", header = F)

#~~ Specify in and out files for structure

infile <- "data/out/structure_sub/rhina_handmade_rad_drop_sort.stru"
outpath <- "data/out/structure_sub/"

#~~ construct job matrix and write to job file

nrep <- 10
up_to_k <- 8
niter <- 100000
burnin <- 50000

#nrep <- 1
#up_to_k <- 6
#niter <- 100
#burnin <- 100

#~~ define variables for job matrix

k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))


#~~ make the job matrix
pop <- "1,2,3,4,5"

jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                 rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(jobs), ncol = length(jobs[1,]), file = "jobs.txt")

#~~ file path to structure

STR_path='/Users/emilyhumble/software/structure/'


#~~ Run Parallel Structure

#setwd("/exports/cmvm/eddie/eb/groups/ogden_grp/structure")

ParallelStructure::parallel_structure(structure_path=STR_path,
                                      joblist="jobs.txt",
                                      n_cpu=8,
                                      infile=infile,
                                      outpath=outpath,
                                      numinds = nrow(gen),
                                      numloci=(ncol(gen)-2)/2,
                                      locdata = 0,
                                      popdata = 1,
                                      popflag = 0,
                                      usepopinfo = 0,
                                      locprior=0,
                                      printqhat=1,
                                      plot_output=0,
                                      onerowperind=1)
