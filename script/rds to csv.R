## CONVERT ALL FILES IN A DIRECTORY

## Load Foreign Library
library(foreign)

## Set working directory in which rds files can be found)
setwd("../data/TCGA_eRNAexp")

## Convert all files in wd from RDS to CSV
### Note: alter the write/read functions for different file types.  
### rds->csv used in this specific example

for (f in Sys.glob('*.rds')) 
  write.csv(read.rds(f), file = gsub('rds$', 'csv', f))