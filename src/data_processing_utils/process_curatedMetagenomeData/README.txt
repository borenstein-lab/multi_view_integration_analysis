# Scripts in this folder are responsible of processing multiple case-control shotgun datasets from the CuratedMetagenomicData package.

--------------------------------------------------------------------
Commands - running without a container:
--------------------------------------------------------------------

# Create dedicated screen (once)
screen -S wgs

# Or: re-enter screen
screen -r wgs

# Install required libraries (once)
/usr/local/lib/R-4.1.2/bin/R
install.packages(c("config","Matrix.utils","fs","BiocManager"),lib="/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/R4_env",repos=c("https://ftp.cc.uoc.gr/mirrors/CRAN/")) 
library(BiocManager)
BiocManager::install("S4Vectors",lib="/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/R4_env")
BiocManager::install("curatedMetagenomicData",lib="/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/R4_env")
quit()

# Update PATH - NOT WORKING (so .libPaths(...) added to main.R)
# setenv PATH $PATH\:/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/R4_env\:/usr/local/stow/R-4.1.2/lib/R-4.1.2/lib/R/library

# Run script
cd /specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/microbiome-project2/src/process_curatedMetagenomeData
/usr/local/lib/R-4.1.2/bin/Rscript main.R >& main.log.03122021

# Exit (detach) screen
# Ctrl+a d

--------------------------------------------------------------------
Commands - running in container with R 4.1 (not working, need to switch to new image):
--------------------------------------------------------------------

# Create dedicated screen (once)
screen -S wgs

# Or: re-enter screen
screen -r wgs

# View running containers: udocker ps
# Create container (once)
udocker create --name=cur_metagenomic_data rocker/tidyverse

# Run container
udocker run --volume=/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/microbiome-project2:/mnt/proj2 cur_metagenomic_data bash

# Install required libraries (once)
<as above>

# Go to mounted directory
cd /mnt/proj2/src/process_curatedMetagenomeData

# Run script
Rscript main.R

# Exit container
exit

# Exit (detach) screen
# Ctrl+a d

# Erase screen
screen -XS <session-id> quit














