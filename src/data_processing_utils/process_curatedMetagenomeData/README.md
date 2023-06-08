Scripts in this folder are responsible of processing multiple case-control shotgun datasets from the CuratedMetagenomicData package.

## Commands

1. Install required libraries (once)
```
R
install.packages(c("config","Matrix.utils","fs","BiocManager"),lib=<...>) 
library(BiocManager)
BiocManager::install("S4Vectors",lib="/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/R4_env")
BiocManager::install("curatedMetagenomicData",lib="/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/R4_env")
quit()
```

2. Update PATH or add `.libPaths(...)` to main.R

3. Run script
`Rscript main.R`
