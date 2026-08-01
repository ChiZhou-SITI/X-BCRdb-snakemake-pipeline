library(alakazam)
library(shazam)
library(dplyr)
library(ggplot2)

args <- commandArgs()

filein=args[6]
fileout=args[7]

dataset  = read.table(filein, sep = '\t', header = 1, stringsAsFactors = FALSE, quote = '', check.names = TRUE)
### 计算SHM
##  regionDefinition = IMGT_V 计算V区域
if (TRUE){
db_obs <- observedMutations(dataset, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            #regionDefinition=NULL,
                            regionDefinition=IMGT_VDJ,
                            frequency=TRUE,
                            combine = TRUE,
                            nproc=10)

#result = db_obs %>%  select(sequence_id, starts_with("mu_freq"))
write.table(db_obs, fileout, sep = '\t', quote = F, row.names = F, col.names = T)
}


















