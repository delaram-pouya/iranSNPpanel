# In this script we will try to 
#   modify the crawled Iranome dataset
#   in the Chr.rds files to gain a
#   clean and easy-to-use data structure


source('Codes/Clean_Iranome_functions.R')
Initialize()


DIR = ("~/iranome/chromosomes/")
chrList <- list.files( DIR, ".rds", full.names = T, include.dirs = T)

## get autosomal chromosomes 
chrList = chrList[1:22]
samples <- gsub("/home/delaram/iranome/chromosomes//" ,"" , sub(".rds", "", chrList ))
chrFiles = mclapply( chrList, readRDS ,mc.cores = detectCores())

## names(chrFiles) <- samples
chrFiles <- do.call(rbind, chrFiles)
FunctionsToApplyList <- list(getChrPos, getDepth, getQual,
                         getACS,getANC,getHets,getHoms,getHetfreq,getHomfreq)
  
chrFilesModified <- do.call(cbind, mclapply(FunctionsToApplyList, function(i) i(chrFiles), mc.cores = 1))


chrFilesMerged <- merge( chrFiles, chrFilesModified , by.x = c('chrom', 'pos'),
                         by.y = c('chrom', 'pos'), all.x = T,sort=F)


ColsToDrop = c("dbnsfp_pred", "dbnsfp_pred_order","exon_num" , "filter" ,
               "orig_alt_alleles","quality_metrics","transcripts" ,
                "vep_annotations" ,"xpos", "xstart" , "xstop", 
               "genotype_depths", "genotype_qualities", "pop_acs",
                "pop_ans", "pop_hetfreq", "pop_hets", "pop_homfreq", "pop_homs")


chrFilesCleaned <- chrFilesMerged[,!colnames(chrFilesMerged) %in% ColsToDrop , drop = F]
temp <- lapply(chrFilesCleaned[,colnames(chrFilesCleaned)[59:106]],  function(i) sapply(i, "[[", 1))
chrFilesCleaned[,59:106] <- temp


saveRDS(h, "totalIranome.rds")
saveRDS(g1, "cleanIranome.rds")






