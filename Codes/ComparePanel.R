# In this script we will check how many of
#    the predicted SNPs are in common with other 
#    published SNP panels

source('Codes/Functions.R')

DIR = 'PaperSNPs/'
SNPfiles <- list.files(DIR, include.dirs = T, full.names = T)
PanelFiles <- lapply(SNPfiles, read.csv)
Names <- c('genotype_RNA_DNA', 'Info_snps', 'univ_iden')
names(PanelFiles) <- Names

PanelFiles_rsID <- list(PanelFiles[['genotype_RNA_DNA']]$refsnp_id, 
                        PanelFiles[['Info_snps']]$Locus_Name, 
                        PanelFiles[['univ_iden']]$rsID)

PanelFiles_rsID <- lapply(PanelFiles_rsID, as.character)
names(PanelFiles_rsID) <- Names



### the predicted target panel
TargetPositions <- readRDS('TargetPositions.rds')
TargetPositionsDF <- data.frame(stringr::str_split_fixed(TargetPositions, '_',2))
colnames(TargetPositionsDF) <- c('chr', 'start')
TargetPositions_ChrPos <- paste0(TargetPositionsDF$chr, sep='_', TargetPositionsDF$start)


### choosing the Biomart dataset
mart.hg19 <-  useMart(biomart="ENSEMBL_MART_SNP",
                      host="grch37.ensembl.org", 
                      path="/biomart/martservice", 
                      dataset="hsapiens_snp")


mart.hg38 <- useMart("ENSEMBL_MART_SNP", 
                     dataset="hsapiens_snp")


snp_attributes <- c("refsnp_id", "chr_name", "chrom_start")
snp_Coordinates <- lapply(PanelFiles_rsID, function(snp_ids){
                           getBM(attributes=snp_attributes, 
                                 filters ="snp_filter", 
                                 values=snp_ids, 
                                 mart=mart.hg38 )} )


saveRDS(snp_Coordinates, 'Data/snp_Coordinates.rds')  
snp_Coordinates <- readRDS('Data/snp_Coordinates.rds')
snp_ChrPos <- lapply(snp_Coordinates, GenerateChrPosID)
  

### WHY ????
lapply(snp_ChrPos, function(x) sum(x %in% TargetPositions_ChrPos) )


# head(ir) >> check if they are present in the iranome dataset
