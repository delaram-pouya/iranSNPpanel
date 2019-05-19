library(data.table)
library(parallel)

setwd("~/iranome/chromosomes/")
l = list.files(".", ".rds")
l = l[-c(23,24,25)]
samples <- sub(".rds", "", l)
f = mclapply(l, readRDS ,mc.cores = detectCores())
## names(f) <- samples
f <- do.call(rbind, f)
h <- do.call(cbind, mclapply(list(getChrPos, getDepth, getQual,
                                  getACS,getANC,getHets,getHoms,getHetfreq,getHomfreq), function(i) i(f), mc.cores = 1))

g <- merge( f, h, by.x = c('chrom', 'pos'), by.y = c('chrom', 'pos'), all.x = T,sort=F)
unselectCol = c("dbnsfp_pred", "dbnsfp_pred_order","exon_num" , "filter" ,"orig_alt_alleles","quality_metrics","transcripts" ,
                "vep_annotations" ,"xpos", "xstart" , "xstop", "genotype_depths", "genotype_qualities", "pop_acs",
                "pop_ans", "pop_hetfreq", "pop_hets", "pop_homfreq", "pop_homs")

g1 <- g[,!colnames(g) %in% unselectCol , drop = F]
g2 <- lapply(g1[,colnames(g1)[59:106]],  function(i) sapply(i, "[[", 1))
g1[,59:106] <- g2
saveRDS(h, "totalIranome.rds")
saveRDS(g1, "cleanIranome.rds")


########################################### functions 

getChrPos <- function(x) data.frame(chrom=x$chrom, pos=as.numeric(x$pos))

getDepth <- function(x){
  # [[i]][[1]]> all_individuals(not just variant carriers), j[2]> (1)value (2)freq
  depth = t(sapply(1:nrow(x), function(i) sapply(x$genotype_depths[[i]][[1]], '[[', 2)))
  colnames(depth) = paste("D", seq(2.5, 97.5, by=5), collapse = NULL, sep="")
  return(depth)}

getQual <- function(x){
  qual = t(sapply(1:nrow(x), function(i) sapply(x$genotype_qualities[[i]][[1]], function(j) j[[2]])))
  colnames(qual) = paste("Q", seq(2.5, 97.5, by=5), collapse = NULL, sep="")
  return(as.data.frame(qual))}

getACS <- function(x){
  acs = t(sapply(1:nrow(x), function(i) x$pop_acs[[i]]))
  colnames(acs) = paste("acs_", colnames(acs), collapse = NULL, sep="")
  return(as.data.frame(acs))}

getANC <- function(x) {
  anc = t(sapply(1:nrow(x), function(i) x$pop_ans[[i]]))
  colnames(anc) = paste("anc_", colnames(anc), collapse = NULL, sep="")
  return(as.data.frame(anc))}

getHets <- function(x) {
  hets = t(sapply(1:nrow(x), function(i) x$pop_hets[[i]]))
  colnames(hets) = paste("hets_", colnames(hets), collapse = NULL, sep="")
  return(as.data.frame(hets))}

getHoms <- function(x) {
  homs = t(sapply(1:nrow(x), function(i) x$pop_homs[[i]]))
  colnames(homs) = paste("homs_", colnames(homs), collapse = NULL, sep="")
  return(as.data.frame(homs)) }


getHomfreq <- function(x) {
  homFreq = t(sapply(1:nrow(x), function(i) x$pop_homfreq[[i]]))
  colnames(homFreq) = paste("homFreq_", colnames(homFreq), collapse = NULL, sep="")
  return(as.data.frame(homFreq)) }

getHetfreq <- function(x) {
  hetFreq = t(sapply(1:nrow(x), function(i) x$pop_hetfreq[[i]]))
  colnames(hetFreq) = paste("hetFreq_", colnames(hetFreq), collapse = NULL, sep="")
  return(as.data.frame(hetFreq)) }


getQualMetr <- function(x){
  qual_metrics =t(sapply(1:nrow(x), function(i) as.numeric(x$quality_metrics[[i]])))
  colnames(qual_metrics) <- names(x$quality_metrics[[1]])
  return(as.data.frame(qual_metrics))}


