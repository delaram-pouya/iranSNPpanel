# This script includes the functions 
#   that are used in the Clean_iranome.R script

###### functions 


Initialize <- function(){
  options(stringsAsFactors = FALSE)
  library(data.table)
  library(parallel)
  library(biomaRt)
  library(UpSetR)
  library(tictoc)
}


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



GenerateChrPosID <- function(x){ paste(x$chr_name,sep='_',x$chrom_start)}

