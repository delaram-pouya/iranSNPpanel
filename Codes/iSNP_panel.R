library(parallel)
setwd("/home/delaram/iranome/chromosomes")
ir= readRDS("cleanIranome.rds")
#
#### Depth > min #Depth>32.5 
ir$Dsum = 0
ir$Dsum = apply(ir[,25:38], 1,sum) # ir[,19:38] > total, should be 800  
hist(ir$Dsum)
d.sorted = sort(ir$Dsum)
d.cut.off = d.sorted[round(length(d.sorted)*0.3)]

#### Quality > min #Qual>37.5
ir$Qsum = 0
ir$Qsum = apply(ir[,46:58], 1,sum) # ir[,39:58] > total, should be 800  
hist(ir$Qsum)
q.sorted = sort(ir$Qsum)
q.cut.off = q.sorted[round(length(q.sorted)*0.3)]

#### Allele_number > min #sequenced individuals
an.sorted = sort(ir$allele_num)
hist(ir$allele_num)
an.cut.off = an.sorted[round(length(an.sorted)*0.1)]

ir$HetI = 0; ir$F.IS = 0; ir$F.ST = 0 ;ir$F.IT = 0
l<- t(apply(ir, 1, calcF))
sapply(1:nrow(ir), function(i) {ir[i,]$F.IS <<- l[i,1]; ir[i,]$F.ST <<- l[i,2]; 
ir[i,]$F.IT <<- l[i,3]; ir[i,]$HetI <<- l[i,4]} )

## filtering > add fst<0.06 , HetI>0.4
ir2 <- ir[  ir$allele_freq>0.4 & ir$allele_freq<0.6 & 
              ir$Dsum> d.cut.off & ir$allele_num>an.cut.off & ir$Qsum>q.cut.off &
              ir$category!= "lof_variant" & 
              ir$major_consequence %in% c("synonymous_variant", "intron_variant", "intergenic_variant" )&
              ir$alt %in% c("A","T","C","G") & ir$ref %in% c("A","T","C","G"),]

#& ir$HetI>0.4 & ir$F.ST<0.06

hist(ir2$site_quality[ir2$site_quality<256174])
summary(ir2$site_quality)
summary(ir$site_quality)
########### distance between them 
spl <- split( ir2 , f = ir2$chrom )
spl <- lapply(spl, function(i) calcDistance(i$pos))
dis <- do.call(c,spl)
sum(dis==0) ### why???
summary(dis)
hist(dis[dis>500000])
length(dis>500000)

ir2$dis=NA
ir2$dis=sapply(1:length(ir2$pos), function(i) ir2$dis<<-calcDistance(i) )
dim(ir2[ir2$dis,])

lt=1
calcDistance <- function(i){
  if(i==1 | i==2) return(FALSE)
  else if(is.na(ir2$pos[i]) |is.na(ir2$pos[i+1]) |is.na(ir2$pos[i-1])) return(FALSE)
  
  else if(ir2$chrom[i]==ir2$chrom[i-1] & ir2$chrom[i]==ir2$chrom[i+1] & 
          (ir2$pos[i]- ir2$pos[i-1])>1e05 & (ir2$pos[i+1]- ir2$pos[i]>1e05)) {
    lt <<- i ; return(TRUE)}
  
  else if(ir2$chrom[i]==ir2$chrom[lt] & (ir2$pos[i]- ir2$pos[lt])>1e05 & lt!=1) {lt <<- i ; return(TRUE)}
  
  else return(FALSE)}



################ 

tribes <- unique(sub("acs_","", grep("acs_", colnames(ir), val=T)))
P.A <- sapply(tribes, function(x) (ir2[,paste0("homFreq_",x)] + ir2[,paste0("hetFreq_",x)])/2)
H.exp <- 2 * P.A * (1-P.A)
H.obs <- sapply(tribes, function(x) 2*ir2[,paste0("hets_",x)]/ir2[,paste0("anc_",x)])
HI <- sapply(1:nrow(H.obs), function(i) sum(sapply(tribes, function(x) ir2[i,paste0("anc_",x)]*H.obs[i,x]))/ir2[i,'allele_num'] ) 
HS <- sapply(1:nrow(H.exp), function(i) sum(sapply(tribes, function(x) ir2[i,paste0("anc_",x)]*H.exp[i,x]))/ir2[i,'allele_num'] ) 

phat  <- (sum(ir2[,83:90])*2 + sum(ir2[,75:82]))/sum(ir2[,67:74]) # (homs*2+hets)/anc
HT <- 2 * phat * (1-phat)
ir2$HI <- HI
ir2$F.IS <- (HS-HI)/HS
ir2$F.ST <- (HT-HS)/HT
ir2$F.IT <- (HT-HI)/HT


#################################
calcF <- function(p){
  
  P.A.Arab = p$homFreq_Arab+ p$hetFreq_Arab/2
  H.exp.Arab = 2*P.A.Arab*(1 - P.A.Arab)
  H.obs.Arab = (2*p$hets_Arab)/p$anc_Arab
  
  P.A.Azeri = p$homFreq_Azeri+ p$hetFreq_Azeri/2
  H.exp.Azeri = 2*P.A.Azeri*(1 - P.A.Azeri)
  H.obs.Azeri = (2*p$hets_Azeri)/p$anc_Azeri
  
  P.A.Baloch = p$homFreq_Baloch+ p$hetFreq_Baloch/2
  H.exp.Baloch = 2*P.A.Baloch*(1 - P.A.Baloch)
  H.obs.Baloch = (2*p$hets_Baloch)/p$anc_Baloch
  
  P.A.Kurd = p$homFreq_Kurd+ p$hetFreq_Kurd/2
  H.exp.Kurd = 2*P.A.Kurd*(1 - P.A.Kurd)
  H.obs.Kurd = (2*p$hets_Kurd)/p$anc_Kurd
  
  P.A.Lur = p$homFreq_Lur+ p$hetFreq_Lur/2
  H.exp.Lur = 2*P.A.Lur*(1 - P.A.Lur)
  H.obs.Lur = (2*p$hets_Lur)/p$anc_Lur
  
  P.A.per = p$homFreq_Persian+ p$hetFreq_Persian/2
  H.exp.per = 2*P.A.per*(1 - P.A.per)
  H.obs.per = (2*p$hets_Persian)/p$anc_Persian
  
  P.A.perG = p$`homFreq_Persian Gulf Islander`+ p$`hetFreq_Persian Gulf Islander`/2
  H.exp.perG = 2*P.A.perG*(1 - P.A.perG)
  H.obs.perG = (2*p$`hets_Persian Gulf Islander`)/p$`anc_Persian Gulf Islander`
  
  P.A.Turk = p$homFreq_Turkmen+ p$hetFreq_Turkmen/2
  H.exp.Turk = 2*P.A.Turk*(1 - P.A.Turk)
  H.obs.Turk = (2*p$hets_Turkmen)/p$anc_Turkmen
  
  ### HI calculation based on observed values
  HI <- (H.obs.Arab*p$anc_Arab+ H.obs.Azeri*p$anc_Azeri+ H.obs.Baloch*p$anc_Baloch+ 
           H.obs.Kurd*p$anc_Kurd+ H.obs.Lur*p$anc_Lur+ H.obs.per*p$anc_Persian+ 
           H.obs.perG*p$`anc_Persian Gulf Islander`+ H.obs.Turk*p$anc_Turkmen)/p$allele_num
  
  ### HS calculation based on expected values
  HS <- (H.exp.Arab*p$anc_Arab+ H.exp.Azeri*p$anc_Azeri+ H.exp.Baloch*p$anc_Baloch+ 
           H.exp.Kurd*p$anc_Kurd+ H.exp.Lur*p$anc_Lur+ H.exp.per*p$anc_Persian+ 
           H.exp.perG*p$`anc_Persian Gulf Islander`+ H.exp.Turk*p$anc_Turkmen)/p$allele_num
  
  ### HT calculation based on expected hetr. overall
  phat  <- (sum(as.data.frame(p)[83:90])*2 + sum(as.data.frame(p)[75:82]))/sum(as.data.frame(p)[67:74]) # (homs*2+hets)/anc
  HT <- 2*phat*(1-phat)
  
  F.IS <- (HS-HI)/HS
  F.ST <- (HT-HS)/HT
  F.IT <- (HT-HI)/HT
  return(c(F.IS,F.ST,F.IT,HI)) 
}
