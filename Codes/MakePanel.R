
## in this script, we will define and apply some filters 
##    to the cleaned iranome dataset, in order to find a set of 
##    optimized SNPs for individual identification in iranian population


library(parallel)
library(UpSetR)
library(tictoc)


ir = readRDS("largeFiles/cleanIranome.rds")

## filtering weird positions in iranome 
##   dataset which have higher than 1 frequency value!

ir$weird = (ir$hetFreq_Arab > 1 | ir$hetFreq_Azeri>1 | ir$hetFreq_Baloch>1 | ir$hetFreq_Kurd>1 |
             ir$hetFreq_Lur > 1| ir$hetFreq_Persian>1| ir$`hetFreq_Persian Gulf Islander`>1 | ir$hetFreq_Turkmen>1)


ir <- subset(ir, !weird)



#### Depth filter

## number of individuals with higher 
##  than 32 depth in that specific snp position(out of 800 people)

ir$NumHighDepth = rowSums(ir[,paste0("D", seq(32.5, 97.5, 5))]) 
hist(ir$NumHighDepth)
MIN_NUMBER_OF_HIGH_DEPTH = quantile(ir$NumHighDepth, 0.3) #103
print(paste0('MIN_NUMBER_OF_HIGH_DEPTH is: ', MIN_NUMBER_OF_HIGH_DEPTH))
# MIN_NUMBER_OF_HIGH_DEPTH = 400




#### Quality score filter

## number of individuals with higher 
##  than 37 quality score in that specific snp position(out of 800 people)

ir$NumHighQuality = rowSums(ir[,paste0("Q", seq(37.5, 97.5, 5))]) 
hist(ir$NumHighQuality)
NumHighQuality_Sorted  = sort(ir$NumHighQuality)
MIN_NUMBER_OF_HIGH_QUALITY_SCORE = NumHighQuality_Sorted[round(length(NumHighQuality_Sorted)*0.3)] #370
print(paste0('MIN_NUMBER_OF_HIGH_QUALITY_SCORE is: ', MIN_NUMBER_OF_HIGH_QUALITY_SCORE))





#### Allele_Number filter 

## this value needs to be 1 for homo positions and 2 for hetero positions
##  theoretically, Summation of this value needs to at least 800(if all positions are homo)
##  if its value is less that 800, it shows that the position have not been sequenced in many individuals
##  we apply a filter on this value to remove positions which can not be(for any reason) sequenced efficiently


AlleleNumber_Sorted = sort(ir$allele_num)
hist(ir$allele_num)
MIN_SUM_ALLELE_NUMBER = AlleleNumber_Sorted[round(length(AlleleNumber_Sorted)*0.05)] #736
print(paste0('MIN_SUM_ALLELE_NUMBER is: ', MIN_SUM_ALLELE_NUMBER))





####  Basic Filters  

ir<- ir[  ir$allele_freq>0.4 & ir$allele_freq<0.6 & 
            ir$category!= "lof_variant" & 
            ir$major_consequence %in% c("synonymous_variant", "intron_variant", "intergenic_variant" )&
            ir$alt %in% c("A","T","C","G") & ir$ref %in% c("A","T","C","G"),]




#### Calculating F-statistics (Fst & HI filter)

## please refer to the link bellow for more information
##  about F-indices and different Heterozygosities 

## http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html


# P.A = P(A)
tribes <- unique(sub("acs_","", grep("acs_", colnames(ir), val=T)))
P.A <- sapply(tribes, function(x) (ir[,paste0("homFreq_",x)] + ir[,paste0("hetFreq_",x)])/2)
H.exp <- 2 * P.A * (1-P.A)
H.obs <- sapply(tribes, function(x) 2*ir[,paste0("hets_",x)]/ir[,paste0("anc_",x)])
HI <- sapply(1:nrow(H.obs), function(i) sum(sapply(tribes, function(x) ir[i,paste0("anc_",x)]*H.obs[i,x]))/ir[i,'allele_num'] ) 
HS <- sapply(1:nrow(H.exp), function(i) sum(sapply(tribes, function(x) ir[i,paste0("anc_",x)]*H.exp[i,x]))/ir[i,'allele_num'] ) 

phat  <- (sum(ir[,83:90])*2 + sum(ir[,75:82]))/sum(ir[,67:74]) # (homs*2+hets)/anc
HT <- 2 * phat * (1-phat)
ir$HI <- HI
ir$F.IS <- (HS-HI)/HS
ir$F.ST <- (HT-HS)/HT
ir$F.IT <- (HT-HI)/HT




## filtering positions with F.st values higher that 1 or less than -1 

ir_Fst_weird <- subset(ir, F.ST >1 | F.ST < -1)
ir <- subset(ir, F.ST <1 & F.ST> (-1) )
ir$F.ST <- abs(ir$F.ST)
  



## Visualizing relationship between HI, Fst, alleleFreq
#### Interpretation ??

plot(ir$F.ST, ir$allele_freq)
plot(ir$F.ST, ir$HI)
plot(ir$HI, ir$allele_freq)
plot(ir$F.IT, ir$F.IS)
plot(ir$F.IT, ir$F.ST)
plot(ir$F.IT, ir$HI)



#### Linkage Filter 

## SNPs need to be at least 500kb apart 
##  from each other to have low linkage

ir$distance = 0
ir_ChrSplit <- split( ir , f = ir$chrom )

sequentialDistance_ChrSplit = lapply(ir_ChrSplit, 
                                     function(x) sapply(1:nrow(x), function(i)  (x[i+1,'pos'] - x[i,'pos']) ))  

sapply(1:length(ir_ChrSplit), function(i){
  ir_ChrSplit[[i]]$distance <<- sequentialDistance_ChrSplit[[i]]
  })

ir <- do.call(rbind, ir_ChrSplit)





#### sensitivity analysis > better ignore distance at first and check it at last 

MIN_NUMBER_OF_HIGH_DEPTH = 103
MIN_NUMBER_OF_HIGH_QUALITY_SCORE = 370
MIN_SUM_ALLELE_NUMBER = 736
MAX_F_INDEX = 0.06  
MIN_INDIVIDUAL_HETEROZYGOSITY_VALUE = 0
MAX_INDIVIDUAL_HETEROZYGOSITY_VALUE = 0.4
MIN_SNP_POSITION_DISTANCE = 1e5


MonoFiltered_ir_list = list(subset(ir,NumHighDepth > MIN_NUMBER_OF_HIGH_DEPTH),
         subset(ir,NumHighQuality > MIN_NUMBER_OF_HIGH_QUALITY_SCORE),
         subset(ir,allele_num > MIN_SUM_ALLELE_NUMBER),
         subset(ir,F.ST < MAX_F_INDEX),
         subset(ir,HI < MAX_INDIVIDUAL_HETEROZYGOSITY_VALUE & HI > MIN_INDIVIDUAL_HETEROZYGOSITY_VALUE),
         subset(ir,distance > MIN_SNP_POSITION_DISTANCE) )




MonoFiltered_ir_listofIDs = lapply(MonoFiltered_ir_list, 
                                   function(x) paste0(x$chrom,sep='_',x$pos))

uniqueIDs = unique(unlist(MonoFiltered_ir_listofIDs))
TableID = as.data.frame(uniqueIDs)

TableID[,2:7] =sapply(1:length(MonoFiltered_ir_listofIDs),function(j) {
  sapply(1:nrow(TableID), function(i) ifelse( TableID[i,'uniqueIDs']%in% MonoFiltered_ir_listofIDs[[j]],1, 0))})

colnames(TableID)=c('snpID','depth','quality','allele.num','Fst','Het','distance')



upset(TableID,sets = c('depth','quality','allele.num','Fst','Het','distance') , 
      matrix.color = "#990000",
      sets.bar.color = c('darkslategrey','cyan4','cyan3','cyan3','cyan2','cyan1'), 
      sets.x.label = c('dep=103,q=370,an=736,Fst=.06,Het=.4,dist=1e5'),
      keep.order = F,nintersects=NA, point.size = 2.6,line.size = 0.7)



ir2<- subset(ir, 
             NumHighDepth > MIN_NUMBER_OF_HIGH_DEPTH &
               NumHighQuality > MIN_NUMBER_OF_HIGH_QUALITY_SCORE &
               allele_num > MIN_SUM_ALLELE_NUMBER &
               F.ST < MAX_F_INDEX & 
               HI > MIN_INDIVIDUAL_HETEROZYGOSITY_VALUE)



TargetPositions <- paste0(ir2$chrom, sep='_', ir2$pos)
length(TargetPositions)






####### wierd detection ! H.obs>1 !! >> check iranome original data for these positions 
weird=sapply(tribes, function(x) ir[(2*ir[,paste0("hets_",x)])>ir[,paste0("anc_",x)],], simplify = F)
mapply(function(i,j) summary(i[,paste0('hetFreq_',names(weird)[j])]), weird, 1:length(weird)) #why het_freq>1 ???


