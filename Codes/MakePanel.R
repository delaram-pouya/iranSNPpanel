library(parallel)
library('UpSetR')
library(tictoc)
setwd("/home/delaram/iranome/chromosomes")
ir= readRDS("cleanIranome.rds")
ir$weird= (ir$hetFreq_Arab > 1 | ir$hetFreq_Azeri>1 | ir$hetFreq_Baloch>1 | ir$hetFreq_Kurd>1 |
             ir$hetFreq_Lur>1| ir$hetFreq_Persian>1| ir$`hetFreq_Persian Gulf Islander`>1 | ir$hetFreq_Turkmen>1)


####### wierd detection ! H.obs>1 !! >> check iranome original data for these positions 
weird=sapply(tribes, function(x) ir[(2*ir[,paste0("hets_",x)])>ir[,paste0("anc_",x)],], simplify = F)
mapply(function(i,j) summary(i[,paste0('hetFreq_',names(weird)[j])]), weird, 1:length(weird)) #why het_freq>1 ???
w = subset(ir,hetFreq_Arab > 1 | hetFreq_Azeri>1 | hetFreq_Baloch>1 | hetFreq_Kurd>1 |
             hetFreq_Lur>1| hetFreq_Persian>1| `hetFreq_Persian Gulf Islander`>1)


ir= subset(ir, !weird)
####### Depth > min #Depth>32.5 
ir$Dsum = rowSums(ir[,paste0("D", seq(32.5, 97.5, 5))]) 
hist(ir$Dsum)
d.cut = quantile(ir$Dsum, 0.3) #103
#d.cut=400

####### Quality > min #Qual>37.5
ir$Qsum = rowSums(ir[,paste0("Q", seq(37.5, 97.5, 5))]) 
hist(ir$Qsum)
q.sorted = sort(ir$Qsum)
q.cut = q.sorted[round(length(q.sorted)*0.3)] #370


####### Allele_number > min #sequenced individuals
an.sorted = sort(ir$allele_num)
hist(ir$allele_num)
an.cut = an.sorted[round(length(an.sorted)*0.05)] #736


####### basic filtering  

ir<- ir[  ir$allele_freq>0.4 & ir$allele_freq<0.6 & 
            ir$category!= "lof_variant" & 
            ir$major_consequence %in% c("synonymous_variant", "intron_variant", "intergenic_variant" )&
            ir$alt %in% c("A","T","C","G") & ir$ref %in% c("A","T","C","G"),]

####### Fst & HI filter

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


####### distance >> this needs to be changed 
ir$dis=0
chr <- split( ir , f = ir$chrom )
dis.l = lapply(chr, function(i) sapply(1:nrow(i), function(j)  (i[j+1,'pos']-i[j,'pos']) ))  
chr.temp = mapply(function(i,x) chr[[i]]$dis <<- x, seq(length(chr)), dis.l)
ir <- do.call(rbind, chr)

####### sesitivity analysis > better ignore distance at first and check it at last 
d.cut = 103
q.cut=370
an.cut=736
f.cut = 0.06  
h.cut.low=0
h.cut.up = 0.4
dis.cut=1e5

l = list(subset(ir,Dsum>d.cut),subset(ir,Qsum>q.cut),subset(ir,allele_num>an.cut),
         subset(ir,F.ST<f.cut),subset(ir,HI<h.cut.up&HI>h.cut.low),subset(ir,dis>dis.cut))

l = lapply(l, function(i) sapply(1:nrow(i), function(j) paste(i[j,1],i[j,2],sep='_',collapse = NULL)))
uni = unique(unlist(l))
tab = as.data.frame(uni)
tab[,2:7] =sapply(1:length(l),function(j) {sapply(1:nrow(tab), function(i) ifelse( tab[i,'uni']%in% l[[j]],1, 0))})
colnames(tab)=c('uni','depth','quality','allele.num','Fst','Het','distance')

upset(tab,sets = c('depth','quality','allele.num','Fst','Het','distance') , 
      matrix.color = "#990000",
      sets.bar.color = c('darkslategrey','cyan4','cyan3','cyan3','cyan2','cyan1'), 
      sets.x.label = c('dep=103,q=370,an=736,Fst=.06,Het=.4,dist=1e5'),
      keep.order = F,nintersects=NA, point.size = 2.6,line.size = 0.7)


ir2<- subset(ir,Dsum>d.cut&Qsum>q.cut&allele_num>an.cut&F.ST<f.cut&HI>h.cut)
length(unique(sapply(1:nrow(ir2), function(i) paste(ir2[i,]$chrom,ir2[i,]$pos,sep = '_',collapse = NULL))))

####################   check Fst calculation !!!!!!!
ir[which.min(ir$HI),]
tmp=subset(ir,F.ST>=0)
plot(tmp$F.ST,tmp$HI)
nrow(subset(ir,F.ST<0)) ##why too much?



################
ir3<-subset(ir2, F.ST<f.cut.off & HI>h.cut.off & dis)
summary(ir3$site_quality)
dim(ir3)
View(head(ir3))


###################
