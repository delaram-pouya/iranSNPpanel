setwd('~/delaram/data/SNPpanel')
library(data.table)
library(GenomicRanges)
library(IRanges)
library(parallel)
#
targets=fread('targeted_405.txt')
annot=fread('~/genecode.annot.gtf')
repeats=fread('~/iranome/hg38repeat.txt')
panel=read.csv('~/delaram/data/SNPpanel/snpPanel.csv',header=T)

colnames(targets)[1]='chr'
targets$chr=paste0('chr',targets$chr)
targets.GR = GRanges(seqnames=targets$chr,ranges=IRanges(start=targets$start,end=targets$end),AF=targets$AF.exome.All)
targets=unique(targets)#396 positions
targets$id=paste0(targets$chr,sep='_',targets$start)
panel$id=paste0(panel$chr,sep='_',panel$pos)

############################ Repeats
#### raw regions
colnames(repeats)=c('chr','start','end','type','V5','strand')
repeats$width=repeats$end-repeats$start
summary(repeats$width)
hist(repeats$width,col = 'blue')
repeats=subset(repeats,repeats$chr%in%paste0('chr',seq(1:22)))
repeats.GR = GRanges(seqnames=repeats$chr, IRanges(start=repeats$start,end=repeats$end),type=repeats$type)
Reps=as.data.frame(mergeByOverlaps(targets.GR,repeats.GR))
Reps=Reps[,colnames(Reps)%in%c('targets.GR.seqnames','targets.GR.start','repeats.GR.start','repeats.GR.end','AF','type','repeats.GR.width')]
nrow(Reps)
summary(Reps$repeats.GR.width)

#### edited regions
rep=repeats
rep$start=ifelse(rep$width>50,rep$start+20,rep$start)
rep$end=ifelse(rep$width>50,rep$end-20,rep$end)
rep$width=rep$end-rep$start
rep.GR = GRanges(seqnames=rep$chr, IRanges(start=rep$start,end=rep$end),type=rep$type)
rep.ed=as.data.frame(mergeByOverlaps(targets.GR,rep.GR))
rep.ed=rep.ed[,colnames(rep.ed)%in%c('targets.GR.seqnames','targets.GR.start','rep.GR.start','rep.GR.end','AF','type','rep.GR.width')]
nrow(rep.ed)
summary(rep.ed$rep.GR.width)
hist(rep.ed$rep.GR.width,col='yellow')

############################# gene_Body
colnames(annot)[1]<-'chr'
annot.GR = GRanges(seqnames=annot$chr, IRanges(start=annot$start,end=annot$end),gene_type=annot$gene_type) 
geneBody = as.data.frame(mergeByOverlaps(targets.GR,annot.GR))
geneBody=geneBody[,colnames(geneBody)%in%c('targets.GR.seqnames','targets.GR.start',
                                           'annot.GR.start','annot.GR.end','AF','gene_type')]
colnames(geneBody)=c('chr','pos','AF','annot.s','annot.e','gene_type')
geneBody$id=paste0(geneBody$chr,sep='_',geneBody$pos)
geneBody=(geneBody[!duplicated(geneBody[,'id']),])
table(geneBody$gene_type)
nrow(geneBody)*100/nrow(targets)
x=subset(targets,!id%in%geneBody$id)
table(x$id%in%panel$id)

