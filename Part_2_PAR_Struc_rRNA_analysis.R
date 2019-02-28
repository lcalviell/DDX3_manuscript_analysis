library("GenomicFeatures")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")
library(ggplot2)
library("gplots")
library("reshape2")
library("ggplot2")
library("cowplot")
library(gridExtra)


#IMPORTANT: SETWD TO DIRECTORY WITH FILES IN IT
#setwd("/bfd/lcalviel/data/riboseq/DDX3_paper_Markus/To_github/")

#LOAD DATA

load("data/gencode.v25.annotation.gtf_Rannot")
load("data/res_DE")
load("data/rRNA_clips")
load("data/deltaG_average_Robj")
load("data/POSTAR2_binding_human_Robj")
peaks<-read.table("data/PARAlyzer_clusters.csv",stringsAsFactors = F,sep=",",header = T)
des_all<-read.table("data/DE_table_ok.csv",sep="\t",header = T,stringsAsFactors = F)


#pick representative transcripts
txok<-GTF_annotation$cds_txs_coords
txok<-txok[txok$reprentative_boundaries]

des_all$tx_type<-des_all$tx_type2
des_all<-des_all[des_all$tx_type2!="Mixed",]

peak_gr<-GRanges(seqnames = peaks$Chromosome,strand = peaks$Strand,ranges = IRanges(start=peaks$ClusterStart,end=peaks$ClusterEnd))
mcols(peak_gr)<-DataFrame(peaks[,colnames(peaks)[!colnames(peaks)%in%c("Chromosome","Strand","ClusterStart","ClusterEnd")]])
peak_gr$ClusterSequence<-DNAStringSet(peak_gr$ClusterSequence)

regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, GTF_annotation$threeutrs,
                GTF_annotation$ncIsof, GTF_annotation$ncRNAs, GTF_annotation$introns, GTF_annotation$intergenicRegions)
names(regions) <- c("cds","fiveutrs","threeutrs",
                    "ncIsof","ncRNAs","introns","intergenic")

regions<-GRangesList(lapply(regions,function(x){mcols(x)<-NULL;x}))

peak_gr$conv_spec<-peak_gr$ConversionEventCount/(peak_gr$ConversionEventCount+peak_gr$NonConversionEventCount)


des_all<-des_all[grep(des_all$biotype,pattern = "pseudo|snRNA|snoRNA",invert = T),]

des_all$biotype[des_all$biotype!="protein_coding"]<-"non_coding"


ok1<-which(des_all$baseMean_rna>20 & des_all$baseMean_ribo>20)
des_all<-des_all[ok1,]

list_DE_genes<-split(des_all$gene_id,des_all$tx_type)

peak_gr_ok<-peak_gr[peak_gr$ReadCount>10]

regions <- list(reduce(unlist(GTF_annotation$exons_txs[names(GTF_annotation$cds_txs)])),GTF_annotation$ncIsof, GTF_annotation$ncRNAs, GTF_annotation$introns, GTF_annotation$intergenicRegions,invertStrand(GTF_annotation$genes))
names(regions) <- c("Coding RNAs",
                    "Non-coding\nRNA isoforms","Non-coding\n genes","Introns","Intergenic space","Antisense")

regions<-GRangesList(lapply(regions,function(x){mcols(x)<-NULL;x}))

peak_gr_ok$region<-"undefined"
for(i in rev(names(regions))){
    peak_gr_ok$region[peak_gr_ok%over%regions[[i]]]<-i
}

ddd<-data.frame(mcols(peak_gr_ok))
ddd$regions<-factor(ddd$region,levels = names(regions))
colssi<-c("red","orange","gold","cornflowerblue","blue","dark blue")
barpeaks<-ggplot(ddd,aes(x=regions,fill=regions)) + geom_bar() +
    scale_fill_manual(values = colssi) +
    theme_classic() +
    ylab("number of peaks") +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) + 
    theme(legend.position="none")

pdf(file = "figures/peak_stats_ddx3.pdf",width = 8,height = 4)
print(barpeaks)
dev.off()

#mRNA BINDING

maptx<-mapToTranscripts(peak_gr_ok,transcripts = GTF_annotation$exons_txs[as.vector(seqnames(txok))])
maptx$score<-peak_gr_ok$ModeScore[maptx$xHits]
maptx$conv_spec<-peak_gr_ok$conv_spec[maptx$xHits]
maptx$sequence<-peak_gr_ok$ClusterSequence[maptx$xHits]
maptx$gc<-letterFrequency(maptx$sequence,letters = "GC",as.prob = T)*100
maptx$c_freq<-letterFrequency(maptx$sequence,letters = "C",as.prob = T)*100
maptx$g_freq<-letterFrequency(maptx$sequence,letters = "G",as.prob = T)*100

seqlevels(maptx)<-names(GTF_annotation$exons_txs[as.vector(seqnames(txok))])
seqlengths(maptx)<-sum(width(GTF_annotation$exons_txs[as.vector(seqnames(txok))]))

strand(maptx)<-"+"
tx_st_cds<-start(txok)[match(as.character(seqnames(maptx)),names(txok))]
tx_end_cds<-end(txok)[match(as.character(seqnames(maptx)),names(txok))]
tx_wdt<-width(txok)[match(as.character(seqnames(maptx)),names(txok))]
tx_end_tx<-seqlengths(txok)[match(as.character(seqnames(maptx)),names(seqlengths(txok)))]

maptx$tx_position<-"5_UTR"
maptx$tx_position[mid(maptx)>=tx_st_cds]<-"CDS"
maptx$tx_position[mid(maptx)>tx_end_cds]<-"3_UTR"


maptx$tx_bin<-as.integer((mid(maptx)/(tx_st_cds))*33)+1
cdspos<-mid(maptx)>=tx_st_cds & mid(maptx)<=tx_end_cds
maptx$tx_bin[cdspos]<-as.integer(((mid(maptx)[cdspos]-tx_st_cds[cdspos])/(tx_wdt[cdspos]))*100)+1
threeutpos<-mid(maptx)>tx_end_cds
maptx$tx_bin[threeutpos]<-as.integer(((mid(maptx)[threeutpos]-tx_end_cds[threeutpos])/(tx_end_tx[threeutpos]-tx_wdt[threeutpos]))*66)+1


mlt<-melt(list_DE_genes)
nnnms<-GTF_annotation$trann$gene_name[match(mlt$value,GTF_annotation$trann$gene_id)]
mlt$L1[intersect(which(mlt$L1=="TE_up"),grep(nnnms,pattern = "HIST|H2B"))]<-"TE_up_Hist"
mlt$L1[intersect(which(mlt$L1=="TE_up"),grep(nnnms,pattern = "HIST|H2B",invert = T))]<-"TE_up_noHist"
mlt$tx_id<-as.character(seqnames(txok))[match(mlt$value,txok$gene_id)]
mlt_pos<-mlt[mlt$L1%in%c("five_up","five_down"),]
maptx$tx_class<-mlt$L1[match(as.character(seqnames(maptx)),mlt$tx_id)]
maptx$tx_class[is.na(maptx$tx_class)]<-"Not_significant"
maptx$tx_class[as.character(seqnames(maptx))%in%mlt$tx_id[mlt$L1=="TE_down"]]<-"TE_down"

maptx$tx_class_pos<-mlt_pos$L1[match(as.character(seqnames(maptx)),mlt_pos$tx_id)]
maptx$tx_class_pos[is.na(maptx$tx_class_pos)]<-"Not_significant"

aggg<-aggregate(maptx$conv_spec,list(as.vector(seqnames(maptx))),sum)
agg_sc<-aggg[,2]
names(agg_sc)<-aggg[,1]
maptx$norm_conv_spec<-maptx$conv_spec/agg_sc[as.vector(seqnames(maptx))]


aggg<-aggregate(maptx$score*width(maptx),list(as.vector(seqnames(maptx))),sum)
agg_sc<-aggg[,2]
names(agg_sc)<-aggg[,1]
maptx$norm_cov_sum<-maptx$score/agg_sc[as.vector(seqnames(maptx))]
aggg<-aggregate(maptx$score*width(maptx),list(as.vector(seqnames(maptx))),max)
agg_sc<-aggg[,2]
names(agg_sc)<-aggg[,1]
maptx$norm_cov_01<-maptx$score/agg_sc[as.vector(seqnames(maptx))]

maptx$tx_bin_ok<-maptx$tx_bin

maptx$tx_bin_ok[maptx$tx_position=="CDS"]<-maptx$tx_bin[maptx$tx_position=="CDS"]+33
maptx$tx_bin_ok[maptx$tx_position=="3_UTR"]<-maptx$tx_bin[maptx$tx_position=="3_UTR"]+133

maptx$conversion_specificity<-"low"
maptx$conversion_specificity[maptx$conv_spec>=.15 & maptx$conv_spec<=.2]<-"medium"
maptx$conversion_specificity[maptx$conv_spec>.2]<-"high"
maptx$tx_class2<-maptx$tx_class
maptx$tx_class2[grep(maptx$tx_class2,pattern = "TE_up")]<-"TE_up"



dfs<-as.data.frame(mcols(maptx))
dfs$tx_position<-factor(dfs$tx_position,levels=c("5_UTR","CDS","3_UTR"))
d<-ggplot(dfs,aes(conv_spec,G.C))+ 
    geom_hex(bins = 50) +
    ylab("%GC content") +
    xlab("T>C conversion\nspecificity") +
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
pdf(file = "figures/par_clip_gc.pdf",width = 5,height = 4)
print(d +scale_fill_viridis_c())
dev.off()



dddf<-data.frame(mcols(maptx))
dddf$conversion_specificity<-factor(dddf$conversion_specificity,levels=c("low","medium","high"))

histpl<-ggplot(dddf,aes(x=conv_spec,fill=conversion_specificity)) + geom_histogram(bins = 1000) +
    geom_vline(xintercept =.15,col="black",lty=3)+
    scale_color_manual(values = alpha(c("dark grey","orange","red"),.7),"") +
    scale_fill_manual(values = alpha(c("dark grey","orange","red"),.7),"") +
    geom_vline(xintercept =.2,col="black",lty=3)+ 
    theme_classic() +
    ylab("count") +
    xlab("T>C conversion\nspecificity") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))
pdf(file = "figures/aggregate_plot_conv_spec_hist.pdf",width = 8,height = 4)
histpl + scale_x_continuous(breaks=c(0,.15,.20,.5,1), labels=c(0,.15,.20,.5,1))
dev.off()

agga2<-aggregate(list(maptx$conv_spec,maptx$norm_cov_sum),by=list(maptx$tx_bin_ok,maptx$tx_class2),mean)
colnames(agga2)<-c("tx_bin","tx_class","conv_spec","norm_cov_sum")

conv_specplot<-ggplot(agga2,aes(x=tx_bin,y=conv_spec,group=tx_class,color=tx_class)) + 
    geom_smooth(method="loess",span=.2,se = T) +  
    geom_point(size=.1) + 
    ylab("T>C conversion\nspecificity (average)") +
    xlab("") +
    geom_vline(xintercept =34,col="black",lty=3)+
    scale_color_manual(values = alpha(c("black","cornflowerblue","orange","blue","dark red"),.7),"Tx_class") + 
    geom_vline(xintercept =134,col="black",lty=3)+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,34,134,200), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "figures/conv_spec_tx_class.pdf",width = 8,height = 4)
print(conv_specplot)
dev.off()


agga2<-aggregate(list(maptx$conv_spec,maptx$norm_cov_sum,maptx$gc),by=list(maptx$tx_bin_ok,maptx$tx_class2),mean)
colnames(agga2)<-c("tx_bin","tx_class","conv_spec","norm_cov_sum","GC")

conv_specplot2<-ggplot(agga2,aes(x=tx_bin,y=conv_spec,group=tx_class,color=tx_class,size=GC)) + 
    geom_point() + 
    ylab("T>C conversion\nspecificity (average)") +
    xlab("") +
    geom_vline(xintercept =34,col="black",lty=3)+
    scale_color_manual(values = alpha(c("black","cornflowerblue","orange","blue","dark red"),.7),"Tx_class") + 
    geom_vline(xintercept =134,col="black",lty=3)+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,34,134,200), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "/bfd/lcalviel/data/riboseq/DDX3_paper_Markus/RNA_up_correct/conv_spec_tx_class_size.pdf",width = 8,height = 4)
print(conv_specplot2 + scale_size_continuous(name="%GC content", 
                                             range = c(0,5), 
                                             breaks = c(20,40,60,80,90),labels = as.character(c(20,40,60,80,90))))

dev.off()



conv_specplot2<-ggplot(agga2,aes(x=tx_bin,y=GC,group=tx_class,color=tx_class,size=conv_spec)) + 
    geom_point() + 
    ylab("%GC content\n(average)") +
    xlab("") +
    geom_vline(xintercept =34,col="black",lty=3)+
    scale_color_manual(values = alpha(c("black","cornflowerblue","orange","blue","dark red"),.7),"Tx_class") + 
    geom_vline(xintercept =134,col="black",lty=3)+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,34,134,200), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "figures/conv_spec_GC_tx_class.pdf",width = 8,height = 4)
print(conv_specplot2 +scale_size_continuous(name="T>C conversion\nspecificity", 
                                            range = c(0,6), 
                                            breaks = c(.1,.2,.3,.4,.6,.8,.9),labels = as.character(c(.1,.2,.3,.4,.6,.8,.9))))

dev.off()




agga<-aggregate(list(maptx$norm_cov_sum,maptx$norm_cov_01),by=list(maptx$tx_bin_ok,maptx$tx_class,maptx$conversion_specificity),sum)
colnames(agga)<-c("tx_bin","tx_class","conv_spec","norm_cov_sum","norm_cov_01")

agga_pos<-aggregate(list(maptx$norm_cov_sum,maptx$norm_cov_01),by=list(maptx$tx_bin_ok,maptx$tx_class_pos,maptx$conversion_specificity),sum)
colnames(agga_pos)<-c("tx_bin","tx_class","conv_spec","norm_cov_sum","norm_cov_01")

covvs<-agga_pos
agg_covv<-aggregate(list(covvs$norm_cov_sum,covvs$norm_cov_01),by=list(covvs$tx_class,covvs$conv_spec),sum)
colnames(agg_covv)<-c("tx_class","conv_spec","norm_cov_sum","norm_cov_01")
agg_covv$id<-paste(agg_covv$tx_class,agg_covv$conv_spec,sep="_")
covvs$id<-paste(covvs$tx_class,covvs$conv_spec,sep="_")

covvs$norm_cov_sum_norm<-covvs$norm_cov_sum/agg_covv$norm_cov_sum[match(covvs$id,agg_covv$id)]
covvs$norm_cov_01_norm<-covvs$norm_cov_01/agg_covv$norm_cov_01[match(covvs$id,agg_covv$id)]

covvs$conv_spec<-factor(covvs$conv_spec,levels=c("low","medium","high"))

agga_all<-aggregate(list(maptx$norm_cov_sum,maptx$norm_cov_01),by=list(maptx$tx_bin_ok,maptx$conversion_specificity),sum)
colnames(agga_all)<-c("tx_bin","conv_spec","norm_cov_sum","norm_cov_01")
covvs<-agga_all

agg_covv<-aggregate(list(covvs$norm_cov_sum,covvs$norm_cov_01),by=list(covvs$conv_spec),sum)
colnames(agg_covv)<-c("conv_spec","norm_cov_sum","norm_cov_01")

covvs$norm_cov_sum_norm<-covvs$norm_cov_sum/agg_covv$norm_cov_sum[match(covvs$conv_spec,agg_covv$conv_spec)]
covvs$norm_cov_01_norm<-covvs$norm_cov_01/agg_covv$norm_cov_01[match(covvs$conv_spec,agg_covv$conv_spec)]

covvs$conv_spec<-factor(covvs$conv_spec,levels=c("low","medium","high"))

covpl<-ggplot(covvs,aes(x=tx_bin,y=norm_cov_sum_norm,group=conv_spec,color=conv_spec)) +
    geom_vline(xintercept =34,col="black",lty=3)+
    scale_color_manual(values = alpha(c("dark grey","orange","red"),.7),"T>C conversion\nspecificity") + geom_line(size=1.1) +
    geom_vline(xintercept =134,col="black",lty=3)+ 
    theme_classic() +
    ylab("Aggregate\nPAR-CLIP peak score") +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,34,134,200), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "figures/aggregate_plot_conv_spec.pdf",width = 8,height = 4)
print(covpl)
dev.off()

names(mcols(gr_peaks_rbps))<-c("RBP","method","cell","score")
paralyz_hek<-gr_peaks_rbps[gr_peaks_rbps$method=="PAR-CLIP,PARalyzer" & gr_peaks_rbps$cell%in%c("HEK293","HEK293T")]

paralyz_hek<-split(paralyz_hek,paralyz_hek$RBP)
paralyz_hek<-endoapply(paralyz_hek,sort)
peak_dd<-peak_gr_ok

mcols(peak_dd)<-mcols(peak_dd)[,c("ClusterID","ReadCount","ConversionEventCount","ModeScore")]
colnames(mcols(peak_dd))<-colnames(mcols(paralyz_hek[[1]]))
peak_dd$RBP<-"DDX3X"
peak_dd$cell<-"HEK293"
peak_dd$method<-"PAR-CLIP,PARalyzer"

paralyz_hek[["DDX3X"]]<-peak_dd
res_clips_meta<-sapply(paralyz_hek,function(x){
    maptx<-mapToTranscripts(x,transcripts = GTF_annotation$exons_txs[as.vector(seqnames(txok))])
    maptx$score<-x$score[maptx$xHits]
    seqlevels(maptx)<-names(GTF_annotation$exons_txs[as.vector(seqnames(txok))])
    seqlengths(maptx)<-sum(width(GTF_annotation$exons_txs[as.vector(seqnames(txok))]))
    
    strand(maptx)<-"+"
    tx_st_cds<-start(txok)[match(as.character(seqnames(maptx)),names(txok))]
    tx_end_cds<-end(txok)[match(as.character(seqnames(maptx)),names(txok))]
    tx_wdt<-width(txok)[match(as.character(seqnames(maptx)),names(txok))]
    tx_end_tx<-seqlengths(txok)[match(as.character(seqnames(maptx)),names(seqlengths(txok)))]
    
    maptx$tx_position<-"5_UTR"
    maptx$tx_position[mid(maptx)>=tx_st_cds]<-"CDS"
    maptx$tx_position[mid(maptx)>tx_end_cds]<-"3_UTR"
    
    
    maptx$tx_bin<-as.integer((mid(maptx)/(tx_st_cds))*33)+1
    cdspos<-mid(maptx)>=tx_st_cds & mid(maptx)<=tx_end_cds
    maptx$tx_bin[cdspos]<-as.integer(((mid(maptx)[cdspos]-tx_st_cds[cdspos])/(tx_wdt[cdspos]))*100)+1
    threeutpos<-mid(maptx)>tx_end_cds
    maptx$tx_bin[threeutpos]<-as.integer(((mid(maptx)[threeutpos]-tx_end_cds[threeutpos])/(tx_end_tx[threeutpos]-tx_wdt[threeutpos]))*66)+1
    
    aggg<-aggregate(maptx$score*width(maptx),list(as.vector(seqnames(maptx))),sum)
    agg_sc<-aggg[,2]
    names(agg_sc)<-aggg[,1]
    maptx$norm_cov_sum<-maptx$score/agg_sc[as.vector(seqnames(maptx))]
    aggg<-aggregate(maptx$score*width(maptx),list(as.vector(seqnames(maptx))),max)
    agg_sc<-aggg[,2]
    names(agg_sc)<-aggg[,1]
    maptx$norm_cov_01<-maptx$score/agg_sc[as.vector(seqnames(maptx))]
    
    maptx$tx_bin_ok<-maptx$tx_bin
    
    maptx$tx_bin_ok[maptx$tx_position=="CDS"]<-maptx$tx_bin[maptx$tx_position=="CDS"]+33
    maptx$tx_bin_ok[maptx$tx_position=="3_UTR"]<-maptx$tx_bin[maptx$tx_position=="3_UTR"]+133
    
    agga_all<-aggregate(list(maptx$norm_cov_sum,maptx$norm_cov_01),by=list(maptx$tx_bin_ok),sum)
    
    colnames(agga_all)<-c("tx_bin","norm_cov_sum","norm_cov_01")
    covvs<-agga_all
    covv<-rep(0,200)
    covv[covvs[,1]]<-covvs[,2]
    covv
})

res_clips_meta_ok<-apply(res_clips_meta,2,function(x){x/sum(x)})

res_clips_meta_ok<-apply(res_clips_meta,2,function(x){
    x<-x/sum(x)
    x[x>.03]<-.03
    x<-x/sum(x)
    x
})

big<-which(res_clips_meta_ok>.03,arr.ind = T)
res_clips_meta_ok[big]<-.03

pdf(file = "figures/aggregate_plot_PARCLIPs.pdf",width = 25,height = 20)
heatmap.2(t(res_clips_meta_ok),Colv = NULL,trace='none')
dev.off()

res_clips_meta_ok<-apply(res_clips_meta,2,function(x){x/sum(x)})


prots=c("DDX3X","EIF3B","FMR1","MOV10")
dat<-melt(data.frame(res_clips_meta_ok))
dat$position<-rep(1:200,dim(res_clips_meta_ok)[2])
dat2<-dat[dat$variable%in%prots,]
dat2$variable<-factor(dat2$variable,prots)
covpl<-ggplot(dat2,aes(x=position,y=value,group=variable,color=variable)) +
    geom_vline(xintercept =34,col="black",lty=3)+
    scale_color_manual(values = alpha(c("red","orange","forestgreen","blue"),.7),"RBP") + 
    geom_line(size=1.1) +
    geom_vline(xintercept =134,col="black",lty=3)+ 
    theme_classic() +
    ylab("Aggregate\nPAR-CLIP peak score") +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,34,134,200), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "figures/aggregate_plot_clipDB_scoring.pdf",width = 8,height = 4)
print(covpl)
dev.off()



eclips<-gr_peaks_rbps[gr_peaks_rbps$method=="eCLIP"]

eclips$id<-paste(eclips$RBP,eclips$cell,sep=";")
eclips<-split(eclips,eclips$id)
res_eclips_meta<-sapply(eclips,function(x){
    maptx<-mapToTranscripts(resize(x,1,"center"),transcripts = GTF_annotation$exons_txs[as.vector(seqnames(txok))])
    maptx$score<-x$score[maptx$xHits]
    seqlevels(maptx)<-names(GTF_annotation$exons_txs[as.vector(seqnames(txok))])
    seqlengths(maptx)<-sum(width(GTF_annotation$exons_txs[as.vector(seqnames(txok))]))
    
    strand(maptx)<-"+"
    tx_st_cds<-start(txok)[match(as.character(seqnames(maptx)),names(txok))]
    tx_end_cds<-end(txok)[match(as.character(seqnames(maptx)),names(txok))]
    tx_wdt<-width(txok)[match(as.character(seqnames(maptx)),names(txok))]
    tx_end_tx<-seqlengths(txok)[match(as.character(seqnames(maptx)),names(seqlengths(txok)))]
    
    maptx$tx_position<-"5_UTR"
    maptx$tx_position[mid(maptx)>=tx_st_cds]<-"CDS"
    maptx$tx_position[mid(maptx)>tx_end_cds]<-"3_UTR"
    
    
    maptx$tx_bin<-as.integer((mid(maptx)/(tx_st_cds))*33)+1
    cdspos<-mid(maptx)>=tx_st_cds & mid(maptx)<=tx_end_cds
    maptx$tx_bin[cdspos]<-as.integer(((mid(maptx)[cdspos]-tx_st_cds[cdspos])/(tx_wdt[cdspos]))*100)+1
    threeutpos<-mid(maptx)>tx_end_cds
    maptx$tx_bin[threeutpos]<-as.integer(((mid(maptx)[threeutpos]-tx_end_cds[threeutpos])/(tx_end_tx[threeutpos]-tx_wdt[threeutpos]))*66)+1
    
    aggg<-aggregate(maptx$score*width(maptx),list(as.vector(seqnames(maptx))),sum)
    agg_sc<-aggg[,2]
    names(agg_sc)<-aggg[,1]
    maptx$norm_cov_sum<-maptx$score/agg_sc[as.vector(seqnames(maptx))]
    aggg<-aggregate(maptx$score*width(maptx),list(as.vector(seqnames(maptx))),max)
    agg_sc<-aggg[,2]
    names(agg_sc)<-aggg[,1]
    maptx$norm_cov_01<-maptx$score/agg_sc[as.vector(seqnames(maptx))]
    
    maptx$tx_bin_ok<-maptx$tx_bin
    
    maptx$tx_bin_ok[maptx$tx_position=="CDS"]<-maptx$tx_bin[maptx$tx_position=="CDS"]+33
    maptx$tx_bin_ok[maptx$tx_position=="3_UTR"]<-maptx$tx_bin[maptx$tx_position=="3_UTR"]+133
    
    agga_all<-aggregate(list(maptx$norm_cov_sum,maptx$norm_cov_01),by=list(maptx$tx_bin_ok),sum)
    
    colnames(agga_all)<-c("tx_bin","norm_cov_sum","norm_cov_01")
    covvs<-agga_all
    covv<-rep(0,200)
    covv[covvs[,1]]<-covvs[,2]
    covv
})


res_eclips_meta_ok<-apply(res_eclips_meta,2,function(x){x/sum(x)})
res_eclips_meta_ok_heat<-res_eclips_meta_ok
big<-which(res_eclips_meta_ok_heat>.03,arr.ind = T)
res_eclips_meta_ok_heat[big]<-.03
pdf(file = "figures/aggregate_plot_eCLIPs.pdf",width = 25,height = 20)
heatmap.2(t(res_eclips_meta_ok_heat),Colv = NULL,trace='none')
dev.off()


prots=c("DDX3X.K562","GEMIN5.K562","RPS3.K562","DDX6.K562")
dat<-melt(data.frame(res_eclips_meta_ok))
dat$position<-rep(1:200,dim(res_eclips_meta_ok)[2])
dat2<-dat[dat$variable%in%prots,]
dat2$variable<-factor(dat2$variable,prots)
#34 vs 35?
covpl<-ggplot(dat2,aes(x=position,y=value,group=variable,color=variable)) +
    geom_vline(xintercept =34,col="black",lty=3)+
    scale_color_manual(values = alpha(c("red","orange","forestgreen","blue"),.7),"RBP") + 
    geom_line(size=1.1) +
    geom_vline(xintercept =134,col="black",lty=3)+ 
    theme_classic() +
    ylab("Aggregate\neCLIP peak score") +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,34,134,200), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "figures/aggregate_plot_clipDB_scoring_eclip.pdf",width = 8,height = 4)
print(covpl)
dev.off()



maptx2<-maptx

tabbo<-table(des_all$tx_type)
tabbo["Not_significant"]<-sum(tabbo)
names(tabbo)<-c("all_genes","RNA_down","RNA_up","TE_down","TE_up")
tabbo<-tabbo[c("all_genes","RNA_up","RNA_down","TE_up","TE_down")]

matta<-unique(cbind(as.character(seqnames(maptx2)),maptx2$conversion_specificity,maptx2$tx_class2))

all<-as.numeric(tabbo["all_genes"])
all_par<-length(unique(matta[,1]))
all_par_pct<-round((all_par/all)*100,digits = 1)
all_par_tc<-length(unique(matta[matta[,2]=="high",1]))
all_par_tc_pct<-round((all_par_tc/all_par)*100,digits = 1)
all_par_tc_pct_all<-round((all_par_tc/all)*100,digits = 1)


list_matta<-list()
for(i in unique(matta[,3])){
    
    class<-as.numeric(tabbo[i])
    class_par<-length(unique(matta[matta[,3]==i,1]))
    class_par_pct<-round((class_par/class)*100,digits = 1)
    class_par_tc<-length(unique(matta[matta[,3]==i & matta[,2]=="high",1]))
    class_par_tc_pct<-round((class_par_tc/class_par)*100,digits = 1)
    class_par_tc_pct_class<-round((class_par_tc/class)*100,digits = 1)
    list_matta[[i]]<-as.character(c(class,class_par,class_par_pct,class_par_tc,class_par_tc_pct,class_par_tc_pct_class))
}
list_matta[["all_genes"]]<-as.character(c(all,all_par,all_par_pct,all_par_tc,all_par_tc_pct,all_par_tc_pct_all))

list_matta<-list_matta[c("all_genes","RNA_up","RNA_down","TE_up","TE_down")]

list_matta<-lapply(list_matta,function(x){
    c(x[1],paste(x[2]," (",x[3],"%)",sep=""),paste(x[4]," (",x[5],"%)",sep=""))
})

tabla<-t(do.call(args = list_matta,what = rbind.data.frame))
colnames(tabla)<-names(list_matta)
rownames(tabla)<-c("Detected genes","DDX3-bound\n(% of Detected genes)","DDX3-bound, high T>C\n(% of DDX3-bound)")
pdf(file = "figures/table_bound_genes.pdf",width = 8,height = 3)
temo2<-ttheme_default(base_family = "Helvetica")
taboo<-tableGrob(tabla,theme=temo2)
grid.arrange(taboo)
dev.off()

tabla<-as.data.frame(tabla)
tabla$category<-gsub(rownames(tabla),pattern = "\n",replacement = " ")
tabla<-tabla[,c("category","all_genes","RNA_up","RNA_down","TE_up","TE_down")]
rownames(tabla)<-NULL
write.table(tabla,file = "data/table_bound_genes.csv",sep="\t",quote = FALSE,row.names = FALSE)


#rRNA



minpos_rdna<-3000
maxpos_rdna<-13000
par1<-melt(lapply(list_clips_rRNA$PAR[[1]],function(x){as.vector(x[[1]])[minpos_rdna:maxpos_rdna]}))
par2<-melt(lapply(list_clips_rRNA$PAR[[2]],function(x){as.vector(x[[1]])[minpos_rdna:maxpos_rdna]}))
par1$replicate<-"PAR-CLIP replicate_1"
par2$replicate<-"PAR-CLIP replicate_2"
parna<-rbind(par1,par2)
parna$position<-rep(minpos_rdna:maxpos_rdna,length(table(par1$L1))+length(table(par2$L1)))
maxxa<-max(parna$value)
parna<-parna[parna$L1!="all_reads",]
parna$mutation<-parna$L1

parna$mutation[!parna$mutation%in%c("T>C","None")]<-"other"
parna$mutation<-factor(parna$mutation,levels=c("None","other","T>C"))
colrna<-c("grey","orange","navyblue")

pos_18s<-c(3657,5527)
pos_5s<-c(6623,6779)
pos_28s<-c(7935,12969)
parnagg<-aggregate(parna$value,by=list(parna$position,parna$replicate),sum)
colnames(parnagg)<-c("position","replicate","value")
parnagg_id<-paste(parnagg$position,parnagg$replicate,sep="_")
parnaid<-paste(parna$position,parna$replicate,sep="_")
parna$norm_value<-parna$value/parnagg$value[match(parnaid,parnagg_id)]
parna$norm_value_rep<-parna$value/parnagg$value[match(parnaid,parnagg_id)]


parrnapl_reps<-ggplot(parna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    facet_wrap(replicate~.,nrow = 2) +
    ylab("Read coverage") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))

parrnapl_noreps<-ggplot(parna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    ylab("PAR-CLIP\nRead coverage") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))


parnatc<-parna[parna$mutation=="T>C",]

parrnapl_norepstc<-ggplot(parnatc,aes(x=position,y=value)) + geom_bar(stat="identity",width=1) +
    ylab("PAR-CLIP\n T>C coverage") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))


par1<-melt(lapply(list_clips_rRNA$PAR[[1]],function(x){as.vector(x[[1]])[minpos_rdna:maxpos_rdna]}))
par2<-melt(lapply(list_clips_rRNA$PAR[[2]],function(x){as.vector(x[[1]])[minpos_rdna:maxpos_rdna]}))

iclip1<-melt(lapply(list_clips_rRNA$iCL[[1]],function(x){as.vector(x[[1]])[minpos_rdna:maxpos_rdna]}))
iclip2<-melt(lapply(list_clips_rRNA$iCL[[2]],function(x){as.vector(x[[1]])[minpos_rdna:maxpos_rdna]}))
iclip1$replicate<-"iCLIP replicate_1"
iclip2$replicate<-"iCLIP replicate_2"
icliprna<-rbind(iclip1,iclip2)
icliprna$position<-rep(minpos_rdna:maxpos_rdna,length(table(iclip1$L1))+length(table(iclip2$L1)))
maxxa<-max(icliprna$value)
icliprna<-icliprna[icliprna$L1!="all_reads",]
icliprna$mutation<-icliprna$L1


icliprna$mutation[!icliprna$mutation%in%c("T>C","None")]<-"other"
icliprna$mutation<-factor(icliprna$mutation,levels=c("None","other","T>C"))
colrna<-c("grey","orange","navyblue")

pos_18s<-c(3657,5527)
pos_5s<-c(6623,6779)
pos_28s<-c(7935,12969)
icliprna$replicate<-factor(icliprna$replicate)

icliprnapl_reps<-ggplot(icliprna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    facet_wrap(replicate~.,nrow = 2) +
    ylab("Read coverage") +
    xlab("position along rDNA") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))


icliprnapl_noreps<-ggplot(icliprna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    ylab("iCLIP\nRead coverage") +
    xlab("position along rDNA") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))


icliprnapl_norepsnomut<-ggplot(icliprna,aes(x=position,y=value)) + geom_bar(stat="identity",width=1) +
    ylab("iCLIP\nRead coverage") +
    xlab("position along rDNA") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))

pdf(file = "figures/new_rDNA_othermuts.pdf",width = 8,height = 4)
plot_grid(parrnapl_noreps,icliprnapl_noreps,ncol = 1,align = "v")
dev.off()
png(filename = "figures/new_rDNA_othermuts.png",width = 8,height = 4,units = "in",res=300)
plot_grid(parrnapl_noreps,icliprnapl_noreps,ncol = 1,align = "v")
dev.off()



pdf(file = "figures/new_rDNA_tconly.pdf",width = 8,height = 4)
plot_grid(parrnapl_norepstc,icliprnapl_norepsnomut,ncol = 1,align = "v")
dev.off()

png(filename = "figures/new_rDNA_tconly.png",width = 8,height = 4,units = "in",res=300)
plot_grid(parrnapl_norepstc,icliprnapl_norepsnomut,ncol = 1,align = "v")
dev.off()


pdf(file = "figures/new_rDNA_replicates_othermuts.pdf",width = 8,height = 7)
plot_grid(parrnapl_reps,icliprnapl_reps,ncol = 1,align = "v")
dev.off()

png(filename = "figures/new_rDNA_replicates_othermuts.png",width = 8,height = 7,units = "in",res=300)
plot_grid(parrnapl_reps,icliprnapl_reps,ncol = 1,align = "v")
dev.off()

icliprna$mutation<-icliprna$L1
icliprna$mutation[!icliprna$mutation%in%c("T>C")]<-"No_T>C"
icliprna$mutation<-factor(icliprna$mutation,levels=c("No_T>C","T>C"))
parna$mutation<-parna$L1
parna$mutation[!parna$mutation%in%c("T>C")]<-"No_T>C"
parna$mutation<-factor(parna$mutation,levels=c("No_T>C","T>C"))
colrna<-c("grey","navyblue")


icliprnapl_reps<-ggplot(icliprna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    facet_wrap(replicate~.,nrow = 2) +
    ylab("Read coverage") +
    xlab("position along rDNA") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))

parrnapl_reps<-ggplot(parna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    facet_wrap(replicate~.,nrow = 2) +
    ylab("Read coverage") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))

icliprnapl_noreps<-ggplot(icliprna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    ylab("iCLIP\nRead coverage") +
    xlab("position along rDNA") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))

parrnapl_noreps<-ggplot(parna,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    ylab("PAR-CLIP\nRead coverage") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))

pdf(file = "figures/new_rDNA_replicates_tcnotc.pdf",width = 8,height = 7)
plot_grid(parrnapl_reps,icliprnapl_reps,ncol = 1,align = "v")
dev.off()

png(filename = "figures/new_rDNA_replicates_tcnotc.png",width = 8,height = 7,units = "in",res=300)
plot_grid(parrnapl_reps,icliprnapl_reps,ncol = 1,align = "v")
dev.off()

pdf(file = "figures/new_rDNA_tcnotc.pdf",width = 8,height = 4)
plot_grid(parrnapl_noreps,icliprnapl_noreps,ncol = 1,align = "v")
dev.off()

png(filename = "figures/new_rDNA_tcnotc.png",width = 8,height = 4,units = "in",res=300)
plot_grid(parrnapl_noreps,icliprnapl_noreps,ncol = 1,align = "v")
dev.off()
#RNA STRUCTURE



txs<-names(list_mean_dg)

ggs<-GTF_annotation$trann$gene_id[match(txs,GTF_annotation$trann$transcript_id)]

str_nog4<-list_mean_dg[ggs%in%unlist(list_DE_genes)]

first_100<-t(as.data.frame(sapply(str_nog4[elementNROWS(str_nog4)>100],FUN = function(x){x[1:100]})))
mid_100<-t(as.data.frame(sapply(str_nog4[elementNROWS(str_nog4)>100],FUN = function(x){
    midp<-as.integer(length(x)/2)
    x[(midp-49):(midp+50)]
})))
end_100<-t(as.data.frame(sapply(str_nog4[elementNROWS(str_nog4)>100],FUN = function(x){x[(length(x)-99):length(x)]})))


txs<-rownames(first_100)
ggs<-GTF_annotation$trann$gene_id[match(txs,GTF_annotation$trann$transcript_id)]

aggs_f<-list()

for(i in names(list_DE_genes)){
    
    aggs<-apply(first_100[ggs%in%list_DE_genes[[i]],],2,function(x){
        x<-x[!is.na(x) & !is.nan(x) & !is.infinite(x)]
        mean(x,na.rm = T)
    })
    aggs_f[[i]]<-aggs
}



aggs_m<-list()

for(i in names(list_DE_genes)){
    
    aggs<-apply(mid_100[ggs%in%list_DE_genes[[i]],],2,function(x){
        x<-x[!is.na(x) & !is.nan(x) & !is.infinite(x)]
        mean(x,na.rm = T)
    })
    aggs_m[[i]]<-aggs
}

aggs_e<-list()

for(i in names(list_DE_genes)){
    
    aggs<-apply(end_100[ggs%in%list_DE_genes[[i]],],2,function(x){
        x<-x[!is.na(x) & !is.nan(x) & !is.infinite(x)]
        mean(x,na.rm = T)
    })
    aggs_e[[i]]<-aggs
}

df<-rbind(melt(aggs_f),melt(aggs_m),melt(aggs_e))
df$position<-rep(1:100,dim(df)[1]/100)
df$section<-c(rep("first100_nt",500),rep("mid_100nt",500),rep("last100_nt",500))
df$section<-factor(df$section,levels=c("first100_nt","mid_100nt","last100_nt"))
gstr<-ggplot(df,aes(x=position,y=value,group=L1,color=L1)) +
    geom_line(size=3) +
    facet_wrap(~section) +
    ylab(expression( "average " * Delta * "G")) +
    xlab("5'UTR position (nt)") +
    scale_color_manual(values = alpha(c("black","cornflowerblue","orange","blue","dark red"),.7),"Tx_class") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) 
cairo_pdf(file = "figures/avg_deltaG_5utrs.pdf",width = 11,height = 5)
print(gstr)
dev.off()


bins_100<-t(as.data.frame(sapply(str_nog4[elementNROWS(str_nog4)>100],FUN = function(x){
    clos<-100*(round(length(x)/100,digits = 0)+1)
    idx<-as.integer(seq(1,length(x),length.out = clos))
    
    apply(matrix(x[idx],ncol = 100),2,function(x){
        x<-x[!is.na(x) & !is.nan(x) & !is.infinite(x)]
        mean(x,na.rm = T)})
    
})))

aggs_bins<-list()
for(i in names(list_DE_genes)){
    
    aggs<-apply(bins_100[ggs%in%list_DE_genes[[i]],],2,function(x){
        x<-x[!is.na(x) & !is.nan(x) & !is.infinite(x)]
        mean(x,na.rm = T)
    })
    aggs_bins[[i]]<-aggs
}

df<-melt(aggs_bins)
df$position<-rep(1:100,dim(df)[1]/100)
gstr_bins<-ggplot(df,aes(x=position,y=value,group=L1,color=L1)) +
    geom_line(size=3) +
    ylab(expression( "average " * Delta * "G")) +
    xlab("5'UTR position (bins)") +
    scale_x_continuous(breaks=c(1,20,40,60,80,100), labels=c("TSS","20","40","60","80","start codon"))+
    scale_color_manual(values = alpha(c("black","cornflowerblue","orange","blue","dark red"),.7),"Tx_class") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) 
cairo_pdf(file = "figures/avg_deltaG_5utrs_bins.pdf",width = 11,height = 5)
print(gstr_bins)
dev.off()

