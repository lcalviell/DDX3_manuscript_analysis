library("GenomicFeatures")
library("rtracklayer")
library(DEXSeq)
library(riborex)
library(ggplot2)
library("gplots")
library("reshape2")
library("ggrepel")
library("topGO")
library("org.Hs.eg.db")
library(gridExtra)
library("cowplot")
library("ggridges")
library(mixtools)
library(ggplot2)
library(corrgram)
library(Biostrings)


#IMPORTANT: SETWD TO DIRECTORY WITH FILES IN IT
#setwd("/bfd/lcalviel/data/riboseq/DDX3_paper_Markus/To_github/")

#LOAD DATA


load("data/gencode.v25.annotation.gtf_Rannot")
load("data/res_DE")
load("data/gene_level_tpms")
load("data/sequences_Robj")

results_DE<-res_DE_markus
DE_type="protein coding genes"

riborex<-results_DE$Riborex_DESeq_cds
ttl<-colnames(results_DE$DEXSeq_ribo_5cds3)[grep(x=colnames(results_DE$DEXSeq_ribo_5cds3),pattern="log2")]
ttl<-gsub(ttl,pattern = "log2fold_",replacement = "")
ttl<-paste(ttl,"\n",DE_type)
if(DE_type=="protein coding genes"){
    des_rna<-results_DE$DESeq_rna_cds
    des_rna<-des_rna[order(des_rna$log2FoldChange,decreasing = T),]
    des_ribo<-results_DE$DESeq_ribo_cds
    des_ribo<-des_ribo[match(rownames(des_rna),rownames(des_ribo)),]
}

if(DE_type=="all genes"){
    des_rna<-results_DE$DESeq_rna_ex
    des_rna<-des_rna[order(des_rna$log2FoldChange,decreasing = T),]
    des_ribo<-results_DE$DESeq_ribo_ex
    des_ribo<-des_ribo[match(rownames(des_rna),rownames(des_ribo)),]
}

#PREPARE DE_PLOT

des_all<-des_rna
des_all$pvalue<-NULL
des_all$baseMean_rna<-des_rna$baseMean
des_all$baseMean_ribo<-des_ribo$baseMean
des_all$log2FoldChange_rna<-des_rna$log2FoldChange
des_all$padj_rna<-des_rna$padj
des_all$log2FoldChange_ribo<-des_ribo$log2FoldChange
des_all$padj_ribo<-des_ribo$padj
des_all$log2FoldChange_riborna<-riborex$log2FoldChange[match(rownames(des_all),rownames(riborex))]
des_all$padj_riborna<-riborex$padj[match(rownames(des_all),rownames(riborex))]

des_ribo_5cds3<-results_DE$DEXSeq_ribo_5cds3[,c(colnames(results_DE$DEXSeq_ribo_5cds3)[grep(x=colnames(results_DE$DEXSeq_ribo_5cds3),pattern="log2")],"pvalue","groupID","featureID","exonBaseMean")]
colnames(des_ribo_5cds3)<-c("log2FoldChange","pval_5ut","groupID","featureID","fiveUTR_BaseMean")
des_ribo_5cds3<-des_ribo_5cds3[des_ribo_5cds3$featureID=="E001",]

des_all$log2FoldChange_5ut<-des_ribo_5cds3$log2FoldChange[match(rownames(des_all),des_ribo_5cds3$groupID)]
des_all$pval_5ut<-des_ribo_5cds3$pval_5ut[match(rownames(des_all),des_ribo_5cds3$groupID)]
des_all$baseMean_5ut<-des_ribo_5cds3$fiveUTR_BaseMean[match(rownames(des_all),des_ribo_5cds3$groupID)]
des_all$chrM<-ifelse(test = seqnames(GTF_annotation$genes[rownames(des_all)])=="chrM",yes = T,no = F)
des_all$gene_id<-rownames(des_all)
des_all<-des_all[,c("gene_id","baseMean_rna","baseMean_ribo","baseMean_5ut",colnames(des_all)[grep(colnames(des_all),pattern = "log2FoldChange_|padj_|pval|exonBaseM")])]
ggmnn<-unique(GTF_annotation$trann[,c("gene_id","gene_name")])
des_all$symb<-GTF_annotation$trann$gene_name[match(rownames(des_all),GTF_annotation$trann$gene_id)]
des_all$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(des_all),GTF_annotation$trann$gene_id)]

oklow<-which(des_all$baseMean_rna>0 | des_all$baseMean_ribo>0)
des_all<-des_all[oklow,]

#CALCULATE 5'UTR OCCUPANCY CHANGE

cnts<-res_DE_markus$DEXSeq_ribo_5cds3$countData
cnts<-apply(cnts,2,FUN = function(x){x/(sum(x)/1e06)})

cnts_5<-cnts[grep(rownames(cnts),pattern = "E001"),]
cnts_cds<-cnts[grep(rownames(cnts),pattern = "E002"),]
rownames(cnts_5)<-sapply(strsplit(rownames(cnts_5),":"),"[[",1)
rownames(cnts_cds)<-sapply(strsplit(rownames(cnts_cds),":"),"[[",1)
cnts_all<-cbind(cnts_5,cnts_cds)
ratios<-t(apply(cnts_all,1,function(x){
    fv_si<-mean(x[1],x[2])
    fv_wt<-mean(x[3],x[4])
    cds_si<-mean(x[5],x[6])
    cds_wt<-mean(x[7],x[8])
    rts_si<-fv_si/cds_si
    rts_wt<-fv_wt/cds_wt
    c(rts_si,rts_wt)
}))
ratios<-ratios[complete.cases(ratios) & ratios[,1]!=0 & ratios[,2]!=0 & ratios[,1]!=Inf & ratios[,2]!=Inf,]
gnms<-GTF_annotation$trann$gene_name[match(rownames(ratios),GTF_annotation$trann$gene_id)]
log2FC_5utrcds<-log2(ratios[,1]/ratios[,2])
des_all$log2FC_5utrcds<-log2FC_5utrcds[des_all$gene_id]


write.table(des_all,file = "data/DE_table.csv",row.names = FALSE,sep="\t",quote = FALSE)

ok1<-which(des_all$baseMean_rna>20 & des_all$baseMean_ribo>20)
des_all<-des_all[ok1,]

des_all$tx_type<-"Not_significant"
des_all$tx_type[des_all$padj_rna<.05 & des_all$log2FoldChange_rna<0]<-"RNA_down"
des_all$tx_type[des_all$padj_rna<.05 & des_all$log2FoldChange_rna>0]<-"RNA_up"

des_all$tx_type[des_all$padj_riborna<.05 & des_all$log2FoldChange_riborna<0]<-"TE_down"
des_all$tx_type[des_all$padj_riborna<.05 & des_all$log2FoldChange_riborna>0]<-"TE_up"
euc_dist<-des_all$log2FoldChange_rna-des_all$log2FoldChange_ribo

#CORRECT FOR MIXED DE GENES

eucrna<-euc_dist[des_all$tx_type=="RNA_up"]

mixt<-normalmixEM(eucrna,mu = c(0,.4),k = 2,mean.constr = c(0,.4))

clstss<-apply(mixt$posterior,1,which.max)
names(clstss)<-des_all$gene_id[des_all$tx_type=="RNA_up"]
des_all$clustr<-3
des_all$clustr[des_all$tx_type=="RNA_up"]<-clstss
des_all$tx_type2<-des_all$tx_type
des_all$tx_type2[des_all$clustr==1]<-"RNA_up_clean"
des_all$tx_type2[des_all$clustr==2]<-"Mixed"
ddf<-rbind(cbind(des_all$log2FoldChange_rna-des_all$log2FoldChange_ribo,des_all$tx_type),cbind(des_all$log2FoldChange_rna[des_all$clustr!=3]-des_all$log2FoldChange_ribo[des_all$clustr!=3],des_all$tx_type2[des_all$clustr!=3]))
ddf<-as.data.frame(ddf)
colnames(ddf)<-c("RNA_Ribo","tx_class")
ddf$RNA_Ribo<-as.numeric(as.character(ddf$RNA_Ribo))
ddf$tx_class<-as.character(ddf$tx_class)
ddf$tx_class[ddf$tx_class=="RNA_up"]<-"RNA_up_original"
ddf$tx_class<-factor(ddf$tx_class,levels=rev(c("Not_significant","RNA_down","TE_up","TE_down","RNA_up_original","RNA_up_clean","Mixed")))
colli<-alpha(c("dark grey","cornflowerblue","dark red","blue","black","orange","dark blue"),.7)
correct_rna<-ggplot(ddf,aes(RNA_Ribo,y=tx_class,fill=tx_class))+
    stat_density_ridges(bandwidth = .04) +
    xlab("RNA - Ribo\nFold Change") +
    ylab("Density") +
    scale_fill_manual(values = rev(colli),"Tx_class") + 
    scale_color_manual(values = rev(colli),"Tx_class") +
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=18))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) 
des_all$tx_type2[des_all$clustr==1]<-"RNA_up"

pdf(file = "figures/correct_DE.pdf",width = 7,height = 4)
print(correct_rna + theme(legend.position="none"))

dev.off()

write.table(des_all,file = "data/DE_table_ok.csv",row.names = FALSE,sep="\t",quote = FALSE)



ok1<-which(des_all$baseMean_rna>20 & des_all$baseMean_ribo>20)
des_all<-des_all[ok1,]

dfa=des_all

dfa<-as.data.frame(dfa)

dfa<-dfa[grep(dfa$biotype,pattern = "pseudo|snRNA|snoRNA",invert = T),]
dfa$biotype[dfa$biotype!="protein_coding"]<-"non_coding"


dfa$tx_type<-"Not_significant"
dfa$tx_type[dfa$padj_rna<.05 & dfa$log2FoldChange_rna<0]<-"RNA_down"
dfa$tx_type[dfa$padj_rna<.05 & dfa$log2FoldChange_rna>0]<-"RNA_up"

dfa$tx_type[dfa$padj_riborna<.05 & dfa$log2FoldChange_riborna<0]<-"TE_down"
dfa$tx_type[dfa$padj_riborna<.05 & dfa$log2FoldChange_riborna>0]<-"TE_up"

dfa2<-dfa
dfa2<-dfa2[dfa2$tx_type2!="Mixed",]
dfa2$Tx_class<-factor(dfa2$tx_type2,levels=rev(c("Not_significant","RNA_up","RNA_down","TE_up","TE_down")))
grd<-ggplot(dfa2,aes(log2FC_5utrcds,group=Tx_class,y=Tx_class,fill=Tx_class))  +theme_ridges()
grd<-grd + scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"))
grd<-grd +xlim(-4,4) + xlab("log2FC_5'UTR_skew")       
grd<-grd + stat_density_ridges(quantile_lines = TRUE, quantiles = 2) + geom_vline(xintercept = 0,lty=2,color="black")
grd<-grd + theme(legend.position="none")
pdf(file = "figures/DE_5riboskew.pdf",width = 7,height = 4)
print(grd)
dev.off()


dfa$log2FoldChange_riborna[dfa$padj_riborna>.05]<-NA
dfa$log10_pval_riborex<--log10(dfa$padj_riborna)

dfa$log10_pval_riborex[is.na(dfa$log10_pval_riborex)]<-0

dfa$tx_type<-factor(dfa$tx_type,levels=c("TE_up","TE_down","RNA_up","RNA_down","Not_significant"))
dfa$biotype[grep(dfa$symb,pattern = "HIST|H2B")]<-"Histone genes"
dfa$biotype[grep(dfa$symb,pattern = "HIST|H2B",invert =  T)]<-"all genes"
dfa$biotype<-factor(dfa$biotype,levels=c("all genes","Histone genes"))
dfa$symb2<-as.character(dfa$symb)
dfa$symb2[grep(dfa$symb2,pattern = "DDX3|DHX36|DHX9|CHGB|CHGA|SCG3|EGR1|FOS|CYP4F11|ANKRD9|COL2A1|COL4A1|C19orf48|MAN1C1|ATF4|NT5DC2",invert = T)]<-NA
dfa$symb2[dfa$tx_type=="Not_significant"]<-NA

dfa$symb3<-NA
toprna<-dfa[dfa$tx_type=="RNA_up",]
toprna<-toprna[toprna$biotype!="Histone genes",]
toprna<-toprna[order(toprna$log2FoldChange_rna,toprna$log2FoldChange_ribo,decreasing = T),"symb"]
if(length(toprna)>6){toprna=toprna[1:5]}
topribo<-dfa[dfa$tx_type=="TE_up",]
topribo<-topribo[topribo$biotype!="Histone genes",]
topribo<-topribo[order(topribo$log2FoldChange_riborna,decreasing = T),"symb"]
if(length(topribo)>6){topribo=topribo[1:4]}
bottomrna<-dfa[dfa$tx_type=="RNA_down",]
bottomrna<-bottomrna[bottomrna$biotype!="Histone genes",]
bottomrna<-bottomrna[order(bottomrna$log2FoldChange_rna,bottomrna$log2FoldChange_ribo,decreasing = F),"symb"]
if(length(bottomrna)>6){bottomrna=bottomrna[1:4]}
bottomribo<-dfa[dfa$tx_type=="TE_down",]
bottomribo<-bottomribo[bottomribo$biotype!="Histone genes",]
bottomribo<-bottomribo[order(bottomribo$log2FoldChange_riborna,decreasing = F),"symb"]
if(length(bottomribo)>10){bottomribo=c(bottomribo[1:4],bottomribo[bottomribo%in%c("ODC1","DVL1","DVL2","HMBS")])}
dfa$symb3[dfa$symb%in%c(toprna,topribo,bottomrna,bottomribo)]<-dfa$symb[dfa$symb%in%c(toprna,topribo,bottomrna,bottomribo)]

dfa$log2FoldChange_5ut[which(dfa$pval_5ut>.05)]<-NA
dfa$log2FoldChange_5ut[which(dfa$log2FoldChange_5ut>1.5)]<-1.5
dfa$log2FoldChange_5ut[which(dfa$log2FoldChange_5ut<(-1.5))]<--1.5


dfa$symb4<-NA
dfa$log10_pval_5utr<--log10(dfa$pval_5ut)
dfa$log10_pval_5utr[is.na(dfa$log10_pval_5utr) | dfa$pval_5ut>.05 ]<-0
top5ut<-dfa[which(dfa$log2FoldChange_5ut>0),]
top5ut<-top5ut[order(top5ut$log10_pval_5utr,decreasing = T),"symb"]

if(length(top5ut)>6){top5ut=top5ut[1:6]}
bottom5ut<-dfa[which(dfa$log2FoldChange_5ut<0),]
bottom5ut<-bottom5ut[order(bottom5ut$log10_pval_5utr,decreasing = T),"symb"]
if(length(bottom5ut)>6){bottom5ut=bottom5ut[1:6]}
dfa$symb4[dfa$symb%in%c(top5ut,bottom5ut)]<-dfa$symb[dfa$symb%in%c(top5ut,bottom5ut)]
dfa$nudgey<-NA
dfa$nudgey[which(nchar(dfa$symb3)>2)]<-ifelse(dfa$log2FoldChange_ribo[which(nchar(dfa$symb3)>2)]>0, .4, -.5)
dfa$nudgex<-NA
dfa$nudgex[which(nchar(dfa$symb3)>2)]<-ifelse(dfa$log2FoldChange_rna[which(nchar(dfa$symb3)>2)]>0, .3, -.4)
dfa$nudgex[dfa$tx_type=="TE_up"]<-0
dfa$nudgey[dfa$tx_type=="RNA_up"]<-0
dfa$nudgey[dfa$tx_type=="RNA_down"]<-0
dfa$nudgex[dfa$tx_type=="RNA_up"]<-0
dfa$nudgex[dfa$tx_type=="RNA_down"]<-0
dfa$symb3[which(dfa$symb3=="C6orf1")]<-NA
dfa$nudgex[which(dfa$symb3=="HMBS")]<-.7
dfa$nudgey[which(dfa$symb3=="DVL2")]<-dfa$nudgey[which(dfa$symb3=="DVL2")]-.1

dfa$nudgey[which(dfa$symb3=="TCFL5")]<-dfa$nudgey[which(dfa$symb3=="TCFL5")]+.1
dfa$nudgex[which(dfa$symb3=="TCFL5")]<-dfa$nudgex[which(dfa$symb3=="TCFL5")]-.3

a<-ggplot(dfa,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=log10_pval_riborex,shape=biotype,label=symb3)) + geom_point()
a<-a + theme_bw() +
    ylab("Ribo change") +
    xlab("RNA change") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    scale_color_manual(values = alpha(c("dark red","blue","firebrick1","cornflowerblue","gray24"),c(.8,.8,.8,.8,.5)),"Tx_class") +
    scale_size_continuous(name = "-log10 adj.pval\nRiborex")

pdf(file = "figures/DE_plot_nolab.pdf",width = 9,height = 5)
grid.arrange(a + ylim(-3,3)+ xlim(-3,3))
dev.off()

png(filename = "figures/DE_plot_nolab.png",width = 9,height = 5,units = "in",res=300)
grid.arrange(a + ylim(-3,3)+ xlim(-3,3))
dev.off()

a<-a + geom_text_repel(size=5,nudge_x = dfa$nudgex[which(nchar(dfa$symb3)>2)],
                       nudge_y = dfa$nudgey[which(nchar(dfa$symb3)>2)],force=2)
pdf(file = "figures/DE_plot.pdf",width = 9,height = 5)
grid.arrange(a + ylim(-3,3)+ xlim(-3,3))
dev.off()


png(filename = "figures/DE_plot.png",width = 9,height = 5,units = "in",res=300)
grid.arrange(a + ylim(-3,3)+ xlim(-3,3))
dev.off()

#CALCULATE SEQUENCE AND TX FEATS

txok<-GTF_annotation$cds_txs_coords
txok<-txok[txok$reprentative_boundaries]
txok2<-txok
txok2$tx_type<-NA
spl_des_all<-split(des_all$gene_id,f = des_all$tx_type2)
spl_des_all<-spl_des_all[-which(names(spl_des_all)=="Mixed")]
for(i in names(spl_des_all)[1:5]){
    txok2$tx_type[txok2$gene_id%in%spl_des_all[[i]]]<-i
}

seqqe<-narrow(seqqqs,start = 1,end = start(txok2))
txok2$GC_content<-letterFrequency(seqqe,letters = "GC",as.prob = T)*100
txok2<-txok2[txok2$gene_id%in%unlist(spl_des_all)]
txok2$len_5ut<-start(txok2)-1
dfs<-as.data.frame(mcols(txok2))
dfs$tx_type<-factor(dfs$tx_type,levels=rev(c("Not_significant","RNA_up","RNA_down","TE_up","TE_down")))

lng<-ggplot(dfs,aes(len_5ut+1,group=tx_type,y=tx_type,fill=tx_type))  +theme_ridges()
lng<-lng + scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"))
lng<-lng + xlab("5'UTR length")       
lng<-lng + stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) 
lng<-lng + theme(legend.position="none")
lng<-lng + scale_x_log10(breaks=c(1,11,101,1001,10001),labels=c(0,10,100,1000,10000))

dfs<-dfs[dfs$len_5ut>0,]
gc5<-ggplot(dfs,aes(G.C,group=tx_type,y=tx_type,fill=tx_type))  +theme_ridges()
gc5<-gc5 + scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"))
gc5<-gc5 + xlab("%GC content in 5'UTR")       
gc5<-gc5 + stat_density_ridges(quantile_lines = TRUE, quantiles = 2) 
gc5<-gc5 + theme(legend.position="none")
pdf(file = "figures/DE_5utr_feats.pdf",width = 8,height = 5)
plot_grid(lng,gc5,align = "v",ncol = 1)
dev.off()


txok2<-txok

txok2$tx_type<-NA
spl_des_all<-split(des_all$gene_id,f = des_all$tx_type2)
spl_des_all<-spl_des_all[-which(names(spl_des_all)=="Mixed")]
for(i in names(spl_des_all)[1:5]){
    txok2$tx_type[txok2$gene_id%in%spl_des_all[[i]]]<-i
}
txok2<-txok2[!is.na(txok2$tx_type)]

seqqqs<-seqqqs[names(seqqqs)%in%as.character(seqnames(txok2))]

seqq1<-narrow(seqqqs,start = 1,end = start(txok2))
utr5_len<-nchar(seqq1)
utr5_GC<-letterFrequency(seqq1,letters = "GC",as.prob = T)*100

seqq2<-narrow(seqqqs,start = start(txok2),end = end(txok2))
cds_len<-nchar(seqq2)
cds_GC<-letterFrequency(seqq2,letters = "GC",as.prob = T)*100

seqq3<-narrow(seqqqs,start = end(txok2),end = txok2$lentx)
utr3_len<-nchar(seqq3)
utr3_GC<-letterFrequency(seqq3,letters = "GC",as.prob = T)*100

df1<-data.frame(cbind(utr5_len,utr5_GC),stringsAsFactors = F)
colnames(df1)<-c("len","GC_content")
df1$region="5_UTR"


df2<-data.frame(cbind(cds_len,cds_GC),stringsAsFactors = F)
colnames(df2)<-c("len","GC_content")
df2$region="CDS"


df3<-data.frame(cbind(utr3_len,utr3_GC),stringsAsFactors = F)
colnames(df3)<-c("len","GC_content")
df3$region="3_UTR"

dfall<-rbind.data.frame(df1,df2,df3)
dfall$Tx_class<-rep(txok2$tx_type,3)

dfall$Tx_class<-factor(dfall$Tx_class,levels=rev(c("Not_significant","RNA_up","RNA_down","TE_up","TE_down")))
dfall$region<-factor(dfall$region,levels=c("5_UTR","CDS","3_UTR"))


lng_all<-ggplot(dfall,aes(len+1,group=Tx_class,color=Tx_class))  +theme_ridges()
lng_all<-lng_all + scale_color_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"))
lng_all<-lng_all + xlab("length (nt)")       
lng_all<-lng_all + stat_density(geom = "line", position = "identity",from = 0) 
lng_all<-lng_all + theme(legend.position="none")
lng_all<-lng_all + scale_x_log10(breaks=c(1,11,101,1001,10001),labels=c(0,10,100,1000,10000))
lng_all<-lng_all + facet_wrap(~region,ncol = 1)



dfs<-dfall[dfall$len>0,]

gc_all<-ggplot(dfs,aes(GC_content,group=Tx_class,color=Tx_class))  +theme_ridges()
gc_all<-gc_all + scale_color_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"))
gc_all<-gc_all + xlab("%GC content")       
gc_all<-gc_all + stat_density(geom = "line", position = "identity",from = 0) 
gc_all<-gc_all + facet_wrap(~region,ncol = 1)
pdf(file = "figures/DE_other_feats.pdf",width = 8,height = 5)
plot_grid(lng_all,gc_all,align = "h",ncol = 2,rel_widths = c(1,1.75))
dev.off()

#GO ANALYSIS

list_DE_genes<-split(des_all$gene_id,f = des_all$tx_type2)
list_DE_genes<-list_DE_genes[-which(names(list_DE_genes)=="Mixed")]

list_DE_gene_symb<-lapply(list_DE_genes,function(x){GTF_annotation$trann$gene_name[match(x,GTF_annotation$trann$gene_id)]})

results_go<-list()
for(i in names(list_DE_gene_symb)){
    
    allG <- unlist(list_DE_gene_symb)
    myG <- list_DE_gene_symb[[i]]
    mapp <- factor(as.integer(allG %in% myG))
    names(mapp) <- allG
    res_go<-list()
    for(ont in c("BP","MF","CC")){
        suppressMessages(GOdata <- new("topGOdata", ontology = ont, # Here one can choose BP/MF/CC
                                       allGenes = mapp,
                                       annot = annFUN.org, 
                                       mapping="org.Hs.eg.db", 
                                       ID="symbol"))
        
       
        suppressMessages(resultFisher <- runTest(GOdata, statistic = "fisher"))
        tmp <- DataFrame(GenTable(GOdata, classicFisher = resultFisher,topNodes = length(score(resultFisher))))
        tmp$classicFisher[tmp$classicFisher=="< 1e-30"]<-min(as.numeric(tmp$classicFisher))
        tmp$gene<-CharacterList(genesInTerm(GOdata,whichGO = tmp$GO.ID))
        res_go[[ont]]<-tmp
    }
    results_go[[i]]<-res_go
}
save(results_go,file = "data/res_go_analysis_ok05_Robj")




tpms_symb<-tpms
rownames(tpms_symb)<-GTF_annotation$trann$gene_name[match(rownames(tpms_symb),GTF_annotation$trann$gene_id)]


allgplots<-list()
allgplots_sel<-list()
for(i in names(results_go)){
    if(i=="Not_significant"){next}
    x<-results_go[[i]]
    list_allgg<-list()
    list_allgg_sel<-list()
    list_ontgo<-list()
    xontlist<-list()
    for(j in names(x)){
        x_ont<-x[[j]]
        
        x_ont$classicFisher<-as.numeric(x_ont$classicFisher)
        x_ont<-x_ont[x_ont$classicFisher<.01 & x_ont$Annotated>2 & x_ont$Annotated<3000,]
        x_ont$Termall<-x_ont$Term
        x_ont$Term<-as.character(x_ont$Term)
        tropp<-which(nchar(x_ont$Term)>30)
        if(length(tropp)>0){
            for(troppp in tropp){
                x_ont$Term[troppp]<-paste(as.character(BString(x_ont$Term[troppp])[1:28]),"...",sep="")
            }
        }
        
        names(x_ont$gene)<-x_ont$Term
        x_ont$gene_diff<-CharacterList(lapply(x_ont$gene,function(x){
            x<-list_DE_gene_symb[[i]][sort(match(x,list_DE_gene_symb[[i]]))]
            x
        }))
        x_ont$TPM_rna_ctrl<-NumericList(lapply(x_ont$gene_diff,function(x){
            as.numeric((tpms_symb[,1][x]+tpms_symb[,2][x])/2)
        }))
        x_ont$TPM_rna_siDDX3<-NumericList(lapply(x_ont$gene_diff,function(x){
            as.numeric((tpms_symb[,3][x]+tpms_symb[,4][x])/2)
        }))
        x_ont$TPM_ribo_ctrl<-NumericList(lapply(x_ont$gene_diff,function(x){
            as.numeric((tpms_symb[,5][x]+tpms_symb[,6][x])/2)
        }))
        x_ont$TPM_ribo_siDDX3<-NumericList(lapply(x_ont$gene_diff,function(x){
            as.numeric((tpms_symb[,7][x]+tpms_symb[,8][x])/2)
        }))
        
        x_ont$gene_diff<-endoapply(x_ont$gene_diff,function(x){
            if(length(x)>2){x<-x[1:2]}
            return(x)
        })
        x_ont<-x_ont[elementNROWS(x_ont$gene_diff)>1,]
        if(dim(x_ont)[1]>10){x_ont<-x_ont[1:10,]}
        
        x_ont$genes<-as.character(paste(x_ont$gene_diff,collapse = "; "))
        x_ont<-as.data.frame(x_ont)
        x_ont$logpadj<--log10(as.numeric(x_ont$classicFisher))
        
        
        x_ont$Term<-factor(x_ont$Term,levels = rev(x_ont$Term))
        x_ont$GO.ID<-factor(x_ont$GO.ID,levels = rev(x_ont$GO.ID))
        xontlist[[j]]<-x_ont
        labba<-paste(x_ont$Term,": (",x_ont$genes,")",sep = "")
        labball<-paste(x_ont$GO.ID,": ",x_ont$Term," (",x_ont$genes,") - p=",x_ont$classicFisher,sep="")
        
        ggont<-ggplot(x_ont, aes(y = logpadj, x=GO.ID, label=genes)) + geom_bar(stat="identity")  + 
            theme(axis.text.x=element_text(angle=-40, hjust=0)) + coord_flip() + theme_classic() + xlab("") + ylab(label = "-log10 adjusted p-value")
        ggont <- ggont + theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
            theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=0, vjust=0.5, size=18))  
        ggont<-ggont + ggtitle(i,subtitle = j)  + geom_hline(yintercept=-log10(.01),col="red") + theme(plot.title = element_text(size = 30, face = "bold")) + theme(plot.subtitle = element_text(size = 20, face = "italic"))
        list_ontgo[[j]]<-ggont
        
        df1<-data.frame(value=unlist(x_ont$TPM_rna_ctrl),stringsAsFactors = F)
        rownames(df1)<-NULL
        df1$tpm_type<-"RNA_ctrl"
        gogo<-unlist(apply(cbind.data.frame(names(elementNROWS(x_ont$TPM_rna_ctrl)),as.numeric(elementNROWS(x_ont$TPM_rna_ctrl))),1,function(x){rep(x[1],as.numeric(x[2]))}))
        names(gogo)<-NULL
        df1$go<-gogo
        
        df2<-data.frame(value=unlist(x_ont$TPM_rna_siDDX3),stringsAsFactors = F)
        rownames(df2)<-NULL
        df2$tpm_type<-"RNA_siDDX3"
        gogo<-unlist(apply(cbind.data.frame(names(elementNROWS(x_ont$TPM_rna_siDDX3)),as.numeric(elementNROWS(x_ont$TPM_rna_siDDX3))),1,function(x){rep(x[1],as.numeric(x[2]))}))
        names(gogo)<-NULL
        df2$go<-gogo
        
        df3<-data.frame(value=unlist(x_ont$TPM_ribo_ctrl),stringsAsFactors = F)
        rownames(df3)<-NULL
        df3$tpm_type<-"Ribo_ctrl"
        gogo<-unlist(apply(cbind.data.frame(names(elementNROWS(x_ont$TPM_ribo_ctrl)),as.numeric(elementNROWS(x_ont$TPM_ribo_ctrl))),1,function(x){rep(x[1],as.numeric(x[2]))}))
        names(gogo)<-NULL
        df3$go<-gogo
        
        df4<-data.frame(value=unlist(x_ont$TPM_ribo_siDDX3),stringsAsFactors = F)
        rownames(df4)<-NULL
        df4$tpm_type<-"Ribo_siDDX3"
        gogo<-unlist(apply(cbind.data.frame(names(elementNROWS(x_ont$TPM_ribo_siDDX3)),as.numeric(elementNROWS(x_ont$TPM_ribo_siDDX3))),1,function(x){rep(x[1],as.numeric(x[2]))}))
        names(gogo)<-NULL
        df4$go<-gogo
        
        
        dfall<-rbind.data.frame(df1,df2,df3,df4)
        dfall$go<-factor(dfall$go,levels = rev(unique(dfall$go)))
        dfall$tpm_type<-factor(dfall$tpm_type,levels = rev(unique(dfall$tpm_type)))
        dfall$experiment_type<-"RNA-seq"
        dfall$experiment_type[dfall$tpm_type%in%c("Ribo_siDDX3","Ribo_ctrl")]<-"Ribo-seq"
        dfall$condition<-"control"
        dfall$condition[dfall$tpm_type%in%c("RNA_siDDX3","Ribo_siDDX3")]<-"siDDX3"
        dfall$condition<-factor(dfall$condition,levels=c("control","siDDX3"))
        dfall$experiment_type<-factor(dfall$experiment_type,levels=c("RNA-seq","Ribo-seq"))
        
        collego<-c("darkorchid4","mediumorchid1","dark green","olivedrab3")
        
        gotpm1<-ggplot(dfall,aes(y=value+1,x=go,color=tpm_type)) + 
            #    geom_point(pch = 19, position = position_jitterdodge()) +
            stat_summary(position = position_dodge(width = .9)) +
            scale_y_log10() + coord_flip() +
            scale_color_manual(values = collego,"experiment") + 
            xlab("") +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) +
            ylab("TPM") + guides(colour = guide_legend(reverse=T)) +
            theme(axis.title.x = element_text(size=30),axis.text.x = element_text(size=25)) 
        
        labba1<-sapply(strsplit(labba,": "),"[[",1)
        labba2<-sapply(strsplit(labba,": "),"[[",2)
        glabba<-ggplot(x_ont, aes(y = 0, x=GO.ID, label=labba1)) +
            geom_text(hjust = 0) + 
            coord_flip() +
            theme_void()+ 
            ylim(0, .5) +
            geom_text(aes(y=.5,label=labba2),hjust = 1) 
        
        
        labball1<-paste(x_ont$GO.ID,": ",x_ont$Term," (",x_ont$genes,")",sep="")
        labball2<-paste("p=",x_ont$classicFisher,sep="")
        sizs<-6
        glabba2<-ggplot(x_ont, aes(y = 0, x=GO.ID, label=labball1)) +
            geom_text(hjust = 0,size=sizs) + 
            coord_flip() +
            theme_void()+ 
            ylim(0, .5) +
            geom_text(aes(y=.5,label=labball2),hjust = 1,size=sizs) + theme(legend.position="none")+
            ggtitle(i,subtitle = j)  + theme(plot.title = element_text(size = 30, face = "bold")) + theme(plot.subtitle = element_text(size = 20, face = "italic"))
        
        
        allgos<-plot_grid(ggont,glabba,gotpm1,align = "h",nrow = 1,rel_widths = c(3,3,3))
        allgos2<-plot_grid(glabba2,gotpm1,align = "h",nrow = 1,rel_widths = c(3,1.5))
        sels<-x_ont[1:4,]
        list_allgg[[j]]<-allgos2
        
        allgos_sel<-plot_grid(ggont+ xlim(rev(as.character(sels$GO.ID))),glabba+ xlim(rev(as.character(sels$GO.ID))),gotpm1 + xlim(rev(as.character(sels$Term)))+ scale_y_log10(limits = c(1,5000)),align = "h",nrow = 1,rel_widths = c(3,3,3))
        allgos_sel2<-plot_grid(glabba2+ xlim(rev(as.character(sels$GO.ID))),gotpm1 + xlim(rev(as.character(sels$Term)))+ scale_y_log10(limits = c(1,5000)),align = "h",nrow = 1,rel_widths = c(3,1.5))
        list_allgg_sel[[j]]<-allgos_sel2
        
        
    }
    allgplots[[i]]<-list_allgg
    allgplots_sel[[i]]<-list_allgg_sel
}





pdf(file = "figures/res_go_analysis.pdf",width =30,height = 6)

plot_grid(allgplots_sel[[1]][[1]],allgplots_sel[[2]][[1]],allgplots_sel[[3]][[1]],allgplots_sel[[4]][[1]],ncol=2,align = "v")

dev.off()



pdf(file = "figures/res_go_analysis_all.pdf",width = 30,height = 35)
plot_grid(allgplots[[1]][[1]],allgplots[[2]][[1]],allgplots[[3]][[1]],allgplots[[4]][[1]],allgplots[[1]][[2]],allgplots[[2]][[2]],
          allgplots[[3]][[2]],allgplots[[4]][[2]],allgplots[[1]][[3]],allgplots[[2]][[3]],allgplots[[3]][[3]],allgplots[[4]][[3]],ncol=2,align = "v")

dev.off()




des_all<-read.table("data/DE_table_ok.csv",sep="\t",header = T,stringsAsFactors = F)

temo2<-ttheme_default(base_family = "Helvetica")
matta<-as.matrix(table(des_all$tx_type))
matta<-cbind(matta,round(matta/sum(matta)*100,digits = 2))
mode(matta)<-"character"
matta[3,1]<-paste(matta[3,1]," (",sum(des_all$tx_type2=="RNA_up"),")",sep = "")
matta[3,2]<-paste(matta[3,2]," (",round(sum(des_all$tx_type2=="RNA_up")/length(des_all$tx_type)*100,digits = 2),")",sep = "")
matta<-t(matta)
rownames(matta)<-c("number of genes","% of total")
temo2<-ttheme_default(base_family = "Helvetica")
taboo<-tableGrob(matta,theme=temo2)
pdf(file = "figures/table_DE_genes.pdf",width = 7,height = 2)
grid.arrange(taboo)
dev.off()


tabla<-as.data.frame(matta)
tabla$category<-gsub(rownames(tabla),pattern = "\n",replacement = " ")
tabla<-tabla[,c("category",colnames(tabla))]
rownames(tabla)<-NULL
write.table(tabla,file = "data/table_DE_genes.csv",sep="\t",quote = FALSE,row.names = FALSE)


tpms_ok<-tpms[rowSums(tpms)>0,]
tpms_ok<-log(tpms_ok+1)

dd<-which(rownames(tpms_ok)==GTF_annotation$trann$gene_id[which(GTF_annotation$trann$gene_name=="DDX3X")][1])
plotta<-function(x,y,corr = NULL,col.regions, cor.method, ...){
    if (!is.null(corr)) 
        return()
    
    plot.xy(xy.coords(x, y), type = "p" ,pch=19,cex=.3, ...);points(x[dd],y[dd],col="red",pch=19, ...)
    box(col = "lightgray")
    
}

panel.cor_mod<-function (x, y, corr = NULL, col.regions, cor.method, digits = 2, 
                         cex.cor, ...) 
{
    if (is.null(corr)) {
        if (sum(complete.cases(x, y)) < 2) {
            warning("Need at least 2 complete cases for cor()")
            return()
        }
        else {
            corr <- cor(x, y, use = "pair", method = cor.method)
        }
    }
    auto <- missing(cex.cor)
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    ncol <- 500
    pal <- col.regions(ncol)
    col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, 
                                                 length.out = ncol + 1), include.lowest = TRUE))
    abscorr <- formatC(abs(corr), digits = digits, format = "f")
    corr <- formatC(corr, digits = digits, format = "f")
    if (auto) 
        cex.cor <- 0.7/strwidth(abscorr)
    text(0.5, 0.5, corr, cex = cex.cor, col = pal[col.ind])
}


colnames(tpms_ok)<-gsub(colnames(tpms_ok),pattern = "_rep",replacement = "\nrep")
png(filename = "figures/replicates_plot.png",width = 7,height = 6,units = "in",res=300)
corrgram(tpms_ok, order=FALSE,
         upper.panel=panel.cor,lower.panel = plotta,digits=3)
dev.off()

pdf(file = "figures/replicates_plot.pdf",width = 7,height = 6)
corrgram(tpms_ok, order=FALSE,
         upper.panel=panel.cor,lower.panel = plotta,digits=3)
dev.off()

png(filename = "figures/replicates_plot_sp.png",width = 7,height = 6,units = "in",res=300)
corrgram(tpms_ok, order=FALSE,
         upper.panel=panel.cor,lower.panel = plotta,digits=3,cor.method = "spearman")
dev.off()

pdf(file = "figures/replicates_plot_sp.pdf",width = 7,height = 6)
corrgram(tpms_ok, order=FALSE,
         upper.panel=panel.cor,lower.panel = plotta,digits=3,cor.method = "spearman")
dev.off()