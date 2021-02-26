library(utils)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(devtools)
library(Biobase)
library(preprocessCore)
library(data.table)

set.seed(2017)

register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))


path = '/data/kuw/biocore/wlku/pipeline/iscChIC_seq_test/'

K27_peaks <- getPeaks(paste0(path,'predata/K4_K27_bivalent_domain.txt'),sort_peaks = TRUE)
K27_peaks<-resize(K27_peaks, width = 5000, fix = "center")
K27_peaks <-unique(K27_peaks)
K27_peaks<-sort(K27_peaks)


sc_bedfile<-read.table(paste0(path,'predata/wc_K27_uniq.txt'));
sc_bedfile2<-paste0(path,'selbed_k27/', sc_bedfile[,2])

nx<-c("K27")
for (i in 1:9207)
{
	nx[i]<-"K27"
}
                        
K27_counts<- getCounts(sc_bedfile2, 
                             K27_peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype =nx ))
                                                          
aa1 <-as.matrix(assays(K27_counts)$counts)
ad1<-colData(K27_counts)$depth            
fwrite(aa1,file=paste0(path,'predata/scK27_4lanes_rc_at_25951_bipeaks.txt'),col.names=F,row.names=F,quote=F,sep="\t")                                        

bedfile<-read.table(paste0(path,'predata/wc_K4_uniq.txt'));
bedfile2<-paste0(path,'selbed_k4/', sc_bedfile[,2])

nx<-c("wbc")
for (i in 1:7798)
{
	nx[i]<-"wbc"
}
fragment_counts <- getCounts(bedfile2, 
                            K27_peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype =nx))
                                          
aa1 <-as.matrix(assays(fragment_counts)$counts)
fwrite(aa1,file=paste0(path,"predata/scK4_rc_at_25951_bipeaks.txt"),col.names=F,row.names=F,quote=F,sep="\t")                          


peaks <- getPeaks(paste0(path,'predata/peakname3.txt'),sort_peaks = TRUE)
peaks<-resize(peaks, width = 3000, fix = "center")
peaks <-unique(peaks)
peaks<-sort(peaks)

bedfile<-read.table(paste0(path,'predata/wc_K4_uniq.txt'));
bedfile2<-paste0(path,'selbed_k4/', sc_bedfile[,2])

nx<-c("wbc")
for (i in 1:7798)
{
	nx[i]<-"wbc"
}
fragment_counts <- getCounts(bedfile2, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype =nx))
                                          
aa1 <-as.matrix(assays(fragment_counts)$counts)
fwrite(aa1,file=paste0(path,"predata/scH3K4me3_7798_cell_counts_at_WBCpeaks.txt"),col.names=F,row.names=F,quote=F,sep="\t")                          



peaks <- getPeaks(paste0(path,'predata/H3K4me3_peakname3.txt'),sort_peaks = TRUE)
peaks<-resize(peaks, width = 3000, fix = "center")
peaks <-unique(peaks)
peaks<-sort(peaks)

bedfile<-read.table(paste0(path,'predata/wc_K4_uniq.txt'));
bedfile2<-paste0(path,'selbed_k4/', sc_bedfile[,2])

nx<-c("wbc")
for (i in 1:7798)
{
	nx[i]<-"wbc"
}
fragment_counts <- getCounts(bedfile2, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype =nx))
                                          
aa1 <-as.matrix(assays(fragment_counts)$counts)
fwrite(aa1,file=paste0(path,"predata/scH3K4me3_7798_cell_counts_at_WBCpeaks.txt"),col.names=F,row.names=F,quote=F,sep="\t")                          



K27_peaks <- getPeaks(paste0(path,'predata/combined_sicer_E100_peaks.sort.merge.txt'),sort_peaks = TRUE)
K27_peaks<-resize(K27_peaks, width = 10000, fix = "center")
K27_peaks <-unique(K27_peaks)
K27_peaks<-sort(K27_peaks)

sc_bedfile<-read.table(paste0(path,'predata/wc_K27_uniq.txt'));
sc_bedfile2<-paste0(path,'selbed_k27/', sc_bedfile[,2])

nx<-c("K27")
for (i in 1:9207)
{
	nx[i]<-"K27"
}
                        
K27_counts<- getCounts(sc_bedfile2, 
                             K27_peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype =nx ))
                                                          
aa1 <-as.matrix(assays(K27_counts)$counts)
ad1<-colData(K27_counts)$depth            
fwrite(aa1,file=paste0(path,'predata/scK27_4lanes_rc_at_79100peaks_width10k.txt'),col.names=F,row.names=F,quote=F,sep="\t")                                        


