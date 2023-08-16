source("Blood_PBMC2.R")
libr <- c("limma", "edgeR", "RColorBrewer", "ggplot2")
invisible(lapply(libr, library, character.only = TRUE))


dir_graphs <- "/home/beckum/Dropbox/0MSc_Bioinfo_Uppsala/sars/3p_Hawman_RockyMountainLabs_sars_cov2_Zika_RM/3p_Hawman_RockyMountainLabs_sars_cov2_Zika_RM_images_graphs_presentations"
list_DEgenes_whole_series <- c('DEgenesCov', 'DEgenesCov_B', 
                               'DEgenesCov_noDPC1', "DEgenesCov_noDPC1_B", 
                               'DEgenesZika', 'DEgenesZika_B',
                               "DEgenesZika_noDPC1", "DEgenesZika_noDPC1_B")    

for (i in length(list_DEgenes_whole_series)) {
    png(paste0(dir_graphs, list_DEgenes_whole_series[i],".png"), width=2160, height=1440, units = "px", bg = "white")
    par(oma=c(1,1,1,1))
    lwid=c(5,8)
    f.heatmap(cbind(0,list_DEgenes_whole_series[i]), keysize=0.4, key.par = list(cex=1.25), key.title=NA, key.ylab=NA,
              labCol=as.factor(c('zero',gsub("DE_", "", colnames(list_DEgenes_whole_series[i])))), # rename columns
              srtCol=0, cexCol=2.5, adjCol = c(0.5,1), offsetCol=1,  # arrange col names
              margins = c(3, 8), # adds to the bottom and right margins
              lwid=lwid
    )
    title(main=paste0("DE genes in ", list_DEgenes_whole_series[i]),line=-0.5,oma=T,adj=0.75,
          cex.main = 2.5)
    dev.off()
}   

 #~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~
for (j in 1:length(l.subfolder_names)) {
  png(paste0("./images/2.11_heatmap_DEgenesAll_B_thr", LETTERS[j], ".png"), width=2160, height=1440, units = "px", bg = "white")
  par(oma=c(1,1,1,1))
  lwid=c(5,8)
  f.heatmap(cbind(0,get(paste0('DEgenesAll_B_thr',LETTERS[j]))), keysize=0.4, key.par = list(cex=1.25), key.title=NA, key.ylab=NA,
            labCol=as.factor(c('zero',gsub("DE_", "", colnames(get(paste0('DEgenesAll_B_thr',LETTERS[j])))))), 
            srtCol=0, cexCol=1.5, adjCol = c(0.5,1), offsetCol=1,  
            margins = c(3, 8), 
            colsep=c(6, 11, 15),
            lwid=lwid
  )
  title(main=paste0("DE genes_baseline BSLN_", l.subfolder_names[j]),line=-0.5,oma=T,adj=0.75,
        cex.main = 2.5)
  dev.off() }
  
  
#### <---------- 2.12 Extract info from heatmaps ---------> ####
dir_clust <- '/home/beckum/Dropbox/0MSc_Bioinfo_Uppsala/sars/3p_Hawman_RockyMountainLabs_sars_cov2_Zika_RM/output_files/C_filtering_threshold/clusters'
dtsets <- c('DEgenesAll_B_noDPC1_thrC', 'DEgenesAll_noDPC1_thrC',
            'DEgenesAll_WPV1_4_thrC', 'DEgenesAll_WPV5_thrC',
            'DEgenesAll_B_thrC', 'DEgenesAll_thrC')

for (k in 1:length(dtsets)){
  folder <- paste0(dir_clust, '/', dtsets[k], '/')
  if (!dir.exists(folder)){
    folder2 <- dir.create(folder, recursive = TRUE)
    }
  else {
    print("Dir already exists!")
  }
}


for (j in 1:length(dtsets)) {
  deins <- f.heatmap(cbind(0,get(dtsets[j])))
  print(paste0('nr of clusters: ', length(unique(deins))))
  for (i in unique(deins)) {
    cluster_i <- deins[deins==i]
    write.table(names(cluster_i), file = paste0(dir_clust, '/', dtsets[j], '/', 'cluster_', i), row.names = F, sep = "\t", quote = F)
    print(length(cluster_i))
  }
}

out <- f.heatmap(cbind(0,get('DEgenesAll_WPV5_thrC')), keysize=0.4, key.par = list(cex=1.25), key.title=NA, key.ylab=NA,
                labCol=as.factor(c('zero',gsub("DE_", "", colnames(get('DEgenesAll_WPV5_thrC'))))), 
                srtCol=0, cexCol=0.8, adjCol = c(0.5,1), offsetCol=1,  
                margins = c(3, 8),
                colsep=c(5, 11, 14),
                lwid=lwid)

out1 <- f.heatmap(get('DEgenesAll_WPV5_thrC'), keysize=0.4, key.par = list(cex=1.25), key.title=NA, key.ylab=NA,
                 labCol=as.factor(gsub("DE_", "", colnames(get('DEgenesAll_WPV5_thrC')))),
                 srtCol=0, cexCol=1, adjCol = c(0.5,1), offsetCol=1, 
                 margins = c(3, 8), 
                 colsep=c(4, 11, 14),
                 lwid=lwid)





#### <---------- 2.13 Functional analyses ---------> ####
library("clusterProfiler")
library("RDAVIDWebService")

clust_1 <- read.table(paste0(dir_clust, '/DEgenesAll_B_noDPC1_thrC/cluster_1'), header = T)


### Plot the clustering
dir_plot <- '/home/beckum/Dropbox/0MSc_Bioinfo_Uppsala/sars/3p_Hawman_RockyMountainLabs_sars_cov2_Zika_RM/output_files/C_filtering_threshold/'
clust_2_plot <- read.table(paste0(dir_plot, 'clusters_nozero/DEgenesAll_WPV1_4_thrC/cluster_2_GO_enrichment'), 
                           sep='\t', header = F,
                           na.strings = "NA")
clust_2_plot <- clust_2_plot[, colSums(is.na(clust_2_plot)) == 0] 

heatmap.2(as.matrix(clust_2_plot), Rowv=FALSE, Colv=FALSE, dendrogram="none", 
          main="Cluster 2_nozero", col=c('black', 'green'), 
          tracecol="#303030", trace="both", 
          hline=F, vline=F,
          notecol="black", notecex=0.8, 
          key.title=NA, key.ylab=NA,
          srtCol=45, cexCol=1, #adjCol = c(0.5,1), offsetCol=1,  
          keysize = 0.6, margins=c(20, 5),
          par(oma=c(1,1,1,1)),
          lwid=c(5,8))


#### <------2.14 Check for 1st degree and 2nd degree neighbors of IRF3/7 -----> ####

#translate the ENSMBL to huSymbols
dir_ITA_files <- '/home/beckum/Dropbox/0MSc_Bioinfo_Uppsala/sars/ITA_Blood_Package/'
dir_files_saved <- "/home/beckum/Desktop/0MSc_Bioinfo_Uppsala/sars/3p_Hawman_RockyMountainLabs_sars_cov2_Zika_RM/output_files/C_filtering_threshold/husymb_WPV1-4_WPV5/"


#save out the gene names to be translated:

write.table(rownames(DEgenesAll_WPV5_thrC), file = paste0(dir_files_saved, 'DEnames_WPV5_thrC'), row.names = F, sep = "\t", quote = F)
write.table(rownames(DEgenesAll_WPV1_4_thrC), file = paste0(dir_files_saved, 'DEnames_WPV1_4_thrC'), row.names = F, sep = "\t", quote = F)

DEnames_WPV5_thrC <- read.table(paste0(dir_files_saved, 'DEnames_WPV5_thrC'), sep="\t", header=T)
DEnames_WPV1_4_thrC <- read.table(paste0(dir_files_saved, 'DEnames_WPV1_4_thrC'), sep="\t", header=T)

mkSymb_ENSMM_file = read.table(paste0(dir_ITA_files, "mmulatta.SYMBOL.txt"), sep="\t", header=F)
ENSMM_ENSG_file = read.table(paste0(dir_ITA_files, "mmulatta.hsapiens_orth.txt"), sep="\t", header=F)
ENSG_hSymb_file = read.table(paste0(dir_ITA_files, "hsapiens.SYMBOL.txt"), sep="\t", header=F)

# translate the data
DE_symb_ENSMM <- merge(DEnames_WPV5_thrC, mkSymb_ENSMM_file, by.x=c('x'), by.y =c('V1'))
write.table(DE_symb_ENSMM, file = paste0(dir_files_saved, 'DEnames_WPV5_thrC', "_hSymb"), row.names = F, sep = "\t", quote = F)
  
DE_symb_ENSMM2 <- merge(DEnames_WPV1_4_thrC, mkSymb_ENSMM_file, by.x=c('x'), by.y =c('V1'))
write.table(DE_symb_ENSMM2, file = paste0(dir_files_saved, 'DEnames_WPV1_4_thrC', "_hSymb"), row.names = F, sep = "\t", quote = F)
######################################

#set of genes that I have to search for:
search_set_IRF3 <- c('DHPS', 'TMC6', 'CENPT', 'IL11RA', 'SLC25A45', 
              'UNC45A', 'CIAO3', 'TM7SF2', 'DAXX', 'TRIM28',
              'NCDN', 'ZNF692', 'APG2', 'DXO', 'RHOT2',
              'INO80E', 'HDAC6')
              
search_set_IRF7 <- c('IFI6', 'MX1', 'ISG15', 'IFI35',
                     'TREX1', 'PML', 'OASL', 'REC8',
                     'MOB3C', 'RHBDF2')

# search for the IRF3 and IRF 7- 1st degree neighbours into WPV1-4
res1 <- DE_symb_ENSMM[DE_symb_ENSMM$V2 %in% search_set_IRF3, ]
res2 <- DE_symb_ENSMM[DE_symb_ENSMM$V2 %in% search_set_IRF7, ]

# search for the IRF3 and IRF 7- 1st degree neighbours into WPV5
res3 <- DE_symb_ENSMM2[DE_symb_ENSMM2$V2 %in% search_set_IRF3, ]
res4 <- DE_symb_ENSMM2[DE_symb_ENSMM2$V2 %in% search_set_IRF7, ]
