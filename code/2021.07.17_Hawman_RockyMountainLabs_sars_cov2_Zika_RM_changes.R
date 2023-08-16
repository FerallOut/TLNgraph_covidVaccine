# All data from whole blood
#### <---------- 1. Data location and ID ---------> ####
#
# Animal species: Macaca mulatta (= rhesus monkey = rhesus macaque); primate 
#
# Animal Groups
# *1 group received RNA vaccine for SARS CoV -2 (6 animals)
# *1 group receive RNA vaccine for Zika virus (6 animals)
#
# Vaccination regime
#
# Challenged
#
# Necropsies
#
# Sample Types and Time points
#     * Whole blood
#        *BSLN,
#        *1wk, 2wk, 3wk,and 4wk (after vaccination)
#        *0d (wk5, just before challenge), 1d, 3d, 5d, 7d (after challenge)
#
# * Experiment type 	Expression profiling by high throughput sequencing 
# * Organization name 	Washington University 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# - 120 samples
# - 6 Mk/tp/condition
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#### <---------- 2. Obtaining the differentially expressed genes ---------> ####
source("Blood_PBMC2.R")
libr <- c("limma", "edgeR", "RColorBrewer", "ggplot2")
invisible(lapply(libr, library, character.only = TRUE))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### <---------- 2.1. Read the counts and the sample information into R ---------> ####
setwd(paste0("/home/", Sys.getenv("SUDO_USER"),"/Dropbox/0MSc_Bioinfo_Uppsala/sars/3p_Hawman_RockyMountainLabs_sars_cov2_Zika_RM"))

# read in the counts 
counts_file <-  "/home/beckum/Dropbox/0MSc_Bioinfo_Uppsala/sars/3p_Hawman_RockyMountainLabs_sars_cov2_Zika_RM/orig_files/SARS_CoV2_01_RM.txt"
counts_gene <- read.delim(counts_file, header = TRUE, sep='\t', na.strings="NA" , row.names=1)
geneid <- row.names(counts_gene)   # get gene names
#summary(counts_gene)

# make the target dataframe by extracting info from file names
sampleID <- sapply(strsplit(colnames(counts_gene),"_"), getElement, 1)
condition <- factor(sapply(strsplit(colnames(counts_gene),"_"), getElement, 4) )
levels(condition) <- c("cov","zika")     

animalID <- factor(sapply(strsplit(colnames(counts_gene),"_"), getElement, 5) )

ord_levels = c("BSLN", paste0("WPV",c(1:5)), paste0("DPC",seq(1,7,2) ))
timepoint <- factor(sapply(strsplit(colnames(counts_gene),"_"), getElement, 6) )
timepoint  <- factor(timepoint, levels = c("BSLN", paste0("w",c(1:5)), paste0("Cd",seq(1,7,2) ))  ) 
levels(timepoint) <- ord_levels


RNAc <- factor(sapply(strsplit(colnames(counts_gene),"_"), getElement, 7) )
libr <- factor(sapply(strsplit(colnames(counts_gene),"_"), getElement, 8) )
GL <- factor(sapply(strsplit(colnames(counts_gene),"_"), getElement, 9) )

# make a 'group' column
group <- factor(paste(timepoint, condition, sep='_'))
group  <- factor(group, levels = paste0(rep(ord_levels, each=2), c("_cov", "_zika") )

target <- data.frame(group, sampleID, condition, animalID, timepoint, RNAc, libr, GL)

# rename columns of 'counts_gene' and rows of 'target' file
colnames(counts_gene) <- sampleID
rownames(target) <- colnames(counts_gene)

#Make the dataframe into a DGEList and attach the information regarding the samples.
counts_list <- DGEList(counts_gene, samples=target, genes=geneid)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~


f.scatter_plot <- function(data, xvar="sampleID", yvar="lib.size", y_annot, fillcol="timepoint", titleY, title_plot) {
    ggplot(data=data, aes(x={{xvar}}, y={{yvar}}, color=fillcol),) +
    geom_point() +
    geom_point(data=data[data$fillcol == "DPC1",],color="red",size=3) +
    #geom_text(aes(label=ifelse(data$fillcol=='Cd1',as.character(data$fillcol),'')),color="black", hjust=0,vjust=1)+   # mark the Cd1 samples in red
    geom_vline(xintercept=data[ !duplicated(data$fillcol), ]$xvar, linetype="dashed", size = 0.3)+ 
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, face="bold", size=10, vjust=0.4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
          plot.title = element_text(hjust = 0.5),                     
          axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0))) +    
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+    
    scale_x_discrete(breaks=data$xvar[seq(1,length(data$xvar),by=3)]) +    
    labs(x = "SampleID", y = title_Y)+
    ggtitle(title_plot)+
    annotate("text",
             x= data[ !duplicated(data$fillcol), ]$xvar,
             y = y_annot,
             label = data[ !duplicated(data$fillcol), ]$fillcol,
             color = "darkred",
             hjust=-0.5) +
    annotate(geom = "rect", xmin = "G062", xmax = "G092",
             ymin = 67000000, ymax = 78000000, fill = "red", alpha = 0.2)
}


# take a peak at the data
png("./images/2.1_Library_sizes.png", width=800, height=400, units = "px", bg = "white")
f.scatter_plot(counts_list$samples, sampleID, lib.size, 60000000, timepoint, "Total number of counts", "Sample raw counts in each sample")
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# check how many genes are expressed (count > 100)
counts_list$samples$nr_expr_genes <- colSums(counts_list$counts > 100)

png("./images/2.1_Gene_count_above_100.png", width=800, height=400, units = "px", bg = "white")
f.scatter_plot(counts_list$samples, sampleID, nr_expr_genes, 11500, timepoint, "Nr of genes expressed (count > 100)", "Number of expressed genes in each sample")
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# how many are 0 expressed genes?
counts_list$samples$zero_counts <- colSums(counts_list$counts == 0)

png("./images/2.1_Gene_count_zero.png", width=800, height=400, units = "px", bg = "white")
f.scatter_plot(counts_list$samples, sampleID, zero_counts, 15500, timepoint, "Nr of genes with count = 0", "Number of expressed genes that have count = 0 in each sample")
dev.off()



#### <---------- 2.2. Data pre-processing ---------> ####
len_time_points <- aggregate(counts_list$samples$sampleID ~ counts_list$samples$timepoint, data=counts_list$samples, FUN=length)

filters_used <-c(100,30,10)
filtered_count_sets <- vector('character', length(filters_used))

for (a in filters_used){
  for (b in (1:nrow(len_time_points))){
    samples_time_point <- counts_list$counts[,row.names(counts_list$samples[counts_list$samples$timepoint==len_time_points[b,1],])]
    assign ( paste0("keep", a), rowSums(samples_time_point>=a)>=round((ncol(samples_time_point)/4), digits =0))			# change 'assign' later!
    assign ( paste0('filtered_counts',a), counts_list[eval(parse(text = paste0("keep", a) )),])
  }
  filtered_count_sets[which(filters_used==a)] <- paste0('filtered_counts',a)       
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# - filterByExpr
keepFBE <- filterByExpr(counts_list, group=group, min.count=1)
filtered_countsFBE <- counts_list[keepFBE, ,]

# the complete list with all the sets that were filtered with different filter values
filtered_count_sets <- append(filtered_count_sets, 'filtered_countsFBE')

# the complete list with all the sets, the initial one and all the filtered ones
count_sets <- append('counts_list', filtered_count_sets)



#### <---------- 2.3 Plot the filtered datasets vs unfiltered ---------> ####
samplenames <- colnames(counts_list)
nsamples <- ncol(counts_list)

# Make a color palette to accommodate 48 colors
pal = c(hcl(0,100, seq(20,100, length.out=8)),
        hcl(60,100, seq(20,100, length.out=8)),
        hcl(120,100, seq(20,100, length.out=8)),
        hcl(180,100, seq(20,100, length.out=8)),
        hcl(240,100, seq(20,100, length.out=8)),
        hcl(300,100, seq(20,100, length.out=8))
)


# calculate cutoff
L <- mean(counts_list$samples$lib.size) * 1e-6
M <- median(counts_list$samples$lib.size) * 1e-6
c(L, M)
lcpm.cutoff <- log2(100/M + 2/L)  # transform 100 counts into cpm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot each of the filtered dataset I obtained
png("./images/2.3_Filtered_unfiltered_datasets.png", width=800, height=680, units = "px", bg = "white")
par(mfrow=c(3,2))

# not all the timepoints exist for all monkeys 
if ((ceiling((length(filtered_count_sets)+1)/2)) %% 2 == 1) {
  par(mfrow=c(ceiling((length(filtered_count_sets)+1)/2),2))
} else {
  par(mfrow=c(ceiling((length(filtered_count_sets)+1)/2)+1,2))
}


for (a in 1:length(count_sets)) {
  plot(density(cpm(get(count_sets[a]), log=TRUE)[,1]), col=pal[1], lwd=2, ylim=c(0,0.65),
       las=1, main="", xlab="")
  title(main=paste0(LETTERS[a-1],'. ',count_sets[a]))
  for (i in 2:nsamples){
    den <- density(cpm(get(count_sets[a]), log=TRUE)[,i])
    lines(den$x, den$y, col=pal[i], lwd=2)
  }
}

par(mfg=c(1,1,3,2))

abline(v=log2(100/M + 2/L), lty=3)
text(x=log2(100/M + 2/L), y=0.4, "A",srt=0,pos=2)

abline(v=log2(30/M + 2/L), lty=3)
text(x=log2(30/M + 2/L), y=0.4, "B",srt=0,pos=2)

abline(v=log2(10/M + 2/L), lty=3)
text(x=log2(10/M + 2/L), y=0.4, "C",srt=0,pos=2)


# plot legend in a spare space in the image
par(mfg=c(3,2,3,2))
plot.new()
legend(x = "center",
       legend = samplenames, ncol=6, bty='n', xpd=TRUE,
       text.col=pal, fill=pal, cex=.9, horiz = FALSE)

dev.off()



#### <---------- 2.4 Normalizing gene expression distributions ---------> ####

#### plot the data
png("./images/2.4_Normalized_datasets.png", width=1440, height=800, units = "px", bg = "white")

counts_list <- calcNormFactors(counts_list, method = "TMM")
par(mfrow = c(2, 1),
    oma = c(4, 2, 2, 2),
    mar = c(1, 1, 0, 0), 
    mgp = c(2, 1, 0),    
    xpd = NA)           

boxplot(cpm(counts_list, log=TRUE), axes=F, col=pal, main="", cex.axis=0.8 )
title(main="A. Unnormalised data", ylab="Log-cpm")
axis(side=2, at=seq(0, 15, by=5), labels=seq(0, 15, by=5), las=1, cex.axis=0.8)


filtered_counts100 <- calcNormFactors(filtered_counts100, method = 'TMM')
boxplot(cpm(filtered_counts100, log=TRUE), axes=F, col=pal, main="", cex.axis=0.8)
title(main="B. Normalised data100", ylab="Log-cpm")
axis(side=2, at=seq(0, 15, by=5), labels=seq(0, 15, by=5), las=1, cex.axis=0.8)

axis(side=1,
     at=1: length(counts_list$samples$sampleID),
     labels= counts_list$samples$sampleID ,
     las=2, cex.axis=0.8)

#title(xlab = "SampleID",
#      outer = TRUE, line = 2)
dev.off()




#### <---------- 2.5 Unsupervised clustering of samples ---------> ####
# only work with > 100 counts filters!

#get the values of the columns you are interested in to later assign colors to each level
tp1 <- counts_list$samples$timepoint
anim1 <- counts_list$samples$animalID
cond1 <- counts_list$samples$condition

col.tp1 <- c("#CC0000", "#006633", '#FFFF33', '#99CCFF', '#FF33CC', '#000000', '#d703fc','#0000FF', '#FF9933', '#00FF00', '#B0B0B0', '#c6fc03')[tp1]
col.anim1<- c("#CC0000", "#006633", '#FFFF33', '#99CCFF', '#FF33CC', '#000000', '#0000FF', '#FF9933', '#00FF00', '#B0B0B0', '#d703fc', '#c6fc03')[anim1]
points <- c(21,24)  # pch=21 circle; 24=triangle
pch.cond1 <- points[cond1]


png("./images/2.5_Unsupervised_clustering.png", width=800, height=800, units = "px", bg = "white")

par(mfrow = c(1, 2),
    oma = c(2, 4, 1, 0),
    mar = c(2, 1, 1, 0),
    mgp = c(2, 1, 0),   
    xpd = NA,           
    pty='s')            

plotMDS(cpm(filtered_counts100, log=TRUE),
        bg = col.anim1, col=col.anim1, pch=pch.cond1, dim=c(1,2),
        xaxt = "n", yaxt = "n", xlab = '')
axis(side=1, at=seq(-1, 3, by=1), labels=seq(-1, 3, by=1), las=1, cex.axis=0.8)
axis(side=2, at=seq(-1.0, 1.5, by=0.5), labels=seq(-1.0, 1.5, by=0.5), las=1, cex.axis=0.8)
title(main="AnimalID colorcode: all data (dim_1,2)", cex.main=1)
legend("topright", legend=levels(anim1), fill=unique(col.anim1),bty="n", cex=1)
legend("bottomright", legend=levels(cond1), pch=points, bty="n", cex=1.1)


plotMDS(cpm(filtered_counts100, log=TRUE),
        bg = col.tp1, col=col.tp1, pch=pch.cond1, dim=c(1,2),
        xaxt = "n", yaxt = "n", xlab = '', ylab = '')
axis(side=1, at=seq(-1, 3, by=1), labels=seq(-1, 3, by=1), las=1, cex.axis=0.8)
axis(side=2, at=seq(-1.0, 1.5, by=0.5), labels=seq(-1.0, 1.5, by=0.5), las=1, cex.axis=0.8)
title(main="Timepoint colorcode: all data (dim_1,2)", cex.main=1)
legend("topright", legend=levels(tp1), fill=unique(col.tp1), bty="n", cex=1)
legend("bottomright", legend=levels(cond1), pch=points, bty="n", cex=1.1)

dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#subset out the DPC1 timepoint
filtered_noCd1 <- filtered_counts100[,!filtered_counts100$samples$timepoint %in%'DPC1']
filtered_noCd1$samples <- droplevels(filtered_noCd1$samples)

tp <- filtered_noCd1$samples$timepoint
anim <- filtered_noCd1$samples$animalID
cond <- filtered_noCd1$samples$condition
col.tp <- c("#CC0000", "#006633", '#FFFF33', '#99CCFF', '#FF33CC', '#000000', '#0000FF', '#FF9933', '#00FF00')[tp]
col.anim <- c("#CC0000", "#006633", '#FFFF33', '#99CCFF', '#FF33CC', '#000000', '#0000FF', '#FF9933', '#00FF00', '#B0B0B0', '#d703fc', '#c6fc03')[anim]
points <- c(21,24)
pch.cond <- points[cond]


#### <---------- 2.6 Creating a design matrix and contrasts on all data points ---------> ####
gr <- factor(filtered_counts100$samples$group)
mk <- factor(filtered_counts100$samples$animalID)
tp <- factor(filtered_counts100$samples$timepoint)
cond <- factor(filtered_counts100$samples$condition)

design1_100 <- model.matrix(~0+gr+mk, data=filtered_counts100$samples)

colnames(design1_100) <- as.factor(gsub("gr", "", colnames(design1_100)))

design1_100 <- design1_100[,!colnames(design1_100)==nonEstimable(design1_100)]

# mean of all D0 columns/ row
cm1_100 <- makeContrasts(
    DE_cov_WPV1 = WPV1_cov - BSLN_cov,
    DE_cov_WPV2 = WPV2_cov - BSLN_cov,
    DE_cov_WPV3 = WPV3_cov - BSLN_cov,
    DE_cov_WPV4 = WPV4_cov - BSLN_cov,
    DE_cov_WPV5 = WPV5_cov - BSLN_cov,
    
    DE_zika_WPV1 = WPV1_zika - BSLN_zika,
    DE_zika_WPV2 = WPV2_zika - BSLN_zika,
    DE_zika_WPV3 = WPV3_zika - BSLN_zika,
    DE_zika_WPV4 = WPV4_zika - BSLN_zika,
    DE_zika_WPV5 = WPV5_zika - BSLN_zika,
   
    levels=colnames(design1_100)
)


# Removing heteroscedascity from count data
png("./images/2.6_voom_orig.png", width=800, height=800, units = "px", bg = "white")
v1_100 <- voom(filtered_counts100, design1_100, plot=TRUE) 
dev.off()

# Fitting linear models for comparisons of interest
vfit1_100 <- lmFit(v1_100, design1_100)
vfit1_100 <- contrasts.fit(vfit1_100, contrasts=cm1_100)

# Calculate t-stats
efit1_100 <- eBayes(vfit1_100)
plotSA(efit1_100, main="Final model: Mean-variance trend")



f.decideTests_thresholds <- function(efit_value, dt_name) {
  #change here the values as you need them!
  l.lfc = list(0, 0.58, 1)
  l.pval = list(0.05, 0.01)
  df.lfc.pval <- expand.grid(l.pval, l.lfc)
  v.approaches <- LETTERS[1:nrow(df.lfc.pval)]
  colnames(df.lfc.pval) <- c("p_val", "lfc_val")
  rownames(df.lfc.pval) <- v.approaches
  print(df.lfc.pval)
  dt_name <- deparse(substitute(dt_name)) 
  l.dt_s <- vector(mode = "list", length = nrow(df.lfc.pval))
  
  for (i in 1: nrow(df.lfc.pval)) {
    assign(paste0(dt_name, '_', v.approaches[i]), decideTests(efit_value, lfc = df.lfc.pval[i,2], p.value = df.lfc.pval[i,1]))
    message("summary for: ", paste0(dt_name, '_', v.approaches[i]))
    print(summary(decideTests(efit_value, lfc = df.lfc.pval[i,2], p.value = df.lfc.pval[i,1])))
    message("table for: ", paste0(dt_name, '_', v.approaches[i]))
    print(table(decideTests(efit_value, lfc = df.lfc.pval[i,2], p.value = df.lfc.pval[i,1])))
    
    l.dt_s[i] <- list(get(paste0(dt_name, '_', v.approaches[i])))
    names(l.dt_s)[i] <- paste0(dt_name, '_', v.approaches[i])
  }
  return(l.dt_s)
} 

l.dt_all <- f.decideTests_thresholds(efit1_100, dt1_100)

  
  
#### <---------- 2.7 Extract the DE genes and whole df for each time point---------> ####

extractDE_write.files <- function(name_decideTests_df, name_efit){
  names_files_DE <- colnames(name_decideTests_df)
  names_files_save <- paste0(colnames(name_decideTests_df), "_logFC_list")
  
  names_files_save <- topTable(name_efit, coef=names_files_DE, n=Inf, p.value=0.05, adjust.method = 'fdr', sort.by = 'P')
  dim(names_files_save)
  write.table(names_files_save, file = paste0(output_dir, names_files_DE), row.names = F, sep = "\t", quote = F)
}

extractDE_write.files(dt1_100, efit1_100)


# because I am testing several cut-off thresholds (adj and logFC), each has one folder:
l.subfolder_names <- paste0(LETTERS[1:length(l.dt_all)], "_filtering_threshold")
l.subsubfolder_names <- c("ENSMM_data", "ENSG_data") 
l.extraction_directory_names <- c("DE_genes_logFC", "DE_genes_name", "all_genes_logFC")  # COME BACK AND REMOVE FOLDERS IF NOT NEEDED


##for just one level of recursiveness:
# for (j in seq_along(subfolder_names)){
#   folder <- dir.create(paste0(getwd(),'/', "output_files", '/', subfolder_names[j]))}

for (j in seq_along(l.subfolder_names)){
  for (k in seq(1:length(l.subsubfolder_names))){
    for (m in seq(1:length(l.extraction_directory_names))) {
      folder <- paste0(getwd(),'/', "output_files", '/', l.subfolder_names[j], '/', l.subsubfolder_names[k], '/', l.extraction_directory_names[m])
      if (!dir.exists(folder)){
      folder2 <- dir.create(folder, recursive = TRUE)
      }
      else {
        print("Dir already exists!")
      }
    }
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# extract the DE genes and write the files to disk: 
f.extractDE_write.files <- function(efit_name, l.subfolder_names, l.subsubfolder_names, l.extraction_directory_names) {
  names_files_DE <- colnames(efit1_100$contrasts)
  names_files_save_all <- paste0(names_files_DE, "_allgenes_logFC")
  names_files_save_DE <- paste0(names_files_DE, "_DEgenes_logFC")
  names_files_save_DEname <- paste0(names_files_DE, "_DEgenes_name")

  l.DEgenes <- vector(mode = "list", length = nrow(df.lfc.pval))
  
  l.lfc = list(0, 0.58, 1)
  l.pval = list(0.05, 0.01)
  df.lfc.pval <- expand.grid(l.pval, l.lfc)
  colnames(df.lfc.pval) <- c("p_val", "lfc_val")
  rownames(df.lfc.pval) <- LETTERS[1:dim(df.lfc.pval)[1]]
  print(df.lfc.pval)
  
  for (k in c(2,1,3)) {
     assign(paste0("target_dir",k), paste0(getwd(),'/', "output_files", '/', l.subfolder_names[j], '/', l.subsubfolder_names[1], '/', l.extraction_directory_names[k], '/') )
  
  #message("Printing combination of p.val and lfc that are tested:")
  for (j in seq_along(l.subfolder_names)) {
    print(df.lfc.pval[j,])
    for (i in 1: length(names_files_DE)) {
      #extract files containing all the genes, together with logFC and p.val
      assign(names_files_save_all[i], topTable(efit1_100, coef=names_files_DE[i], n=Inf))
      
      #write files containing all the genes, together with logFC and p.val
      write.table(get(names_files_save_all[i]), file = paste0(target_dir, names_files_save_all[i]), row.names = F, sep = "\t", quote = F)
     
      #extract files containing just the DE genes, together with logFC and p.val
      assign(names_files_save_DE[i], topTable(efit1_100, coef=names_files_DE[i], n=Inf, lfc = df.lfc.pval[j,2], p.value = df.lfc.pval[j,1], adjust.method = 'fdr', sort.by = 'P'))
      
      #message("dimensions for: ", paste0(names_files_save_DE[i], '_', LETTERS[j]))
      #print(paste0(dim(get(names_files_save_DE[i]))[1], 'rows'))
      
      #write files containing only DE genes, logFC and p.val, skip empty files
      if (nrow(get(names_files_save_DE[i])) >= 1){
        write.table(get(names_files_save_DE[i]), file = paste0(target_dir2, names_files_save_DE[i]), row.names = F, sep = "\t", quote = F)}
      else {next}
      
      #extract just the names of the DE genes
      assign(names_files_save_DEname[i], get(names_files_save_DE[i])$genes)
  
      #write files containing only the names of the DE genes, skip empty files
      if (length(get(names_files_save_DEname[i])) >= 1){
           write.table(get(names_files_save_DEname[i]), file = paste0(target_dir3, names_files_save_DEname[i]), row.names = F, sep = "\t", quote = F)}
      else {next}
    }      
  }
 }
}

f.extractDE_write.files(efit1_100, l.subfolder_names, l.subsubfolder_names, l.extraction_directory_names)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# remove the name of the variables containing the DE genes or names from the environment to be able to load them later from files and to loop through all of them
rm(list = ls(pattern = "^DE_cov_|DE_zika_"))


#### <---------- 2.8 Graphical representations of DE results ---------> ####

# mean-difference plot displaying log-FCs from the linear model fit against the average log-CPM values
for (thresholds in 1:length(l.dt_all)) {
  png(paste0("./images/2.8_DEgenes_plotMD", "_", names(l.dt_all)[thresholds], ".png"), width=800, height=2400, units = "px", bg = "white")
  par(mfrow = c(9, 3),     
      oma = c(0, 2, 0, 0), 
      mar = c(4, 3, 1, 0), 
      mgp = c(2, 1, 0),    
      xpd = NA,            
      pty='m')            
  for (i in 1: length(colnames(vfit1_100))) {
    plotMD(vfit1_100, column=i, status=l.dt_all[[thresholds]][,1], col=c('red', 'blue'), main=colnames(vfit1_100)[i])
  }
  dev.off()
}  


#### <---------- 2.9 Extract gene lists DE in each condition ---------> ####

## only the genes for each time point are in 'l.DE_genes' list from above

# load the necessary files and variables into the memory
library(stringr)

f.subset_safely <- function(x, index) {
  if (length(x) < index) {
    return(NA_character_) }
  x[[index]] }

f.str_split_n <- function(string, pattern, n) {
  out <- str_split(string, pattern)
  vapply(out, f.subset_safely, character(1L), index = n) }


for (j in seq_along(l.subfolder_names)) {
  target_dir <- paste0(getwd(),'/', "output_files", '/', l.subfolder_names[j], '/',
                       l.subsubfolder_names[1], '/', l.extraction_directory_names[2])
  v.filenames <- list.files(target_dir, full.names=TRUE)
  
  assign(paste0('l.DEgene_names_threshold_', LETTERS[j]), sapply(v.filenames, read.table, header=TRUE))
  l.threshold_DE_names <- setNames(eval(parse(text = paste0('l.DEgene_names_threshold_', LETTERS[j]))), 
                  f.str_split_n(substr(v.filenames, 154, 190), '_DEgenes', 1))
  assign(  paste0('l.DEgene_names_threshold_', LETTERS[j]), l.threshold_DE_names) }


#### <---------- 2.10 Extract the logFC for the DE genes in each condition ---------> ####

# load the necessary files and variables into the memory - all data has been written into files after 2.7

for (j in seq_along(l.subfolder_names)) {
  target_dir3 <- paste0(getwd(),'/', "output_files", '/', l.subfolder_names[j], '/',
                       l.subsubfolder_names[1], '/', l.extraction_directory_names[3])
  v.filenames3 <- list.files(target_dir3, full.names=TRUE)
    
  assign(paste0('l.gene_logFC_threshold_', LETTERS[j]), lapply(v.filenames3, read.table, header=TRUE))
  l.genes_logFC_theshold <- setNames(eval(parse(text = paste0('l.gene_logFC_threshold_', LETTERS[j]))), 
                                   f.str_split_n(substr(v.filenames3, 156, 190), '_allgenes', 1))
  assign(  paste0('l.gene_logFC_threshold_', LETTERS[j]), l.genes_logFC_theshold) }


# extract needed datasets
for (j in seq_along(l.subfolder_names)) {
    assign( paste0('c_thr', LETTERS[j]) , topTable(efit1_100, n=Inf))
    print(paste0('Extract logFC values for threshold ', LETTERS[j]))
    target_dir2 <- paste0(getwd(),'/', "output_files", '/', l.subfolder_names[j], '/', l.subsubfolder_names[1], '/', l.extraction_directory_names[1], '/')
  
    
# all genes in cov_zika baseline BSLN with DPC1 data                                                          
    choice_genes <- paste0('DEgenesAll_B_thr', LETTERS[j])    
    assign( choice_genes, 
            get(paste0('c_thr', LETTERS[j]))[which(rownames(get(paste0('c_thr', LETTERS[j]))) 
            %in% get(paste0('v.DEnames_cov_zika_bslBSLN_thr', LETTERS[j]))),])


    assign( paste0('DEgenesAll_B_thr', LETTERS[j]), as.matrix(get(paste0('DEgenes_all_B_subset_thr', LETTERS[j]))[,which(colnames(get(paste0('DEgenes_all_B_subset_thr', LETTERS[j]))) %in% 
            paste0("DE_", rep(c("cov_", "zika_"), each=5), c(paste0("WPV",c(1:5), "_cov") )) )]))
    print(dim(get(paste0('DEgenesAll_B_thr', LETTERS[j])))[1])
    write.table(get(paste0('DEgenesAll_B_thr', LETTERS[j])), file = paste0(target_dir2, 'DEgenesAll_B_thr', LETTERS[j]), row.names = T, sep = "\t", quote = F)

common_genes_cov_DPC <- as.integer(as.logical(DEgenes_cov_DPC1 %in% DEgenes_cov_DPC_noDPC1))
common_genes_zika_DPC <- as.integer(as.logical(DEgenes_zika_DPC1 %in% DEgenes_zika_DPC_noDPC1))

all_counts <- matrix(0, nrow=length(DEgenes_cov_DPC), ncol=2)
for (i in 1: length(DEgenes_cov_DPC)){
    all_counts[i, 1] <- DEgenes_cov_DPC[i] %in% DEgenes_cov_DPC1
    all_counts[i, 2] <- DEgenes_cov_DPC[i] %in% DEgenes_cov_DPC_noDPC1
}
colnames(all_counts) <- c('DE genes in DPC1', 'DE genes in DPC3, DPC5, DPC7')
vennDiagram(vennCounts(all_counts), 
            counts.col=c("black"), 
            circle.col = c("red", "blue"), 
            cex = 1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

