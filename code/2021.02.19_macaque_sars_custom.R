#### <---------- 1. Data location and ID ---------> ####
#     * data location:	 
#     * <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155363> 
#     * <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP274118&o=acc_s%3Aa> 
#     
# * Title		Transcriptomic profiles of whole blood from SARS-CoV-2 rhesus macaques 
# * Organism 		Macaca mulatta (= rhesus monkey = rhesus macaque); primate 
# * Experiment type 	Expression profiling by high throughput sequencing 
# * Summary 		Total RNAseq profiling of whole blood collected longitudinally from 8 macaques infected with SARS-CoV-2 
# * Overall design 	Longitudinal host response profiling was performed to assess host responses associated with moderate COVID-19 pneumonia in a rhesus macaque model. 
# * tissue:		whole blood 
# * experimental conditions:	control & SARS 
# 
# * Supplementary files: 
# *  GSE155363_COM_master_counts_refseq.txt.gz   - transcript counts (quantified expression for each associated sample )
# * GSE155363_COM_master_targets.txt.gz           - data about animals and time points and experimental conditions (associated sample metadata)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### <---------- Data ID ---------> ####
# - 48 samples:	
#     - BioSample accession numbers:		SAMN156676(02) - (50)
# - GEO sample accession number:		GSM47001(01) - (48)
# - SRA sample identifier:		    SRS7106(789) - (836)
# 
# 
# - 1 Project:
#     - BioProject accession number:		PRJNA649525
# - Series (same as project?):
#     - GEO series accession number:		GSE155363
# 
# 
# - 1 platform: Platform ID -  a summary description of the array or sequencer 
# - GEO platform accession number: 	GPL23949 	Illumina HiSeq 4000 (Macaca mulatta)
# - SRA (NCBI accession numbers)
# - SRA experiment:	                SRX88438(21) - (68)		an object that contains the metadata describing the library, platform selection, and processing parameters involved in a particular sequencing experiment.
# - SRA study:                        SRP274118 		A Study is an object that contains the project metadata describing a sequencing study or project. Imported from BioProject.
# - SRA runs:                         SRR123441(31) - (78)		an object that contains actual sequencing data for a particular sequencing experiment. Experiments may contain many Runs depending on the number of sequencing instrument runs that were needed. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("Blood_PBMC2.R")
library(limma)
library(edgeR)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### <---------- 2.1. Read the counts and the sample information into R ---------> ####
# Set your working directory to the folder where you have downloaded the datasets
setwd(paste0("/home/", Sys.getenv("SUDO_USER"),"/Dropbox/0MSc_Bioinfo_Uppsala/sars/1p_GSE155363_Rasmussen_ColumbiaU_sars/GSE155363_Rasmussen_ColumbiaU_code"))


counts_file <- "../GSE155363_Rasmussen_ColumbiaU_COVID_orig_files/GSE155363_COM_master_counts_refseq.txt"
sample_file <- "../GSE155363_Rasmussen_ColumbiaU_COVID_orig_files/GSE155363_COM_master_targets.txt"
counts_gene <- read.delim(counts_file, header = TRUE, sep='\t', na.strings="NA", row.names=1)
sample_info <- read.delim(sample_file, header = TRUE, sep='\t', na.strings="NA", row.names=1)

geneid <- row.names(counts_gene)   # get gene names
colnames(counts_gene) <- gsub('00','', colnames(counts_gene))
rownames(sample_info) <- gsub('00','', rownames(sample_info))

# shorten the names used in describing the experimental conditions
sample_info$MonkeyID <- as.factor(gsub("Monkey", "Mk", sample_info$MonkeyID))
sample_info$SampleID <- as.factor(gsub("00", "", sample_info$SampleID))
sample_info$Condition <- as.factor(gsub("Control", "Ctr", sample_info$Condition))
sample_info$Time_Condition <- as.factor(gsub("Control", "Ctr", sample_info$Time_Condition))
sample_info$Monkey_Time_Condition <- as.factor(gsub("Control", "Ctr", gsub("Monkey", "Mk", sample_info$Monkey_Time_Condition)))  # nested replacement!
#levels(sample_info$Monkey_Time_Condition)
sample_info$Monkey_Time <- as.factor(gsub("Monkey", "Mk", sample_info$Monkey_Time))


#re-arrange the levels of the Time column, will need it for correct plotting later
# First introduce the extra time points that we need to sort the data: D0x, D1x and D3x
sample_info$Time <- ifelse(sample_info$MonkeyID %in% c('Mk1','Mk2','Mk3', 'Mk4'), 
                       paste(sample_info$Time, 'x', sep=''), sample_info$Time)
sample_info$Time <- as.factor(sample_info$Time)

# order the levels
sample_info$Time  <- factor(sample_info$Time, levels = c("D0", "D0x","D1","D1x", "D3","D3x", "D5", "D7", "D10", "D12", "D15", "D17"))
sample_info$Condition  <- factor(sample_info$Condition, levels = c("Ctr", "Case"))

# change the rest of the columns and order 
sample_info$Time_Condition <- as.factor(paste(sample_info$Time, sample_info$Condition, sep='_'))
sample_info$Time_Condition <- factor(sample_info$Time_Condition, 
	levels = c("D0_Ctr", "D0x_Ctr" ,"D1_Case", "D1x_Case", "D3_Case", "D3x_Case", "D5_Case", "D7_Case", "D10_Case", "D12_Case", "D15_Case", "D17_Case")) 

sample_info$Monkey_Time_Condition <- as.factor(paste(sample_info$MonkeyID, sample_info$Time_Condition, sep='_'))
sample_info$Monkey_Time <- as.factor(paste(sample_info$MonkeyID, sample_info$Time, sep='_'))
head(sample_info)


#Make the dataframe into a DGEList and attach the information regarding the samples.
counts_list <- DGEList(counts_gene, samples=sample_info, genes=geneid)
counts_list$samples$group <-counts_list$samples$Time_Condition 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### <---------- 2.2. Data pre-processing ---------> ####
# - manual filtering: threshold > 100 counts/ 25% of data of timepoint

# could make this easier since all tp have the same length ...
len_time_points <- aggregate(counts_list$samples$SampleID ~ counts_list$samples$Time, data=counts_list$samples, FUN=length)

for (k in (1:nrow(len_time_points))){
    samples_time_point <- counts_list$counts[,row.names(counts_list$samples[counts_list$samples$Time==len_time_points[k,1],])]
    dim(counts_list$counts)
    keep1 <- rowSums(samples_time_point>=100)>=round((ncol(samples_time_point)/4), digits =0)
    filtered_list_manual <- counts_list[keep1,]
}

lcpm_init <- cpm(counts_list, log=TRUE)  
lcpm_fl <- cpm(filtered_list_manual, log=TRUE)  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### <---------- 2.3 Plot the 4 filtered datasets ---------> ####
samplenames <- colnames(filtered_list_manual)
nsamples <- ncol(filtered_list_manual)

# Make a color palette to accommodate 48 colors: 3 concatenated palettes
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
par(mar=c(5.1,4.1,4.1,2.1))

# plot each of the filtered dataset I obtained:
par(mfrow=c(2,2))
plot(density(lcpm_init[,1]), col=pal[1], lwd=2, ylim=c(0,0.65),
     las=1, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm_init[,i])
    lines(den$x, den$y, col=pal[i], lwd=2)
}
#~~~~~~~~~~~~~~~~~
plot(density(lcpm_fl[,1]), col=pal[1], lwd=2, ylim=c(0,0.65),
     las=1, main="", xlab="")
title(main="B. Filtered data manual, >100/25%", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm_fl[,i])
    lines(den$x, den$y, col=pal[i], lwd=2)
}
#~~~~~~~~~~~~~~~
# plot legend in a spare space in the image
plot.new()
legend(x = "center",
       legend = samplenames, ncol=6, bty='n', xpd=TRUE,
       text.col=pal, fill=pal, cex=.9, horiz = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### <---------- 2.4 Normalizing gene expression distributions ---------> ####
filtered_list_manual <- calcNormFactors(filtered_list_manual, method = "TMM")

par(mfrow=c(2,1))
boxplot(lcpm_init, las=2,  col=pal, main="", cex.axis=0.8 )
title(main="A. Unnormalised data", ylab="Log-cpm")

lcpm_fl <- cpm(filtered_list_manual, log=TRUE)  # transform for plotting
boxplot(lcpm_fl, las=2, col=pal, main="", cex.axis=0.8)
title(main="B. Normalised data", ylab="Log-cpm")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### <---------- 2.5 Unsupervised clustering of samples ---------> ####
par(mfrow=c(3,1))
#get the values of the columns you are interested in to later assign colors to each level
col.time <-as.factor(filtered_list_manual$samples$Time)
col.condition <- as.factor(filtered_list_manual$samples$Condition)
col.monkey <- as.factor(filtered_list_manual$samples$MonkeyID)

levels(col.time) <-   c("#CC0000", "#006633", '#FFFF33', '#99CCFF',  '#FF33CC', '#000000', '#0000FF', '#FF9933', '#00FF00', '#03f0fc', '#d703fc', '#c6fc03')  
levels(col.condition) <- c("#CC0000", "#006633") 
levels(col.monkey) <-  c("#CC0000", "#006633", '#FFFF33', '#0000FF', '#FF33CC', '#000000', '#99CCFF', '#FF9933')

col.time <- as.character(col.time)
col.condition <- as.character(col.condition)
col.monkey <- as.character(col.monkey)

plotMDS(lcpm_fl, labels=filtered_list_manual$samples$Monkey_Time_Condition, col = col.time, pch=col.condition, dim=c(1,2), cex.axis=1.5)
title(main="A. Time colorcoded: dimensions 1, 2")
legend("bottomright", legend=levels(filtered_list_manual$samples$Time), text.col=unique(col.time),  bty="n", cex=1.5)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### <---------- 2.6 Creating a design matrix and contrasts on all data points ---------> ####
tc <- factor(filtered_list_manual$samples$Time_Condition)
mk <- factor(filtered_list_manual$samples$MonkeyID)

design <- model.matrix(~0+tc+mk, data=filtered_list_manual$samples)

#rename columns
colnames(design) <- as.factor(gsub("tc", "", gsub("mk", "", colnames(design))))

nonEstimable(design)
design <- design[,!colnames(design)==nonEstimable(design)]

# mean of all D0 columns/ row
cm <- makeContrasts(
                    DE_tp1 = (D1_Case + D1x_Case) - (D0_Ctr + D0x_Ctr),
                    DE_tp3 = (D3_Case + D3x_Case) - (D0_Ctr + D0x_Ctr),
                    DE_tp5 = D5_Case - D0_Ctr,
                    DE_tp7 = D7_Case - D0_Ctr,
                    DE_tp10 = D10_Case - D0_Ctr,
                    DE_tp12 = D12_Case - D0_Ctr,
                    DE_tp15 = D15_Case - D0_Ctr,
                    DE_tp17 = D17_Case - D0_Ctr,
                    levels=colnames(design)
                    )

# Removing heteroscedascity from count data
v <- voom(filtered_list_manual, design, plot=TRUE) # convert raw counts to log-CPM and calc precision weights

# Fitting linear models for comparisons of interest
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cm)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

#Examining the number of DE genes
dt2 <-decideTests(efit)
summary(decideTests(efit))
#----------------------------------------------------

# get the DE genes and write each per time point to a file
f.makeDir <- function(myDir, subDir) {
    if (!dir.exists(paste0(myDir,subDir))) {
        dir.create(file.path(paste0(myDir,subDir), showWarnings = FALSE)
    } else {
        print("Dir already exists!")
}}


initial_loc <- getwd()

output_dir <- paste0(basename(dir_files_DEgenes), "_hSymb_transl_All_genes/")
makeDir(initial_loc, output_dir)


coef_range <- c(1,3,5,7,10,12,15,17)

for (i in coef_range) {      
	DE_genes <- topTable(efit, coef=which(coef_range == i), n=Inf, p.value=0.05, adjust.method = 'fdr', sort.by = 'P')                                  
	assign(paste0("DE_tp_",i), DE_genes)  # replace 'assign' this when you have time!!
	write.table(DE_genes, file = paste0(output_folder, paste0("C_8mk_12tp_D0_vs_D",i,".txt")), row.names = F, sep = "\t", quote = F)
}
#-----------------------------------------------------------------

# graphical representations of DE results
output_plots <- paste0("GSE155363_Rasmussen_ColumbiaU_COVID_output_files/", Sys.Date(), '_8mk_12tp', '/')
makeDir(initial_loc,output_plots)

range_cols <- 1:8

for (i in range_cols) {      
    jpeg(paste0(initial_loc,output_plots,"plotMD_col",i)
    plotMD(vfit, column=i, status=dt2[,i], col=c('red', 'blue'), main=colnames(vfit)[i])                                  
    dev.off()
}
#-----------------------------------------------------------------


# data (logFC) from TopTable
f.extractCol <- function(col_range, col_name){
    for (i in col_range) {                                       
    assign(paste0("DE_tp_",i), get(paste0("DE_tp_",i)$col_name))              # replace 'assign' this when you have time!!
    }
}

range_12tp <- c(1,3,5,7,12)
extractCol(range_12tp, genes)
DE_12tp_only <- unique(c(DEgenes_tp1, DEgenes_tp3, DEgenes_tp5, DEgenes_tp7, DEgenes_tp12))

c <- topTable(efit, n=Inf)
counts_12tp <- c[which(rownames(c) %in% DE_12tp),]
count_mx3 <- as.matrix(counts_12tp[,which(colnames(counts_12tp) %in% c('DEgenes_tp1', 'DEgenes_tp3', 'DEgenes_tp5', 'DEgenes_tp7', 'DEgenes_tp12'))])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

jpeg(paste0(initial_loc,output_plots,"Heatmap_counts_mx3"))
f.heatmap(count_mx3)
dev.off()
