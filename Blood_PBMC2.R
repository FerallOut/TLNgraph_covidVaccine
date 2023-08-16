library(igraph)
library(dynamicTreeCut)
library(topGO)
library(GOSim)
setwd("/home/ferall/Downloads/ITA-master/ITA_Blood_Package")


######################## define main functions

## heatmap function adapted to various parameters

f.heatmap = function(dataM, 
                     colors=c('darkblue', 'mediumblue', 'dodgerblue', 'white', 'orange', 'red', 'darkred'),
                     colorsaturation=0.25,distmethod='bicor',clustermethod='ward.D2',
                     clusterdim='both',exclude=NULL,cutoffmethod='depth',cutoff=1,
                     RowColors=NA,main=paste(c(distmethod, clustermethod), collapse=", "),
                     ...){
  if(length(colnames(dataM))==0){colnames(dataM) = as.character(1:ncol(dataM))}
  
  
  fx = function(v){return(!all(is.na(v)))}
  dataM = dataM[apply(FUN=fx, X=dataM, MARGIN=1),apply(FUN=fx, X=dataM, MARGIN=2)]
  
  distmethod = match.arg(arg=distmethod, choices=c('bicor', 'pearson', 'spearman', 'direct', 'euclidean', 'binary'))
  clusterdim = match.arg(arg=clusterdim, choices=c('both', 'row', 'column', 'none'))
  cutoffmethod = match.arg(arg=cutoffmethod, choices=c('depth', 'height', 'number', 'none'))
  
  Rowv=T;Colv=T;dendrogram='both'
  if(clusterdim=='row'){Colv=F; dendrogram='row'}
  if(clusterdim=='column'){Rowv=F; dendrogram='column'}
  if(clusterdim=='none'){Rowv=F; Colv=F; dendrogram='none'} 
  
  # Removing excluded samples:
  dataM = dataM[,which(colnames(dataM)%in%exclude==F)]
  
  # Color scale:  
  neg.col = rev(colors[1:ceiling(length(colors)/2)])
  pos.col = colors[ceiling(length(colors)/2):length(colors)]
  
  bins = diff(round(quantile(abs(dataM), seq(from=0, to=1, length.out=length(pos.col))^colorsaturation, na.rm=T)/max(abs(range(dataM, na.rm=T)), na.rm=T), digits=3)*1000)
  a1 = list()
  a2 = list()
  for(i in 1:length(bins)){
    a1[[i]] = t(colorRamp(neg.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
    a2[[i]] = t(colorRamp(pos.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
  }
  a1 = matrix(unlist(a1), ncol=3, byrow=T); a2 = matrix(unlist(a2), ncol=3, byrow=T)
  
  b1 = rgb(red=a1[,1], green=a1[,2], blue=a1[,3], max=255)
  b2 = rgb(red=a2[,1], green=a2[,2], blue=a2[,3], max=255)
  
  #c = abs(range(dataM))
  #if(c[1]>c[2]){
  # b2 = b2[1:round(length(b2)*c[2]/c[1])]
  #} else {
  # b1 = b1[1:round(length(b1)*c[1]/c[2])]
  #}
  color.vector = c(rev(b1), b2)
  
  # Distance metric
  if(distmethod=='bicor'){fd = function(x){return(as.dist(1-bicor(t(x), use = 'pairwise.complete.obs')))}}
  if(distmethod=='pearson'){fd = function(x){return(as.dist(1-cor(t(x), method='pearson', use = 'pairwise.complete.obs')))}}
  if(distmethod=='spearman'){fd = function(x){return(as.dist(1-cor(t(x), method='spearman', use = 'pairwise.complete.obs')))}}
  if(distmethod=='direct'){fd = function(x){return(as.dist(x))}}
  if(distmethod=='euclidean'){fd = function(x){return(dist(x))}}
  if(distmethod=='binary'){fd = function(x){return(dist(x, method="binary"))}}
  
  # Clustering method
  fh = function(x){return(stats::hclust(x,method=clustermethod))}
  
  # Rowside colors
  
  if(cutoffmethod=='depth'){fc = function(M){return(cutreeHybrid(dendro=fh(fd(M)), distM=as.matrix(fd(M)), deepSplit=cutoff, verbose=0, minClusterSize=1)$labels)}}
  if(cutoffmethod=='height'){fc = function(M){return(cutree(fh(fd(M)), h=cutoff))}}
  if(cutoffmethod=='number'){fc = function(M){return(cutree(fh(fd(M)), k=cutoff))}}
  if(cutoffmethod=='none'){fc = function(M){return(rep(1, nrow(M)))}}
  
  if(dendrogram%in%c('none','column')){rowcol=rep('grey70', nrow(dataM))} else {
    rowcol = c('blue', 'red', 'orange', 'skyblue', 'yellow', 'black', 'darkblue', 'cyan', 'darkred', 'darkgreen', 'pink', 'purple', 'gray10')[suppressWarnings(fc(dataM))]
    if(length(rowcol)!=nrow(dataM)){rowcol=rep("gray", nrow(dataM))}
  }
  if(!all(is.na(RowColors))){rowcol=RowColors}
  
  hm = heatmap.2(dataM,
                 col=color.vector,
                 hclustfun=fh,
                 distfun=fd,
                 trace='none',
                 Rowv=Rowv,
                 Colv=Colv,
                 dendrogram=dendrogram,
                 RowSideColors=rowcol,
                 symbreaks=T, main=main, ...)
  
  names(rowcol) = rownames(dataM)
  #return(rowcol)
  return(rev(rowcol[hm$rowInd]))
}
#------------------------------------------------------------


## GO term enrichment
f.GOenrich = function(fg, cut=0.01){
  fg.x = Ann.Entrez[fg]; fg.x = fg.x[!is.na(fg.x)]
  bg.x = Ann.Entrez[rownames(TopCorrGene)]; bg.x = bg.x[!is.na(bg.x)]
  A = GOenrichment(genesOfInterest=fg.x, allgenes=bg.x, cutoff=cut)
  a = sort(A$p.values)
  b = A$GOTerms$Term
  names(b) = A$GOTerms$go_id
  return(cbind(a, b[names(a)]))
}
#------------------------------------------------------------


## TLN construction and plotting
f.tree.leaf.network = function(hc, nodesize=7){
  ftM1 = -hc$merge
  fx = function(i){
    v = ftM1[i,]
    #if(v[1]<0 & v[2]>0){return(cbind(v,c(v[1],-i)))} else {
    return(rbind(v,rep(-i,2)))
  }
  ftM2 = t(matrix(as.vector(sapply(FUN=fx,X=1:nrow(ftM1))),nrow=2))
  g = graph.data.frame(ftM2,directed=F)
  
  aa = hc$height; names(aa) = -(1:length(aa))
  col.x = sapply(FUN=function(s){return(paste(c("grey",s),collapse=""))},X=40:95) 
  E(g)$color = rep(col.x[floor((aa-min(aa))*55/max((aa-min(aa))))+1], each=2)
  E(g)$width = 2
  
  label.index = as.numeric(V(g)$name)
  label = V(g)$name
  label[label.index>0] = hc$label[label.index[label.index>0]]
  label[as.numeric(label.index)<0] = make.names(signif(aa[label[as.numeric(label.index)<0]],3))
  V(g)$label = label
  V(g)$size = nodesize; V(g)$size[igraph::degree(g)>1] = 0
  V(g)$label.family = "sans"
  #g$layout = layout_with_sugiyama(g)$layout
  
  V(g)$size[igraph::degree(g)==1]=nodesize
  V(g)$size[igraph::degree(g)>1]=0.1
  V(g)$label.cex=0.6
  #plot(g)
  return(g)
}
#------------------------------------------------------------

                              
f.TLNplot = function(TLN, ValueM, colorsat=1, colors=NA, title="NoTitle", markedNodes=NA, valueCap=NA, writeFile=F){
  if(all(is.na(colors))){
    colors=c('darkblue', 'mediumblue', 'dodgerblue', 'white', 'orange', 'red', 'darkred')
  	}
  
  if(is.vector(ValueM)){ValueM = t(as.matrix(ValueM))}
  neg.col = rev(colors[1:ceiling(length(colors)/2)])
  pos.col = colors[ceiling(length(colors)/2):length(colors)]
  bins = rep(100,(length(neg.col)-1))
  a1 = list()
  a2 = list()
  for(i in 1:length(bins)){
    a1[[i]] = t(colorRamp(neg.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
    a2[[i]] = t(colorRamp(pos.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
  }
  a1 = matrix(unlist(a1), ncol=3, byrow=T); a2 = matrix(unlist(a2), ncol=3, byrow=T)
  b1 = rgb(red=a1[,1], green=a1[,2], blue=a1[,3], max=255)
  b2 = rgb(red=a2[,1], green=a2[,2], blue=a2[,3], max=255)
  color.vector = c(rev(b1), b2, b2[length(b2)])
  nn = V(TLN)$name[igraph::degree(TLN)==1]
  ps = lapply(FUN=function(n, x){return(rep(1,n))}, X=1:ncol(ValueM), n=nrow(ValueM))
  if(!is.na(valueCap)){ValueM[ValueM<(-valueCap)] = -valueCap; ValueM[ValueM > valueCap] = valueCap}
  if(is.na(colorsat)){mm = max(abs(range(ValueM,na.rm=T)))} else {mm = colorsat}
  if(!is.na(colorsat) & max(abs(range(ValueM,na.rm=T)))>colorsat){print("WARNING: Color saturation too low")}

  pc = color.vector[round(ValueM[,nn]*((length(neg.col)-1)*100)/mm)+(((length(neg.col)-1)*100)+1)]
  pc = matrix(pc, nrow=nrow(ValueM))
  pc = lapply(FUN=function(i){return(rev(pc[,i]))}, X=1:ncol(pc))

  TLNx = TLN
  V(TLNx)$label[is.na(suppressWarnings(as.numeric(V(TLNx)$label)))] = ""
  V(TLNx)$shape[igraph::degree(TLNx)==1] = "pie"
  #if(nrow(ValueM)==1){V(TLNx)$shape[igraph::degree(TLNx)==1] = "circle"}
  V(TLNx)$shape[is.na(V(TLNx)$shape)] = "circle"
  V(TLNx)$pie[igraph::degree(TLNx)==1] = ps
  V(TLNx)$pie.color[igraph::degree(TLNx)==1] = pc
  if(nrow(ValueM)==1){V(TLNx)$color[igraph::degree(TLNx)==1] = unlist(pc)}
  coord.v = signif(seq(from=-0.5, to=0.5, length.out=8),3)
  cv = color.vector[round(seq(from=1, to=length(color.vector), length.out=7))]
  plotColorScale = function(){
    for(i in 1:length(cv)){rect(xleft=coord.v[i], xright=coord.v[i+1], ytop=-1.2, ybottom=-1.3, col=cv[i], border=NA)}}
  nv = signif(seq(from=-mm, to=mm, length.out=length(coord.v)),2)
  if(writeFile==F){plot(TLNx, main=title); plotColorScale(); text(x=coord.v, y=-1.34, labels=nv, adj=0.5)}
  if(writeFile==T){
    pdf(paste(c(title, "pdf"),collapse="."))
    plot(TLNx)
    plotColorScale()
    text(x=coord.v, y=-1.34, labels=nv, adj=0.5)
    title(title, cex.main=2)
    dev.off()
  }
  #return(TLNx)
}
#------------------------------------------------------------

f.TLNaddSpiderPlot = function(TLN, ValueM, minratio=0.25, sizefactor=2.7, minvalue.override=NULL, maxvalue.override=NULL, spokes=T){
  ## Will add spider plots to nodes, translating node names to ValueM colnames.
  if(is.null(colnames(ValueM))){colnames(ValueM) = as.character(1:ncol(ValueM))}
  ValueM[is.na(ValueM)] = min(ValueM, na.rm=T)
  
  coord.x = layout.norm(TLN$layout)
  angle.x = seq(from=0, to=(2*pi), length.out=nrow(ValueM)+1)[1:nrow(ValueM)]
  angle.x = angle.x+angle.x[2]/2
  
  if(!is.null(minvalue.override)){minv = minvalue.override} else {minv = min(ValueM)}
  if(!is.null(maxvalue.override)){maxv = maxvalue.override} else {maxv = max(ValueM)}
  
  ValueM.norm = (ValueM-minv)/(maxv-minv)*(1-minratio)+minratio
  points(x=coord.x[V(TLN)$size>1,1], y=coord.x[V(TLN)$size>1,2], cex=(V(TLN)$size[V(TLN)$size>1])*minratio/2.35, col="gray60")
  if(spokes==T){
    for(ii in 1:ncol(ValueM)){
      node.n=colnames(ValueM)[ii]
      radius.x = V(TLN)$size[V(TLN)$name==node.n]*sizefactor/diff(range(TLN$layout))
      coord.xx = c()
      for(i in 1:length(angle.x)){coord.xx[i] = coord.x[node.n,1]+sin(angle.x[i])*radius.x}
      coord.xy = c()
      for(i in 1:length(angle.x)){coord.xy[i] = coord.x[node.n,2]+cos(angle.x[i])*radius.x}
      for(i in 1:length(angle.x)){lines(x=c(coord.x[node.n,1],coord.xx[i]),y=c(coord.x[node.n,2],coord.xy[i]),col="gray60")}
    }
  }
  for(ii in 1:ncol(ValueM)){
    node.n=colnames(ValueM)[ii]
    radius.x = V(TLN)$size[V(TLN)$name==node.n]*sizefactor/diff(range(TLN$layout))
    coord.xx = c()
    for(i in 1:length(angle.x)){coord.xx[i] = coord.x[node.n,1]+sin(angle.x[i])*radius.x*ValueM.norm[i,node.n]}
    coord.xy = c()
    for(i in 1:length(angle.x)){coord.xy[i] = coord.x[node.n,2]+cos(angle.x[i])*radius.x*ValueM.norm[i,node.n]}
    lines(x=c(coord.xx, coord.xx[1]), y=c(coord.xy, coord.xy[1]))
  }
  text(paste(rownames(ValueM), collapse="\n"), x=-1,y=1.5, cex=0.7, adj=c(0,1))
}
#------------------------------------------------------------


f.CalculateModRatios = function(DEfiles){
  fa = function(DEfile){
    de.table = as.matrix(read.table(DEfile,row.names=1,header=T,sep="\t"))
    vx = names(Ann.Symbol); names(vx) = Ann.Symbol
    rownames(de.table) = vx[rownames(de.table)]
    de.up = rownames(de.table)[de.table[,'adj.P.Val']<0.05 & de.table[,'logFC']>log2(1.5)]             
    de.down = rownames(de.table)[de.table[,'adj.P.Val']<0.05 & de.table[,'logFC']<(-log2(1.5))]         # same
    #fb = function(mod){return((length(intersect(de.up,mod))-length(intersect(de.down,mod)))/length(mod))}
    fb = function(mod){return((length(intersect(de.up,mod))-length(intersect(de.down,mod)))/length(intersect(rownames(de.table),mod)))}
    return(as.matrix(sapply(FUN=fb, X=ITA.modules)))
  }
  ModRatios = t(sapply(FUN=fa, X=DEfiles))
  colnames(ModRatios) = as.character(1:ncol(ModRatios))
  return(ModRatios)
}
#----------------------------------------------------


### Calls ###
### 1. Read Annotation files and read functions.

## Annotation files
xx = read.table("hsapiens.ENTREZ.txt", sep="\t", header=F)
Ann.Entrez = as.character(unlist(xx[,2]))
names(Ann.Entrez) = as.vector(unlist(xx[,1]))

xx = read.table("hsapiens.SYMBOL.txt",sep="\t", header=F)
Ann.Symbol = as.character(unlist(xx[,2]))
names(Ann.Symbol) = as.vector(unlist(xx[,1]))


### 2. TopCorr Network

TopCorrGene = read.csv("TopCorrGene.csv", row.names=1)
TopCorrVal = as.matrix(read.csv("TopCorrVal.csv", row.names=1))

TopCorrVal[19305,1] = TopCorrVal[19305,2] + (TopCorrVal[19305,2]-TopCorrVal[19305,3])
TopCorrVal[19312,1] = TopCorrVal[19312,2] + (TopCorrVal[19312,2]-TopCorrVal[19312,3])

### 2.1. Network construction
K=5
ftM1 = cbind(rep(rownames(TopCorrGene), each=K), as.vector(t(TopCorrGene[,1:K])))
ftM1rank = rep(1:K, times=nrow(TopCorrGene))
ftM1tval = as.vector(t(TopCorrVal[,1:K]))
ftM1dist = as.vector(t(TopCorrVal[,1]/TopCorrVal[,1:K]))
ftM1distNorm = rep(mean(TopCorrVal[,1])/colMeans(TopCorrVal[,1:5]),nrow(TopCorrGene))
ftM2 = t(apply(FUN=sort, X=ftM1, MARGIN=1))
ftM3 = unique(cbind(ftM2,ftM1rank,ftM1tval, ftM1dist, ftM1distNorm))

TC.Graph = graph.data.frame(ftM3[,1:2],directed=F)
E(TC.Graph)$rank = as.numeric(ftM3[,3])
E(TC.Graph)$tval = as.numeric(ftM3[,4])
E(TC.Graph)$dist = as.numeric(ftM3[,5])
E(TC.Graph)$distNorm = as.numeric(ftM3[,6])
TC.Graph = simplify(TC.Graph,edge.attr.comb="min")
V(TC.Graph)$label = Ann.Symbol[V(TC.Graph)$name]
V(TC.Graph)$label.family = "sans"


#### 2.2. Edge weights per dataset
XX = readRDS("Blood_PBMC2.TC.Graph.DSweights.rds")
TC.Graph.DSweights = matrix(as.numeric(XX),ncol=ncol(XX))
colnames(TC.Graph.DSweights) = colnames(XX)
rownames(TC.Graph.DSweights) = rownames(XX)
rm(XX)


TC.Gr.Dist = distances(TC.Graph, weights=E(TC.Graph)$distNorm, 1:10)  # modif
#TC.Gr.Dist = readRDS("TopCorrDist.RDS")                              # check why it is not working
#TC.Gr.Dist.HC = hclust(as.dist(TC.Gr.Dist), "average")
TC.Gr.Dist.HC = readRDS("TopCorrHclust.RDS")

#TC.Gr.Dist.HC.cut = cutreeHybrid(dendro=TC.Gr.Dist.HC, distM=TC.Gr.Dist, 
#                                 deepSplit=3, minClusterSize = 100)

TC.Gr.Dist.HC.cut = readRDS("TC.Gr.Dist.HC.cut.rds") 


#table(TC.Gr.Dist.HC.cut$labels)
#ITA.modules = list()
#for(ii in 1:max(TC.Gr.Dist.HC.cut$labels)){
#  ITA.modules[[ii]] = TC.Gr.Dist.HC$labels[TC.Gr.Dist.HC.cut$labels==ii]}
#names(ITA.modules) = as.character(1:length(ITA.modules))
ITA.modules = readRDS("ITA.TLN.Modules.RDS")

## Distances between modules:
table(TC.Gr.Dist.HC.cut$labels) 

# Mod.Dist = matrix(ncol=max(TC.Gr.Dist.HC.cut$labels), nrow=max(TC.Gr.Dist.HC.cut$labels))
# rownames(Mod.Dist) = as.character(1:nrow(Mod.Dist))
# for(ii in 2:nrow(Mod.Dist)){
#   for(jj in 1:(ii-1)){
#     print(c(ii,jj))
#     Mod.Dist[ii,jj] = mean(TC.Gr.Dist[TC.Gr.Dist.HC.cut$labels==ii,TC.Gr.Dist.HC.cut$labels==jj])
#     Mod.Dist[jj,ii] = Mod.Dist[ii,jj]
#   }
# }

Mod.Dist = readRDS("Mod.Dist.rds")    
hist(Mod.Dist[lower.tri(Mod.Dist)],40)
Mod.Dist.HC = hclust(as.dist(Mod.Dist),"average")
plot(Mod.Dist.HC)

ITA.TLN = f.tree.leaf.network(Mod.Dist.HC,nodesize=11)
sv = V(ITA.TLN)$size; for(i in 1:length(ITA.modules)){sv[V(ITA.TLN)$name==i] = length(ITA.modules[[i]])}
V(ITA.TLN)$ModuleSize = sv
V(ITA.TLN)$size = sqrt(V(ITA.TLN)$ModuleSize)*0.95; V(ITA.TLN)$size[V(ITA.TLN)$size==min(V(ITA.TLN)$size)] = 0.1
#write.table(file="/home/ferall/Downloads/ITA-master/ITA_Blood_Package/ITA.TLN.Blood_PBMC2.ftM.txt", get.edgelist(ITA.TLN), sep="\t", col.names=F, row.names=F, quote=F)
#write.table(file="/home/ferall/Downloads/ITA-master/ITA_Blood_Package/ITA.TLN.Blood_PBMC2.nodes.txt", get.data.frame(ITA.TLN,"vertices"), sep="\t", col.names=T, row.names=F, quote=F)

A = readLines("ITA.TLN.Blood_PBMC2.ftM.txt.cyjs")
a1 = A[sapply(FUN=function(s){return(grepl(x=s, '\"x\"'))}, X=A)]
a1 = as.numeric(sapply(FUN=function(s){return(substring(s, first=15, last=(nchar(s)-1)))}, X=a1))

a2 = A[sapply(FUN=function(s){return(grepl(x=s, '\"y\"'))}, X=A)]
a2 = as.numeric(sapply(FUN=function(s){return(substring(s, first=15, last=(nchar(s)-1)))}, X=a2))

a3 = A[sapply(FUN=function(s){return(grepl(x=s, '\"name\"'))}, X=A)][2:(length(a1)+1)]
a3 = unname(sapply(FUN=function(s){return(substring(s, first=19, last=(nchar(s)-2)))}, X=a3))

aa = cbind(a1,-a2); rownames(aa) = a3
ITA.TLN$layout = aa[V(ITA.TLN)$name,]

xx = as.matrix(read.table("ModuleNames.txt",sep="\t",header=T)) # modif!!
vx = as.character(xx[,1]); names(vx) =  as.character(xx[,2])
V(ITA.TLN)$label = vx[V(ITA.TLN)$name]

plot(ITA.TLN)
#saveRDS(ITA.TLN, file="Blood_PBMC2_ITA.TLN.rds")
ITA.TLN = readRDS(file="Blood_PBMC2_ITA.TLN.rds")



### 3. GO term enrichment

ITA.modules.GOenrich = list(); length(ITA.modules.GOenrich) = length(ITA.modules)
for(i in 1:length(ITA.modules)){ITA.modules.GOenrich[[i]] = try(f.GOenrich(ITA.modules[[i]])); print(head(ITA.modules.GOenrich[[i]]))}
for(i in 1:length(ITA.modules)){write.table(ITA.modules.GOenrich[[i]], paste("/home/ferall/Downloads/ITA-master/ITA_Blood_Package/Modules.GOenrich/Module.GOenrich.", i, ".txt", sep=""), row.names=T, col.names=F, quote=F, sep="\t")}
ITA.modules.GOenrich = list()
for(i in 1:length(ITA.modules)){
  ITA.modules.GOenrich[[i]] = read.table(paste("/home/ferall/Downloads/ITA-master/ITA_Blood_Package/Modules.GOenrich/Module.GOenrich.", i, ".txt", sep=""),sep="\t", quote="", row.names=1)
}
sapply(FUN=function(df){return(unname(df[1:10,2]))}, X=ITA.modules.GOenrich[1:5])
#################################################

# read in the file(s) with DE genes
dir_files <- "/home/ferall/Desktop/MSc/sars/output_files/B_4mk_9tp_all_genes/"
gene_file_list = list.files(path=dir_files, pattern="_all.txt")  

mkSymb_ENSMM_file = read.table("mmulatta.SYMBOL.txt", sep="\t", header=F)

ENSMM_ENSG_file = read.table("mmulatta.hsapiens_orth.txt", sep="\t", header=F)

ENSG_hSymb_file = read.table("hsapiens.SYMBOL.txt", sep="\t", header=F)

######################################

translate_data <- function(x){
  DE_file <- read.table(file=paste(dir_files, x, sep=''), sep="\t", header=T)
  DE_symb_ENSMM <- merge(DE_file, mkSymb_ENSMM_file, by.x=c('genes'), by.y =c('V2'))
  DE_symb_ENSMM_ENSG <- merge(DE_symb_ENSMM, ENSMM_ENSG_file[,c('V1','V2')], by=c('V1'))
  DE_symb_ENSMM_ENSG_hsymb <- merge(DE_symb_ENSMM_ENSG[,c('V2','genes','logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')], ENSG_hSymb_file, by.x=c('V2'), by.y =c('V1'))
  df = subset(DE_symb_ENSMM_ENSG_hsymb, select = -c(genes,V2))
  df2 <- df[, c('V2.y','logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')]
  write.table(df2, file = paste("/home/ferall/Downloads/ITA-master/ITA_Blood_Package/hSymb_translated_files_all_genes/", x, "_hSymb",sep='') ,row.names = F, sep = "\t", quote = F)
}

lapply(gene_file_list, translate_data)
###########################


#### <----- 4. Analyzing data -----> ####

dir_tr_files <- "/home/ferall/Downloads/ITA-master/ITA_Blood_Package/hSymb_translated_files_all_genes/"
tr_file_list = list.files(path=dir_tr_files, pattern="_hSymb") 


# read in the file(s) with DE genes
setwd(grep("_output_files", list.dirs(dirname(getwd())), value=TRUE))
#getwd()

dir_files_DEgenes <- paste0(getwd(), '/', '2021.02.28_8mk_12tp', '/')
gene_file_list = list.files(path=dir_files_DEgenes, pattern=".txt")  

output_dir <- paste0(getwd(), '/', basename(dir_files_DEgenes), "_hSymb_transl_All_genes/")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}



# get ratios for graph
ratios <- f.CalculateModRatios(tr_file_list)

# TLN graph with piecharts in each node
f.TLNplot(ITA.TLN, ratios)


## (will work only with DE_D1 and DE_D5 (the others files didn't have enough rows))

#get top 5 and bottom 5 DE values:
bottom_5DE <- head(sort(ratios, decreasing=F), 5)
top_5DE <- head(sort(ratios, decreasing=T), 5)

##get index of the columns from the matrix where the top/bottom values are:
# gives indexes from a matrix-vector -> need to get the correct column values!
index_col_bottom_vect <- match(bottom_5DE, ratios)
index_col_top_vect <- match(top_5DE, ratios)   

### FUNCTION to correct the index
get_correct_index <- function(index_vect){
  index_correct <- c()
  for (i in 1: length(index_vect)){
    #print(i)
    index_correct[i] <- which(ratios==ratios[index_vect[i]], arr.ind = T)[,'col']
    #print(index_correct)
  }
  return(index_correct)
}


index_col_top <- get_correct_index(index_col_top_vect)
index_col_bottom <- get_correct_index(index_col_bottom_vect)

#get index of the columns from the matrix in the V(ITA.TLN)$name:
index_ITAname_top <- match(index_col_top, V(ITA.TLN)$name)
index_ITAname_bottom <- match(index_col_bottom, V(ITA.TLN)$name)

# get the label of the nodes for these indexes
node_label_top <- V(ITA.TLN)$label[index_ITAname_top]
node_label_bottom <- V(ITA.TLN)$label[index_ITAname_bottom]
cat("Top 5 upregulated nodes: ", node_label_top) 
cat("Top 5 downregulated nodes: ", node_label_bottom) 


# get the GO terms associated with each of the 5 up-/ down- regulated nodes
###FUNCTION get all the GO terms associated with the list of nodes
get_all_GO_terms <- function(node_label){
  GO_terms <- matrix(ncol=1)
  col_name <- c('GO_terms')
  colnames(GO_terms) <- col_name
  for (i in 1: length(node_label)){
    GOs <- f.GOenrich((as.numeric(node_label[i])))[,2, drop=FALSE]
    GO_terms <- rbind(GO_terms,'---',node_label[i], GOs)
  }
  return(GO_terms)
}


GO_terms_5upreg <- get_all_GO_terms(node_label_top)
print(GO_terms_5upreg)
GO_terms_5downreg <- get_all_GO_terms(node_label_bottom)
print(GO_terms_5downreg)
#######################################################

# spider plots
f.TLNaddSpiderPlot(ITA.TLN, ratios, minratio=0.25, sizefactor=1, minvalue.override=NULL, maxvalue.override=NULL, spokes=T)

