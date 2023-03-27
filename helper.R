###LOAD LIBRARIES AND DATA, AND DEFINE CERTAIN HELPER FUNCTIONS


library("dada2")
library("devtools")
library("dplyr")
library("ggplot2")
library('ggrepel')
library('ggpubr')
library("microbiome")
library("phyloseq")
library("Rcpp")
library("reshape2")
library("tidyr")
library("vegan")
library("gtools")
library("DECIPHER")
library("ggvenn")
library("decontam")


load('RequiredDat.RData')

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


AbunExplained <- function(ps,taxList,plothead){
  otu.rel = abundances(ps,'compositional')
  plotdf = data.frame(sample_data(ps)[,'Sample_ID'],colSums(otu.rel[taxList,]))
  rownames(plotdf) = rownames(sample_data(ps))
  colnames(plotdf) = c('Sample_ID','Abundance_Explained')
  plotdf[rownames(plotdf),'Sample_ID'] = sample_data(ps)[rownames(plotdf),'Sample_ID']
  plotdf[rownames(plotdf),'Abundance_Explained'] = colSums(otu.rel[taxList,rownames(plotdf)])
  plotdf$is.neg = sample_data(ps)$is.neg
  
  ggplot(plotdf,aes(x=Sample_ID,y=Abundance_Explained,colour = is.neg,label = Sample_ID)) + geom_point() + geom_text_repel() + ggtitle(plothead) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
}

CountPlot <- function(ps,sampName,tax_new){
  otu.count = pseq@otu_table@.Data[sampName,]
  
  plotdf = data.frame(ReadCount = otu.count, ASVs = colnames(pseq@otu_table@.Data), names = tax_new[colnames(pseq@otu_table@.Data),"FinalAssign"])
  
  ggplot(plotdf,aes(x=ASVs,y=ReadCount,label=paste(ASVs,names))) + 
    geom_point() + geom_text_repel() + ggtitle(paste("Count levels in",sample_data(pseq)[sampName,'Sample_ID']))
  
}





scatPlotASVs <- function(pseq,taxList,plothead){
  otu.rel = abundances(ps,'compositional')
  PosSamps = rownames(sample_data(ps)[which(sample_data(ps)$is.neg == FALSE),])
  NegSamps = rownames(sample_data(ps)[which(sample_data(ps)$is.neg == TRUE),])
  abunPos = rowMeans(otu.rel[taxList,PosSamps])
  abunNeg = rowMeans(otu.rel[taxList,NegSamps])
  
  plotdf = data.frame(ASVs = taxList,AvgAbunPos = abunPos,AvgAbunNeg = abunNeg)
  
  ggplot(plotdf,aes(x = log(AvgAbunNeg),y = log(AvgAbunPos),label = ASVs)) + 
    geom_point() + geom_text_repel() + ggtitle(plothead) +
    xlab('Average Relative Abundance in Negative Controls') + 
    ylab('Average Relative Abundance in Positive Samples') 
}


bcDistDayPlotInc_comp <- function(pseq_A,Medium,Inoc){
  
  pseq_Beta = subset_samples(pseq_A,is.neg == FALSE)
  beta_plot <- transform_sample_counts(pseq_Beta, function(otu) otu/sum(otu))
  
  jsdDist = (phyloseq::distance(beta_plot,Type = "samples", "bray"))
  jsdDF <- melt(as.matrix(jsdDist))
  
  colnames(jsdDF) = c('Samp1','Samp2','value')
  
  jsdWide = spread(jsdDF, Samp2, value)
  
  rownames(jsdWide) = jsdWide$Samp1
  jsdWide$Samp1 = NULL
  
  sampInfo = sample_data(pseq_A)
  
  relSampDat = sampInfo[intersect(which(sampInfo$Medium == Medium),which(sampInfo$Community_Type == Inoc)),]
  
  relSampDat$Day = as.numeric(relSampDat$Day)
  
  #EqrelDat = relSampDat[which(relSampDat$Community_Type == 'Equal'),]
  #FecrelDat = relSampDat[which(relSampDat$Community_Type == 'Fecal'),]
  
  distdf = data.frame(seq(3, 27, 2),0,0)
  
  colnames(distdf) = c('Days','Distance_PC','Distance_DC')
  
  rownames(distdf) = as.character(distdf$Days)
  
  for (i in rownames(distdf)) {
    Dayno = as.numeric(i)
    
    dayi_PC = intersect(rownames(relSampDat)[which(relSampDat$Day == Dayno)],rownames(relSampDat)[which(relSampDat$Compartment == 'PC')])
    dayione_PC = intersect(rownames(relSampDat)[which(relSampDat$Day == (Dayno-2))],rownames(relSampDat)[which(relSampDat$Compartment == 'PC')])
    
    dayi_DC = intersect(rownames(relSampDat)[which(relSampDat$Day == Dayno)],rownames(relSampDat)[which(relSampDat$Compartment == 'DC')])
    dayione_DC = intersect(rownames(relSampDat)[which(relSampDat$Day == (Dayno-2))],rownames(relSampDat)[which(relSampDat$Compartment == 'DC')])
    
    distdf[i,"Distance_PC"] = jsdWide[dayi_PC,dayione_PC]
    distdf[i,"Distance_DC"] = jsdWide[dayi_DC,dayione_DC]
    
    #distdf[i,"SSD"] = sum((log(beta_plot@otu_table[dayi,]+ 1.e-06) - log(beta_plot@otu_table[dayione,]+ 1.e-06))**2)
    
  }
  #distdf = rbind(c(0,1041.026),distdf)
  
  plotobj = ggplot(distdf) + geom_point(aes(x = Days, y = Distance_PC,color = 'PC')) + geom_point(aes(x = Days, y = Distance_DC,color = 'DC')) + 
    geom_smooth(aes(x = Days, y = Distance_PC,color = 'PC'),se = FALSE) + geom_smooth(aes(x = Days, y = Distance_DC,color = 'DC'),se = FALSE) + 
    ggtitle(paste(Inoc,'Inoculum in',Medium,'Medium')) + ylim(0,1) + theme(aspect.ratio=1) + scale_x_discrete(limits = seq(3, 27, 2),labels = c('D1-D3','D3-D5','D5-D7','D7-D9','D9-D11','D11-D13','D13-D15','D15-D17','D17-D19','D19-D21','D21-D23','D23-D25','D25-D27'), guide = guide_axis(angle = 45)) + theme_classic() +
    scale_colour_manual(name = 'Compartment',values =c('PC'= PC_col,'DC'=DC_col), labels = c('PC','DC')) + ylab("Bray-Curtis dissimilarity") + xlab('Interval')
  
  return(plotobj)
}

function(pseq_A,Medium){
  
  pseq_Beta = subset_samples(pseq_A,is.neg == FALSE)
  beta_plot <- transform_sample_counts(pseq_Beta, function(otu) otu/sum(otu))
  jsdDist = (phyloseq::distance(beta_plot,Type = "samples", "bray"))
  jsdDF <- melt(as.matrix(jsdDist))
  
  colnames(jsdDF) = c('Samp1','Samp2','value')
  
  jsdWide = spread(jsdDF, Samp2, value)
  
  rownames(jsdWide) = jsdWide$Samp1
  jsdWide$Samp1 = NULL
  
  sampInfo = sample_data(pseq_A)
  
  relSampDat = sampInfo[which(sampInfo$Medium == Medium),]
  
  relSampDat$Day = as.numeric(relSampDat$Day)
  
  EqrelDat = relSampDat[which(relSampDat$Community_Type == 'Equal'),]
  FecrelDat = relSampDat[which(relSampDat$Community_Type == 'Fecal'),]
  
  distdf = data.frame(seq(1, 27, 2),0,0)
  
  colnames(distdf) = c('Days','Distance_PC','Distance_DC')
  
  rownames(distdf) = as.character(distdf$Days)
  
  for (i in rownames(distdf)) {
    Dayno = as.numeric(i)
    
    Eqdayi_PC = intersect(rownames(EqrelDat)[which(EqrelDat$Day == Dayno)],rownames(EqrelDat)[which(EqrelDat$Compartment == 'PC')])
    Fecdayi_PC = intersect(rownames(FecrelDat)[which(FecrelDat$Day == Dayno)],rownames(FecrelDat)[which(FecrelDat$Compartment == 'PC')])
    
    Eqdayi_DC = intersect(rownames(EqrelDat)[which(EqrelDat$Day == Dayno)],rownames(EqrelDat)[which(EqrelDat$Compartment == 'DC')])
    Fecdayi_DC = intersect(rownames(FecrelDat)[which(FecrelDat$Day == Dayno)],rownames(FecrelDat)[which(FecrelDat$Compartment == 'DC')])
    
    
    distdf[i,"Distance_PC"] = jsdWide[Eqdayi_PC,Fecdayi_PC]
    distdf[i,"Distance_DC"] = jsdWide[Eqdayi_DC,Fecdayi_DC]
    
  }
  
  
  plotobj = ggplot(distdf) + geom_point(aes(x = Days, y = Distance_PC,color = 'PC'), size = 0.75)+ geom_point(aes(x = Days, y = Distance_DC,color = 'DC'), size = 0.75) +
    geom_smooth(aes(x = Days, y = Distance_PC,color = 'PC'),se = FALSE, size = 0.5) + geom_smooth(aes(x = Days, y = Distance_DC,color = 'DC'),se = FALSE, size = 0.5) +
    ggtitle(paste(Medium,"Medium")) + geom_hline(yintercept = 0.6045086,color = 'grey', linetype = 'dashed',size = 0.2) + geom_hline(yintercept = 0.4450767,color = 'grey', linetype = 'dashed', size = 0.2) +
    geom_text(aes( 0, 0.6045086, label = 'Theoretical Inoculums', vjust = -1, hjust = -1), size = 2) + geom_text(aes( 0, 0.4450767, label = 'Observed Inoculums', vjust = -1, hjust = -1), size = 2) +
    ylim(0,0.75) + theme(aspect.ratio = 1) + theme_classic() + scale_x_discrete(limits = seq(1,27,2)) +
    scale_colour_manual(name = 'Compartment',values =c('PC'= PC_col,'DC'=DC_col), labels = c('PC','DC')) + ylab("Bray-Curtis dissimilarity") + xlab('Day') +
    theme(plot.title = element_text(size = 7, face = "bold"),axis.text=element_text(size=5),axis.title=element_text(size=6), legend.title = element_text(size = 6), legend.text = element_text(size = 4))
  return(plotobj)
}

ssdDistDayPlotComm_comp <- function(pseq_A,Medium){
  
  pseq_Beta = subset_samples(pseq_A,is.neg == FALSE)
  beta_plot <- transform_sample_counts(pseq_Beta, function(otu) otu/sum(otu))
  
  
  sampInfo = sample_data(pseq_A)
  
  relSampDat = sampInfo[which(sampInfo$Medium == Medium),]
  
  relSampDat$Day = as.numeric(relSampDat$Day)
  
  EqrelDat = relSampDat[which(relSampDat$Community_Type == 'Equal'),]
  FecrelDat = relSampDat[which(relSampDat$Community_Type == 'Fecal'),]
  
  distdf = data.frame(seq(1, 27, 2),0,0)
  
  colnames(distdf) = c('Days','SSD_PC','SSD_DC')
  
  rownames(distdf) = as.character(distdf$Days)
  
  for (i in rownames(distdf)) {
    Dayno = as.numeric(i)
    
    Eqdayi_PC = intersect(rownames(EqrelDat)[which(EqrelDat$Day == Dayno)],rownames(EqrelDat)[which(EqrelDat$Compartment == 'PC')])
    Fecdayi_PC = intersect(rownames(FecrelDat)[which(FecrelDat$Day == Dayno)],rownames(FecrelDat)[which(FecrelDat$Compartment == 'PC')])
    
    Eqdayi_DC = intersect(rownames(EqrelDat)[which(EqrelDat$Day == Dayno)],rownames(EqrelDat)[which(EqrelDat$Compartment == 'DC')])
    Fecdayi_DC = intersect(rownames(FecrelDat)[which(FecrelDat$Day == Dayno)],rownames(FecrelDat)[which(FecrelDat$Compartment == 'DC')])
    
    
    distdf[i,"SSD_PC"] = sum((log(beta_plot@otu_table[Eqdayi_PC,]+ 1.e-06) - log(beta_plot@otu_table[Fecdayi_PC,]+ 1.e-06))**2)
    distdf[i,"SSD_DC"] = sum((log(beta_plot@otu_table[Eqdayi_DC,]+ 1.e-06) - log(beta_plot@otu_table[Fecdayi_DC,]+ 1.e-06))**2)
    
  }
  
  #distdf = rbind(c(0,1041,1041),distdf)
  
  plotobj = ggplot(distdf) + geom_point(aes(x = Days, y = SSD_PC,color = 'PC'), size = 0.75)+ geom_point(aes(x = Days, y = SSD_DC,color = 'DC'), size = 0.75) +
    geom_smooth(aes(x = Days, y = SSD_PC,color = 'PC'),se = FALSE, size = 0.5) + geom_smooth(aes(x = Days, y = SSD_DC,color = 'DC'),se = FALSE, size = 0.5) +
    ggtitle(paste(Medium,"Medium")) + geom_hline(yintercept = 1318.4,color = 'grey', linetype = 'dashed', size = 0.2) +
    geom_text(aes( 0, 1302.4, label = 'Observed Inoculums', vjust = 1, hjust = -1), size = 3,check_overlap = T) +
    geom_hline(yintercept = 1302.4,color = 'grey', linetype = 'dashed', size = 0.2) +
    geom_text(aes( 0, 1318.4, label = 'Theoretical Inoculums', vjust = -1, hjust = -1), size = 3,check_overlap = T) +
    theme(aspect.ratio = 1) + theme_classic() + scale_x_discrete(limits = c(seq(1,27,2))) +
    scale_colour_manual(name = 'Compartment',values =c('PC'= PC_col,'DC'=DC_col), labels = c('PC','DC')) + ylab("SSD") + xlab('Day') +
    theme(plot.title = element_text(size = 14, face = "bold"),axis.text=element_text(size=9),axis.title=element_text(size=10), legend.title = element_text(size = 12), legend.text = element_text(size = 10))
  return(plotobj)
}

ssdDistDayPlotInc_comp <- function(pseq_A,Medium,Inoc){
  
  pseq_Beta = subset_samples(pseq_A,is.neg == FALSE)
  beta_plot <- transform_sample_counts(pseq_Beta, function(otu) otu/sum(otu))
  
  
  sampInfo = sample_data(pseq_A)
  
  relSampDat = sampInfo[intersect(which(sampInfo$Medium == Medium),which(sampInfo$Community_Type == Inoc)),]
  
  relSampDat$Day = as.numeric(relSampDat$Day)
  
  #EqrelDat = relSampDat[which(relSampDat$Community_Type == 'Equal'),]
  #FecrelDat = relSampDat[which(relSampDat$Community_Type == 'Fecal'),]
  
  distdf = data.frame(seq(3, 27, 2),0,0)
  
  colnames(distdf) = c('Days','SSD_PC','SSD_DC')
  
  rownames(distdf) = as.character(distdf$Days)
  
  for (i in rownames(distdf)) {
    Dayno = as.numeric(i)
    
    dayi_PC = intersect(rownames(relSampDat)[which(relSampDat$Day == Dayno)],rownames(relSampDat)[which(relSampDat$Compartment == 'PC')])
    dayione_PC = intersect(rownames(relSampDat)[which(relSampDat$Day == (Dayno-2))],rownames(relSampDat)[which(relSampDat$Compartment == 'PC')])
    
    dayi_DC = intersect(rownames(relSampDat)[which(relSampDat$Day == Dayno)],rownames(relSampDat)[which(relSampDat$Compartment == 'DC')])
    dayione_DC = intersect(rownames(relSampDat)[which(relSampDat$Day == (Dayno-2))],rownames(relSampDat)[which(relSampDat$Compartment == 'DC')])
    
    distdf[i,"SSD_PC"] = sum((log(beta_plot@otu_table[dayi_PC,]+ 1.e-06) - log(beta_plot@otu_table[dayione_PC,]+ 1.e-06))**2)
    distdf[i,"SSD_DC"] = sum((log(beta_plot@otu_table[dayi_DC,]+ 1.e-06) - log(beta_plot@otu_table[dayione_DC,]+ 1.e-06))**2)
    
    #distdf[i,"SSD"] = sum((log(beta_plot@otu_table[dayi,]+ 1.e-06) - log(beta_plot@otu_table[dayione,]+ 1.e-06))**2)
    
  }
  #distdf = rbind(c(0,1041.026),distdf)
  
  plotobj = ggplot(distdf) + geom_point(aes(x = Days, y = SSD_PC,color = 'PC')) + geom_point(aes(x = Days, y = SSD_DC,color = 'DC')) + 
    geom_smooth(aes(x = Days, y = SSD_PC,color = 'PC'),se = FALSE) + geom_smooth(aes(x = Days, y = SSD_DC,color = 'DC'),se = FALSE) + 
    ggtitle(paste(Inoc,'Inoculum in',Medium,'Medium')) + scale_x_discrete(limits = seq(3, 27, 2),labels = c('D1-D3','D3-D5','D5-D7','D7-D9','D9-D11','D11-D13','D13-D15','D15-D17','D17-D19','D19-D21','D21-D23','D23-D25','D25-D27'), guide = guide_axis(angle = 45)) + theme_classic() +
    scale_colour_manual(name = 'Compartment',values =c('PC'= PC_col,'DC'=DC_col), labels = c('PC','DC')) + ylab("SSD") + xlab('Day') + theme(aspect.ratio = 1)
  
  return(plotobj)
}
