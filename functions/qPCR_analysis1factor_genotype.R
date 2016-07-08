###
#qPCR data analysis function - Single factor - genotype
#Peter Crisp
#8/7/16
#Usage
#This function will summarise the data output from LinReg calculating means and SE and outputing the individual and group averages in a list of dataframes
#It expects the LinReg compact table with 3 extra columns at the begining
#Column A is "Sample" which is the list if Sample (individual replicate) names or numbers
#Column B is "Group" which is a replicate group eg Genotype

###NOTE - nedd to copy the strategy of the 2 factor code and output the dataframes.

#### summarise ####
qPCR_analysis <- function(qPCR_data, ref_amplicon, ref_sample, ref_group, qPCR_exp, data_analysis, OutPutFolder){
  
  #qPCR_data <- read.csv("qPCR#64 xrn4 light stress run 1.csv", skip=3)
  #ref_amplicon <- "PP2AA3_3'"
  #ref_sample <- "1"
  #ref_group <- "Col"
  #qPCR_exp <- "qPCR#64.1"
  #data_analysis <- "xrn4_lightstress"
  
  #Calculates mean, sd, n, and se on the value of N0 broken down by group, amplicon, sample and PCR_eff
  Data <- ddply(qPCR_data, c("Group", "Amplicon", "Sample", "mean_PCR_eff"), summarise,
                AverageAbundance = mean(N0, na.rm=TRUE),
                sd = sd(N0, na.rm=TRUE),
                n = sum(!is.na(N0)),
                se = sd/sqrt(n))
  
  #### ref gene normalise ####
  Data$norm_abundance <- "Null"
  
  reference_amplicon <- ref_amplicon
  # reference_amplicon <- "CYP5_plate2"
  # reference_amplicon <- "GAPC2_3'"
  
  for(i in 1:length(Data$Sample)){
    efficiency_value <- Data[i,"mean_PCR_eff"]
    target_value <- Data[i,"AverageAbundance"]
    target_sample <- as.character(Data[i,"Sample"])
    ref_value <- Data[which(Data[,"Sample"]==target_sample &
                              Data[,"Amplicon"]==reference_amplicon),"AverageAbundance"]
    norm_value <- target_value/ref_value
    Data[i,"norm_abundance"] <- as.numeric(norm_value)
  }
  
  Data$norm_abundance <- as.numeric(Data$norm_abundance)
  
  #### ref sample normalise} ####
  Data$norm_abundance2 <- "Null"
  
  ### User to define reference Sample
  reference_sample <- ref_sample
  
  for(i in 1:length(Data$Sample)){
    #i=1
    target_value <- Data[i,"norm_abundance"]
    target_amplicon <- as.character(Data[i,"Amplicon"])
    ref_value <- Data[which(Data[,"Amplicon"]==target_amplicon &
                              Data[,"Sample"]==reference_sample),"norm_abundance"]
    norm_value2 <- target_value/ref_value
    Data[i,"norm_abundance2"] <- as.numeric(norm_value2)
  }
  
  Data$norm_abundance2 <- as.numeric(Data$norm_abundance2)
  
  #### r summarise by group ####
  
  GroupData <- ddply(Data, c("Group", "Amplicon", "mean_PCR_eff"), summarise,
                     AverageAbundance = mean(norm_abundance2, na.rm=TRUE),
                     sd = sd(norm_abundance2, na.rm=TRUE),
                     n = sum(!is.na(norm_abundance2)),
                     se = sd/sqrt(n))
  #### r re-norm to control average ####
  GroupData$norm_abundance2 <- "Null"
  GroupData$se2 <- "Null"
  
  ### User to define reference Sample
  reference_sample <- ref_group
  
  for(i in 1:length(GroupData$AverageAbundance)){
    #i=1
    target_value <- GroupData[i,"AverageAbundance"]
    target_amplicon <- as.character(GroupData[i,"Amplicon"])
    ref_value <- GroupData[which(GroupData[,"Amplicon"]==target_amplicon &
                                   GroupData[,"Group"]==reference_sample),"AverageAbundance"]
    norm_value2 <- target_value/ref_value
    GroupData[i,"norm_abundance2"] <- as.numeric(norm_value2)
    
    target_value <- GroupData[i,"se"]
    target_amplicon <- as.character(GroupData[i,"Amplicon"])
    ref_value <- GroupData[which(GroupData[,"Amplicon"]==target_amplicon &
                                   GroupData[,"Group"]==reference_sample),"AverageAbundance"]
    norm_value2 <- target_value/ref_value
    GroupData[i,"se2"] <- as.numeric(norm_value2)
  }
  
  GroupData$norm_abundance2 <- as.numeric(GroupData$norm_abundance2)
  GroupData$se2 <- as.numeric(GroupData$se2)
  
  #### plot ####
  
  
  plotData <- GroupData
  #grey.pallet <- gray.colors(length(unique(plotData$Amplicon)), start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)
  plot.out <- ggplot(plotData, aes(x=Amplicon, y=norm_abundance2, fill=Group)) +
    geom_bar(stat="identity", 
             #fill="grey", 
             colour="black", 
             position="dodge") +
    geom_errorbar(aes(ymin=norm_abundance2-se2, 
                      ymax=norm_abundance2+se2), 
                  width=.2, 
                  position=position_dodge(width = 0.90)) +
    ylab("Relative abundance") +
    #scale_y_continuous(limits=c(0, 10)) +
    scale_fill_manual(values = c("white", "grey", "darkgrey", "black", "lightgrey", "brown")) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_text(size = rel(0.7)),
          axis.text.x = element_text(angle = 45, hjust = 1, size=7),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(),
          #         legend.position=c(0.85,0.8), 
          legend.text = element_text(size =5),
          legend.title = element_text(size =10))
  
  print(plot.out)
  
  pdf(paste0(OutPutFolder, "/", qPCR_exp, "_", data_analysis, ".pdf"), width=12/2.54, height=8/2.54)
  print(plot.out)
  dev.off()
  
  write.csv(plotData, paste0(OutPutFolder, "/", "Results_", qPCR_exp, "_", data_analysis, ".csv"))
  
  #### r stats, eval=FALSE ####
  ### anova and TukeyHSD
  
  for(gene in as.character(unique(Data$Amplicon))){
    #   gene="ABF3"
    data_stat <- Data[which(Data$Amplicon==gene),]
    a <- with(data_stat, aov(norm_abundance2~Group))
    summary(a)
    hsd <- TukeyHSD(a)
    print(gene)
    print(hsd)
    print(lapply(hsd, function(h) h[which(h[,4]<0.05),]))
    
  }
  
  #try returning the dataframes that this function creates... comment out if this breakes the code...
  #THIS CODE LINE IS UN-TESTED AND COULD BREAK THE CODE...
  return(list(Data, GroupData))
}