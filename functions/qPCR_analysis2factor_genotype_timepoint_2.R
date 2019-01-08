###
#qPCR data analysis function - genotype ~ timepoint
#Peter Crisp
#8/7/16
#Usage
#This function will summarise the data output from LinReg calculating means and SE and outputing the individual and group averages in a list of dataframes
#It expects the LinReg compact table with 3 extra columns at the begining
#Column A is "Sample" which is the list if Sample (individual replicate) names or numbers
#Column B is "Group" which is a replicate group eg Genotype
#Column C is "Timepoint" for a time course, wlthough it might also be a different second variable, although it must still be named "Timepoint"

# "qPCR_data" =input data table
# "ref_amplicon" = housekeeping gene
# "ref_sample" = sample to normalise data to from column A of input table. 
  # Every biological replicate should have individual sample number. 
  # eg Col-0, untreated rep1 generally is sample "1" 
  # If within genotype normalisation is desired then ref_sample = "within_genotype"
# "ref_timepoint"  = The timepoint the reference sample belongs too. 
  # This variable is used to re-normalised reference point to zero after averaging the biological reps
# "Groups_to_plot" = which genotypes to plot. This final results table is filtered for only genotypes listed in this variable.
# "qPCR_exp"; "data_analysis"; "OutPutFolder" = these variable determine the output file prefixes and location

#### summarise ####
qPCR_analysis_2 <- function(qPCR_data, ref_amplicon, ref_sample, ref_timepoint, Groups_to_plot, qPCR_exp, data_analysis, OutPutFolder){
  
 print("qPCR_data")
 print(qPCR_data)
 print("ref_amplicon")
 print(ref_amplicon)
 print("ref_sample")
 print(ref_sample)
 print("ref_timepoint")
 print(ref_timepoint)
 print("Groups_to_plot")
 print(Groups_to_plot)
 print("qPCR_exp")
 print(qPCR_exp)
 print("data_analysis")
 print(data_analysis)
 print("OutPutFolder")
 print(OutPutFolder)
  
  #debugging
 qPCR_data = DataRaw
 ref_amplicon = "PP2A"
 ref_sample = "within_genotype"
 # ref_sample = "1"
 ref_timepoint = "0" 
 genotypes <- unique(DataRaw$Group)
 Groups_to_plot = genotypes
 qPCR_exp = "Exp423"
 data_analysis = "qPCR66_67_68_xrn4_RRGS_PP2A"
 OutPutFolder = "qPCR_results"
  
 
  #### summarise technical reps ####
 
  #qPCR_data$Timepoint <-  as.factor(qPCR_data$Timepoint)
  #Calculates mean, sd, n, and se on the value of N0 broken down by group, amplicon, sample and PCR_eff
  Data <- ddply(qPCR_data, c("Group", "Timepoint", "Amplicon", "Sample", "mean_PCR_eff"), summarise,
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
    # i = 2
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
  
  ### User to define reference Sample
  reference_sample_timepoint <- ref_timepoint
  
  #for single genotype reference (eg Col-0) calculate normalisation values
  for(i in 1:length(Data$Sample)){
    #i=111
    #Data[i,]
    target_value <- Data[i,"norm_abundance"]
    target_amplicon <- as.character(Data[i,"Amplicon"])
    # apply different reference values depending if within genotype normalisation or ref Col-0
    if(reference_sample=="within_genotype"){
      ref_value <- mean(Data[which(Data[,"Amplicon"]==target_amplicon &
                                     Data[,"Group"]==target_genotype &
                                     Data[,"Timepoint"]==reference_sample_timepoint)
                             ,"norm_abundance"])
    } else {
    ref_value <- Data[which(Data[,"Amplicon"]==target_amplicon &
                              Data[,"Sample"]==reference_sample),"norm_abundance"]
    }
    norm_value2 <- target_value/ref_value
    Data[i,"norm_abundance2"] <- as.numeric(norm_value2)
  }
  
  #for within genotype reference calculate normalisation values
  
  for(i in 1:length(Data$Sample)){
    #i=1
    target_value <- Data[i,"norm_abundance"]
    target_amplicon <- as.character(Data[i,"Amplicon"])
    target_genotype <- as.character(Data[i,"Group"])
    ref_value <- mean(Data[which(Data[,"Amplicon"]==target_amplicon &
                              Data[,"Group"]==target_genotype &
                              Data[,"Timepoint"]==reference_sample_timepoint)
                              ,"norm_abundance"])
    norm_value2 <- target_value/ref_value
    Data[i,"norm_abundance2"] <- as.numeric(norm_value2)
  }
  
  
  Data$norm_abundance2 <- as.numeric(Data$norm_abundance2)
  
  #### r summarise by group ####
  
  GroupData <- ddply(Data, c("Group", "Timepoint", "Amplicon", "mean_PCR_eff"), summarise,
                     AverageAbundance = mean(norm_abundance2, na.rm=TRUE),
                     sd = sd(norm_abundance2, na.rm=TRUE),
                     n = sum(!is.na(norm_abundance2)),
                     se = sd/sqrt(n))
 
  
   #### r re-norm to control average ####
  GroupData$norm_abundance2 <- "Null"
  GroupData$se2 <- "Null"
  
  
  for(i in 1:length(GroupData$AverageAbundance)){
    #i=27
    target_value <- GroupData[i,"AverageAbundance"]
    target_amplicon <- as.character(GroupData[i,"Amplicon"])
    
    if(reference_sample=="within_genotype"){
      #add code here...
      ref_value <- GroupData[which(GroupData[,"Amplicon"]==target_amplicon &
                                     GroupData[,"Group"]==reference_sample_group1 & 
                                     GroupData[,"Timepoint"]==reference_sample_timepoint),"AverageAbundance"]
      norm_value2 <- target_value/ref_value
      GroupData[i,"norm_abundance2"] <- as.numeric(norm_value2)
      
      #normalised abundance se
      target_value <- GroupData[i,"se"]
      target_amplicon <- as.character(GroupData[i,"Amplicon"])
      ref_value <- GroupData[which(GroupData[,"Amplicon"]==target_amplicon &
                                     GroupData[,"Group"]==reference_sample_group1 & 
                                     GroupData[,"Timepoint"]==reference_sample_timepoint),"AverageAbundance"]
    } else {
    #normalised abundance
      reference_sample_group1 <- unique(as.character(Data[which(Data[,"Sample"]==reference_sample),"Group"]))
    ref_value <- GroupData[which(GroupData[,"Amplicon"]==target_amplicon &
                                   GroupData[,"Group"]==reference_sample_group1 & GroupData[,"Timepoint"]==reference_sample_timepoint),"AverageAbundance"]
    norm_value2 <- target_value/ref_value
    GroupData[i,"norm_abundance2"] <- as.numeric(norm_value2)
    
    #normalised abundance se
    target_value <- GroupData[i,"se"]
    target_amplicon <- as.character(GroupData[i,"Amplicon"])
    ref_value <- GroupData[which(GroupData[,"Amplicon"]==target_amplicon &
                                   GroupData[,"Group"]==reference_sample_group1 & GroupData[,"Timepoint"]==reference_sample_timepoint),"AverageAbundance"]
    }
    
    norm_value2 <- target_value/ref_value
    GroupData[i,"se2"] <- as.numeric(norm_value2)
  }
  
  GroupData$norm_abundance2 <- as.numeric(GroupData$norm_abundance2)
  GroupData$se2 <- as.numeric(GroupData$se2)
  
  ###time course plot###
  
  # unique(qPCRdata.summary$Genotype)[-(2:3)]
  
  plot_colours = plot.colours <- c("#4E8E9B",
                                   "#8D66D3",
                                   "#6CAA37",
                                   "#C74C32",
                                   "#5D5626",
                                   #                                  "#6A4873",
                                   #                                  "#58A36F",
                                   "#C48F3B")
  # "#C0546E",
  # "#D150B7",
  # "#8A8FCA")
  
  pdf(paste0(OutPutFolder, "/", qPCR_exp, "_", data_analysis, "_log2.pdf"),
      height=3,
      width=4)
  for(i in unique(GroupData$Amplicon)){
    # i="Col"
    geno=factor(Groups_to_plot)
    
    plot.data <- GroupData[GroupData$Amplicon==i 
                           & GroupData$Group %in% geno
                           , ]
    print(
      ggplot(data=plot.data, aes(x=Timepoint, y=log2(norm_abundance2), color=Group, fill=Group,
                                 #                            , group = factor(Genotype)
      )) +
        geom_point() +
        geom_errorbar(aes(ymin=log2(norm_abundance2-se2), 
                          ymax=log2(norm_abundance2+se2)), 
                      width=30, 
                      position=position_dodge(width = 0.90)) +
        geom_smooth(aes(stat="identity"),
                    #               method = "lm", formula = y ~ poly(x, 4)
                    #               method =family = binomial, formula = y ~ poly(x,2)
        ) +
        #   geom_line() +
        #   scale_x_discrete(breaks=c(0,30,60,75,90,120,125), labels=c(0,30,60,75,90,120, 125)) +
        ylab(paste0(i," Relative abundance")) +
        geom_vline(xintercept = 65, linetype = "dotted") + 
        theme_bw() +
        scale_fill_manual(values = plot_colours) +
        scale_colour_manual(values = plot_colours) +
        theme(text = element_text(size=10))
    )
  }
  
  dev.off()
  
  pdf(paste0(OutPutFolder, "/", qPCR_exp, "_", data_analysis, ".pdf"),
      height=3,
      width=4)
  for(i in unique(GroupData$Amplicon)){
    # i="Col"
    geno=factor(Groups_to_plot)
    
    plot.data <- GroupData[GroupData$Amplicon==i 
                           & GroupData$Group %in% geno
                           , ]
    print(
      ggplot(data=plot.data, aes(x=Timepoint, y=norm_abundance2, color=Group, fill=Group,
                                 #                            , group = factor(Genotype)
      )) +
        geom_point() +
        geom_errorbar(aes(ymin=norm_abundance2-se2, 
                          ymax=norm_abundance2+se2), 
                      width=30, 
                      position=position_dodge(width = 0.90)) +
        geom_smooth(aes(stat="identity"),
                    #               method = "lm", formula = y ~ poly(x, 4)
                    #               method =family = binomial, formula = y ~ poly(x,2)
        ) +
        #   geom_line() +
        #   scale_x_discrete(breaks=c(0,30,60,75,90,120,125), labels=c(0,30,60,75,90,120, 125)) +
        ylab(paste0(i," Relative abundance")) +
        geom_vline(xintercept = 65, linetype = "dotted") + 
        theme_bw() +
        scale_fill_manual(values = plot_colours) +
        scale_colour_manual(values = plot_colours) +
        theme(text = element_text(size=10))
    )
  }
  
  dev.off()
  
  #try returning the dataframes that this function creates... comment out if this breakes the code...
  return(list(Data, GroupData))
  
  #### plot ####
  
  
  # plotData <- GroupData
  # #grey.pallet <- gray.colors(length(unique(plotData$Amplicon)), start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)
  # plot.out <- ggplot(plotData, aes(x=Group, y=norm_abundance2, fill=Group2)) +
  #   geom_bar(stat="identity", 
  #            #fill="grey", 
  #            colour="black", 
  #            position="dodge") +
  #   geom_errorbar(aes(ymin=norm_abundance2-se2, 
  #                     ymax=norm_abundance2+se2), 
  #                 width=.2, 
  #                 position=position_dodge(width = 0.90)) +
  #   ylab("Relative abundance") +
  #   #scale_y_continuous(limits=c(0, 10)) +
  #   scale_fill_manual(values = c("white", "grey", "darkgrey", "black", "lightgrey", "brown")) +
  #   theme(axis.title.x=element_blank(),
  #         axis.title.y=element_text(size = rel(0.7)),
  #         axis.text.x = element_text(angle = 45, hjust = 1, size=7),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.border=element_blank(),
  #         axis.line=element_line(),
  # #         legend.position=c(0.85,0.8), 
  #         legend.text = element_text(size =5),
  #         legend.title = element_text(size =10))
  # 
  # print(plot.out)
  # 
  # pdf(paste0(qPCR_exp, "_", data_analysis, ".pdf"), width=12/2.54, height=8/2.54)
  # print(plot.out)
  # dev.off()
  
  #write.csv(plotData, paste0("Results_", qPCR_exp, "_", data_analysis, ".csv"))
  
  #### r stats, eval=FALSE ####
  ### anova and TukeyHSD
  
  # figure out how to do stats on time course later
  # for(gene in as.character(unique(Data$Amplicon))){
  # #   gene="ABF3"
  #   data_stat <- Data[which(Data$Amplicon==gene),]
  # a <- with(data_stat, aov(norm_abundance2~Group))
  # summary(a)
  # hsd <- TukeyHSD(a)
  # print(gene)
  # print(hsd)
  # print(lapply(hsd, function(h) h[which(h[,4]<0.05),]))
  
}
#}