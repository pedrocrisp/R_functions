###

#### NOTE:
## This code does not renormalise the final summarised data to make Col-0 = 1 at time 0; was there a reason for this??? 
## Or perhaps I though that it did do this but actually it doesnt work?


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
# "qPCR_exp"; "OutPutFolder" = these variable determine the output file prefixes and location

#### summarise ####
qPCR_analysis_3 <- function(qPCR_data, ref_amplicon, ref_sample, ref_timepoint, Groups_to_plot, qPCR_exp, OutPutFolder){
  
 print("qPCR_data")
 print(str(qPCR_data))
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
 print("OutPutFolder")
 print(OutPutFolder)
 data_analysis = paste0(ref_amplicon, "_", ref_sample)
 print(data_analysis)
  
  #debugging
#  qPCR_data = DataRaw
#  ref_amplicon = "PP2A"
#  ref_sample = "within_genotype"
#  # ref_sample = "1"
#  ref_timepoint = "0" 
#  genotypes <- unique(DataRaw$Group)
#  Groups_to_plot = genotypes
#  qPCR_exp = "Exp423_qPCR66_75_76_xrn4_RRGS"
#  OutPutFolder = "qPCR_results"
  
 
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
      #using mean of the reference samples means that we dont have to re-norm to 1 later!
      target_genotype <- as.character(Data[i,"Group"])
      ref_value <- mean(Data[which(Data[,"Amplicon"]==target_amplicon &
                                     Data[,"Group"]==target_genotype &
                                     Data[,"Timepoint"]==reference_sample_timepoint)
                             ,"norm_abundance"])
    } else {
      
      #using mean of the reference samples means that we dont have to re-norm to 1 later!
      target_genotype <- unique(as.character(Data[which(Data[,"Sample"]==reference_sample),"Group"]))
      ref_value <- mean(Data[which(Data[,"Amplicon"]==target_amplicon &
                                     Data[,"Group"]==target_genotype &
                                     Data[,"Timepoint"]==reference_sample_timepoint)
                             ,"norm_abundance"])

      # old code that did not renormalise properly to Col-0:
#     ref_value <- Data[which(Data[,"Amplicon"]==target_amplicon &
#                               Data[,"Sample"]==reference_sample),
#                       "norm_abundance"]
    }
    norm_value2 <- target_value/ref_value
    Data[i,"norm_abundance2"] <- as.numeric(norm_value2)
  }

  
  
  Data$norm_abundance2 <- as.numeric(Data$norm_abundance2)
  
  write.csv(Data, paste0(OutPutFolder, "/", "Results_", qPCR_exp, "_", data_analysis, "_eachRep.csv"))
  
  
  #### r summarise by group ####
  
  GroupData <- ddply(Data, c("Group", "Timepoint", "Amplicon", "mean_PCR_eff"), summarise,
                     AverageAbundance = mean(norm_abundance2, na.rm=TRUE),
                     sd = sd(norm_abundance2, na.rm=TRUE),
                     n = sum(!is.na(norm_abundance2)),
                     se = sd/sqrt(n))
 
  
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
      ggplot(data=plot.data, aes(x=Timepoint, y=log2(AverageAbundance), color=Group, fill=Group
                                 #                            , group = factor(Genotype)
      )) +
        geom_point() +
        geom_errorbar(aes(ymin=log2(AverageAbundance-se), 
                          ymax=log2(AverageAbundance+se)), 
                      width=30, 
                      position=position_dodge(width = 0.90)) +
        geom_smooth(aes(stat="identity")
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
      ggplot(data=plot.data, aes(x=Timepoint, y=AverageAbundance, color=Group, fill=Group
                                 #                            , group = factor(Genotype)
      )) +
        geom_point() +
        geom_errorbar(aes(ymin=AverageAbundance-se, 
                          ymax=AverageAbundance+se), 
                      width=30, 
                      position=position_dodge(width = 0.90)) +
        geom_smooth(aes(stat="identity")
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
  
  write.csv(GroupData, paste0(OutPutFolder, "/", "Results_", qPCR_exp, "_", data_analysis, ".csv"))
  
  #try returning the dataframes that this function creates... comment out if this breakes the code...
  return(list(Data, GroupData))
  
}



  