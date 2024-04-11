#install packages once
install.packages("dplyr")
install.packages("tidyr")
install.packages("data.table")
install.packages("lessR")
install.packages("foreach")
install.packages("doParallel")
install.packages("doSNOW")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("reshape2")
install.packages("tidyselect")
install.packages("tidyverse")
install.packages("devtools")
install.packages("githubinstall")
install.packages("survival")
install.packages('lubridate')
install.packages('ggsurvfit')
install.packages('gtsummary')
install.packages('tidycmprsk')
install.packages('Rtools')
install.packages('condSURV')
install.packages('ggplot2')
install.packages('cutpointr')
install.packages('survminer')
install.packages('ROCR')
install.packages('qvalue')
install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")
#library packages
library(dplyr)
library(tidyr)
library(data.table)
library(lessR)
library(foreach)
library(doParallel)
library(doSNOW)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(tidyselect)
library(tidyverse)
library(devtools)
library(githubinstall)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condSURV)
library(ggplot2)
library(cutpointr)
library(survminer)
library(ROCR)
library(qvalue)

#Set Directory
setwd("C:/Users/soldi/Documents/R/CTRP2")

#Selecting LIHC dataset, can pick other TCGA cancers if you swap it
TCGA_Cancer_ID <- "lihc"
setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018"))

# Load the data
PtClinical <- read.table(file = "data_clinical_patient.txt", sep = "\t",quote = "", header = TRUE)
PtMutation <- read.table(file = "data_mutations.txt", sep = "\t",quote = "", header = TRUE)
PtExpression <- read.table(file = "data_mrna_seq_v2_rsem.txt", sep = "\t",quote = "", header = TRUE)
PtExpressionZ <- read.table(file = "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", sep = "\t",quote = "", header = TRUE)


#LIHC E2F1 Figure

j <- grep("E2F1", PtExpressionZ$Hugo_Symbol)

#create function to calculate and generate LOCC and AUC graphs
Loccfunc <- function(j){
  
  Mutation1 <- PtExpressionZ[j,1]
  Gene_set <- Mutation1
  #Pick Gene Mutation Data
  if (all(Mutation1 != "")) {
    PtSpecificMutation <- subset(PtMutation, select = c(1,9:20,39))
    PtSpecificMutation <- filter(PtSpecificMutation, Hugo_Symbol == Mutation1)
    PtSpecificMutation <- filter(PtSpecificMutation, Variant_Classification != "Silent")
    PtSpecificMutation <-  unique(PtSpecificMutation[ , c("Hugo_Symbol", "Tumor_Sample_Barcode")], .keep_all = TRUE)
    PtSpecificMutation$Tumor_Sample_Barcode <- substr(PtSpecificMutation$Tumor_Sample_Barcode,1,nchar(PtSpecificMutation$Tumor_Sample_Barcode)-3)
    colnames(PtSpecificMutation)[2] <- "PATIENT_ID"
    
    # Merge the clinical and mutation data
    PtClinical2 <- merge(PtClinical, PtSpecificMutation, by = "PATIENT_ID", all.x = TRUE)
    PtClinical2$Hugo_Symbol <- replace_na(PtClinical2$Hugo_Symbol, "WT")
    PtClinical2$Mutant <- PtClinical2$Hugo_Symbol
    
    # Process the expression data
    PtExpressionGene <- PtExpressionZ[PtExpressionZ$Hugo_Symbol %in% Mutation1,]
    PtExpressionGene <- subset(PtExpressionGene, select = c(-2))
    PtExpressionGene <- na.omit(PtExpressionGene)
    Averages <- summarize_all(PtExpressionGene[,-c(1)], mean)
    Averages <- data.frame(0, Averages)
    colnames(Averages)[1] <- "Hugo_Symbol"
    PtExpressionGene <- rbind(PtExpressionGene, Averages)
    PtExpressionGene[nrow(PtExpressionGene),1] <- "Mean"
    
    # Compute Gene Z-Score
    PtGeneScore <- PtExpressionGene[nrow(PtExpressionGene),]
    PtGeneScore <- rbind(colnames(PtGeneScore), PtGeneScore)
    PtGeneScore2 <- t(PtGeneScore)
    PtGeneScore2 <- PtGeneScore2[2:nrow(PtGeneScore2),]
    colnames(PtGeneScore2) <- c("PATIENT_ID", "Gene_Score")
    PtGeneScore2[,1] <- substr(PtGeneScore2[,1],1,nchar(PtGeneScore2[,1])-3)
    PtGeneScore2[,1] <- str_replace(PtGeneScore2[,1], "\\.", "-")
    PtGeneScore2[,1] <- str_replace(PtGeneScore2[,1], "\\.", "-")
  
  
  
  # Set the initial cut-off point
  cutter <- 22
  
  # Merge clinical data with gene score data
  PtClinical4 <- merge(PtClinical, PtGeneScore2, by = "PATIENT_ID", all = TRUE)  
  
  # Select relevant columns
  PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
  
  # Further subset the data
  PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
  
  # Rename columns
  colnames(PtSurv) <- c("time", "status", "score")
  
  # Recode 'status' to numerical values and remove NA rows
  PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
  PtSurv <- na.omit(PtSurv)
  
  # Convert 'score' to numeric and replace NA values with 0
  PtSurv$score <- as.numeric(PtSurv$score)
  PtSurv[is.na(PtSurv)] <- 0
  
  # Order by 'score'
  PtSurv <- PtSurv[order(as.numeric(PtSurv$score), decreasing = TRUE), ]
  
  # Save the gene score value at the cut-off point
  GeneScoreValue <- PtSurv[cutter,3]
  
  # Assign high and low gene score groups
  PtSurv[1:cutter,3] <- "High_Gene"
  PtSurv[cutter:(nrow(PtSurv)+1),3] <- "Low_Gene"
  
  # Filter out only the high and low gene score groups
  PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
  
  # Calculate survival difference using the survdiff function
  CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
  
  # Create a dataframe with the results
  CutTrial <- data.frame(CutTrial$pvalue, CutTrial$n, CutTrial$obs, CutTrial$exp, GeneScoreValue)
  
  # Calculate HR
  HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
  CutTrial$HR <- HR
  
  # Initialize AggreCut dataframe with the results from the first cut-off point
  AggreCut <- CutTrial
  
  # Initialize a data frame to store the results
  AggreCut <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c(colnames(CutTrial)))))
  # Initialize a variable to track duplicates
  Dupe <- 0
  
  # Loop over different cut-off points
  for (cutter in 2:(nrow(PtSurv)-1)){  
    # Repeat the same process with different cut-off points
    PtClinical4 <- merge(PtClinical, PtGeneScore2,  by = "PATIENT_ID", all = TRUE)  
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    colnames(PtSurv) <- c("time", "status", "score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv[is.na(PtSurv)] <- 0
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    GeneScoreValue <- PtSurv[cutter,3]
    PtSurv[1:cutter,3] <- "High_Gene"
    PtSurv[cutter:(nrow(PtSurv)),3] <- "Low_Gene"
    PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
    
    # Calculate survival difference for the current cut-off point
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    CutTrial <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, GeneScoreValue)
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
    
    # Add the results to the AggreCut data frame
    AggreCut <- rbind(AggreCut, CutTrial)
  }
  
  # Now we process the AggreCut data frame
  AggreCut <- AggreCut[order(as.numeric(AggreCut[,1]), decreasing = FALSE),]
  AggreCut2 <- AggreCut[AggreCut$groups == "score=High_Gene",]
  AggreCut2 <- AggreCut2[order(as.numeric(AggreCut2[,3]), decreasing = FALSE),]
  AggreCut2$Fraction <- AggreCut2$Freq/nrow(AggreCut2)
  AggreCut2$logp <- -log(AggreCut2$CutTrial.pvalue, base = 10)
  
  #Check for Duplicates
  if (AggreCut2[nrow(AggreCut2), 6] == AggreCut2[nrow(AggreCut2)-1, 6] ){
    Dupe = 1 
  }
  
  
  #Count Dupes
  if (Dupe == 1){
    PreV  <- AggreCut2[1,6]
    AggreCut2[1,10] <- Dupe
    for (n in 2:nrow(AggreCut2)){
      NowV <- AggreCut2[n,6]
      if (PreV==NowV){
        Dupe = Dupe +  1
      } else {
        Dupe = 1
      }
      PreV <- NowV
      AggreCut2[n,10] <- Dupe
    }
  }
  
  DupeStart <- nrow(AggreCut2)-Dupe
  
  #Dupe Penalty
  for (n in DupeStart:nrow(AggreCut2)){
    PreVHR <- AggreCut2[n-1,7] - 0.1
    PreVHR <- max(PreVHR, 1)
    PreVP <- AggreCut2[n-1,9] - 0.1
    PreVP <- max(PreVP, 0)
    AggreCut2[n,7] <- PreVHR
    AggreCut2[n,9] <- PreVP
  }
  
  
  # output data
  write.csv(AggreCut2, file = paste0("AggreCut", Gene_set,cutter, ".csv"))
  
  #Plot LOCC Graph
  {
    df.cut <- subset(AggreCut2, select = c(7,8))
    df.cut <- melt(df.cut, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2 <- subset(AggreCut2, select = c(8,9))
    df.cut2$logp <- (df.cut2$logp/2)
    df.cut2 <- melt(df.cut2, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2$variable <- "-Log (p value)"
    df.cut2 <- rbind(df.cut2, df.cut)
    df.cut2 <- df.cut2[df.cut2$Fraction > 0.05,]
    df.cut2 <- df.cut2[df.cut2$Fraction < 0.95,]
    df.cut2$variable <- factor(df.cut2$variable, levels = c("HR", "-Log (p value)"))
    p <- ggplot(df.cut2, aes(x = Fraction, y = value,color = variable)) +geom_point() + geom_line(aes(color = variable)) + scale_color_manual(values=c("black", "orange")) + ggtitle(paste0(Gene_set, " TCGA ", TCGA_Cancer_ID ," LOCC")) + scale_y_continuous(
      # Features of the first axis
      name = "Hazard Ratio (HR)",
      # Add a second axis and specify its features
      sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
    ) + scale_x_continuous( name = "Fraction in High Activity Group", limits = c(0, 1), expand = c(0, 0)) + geom_hline(aes(yintercept=1), color = "red") + geom_hline(aes(yintercept=0.65),  color = "green") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30))  + labs(colour="Legend",x="xxx",y="yyy")
    
    
    
    ggsave(paste0(Gene_set,"_LOCC.pdf"), plot=p,width = 7, height = 7, dpi = 300,)
    # dev.off
  }
  
  PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
  PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
  colnames(PtSurv) <- c("time", "status", "score")
  PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
  PtSurv <- na.omit(PtSurv)
  PtSurv$score <- as.numeric(PtSurv$score)
  PtSurv$time <- as.numeric(PtSurv$time)
  PtSurv[is.na(PtSurv)] <- 0
  # PtSurv <- LiziSurv2
  PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
  
  {
    #Check if gene is expressed
    if (nrow(PtSurv) > 100 & PtSurv[nrow(PtSurv)/5,3] != PtSurv[nrow(PtSurv)-1,3]){
    
      #Mutation Analysis
      {
        # Subset the dataframe `PtClinical2` to keep only relevant columns
        PtSurvivalCurve <- subset(PtClinical2, select = c(1,30:37,38))
        
        # Further subset `PtSurvivalCurve` to keep only essential columns
        PtSurv <- subset(PtSurvivalCurve, select = c(1,3,2,10))
        
        # Rename the columns of `PtSurv`
        colnames(PtSurv) <- c("ID","time", "status", "mutation")
        
        # Recode the `status` column such that '1:DECEASED' is 1 and '0:LIVING' is 0
        PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
        
        # Remove NA values from `PtSurv`
        PtSurv <- na.omit(PtSurv)
        
        # Replace the mutations with "Mutant"
        PtSurv$mutation[PtSurv$mutation == Mutation1] <- "Mutant"
        
        
        # Order the `PtSurv` dataframe based on `mutation`
        PtSurv <- PtSurv[order(as.character(PtSurv$mutation), decreasing = FALSE), ]
        
        # Keep only distinct rows of `PtSurv` based on `ID`
        PtSurv <- distinct(PtSurv, ID, .keep_all = T)
        
        # Only perform survival analysis if there are more than 5 "Mutant"
        if (length(PtSurv$mutation[PtSurv$mutation == "Mutant"]) > 5){
          # Compute the survival difference
          surv_diff <- survdiff(Surv(time, Status) ~ mutation, data = PtSurv) 
          
          # Compute the survival fit
          surv_fit <- survfit(Surv(time, Status) ~ mutation, data = PtSurv)
          
          # Plot the survival analysis
          r <- survfit2(Surv(time, Status) ~ mutation, data = PtSurv) 
          t <-  ggsurvplot(r, data = PtSurv, risk.table = TRUE, size=1.2, fontsize = 7.5, 
                           risk.table.y.text = FALSE, tables.y.text = FALSE, risk.table.height = 0.35,
                           font.title = c(16, "bold", "black"), font.subtitle = c(16, "bold", "black"),
                           font.caption = c(16, "bold", "black"), font.x = c(16, "bold", "black"),
                           font.y = c(16, "bold", "black"), font.tickslab = c(16, "bold", "black"))
          t$plot <- t$plot + theme(legend.title = element_text(size = 18, color = "black", face = "bold"),
                                   legend.text = element_text(size = 18, color = "black", face = "bold"),
                                   axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                   axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                                   axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                                   axis.title.y = element_text(size = 20, color = "black", face = "bold"))
          t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                                     axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                     axis.title.x = element_text(size = 20, color = "black", face = "bold"))
          
          pdf(paste0(Gene_set,"_mutsurv.pdf"))
          print(t, newpage = FALSE)
          dev.off()
        }
      }
      
      
      
      
    #Cut at best cutpoint
    {
      
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      # PtSurv <- LiziSurv2
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      res.cut <- surv_cutpoint(PtSurv, time = "time", event = "Status",
                               variables = "score", minprop = 0.1)
      summary(res.cut)
      PtSurv$score[as.numeric(PtSurv$score) > res.cut$score[4]] <- as.character("High_Score")
      PtSurv$score[(PtSurv$score) <= res.cut$score[4]] <- as.character("Low_Score")
      
    }
    # Print a summary of `res.cut`
    summary(res.cut)
    
    # Open a new device
    dev.new(width = 200, height = 200, unit = "px")
    
    # Plot the survival curve with risk table
    r <- survfit2(Surv(time, Status) ~ score, data = PtSurv) 
    t <- ggsurvplot(r, data = PtSurv, risk.table = TRUE,
                    size = 1.2,
                    fontsize = 7.5,
                    risk.table.y.text = FALSE,
                    tables.y.text = FALSE,
                    risk.table.height = 0.35,
                    font.title = c(16, "bold", "black"),
                    font.subtitle = c(16, "bold", "black"),
                    font.caption = c(16, "bold", "black"),
                    font.x = c(16, "bold", "black"),
                    font.y = c(16, "bold", "black"),
                    font.tickslab = c(16, "bold", "black")) 
    
    # Add themes to the plot and table
    t$plot <- t$plot + theme(legend.title = element_text(size = 14, color = "black", face = "bold"),
                             legend.text = element_text(size = 14, color = "black", face = "bold"),
                             axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.y = element_text(size = 20, color = "black", face = "bold"))
    t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                               axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                               axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    )
    
    # Save the plot as a pdf
    pdf(paste0(Gene_set, "survplot.pdf"))
    print(t, newpage = FALSE)
    dev.off()
    
    # Perform Cox proportional hazards model
    fit <- coxph(Surv(time,Status) ~ score, data = PtSurv)
    
    test.ph <- cox.zph(fit)
    
    # Compute and save summary of fit, survival difference, and survival fit
    SumFit <- summary(fit)
    Diffsurv <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    Fitsurv <- survfit(Surv(time, Status) ~ score, data = PtSurv)
    
    capture.output(SumFit, file = paste0(Gene_set,"SumFit.csv"))
    capture.output(Diffsurv, file = paste0(Gene_set,"Diffsurv.csv"))
    capture.output(Fitsurv, file = paste0(Gene_set,"Fitsurv.csv"))
    
    library(ROCR)
    #Calculate and Plot ROC Curve
    {
      #Prepare Dataset
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      PtSurv2 <- PtSurv
      PtSurv2$time <- PtSurv2$time + 0.01
      # Filter by time over 1 month
      # PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
      PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
      PtSurv2 <- rbind(PtSurv2, PtSurv3)
      
      #Time-dependent ROC Curve
      midsurv <- 1000000
      
      df.y <- data.frame(PtSurv$time, PtSurv$Status)
      df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
      colnames(df) <- c("predictions", "time", "status")
      for (o in 1:nrow(df)) {
        df[o,4] <- as.numeric(0)
        if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
          df[o,4] <- 1
        }
      }
      
      
      #Plot ROC Curve
      df <- subset(df, select = c(1,4))
      colnames(df) <- c("predictions", "labels")
      
      
      pred <- prediction(df$predictions, df$labels)
      perf <- performance(pred,"tpr","fpr")
      line = data.frame(1:nrow(df),1:nrow(df))
      line = line/nrow(df)
      colnames(line) <- c("False Positive Rate", "True Positive Rate")
      
      
      pdf(paste0(Gene_set, "AUC.pdf"))
      plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(a=0, b= 1, col = "red")
      
      
      dev.off()
      auc_ROCR <- performance(pred, measure = "auc")
      auc_ROCR <- auc_ROCR@y.values[[1]]
    }
    
    
    #Select cutoff for activity
    {
      PtScore3 <- merge(PtGeneScore2, PtSpecificMutation, all.x = TRUE)
      bestcutoff <- res.cut$score
      
      write.csv(PtScore3, file = paste0(Gene_set,"vsMut.csv"))
      
      
      bestcutoff <- res.cut$cutpoint$cutpoint
    }
    
    
    #Prepare Activity data
    {
      
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(1,3,2,10))
      
      colnames(PtSurv) <- c("PATIENT_ID", "time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      
      df.activity <- data.frame(PtSurv)
      df.activity <- df.activity[order(as.numeric(df.activity[,4]), decreasing = TRUE),]
      df.activity <- df.activity %>% 
        mutate(row_id=row_number())
      df.activity <- merge(df.activity, PtSpecificMutation, by = "PATIENT_ID", keep.all = TRUE, all.x = TRUE)
      df.activity[is.na(df.activity)] <- "WT"
      df.activity <- subset(df.activity, select = c(4,6,7))
      df.activity$N_Score <- as.numeric(df.activity$score)
      df.activity <- df.activity[order(df.activity$score, decreasing = TRUE),]
    }
    
    #Plot Activity Graph
    {
      q <- ggplot(df.activity, aes(x = row_id, y = N_Score,color = Hugo_Symbol)) + geom_point(size = 1.8)  + ggtitle(paste0(Gene_set," TCGA ", TCGA_Cancer_ID,  " Curve")) + scale_y_continuous(
        #+ scale_color_manual(values=c("cornflower blue")) 
        # Features of the first axis
        name = paste0(Gene_set, "Score"), expand = c(0, 0)
        # Add a second axis and specify its features
        #sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
      ) + scale_x_continuous( name = paste0("Ranking by ", Gene_set), expand = c(0, 1))  + geom_hline(aes(yintercept=bestcutoff),  color = "gray50") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30), plot.title = element_text(hjust = 0.5))  + labs(colour="Legend",x="xxx",y="yyy")
      
      # Open a new device
      dev.new(width = 200, height = 200, unit = "px")
      # q
      ggsave(paste0(Gene_set," Activityv3.pdf"), plot=q)
      graphics.off()
    }
    
    
    
    #LOCC Calculations
    {
      
      #Calculate Parameters
      length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])
      max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      #Create Dataframe and Save numbers
      DF.qkscore <- data.frame(Mutation1, max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])/nrow(AggreCut2), max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]),res.cut$cutpoint$cutpoint, Dupe, AggreCut2$HR[AggreCut2$logp == (max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]))], auc_ROCR)
      colnames(DF.qkscore) <- c("Gene", "-log (p value)", "Percentage highly significantly", "Highest HR", "Lowest HR", "Cut","dupe", "significant HR", "auc_ROCR")
      write.csv(DF.qkscore, paste0(Gene_set, "_LOCC_Scorev3.csv"))
    }
    }
  }
  }
  
  graphics.off()
  
}

#Locc function for sampling (2-Fold Cross Validation), s for number for times to run sample simulation
Loccfunc2 <- function(j,s){
  
  Mutation1 <- PtExpressionZ[j,1]
  Gene_set <- Mutation1
  #Pick Gene Mutation Data
  if (all(Mutation1 != "")) {
    PtSpecificMutation <- subset(PtMutation, select = c(1,9:20,39))
    PtSpecificMutation <- filter(PtSpecificMutation, Hugo_Symbol == Mutation1)
    PtSpecificMutation <- filter(PtSpecificMutation, Variant_Classification != "Silent")
    PtSpecificMutation <-  unique(PtSpecificMutation[ , c("Hugo_Symbol", "Tumor_Sample_Barcode")], .keep_all = TRUE)
    PtSpecificMutation$Tumor_Sample_Barcode <- substr(PtSpecificMutation$Tumor_Sample_Barcode,1,nchar(PtSpecificMutation$Tumor_Sample_Barcode)-3)
    colnames(PtSpecificMutation)[2] <- "PATIENT_ID"
    
    # Merge the clinical and mutation data
    PtClinical2 <- merge(PtClinical, PtSpecificMutation, by = "PATIENT_ID", all.x = TRUE)
    PtClinical2$Hugo_Symbol <- replace_na(PtClinical2$Hugo_Symbol, "WT")
    PtClinical2$Mutant <- PtClinical2$Hugo_Symbol
    
    # Process the expression data
    PtExpressionGene <- PtExpressionZ[PtExpressionZ$Hugo_Symbol %in% Mutation1,]
    PtExpressionGene <- subset(PtExpressionGene, select = c(-2))
    PtExpressionGene <- na.omit(PtExpressionGene)
    Averages <- summarize_all(PtExpressionGene[,-c(1)], mean)
    Averages <- data.frame(0, Averages)
    colnames(Averages)[1] <- "Hugo_Symbol"
    PtExpressionGene <- rbind(PtExpressionGene, Averages)
    PtExpressionGene[nrow(PtExpressionGene),1] <- "Mean"
    
    # Compute Gene Z-Score
    PtGeneScore <- PtExpressionGene[nrow(PtExpressionGene),]
    PtGeneScore <- rbind(colnames(PtGeneScore), PtGeneScore)
    PtGeneScore2 <- t(PtGeneScore)
    PtGeneScore2 <- PtGeneScore2[2:nrow(PtGeneScore2),]
    colnames(PtGeneScore2) <- c("PATIENT_ID", "Gene_Score")
    PtGeneScore2[,1] <- substr(PtGeneScore2[,1],1,nchar(PtGeneScore2[,1])-3)
    PtGeneScore2[,1] <- str_replace(PtGeneScore2[,1], "\\.", "-")
    PtGeneScore2[,1] <- str_replace(PtGeneScore2[,1], "\\.", "-")
  
  
  
  # Set the initial cut-off point
  cutter <- 22
  
  # Merge clinical data with gene score data
  PtClinical4 <- merge(PtClinical, PtGeneScore2, by = "PATIENT_ID", all = TRUE)  
  
  # Select relevant columns
  PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
  
  # Further subset the data
  PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
  
  # Rename columns
  colnames(PtSurv) <- c("time", "status", "score")
  
  # Recode 'status' to numerical values and remove NA rows
  PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
  PtSurv <- na.omit(PtSurv)
  
  # Convert 'score' to numeric and replace NA values with 0
  PtSurv$score <- as.numeric(PtSurv$score)
  PtSurv[is.na(PtSurv)] <- 0
  
  # Order by 'score'
  PtSurv <- PtSurv[order(as.numeric(PtSurv$score), decreasing = TRUE), ]
  
  # Save the gene score value at the cut-off point
  GeneScoreValue <- PtSurv[cutter,3]
  
  # Assign high and low gene score groups
  PtSurv[1:cutter,3] <- "High_Gene"
  PtSurv[cutter:(nrow(PtSurv)+1),3] <- "Low_Gene"
  
  # Filter out only the high and low gene score groups
  PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
  
  # Calculate survival difference using the survdiff function
  CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
  
  # Create a dataframe with the results
  CutTrial <- data.frame(CutTrial$pvalue, CutTrial$n, CutTrial$obs, CutTrial$exp, GeneScoreValue)
  
  # Calculate HR
  HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
  CutTrial$HR <- HR
  
  # Initialize AggreCut dataframe with the results from the first cut-off point
  AggreCut <- CutTrial
  
  # Initialize a data frame to store the results
  AggreCut <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c(colnames(CutTrial)))))
  # Initialize a variable to track duplicates
  Dupe <- 0
  
  # Loop over different cut-off points
  for (cutter in 2:(nrow(PtSurv)-1)){  
    # Repeat the same process with different cut-off points
    PtClinical4 <- merge(PtClinical, PtGeneScore2,  by = "PATIENT_ID", all = TRUE)  
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    colnames(PtSurv) <- c("time", "status", "score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv[is.na(PtSurv)] <- 0
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    GeneScoreValue <- PtSurv[cutter,3]
    PtSurv[1:cutter,3] <- "High_Gene"
    PtSurv[cutter:(nrow(PtSurv)),3] <- "Low_Gene"
    PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
    
    # Calculate survival difference for the current cut-off point
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    CutTrial <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, GeneScoreValue)
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
    
    # Add the results to the AggreCut data frame
    AggreCut <- rbind(AggreCut, CutTrial)
  }
  
  # Now we process the AggreCut data frame
  AggreCut <- AggreCut[order(as.numeric(AggreCut[,1]), decreasing = FALSE),]
  AggreCut2 <- AggreCut[AggreCut$groups == "score=High_Gene",]
  AggreCut2 <- AggreCut2[order(as.numeric(AggreCut2[,3]), decreasing = FALSE),]
  AggreCut2$Fraction <- AggreCut2$Freq/nrow(AggreCut2)
  AggreCut2$logp <- -log(AggreCut2$CutTrial.pvalue, base = 10)
  
  #Check for Duplicates
  if (AggreCut2[nrow(AggreCut2), 6] == AggreCut2[nrow(AggreCut2)-1, 6] ){
    Dupe = 1 
  }
  
  
  #Count Dupes
  if (Dupe == 1){
    PreV  <- AggreCut2[1,6]
    AggreCut2[1,10] <- Dupe
    for (n in 2:nrow(AggreCut2)){
      NowV <- AggreCut2[n,6]
      if (PreV==NowV){
        Dupe = Dupe +  1
      } else {
        Dupe = 1
      }
      PreV <- NowV
      AggreCut2[n,10] <- Dupe
    }
  }
  
  DupeStart <- nrow(AggreCut2)-Dupe
  
  #Dupe Penalty
  for (n in DupeStart:nrow(AggreCut2)){
    PreVHR <- AggreCut2[n-1,7] - 0.1
    PreVHR <- max(PreVHR, 1)
    PreVP <- AggreCut2[n-1,9] - 0.1
    PreVP <- max(PreVP, 0)
    AggreCut2[n,7] <- PreVHR
    AggreCut2[n,9] <- PreVP
  }
  
  
  # output data
  write.csv(AggreCut2, file = paste0("AggreCut", Gene_set,cutter,s, ".csv"))
  
  #Plot LOCC Graph
  {
    df.cut <- subset(AggreCut2, select = c(7,8))
    df.cut <- melt(df.cut, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2 <- subset(AggreCut2, select = c(8,9))
    df.cut2$logp <- (df.cut2$logp/2)
    df.cut2 <- melt(df.cut2, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2$variable <- "-Log (p value)"
    df.cut2 <- rbind(df.cut2, df.cut)
    df.cut2 <- df.cut2[df.cut2$Fraction > 0.05,]
    df.cut2 <- df.cut2[df.cut2$Fraction < 0.95,]
    df.cut2$variable <- factor(df.cut2$variable, levels = c("HR", "-Log (p value)"))
    p <- ggplot(df.cut2, aes(x = Fraction, y = value,color = variable)) +geom_point() + geom_line(aes(color = variable)) + scale_color_manual(values=c("black", "orange")) + ggtitle(paste0(Gene_set, " TCGA ", TCGA_Cancer_ID ," LOCC")) + scale_y_continuous(
      # Features of the first axis
      name = "Hazard Ratio (HR)",
      # Add a second axis and specify its features
      sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
    ) + scale_x_continuous( name = "Fraction in High Activity Group", limits = c(0, 1), expand = c(0, 0)) + geom_hline(aes(yintercept=1), color = "red") + geom_hline(aes(yintercept=0.65),  color = "green") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30))  + labs(colour="Legend",x="xxx",y="yyy")
    
    
    ggsave(paste0(Gene_set,s,"_LOCC.pdf"), plot=p)
    # dev.off
  }
  
  PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
  PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
  colnames(PtSurv) <- c("time", "status", "score")
  PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
  PtSurv <- na.omit(PtSurv)
  PtSurv$score <- as.numeric(PtSurv$score)
  PtSurv$time <- as.numeric(PtSurv$time)
  PtSurv[is.na(PtSurv)] <- 0
  # PtSurv <- LiziSurv2
  PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
  
  
  {
    #Check if gene is expressed
    if (nrow(PtSurv) > 100 & PtSurv[nrow(PtSurv)/5,3] != PtSurv[nrow(PtSurv)-1,3]){
      
      #Mutation Analysis
      {
        # Subset the dataframe `PtClinical2` to keep only relevant columns
        PtSurvivalCurve <- subset(PtClinical2, select = c(1,30:37,38))
        
        # Further subset `PtSurvivalCurve` to keep only essential columns
        PtSurv <- subset(PtSurvivalCurve, select = c(1,3,2,10))
        
        # Rename the columns of `PtSurv`
        colnames(PtSurv) <- c("ID","time", "status", "mutation")
        
        # Recode the `status` column such that '1:DECEASED' is 1 and '0:LIVING' is 0
        PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
        
        # Remove NA values from `PtSurv`
        PtSurv <- na.omit(PtSurv)
        
        # Replace the mutations with "Mutant"
        PtSurv$mutation[PtSurv$mutation == Mutation1] <- "Mutant"
        
        
        # Order the `PtSurv` dataframe based on `mutation`
        PtSurv <- PtSurv[order(as.character(PtSurv$mutation), decreasing = FALSE), ]
        
        # Keep only distinct rows of `PtSurv` based on `ID`
        PtSurv <- distinct(PtSurv, ID, .keep_all = T)
        
        # Only perform survival analysis if there are more than 5 "Mutant"
        if (length(PtSurv$mutation[PtSurv$mutation == "Mutant"]) > 5){
          # Compute the survival difference
          surv_diff <- survdiff(Surv(time, Status) ~ mutation, data = PtSurv) 
          
          # Compute the survival fit
          surv_fit <- survfit(Surv(time, Status) ~ mutation, data = PtSurv)
          
          # Plot the survival analysis
          r <- survfit2(Surv(time, Status) ~ mutation, data = PtSurv) 
          t <-  ggsurvplot(r, data = PtSurv, risk.table = TRUE, size=1.2, fontsize = 7.5, 
                           risk.table.y.text = FALSE, tables.y.text = FALSE, risk.table.height = 0.35,
                           font.title = c(16, "bold", "black"), font.subtitle = c(16, "bold", "black"),
                           font.caption = c(16, "bold", "black"), font.x = c(16, "bold", "black"),
                           font.y = c(16, "bold", "black"), font.tickslab = c(16, "bold", "black"))
          t$plot <- t$plot + theme(legend.title = element_text(size = 18, color = "black", face = "bold"),
                                   legend.text = element_text(size = 18, color = "black", face = "bold"),
                                   axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                   axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                                   axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                                   axis.title.y = element_text(size = 20, color = "black", face = "bold"))
          t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                                     axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                     axis.title.x = element_text(size = 20, color = "black", face = "bold"))
          
          pdf(paste0(Gene_set,s,"_mutsurv.pdf"))
          print(t, newpage = FALSE)
          dev.off()
        }
      }
      
    #Cut at best cutpoint
    {
      
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      # PtSurv <- LiziSurv2
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      res.cut <- surv_cutpoint(PtSurv, time = "time", event = "Status",
                               variables = "score", minprop = 0.1)
      summary(res.cut)
      PtSurv$score[as.numeric(PtSurv$score) > res.cut$score[4]] <- as.character("High_Score")
      PtSurv$score[(PtSurv$score) <= res.cut$score[4]] <- as.character("Low_Score")
      
    }
    # Print a summary of `res.cut`
    summary(res.cut)
    
    # Open a new device
    dev.new(width = 200, height = 200, unit = "px")
    
    # Plot the survival curve with risk table
    r <- survfit2(Surv(time, Status) ~ score, data = PtSurv) 
    t <- ggsurvplot(r, data = PtSurv, risk.table = TRUE,
                    size = 1.2,
                    fontsize = 7.5,
                    risk.table.y.text = FALSE,
                    tables.y.text = FALSE,
                    risk.table.height = 0.35,
                    font.title = c(16, "bold", "black"),
                    font.subtitle = c(16, "bold", "black"),
                    font.caption = c(16, "bold", "black"),
                    font.x = c(16, "bold", "black"),
                    font.y = c(16, "bold", "black"),
                    font.tickslab = c(16, "bold", "black")) 
    
    # Add themes to the plot and table
    t$plot <- t$plot + theme(legend.title = element_text(size = 14, color = "black", face = "bold"),
                             legend.text = element_text(size = 14, color = "black", face = "bold"),
                             axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.y = element_text(size = 20, color = "black", face = "bold"))
    t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                               axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                               axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    )
    
    # Save the plot as a pdf
    pdf(paste0(Gene_set,s, "survplot.pdf"))
    print(t, newpage = FALSE)
    dev.off()
    
    # Perform Cox proportional hazards model
    fit <- coxph(Surv(time,Status) ~ score, data = PtSurv)
    
    # Compute and save summary of fit, survival difference, and survival fit
    SumFit <- summary(fit)
    Diffsurv <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    Fitsurv <- survfit(Surv(time, Status) ~ score, data = PtSurv)
    
    capture.output(SumFit, file = paste0(Gene_set,s,"SumFit.csv"))
    capture.output(Diffsurv, file = paste0(Gene_set,s,"Diffsurv.csv"))
    capture.output(Fitsurv, file = paste0(Gene_set,s,"Fitsurv.csv"))
    
    library(ROCR)
    #Calculate and Plot ROC Curve
    {
      #Prepare Dataset
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      PtSurv2 <- PtSurv
      PtSurv2$time <- PtSurv2$time + 0.01
      # Filter by time over 1 month
      # PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
      PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
      PtSurv2 <- rbind(PtSurv2, PtSurv3)
      
      #Time-dependent ROC Curve
      midsurv <- 1000000
      
      df.y <- data.frame(PtSurv$time, PtSurv$Status)
      df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
      colnames(df) <- c("predictions", "time", "status")
      for (o in 1:nrow(df)) {
        df[o,4] <- as.numeric(0)
        if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
          df[o,4] <- 1
        }
      }
      
      
      #Plot ROC Curve
      df <- subset(df, select = c(1,4))
      colnames(df) <- c("predictions", "labels")
      
      
      pred <- prediction(df$predictions, df$labels)
      perf <- performance(pred,"tpr","fpr")
      line = data.frame(1:nrow(df),1:nrow(df))
      line = line/nrow(df)
      colnames(line) <- c("False Positive Rate", "True Positive Rate")
      
      
      pdf(paste0(Gene_set,s, "AUC.pdf"))
      plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(a=0, b= 1, col = "red")
      
      
      dev.off()
      auc_ROCR <- performance(pred, measure = "auc")
      auc_ROCR <- auc_ROCR@y.values[[1]]
    }
    
    #Select cutoff for activity
    {
      PtScore3 <- merge(PtGeneScore2, PtSpecificMutation, all.x = TRUE)
      bestcutoff <- res.cut$score
      
      write.csv(PtScore3, file = paste0(Gene_set,s,"vsMut.csv"))
      
      
      bestcutoff <- res.cut$cutpoint$cutpoint
    }
    
    
    #Prepare Activity data
    {
      
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(1,3,2,10))
      
      colnames(PtSurv) <- c("PATIENT_ID", "time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      
      df.activity <- data.frame(PtSurv)
      df.activity <- df.activity[order(as.numeric(df.activity[,4]), decreasing = TRUE),]
      df.activity <- df.activity %>% 
        mutate(row_id=row_number())
      df.activity <- merge(df.activity, PtSpecificMutation, by = "PATIENT_ID", keep.all = TRUE, all.x = TRUE)
      df.activity[is.na(df.activity)] <- "WT"
      df.activity <- subset(df.activity, select = c(4,6,7))
      df.activity$n_Score <- as.numeric(df.activity$score)
      df.activity <- df.activity[order(df.activity$score, decreasing = TRUE),]
    }
    
    #Plot Activity Graph
    {
      q <- ggplot(df.activity, aes(x = row_id, y = n_Score,color = Hugo_Symbol)) + geom_point(size = 1.8)  + ggtitle(paste0(Gene_set," TCGA ", TCGA_Cancer_ID,  " Curve")) + scale_y_continuous(
        #+ scale_color_manual(values=c("cornflower blue")) 
        # Features of the first axis
        name = paste0(Gene_set, "Score"), expand = c(0, 0)
        # Add a second axis and specify its features
        #sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
      ) + scale_x_continuous( name = paste0("Ranking by ", Gene_set), expand = c(0, 1))  + geom_hline(aes(yintercept=bestcutoff),  color = "gray50") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30), plot.title = element_text(hjust = 0.5))  + labs(colour="Legend",x="xxx",y="yyy")
      
      # q
      ggsave(paste0(Gene_set,s," Activityv3.pdf"), plot=q)
      # graphics.off()
    }
    
    
    #LOCC Calculations
    {
      
      
      #Calculate Parameters
      length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])
      max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      #Create Dataframe and Save numbers
      DF.qkscore <- data.frame(Mutation1, max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])/nrow(AggreCut2), max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]),res.cut$cutpoint$cutpoint, Dupe, AggreCut2$HR[AggreCut2$logp == (max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]))], auc_ROCR)
      colnames(DF.qkscore) <- c("Gene", "-log (p value)", "Percentage highly significantly", "Highest HR", "Lowest HR", "Cut","dupe", "significant HR", "auc_ROCR")
      write.csv(DF.qkscore, paste0(Gene_set,s, "_LOCC_Scorev3.csv"))
    }
    }
  }
  }
  
  graphics.off()
  
}

#Locc score of Random ordered data
LoccfuncR <- function(r){
  
  Mutation1 <- PtExpressionZ[9,1]
  Gene_set <- Mutation1
  #Pick Gene Mutation Data
  if (all(Mutation1 != "")) {
    PtSpecificMutation <- subset(PtMutation, select = c(1,9:20,39))
    PtSpecificMutation <- filter(PtSpecificMutation, Hugo_Symbol == Mutation1)
    PtSpecificMutation <- filter(PtSpecificMutation, Variant_Classification != "Silent")
    PtSpecificMutation <-  unique(PtSpecificMutation[ , c("Hugo_Symbol", "Tumor_Sample_Barcode")], .keep_all = TRUE)
    PtSpecificMutation$Tumor_Sample_Barcode <- substr(PtSpecificMutation$Tumor_Sample_Barcode,1,nchar(PtSpecificMutation$Tumor_Sample_Barcode)-3)
    colnames(PtSpecificMutation)[2] <- "PATIENT_ID"
    
    # Merge the clinical and mutation data
    PtClinical2 <- merge(PtClinical, PtSpecificMutation, by = "PATIENT_ID", all.x = TRUE)
    PtClinical2$Hugo_Symbol <- replace_na(PtClinical2$Hugo_Symbol, "WT")
    PtClinical2$Mutant <- PtClinical2$Hugo_Symbol
    
    # Process the expression data
    PtExpressionGene <- PtExpressionZ[PtExpressionZ$Hugo_Symbol %in% Mutation1,]
    PtExpressionGene <- subset(PtExpressionGene, select = c(-2))
    PtExpressionGene <- na.omit(PtExpressionGene)
    Averages <- summarize_all(PtExpressionGene[,-c(1)], mean)
    Averages <- data.frame(0, Averages)
    colnames(Averages)[1] <- "Hugo_Symbol"
    PtExpressionGene <- rbind(PtExpressionGene, Averages)
    PtExpressionGene[nrow(PtExpressionGene),1] <- "Mean"
    
    # Compute Gene Z-Score
    PtGeneScore <- PtExpressionGene[nrow(PtExpressionGene),]
    PtGeneScore <- rbind(colnames(PtGeneScore), PtGeneScore)
    PtGeneScore2 <- t(PtGeneScore)
    PtGeneScore2 <- PtGeneScore2[2:nrow(PtGeneScore2),]
    colnames(PtGeneScore2) <- c("PATIENT_ID", "Gene_Score")
    PtGeneScore2[,1] <- substr(PtGeneScore2[,1],1,nchar(PtGeneScore2[,1])-3)
    PtGeneScore2[,1] <- str_replace(PtGeneScore2[,1], "\\.", "-")
    PtGeneScore2[,1] <- str_replace(PtGeneScore2[,1], "\\.", "-")
    Gene_set <- paste0('random',r)
    Mutation1 <- paste0('random',r)
    PtGeneScore2 <- as.data.frame(PtGeneScore2)
    PtGeneScore2 <- mutate(PtGeneScore2, Gene_Score = sample(c(0:1000), nrow(PtGeneScore2), replace = FALSE))
    
    
    # Set the initial cut-off point
    cutter <- 22
    
    # Merge clinical data with gene score data
    PtClinical4 <- merge(PtClinical, PtGeneScore2, by = "PATIENT_ID", all = TRUE)  
    
    # Select relevant columns
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    
    # Further subset the data
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    
    # Rename columns
    colnames(PtSurv) <- c("time", "status", "score")
    
    # Recode 'status' to numerical values and remove NA rows
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    
    # Convert 'score' to numeric and replace NA values with 0
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv[is.na(PtSurv)] <- 0
    
    # Order by 'score'
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score), decreasing = TRUE), ]
    
    # Save the gene score value at the cut-off point
    GeneScoreValue <- PtSurv[cutter,3]
    
    # Assign high and low gene score groups
    PtSurv[1:cutter,3] <- "High_Gene"
    PtSurv[cutter:(nrow(PtSurv)+1),3] <- "Low_Gene"
    
    # Filter out only the high and low gene score groups
    PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
    
    # Calculate survival difference using the survdiff function
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    
    # Create a dataframe with the results
    CutTrial <- data.frame(CutTrial$pvalue, CutTrial$n, CutTrial$obs, CutTrial$exp, GeneScoreValue)
    
    # Calculate HR
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
    
    # Initialize AggreCut dataframe with the results from the first cut-off point
    AggreCut <- CutTrial
    
    # Initialize a data frame to store the results
    AggreCut <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c(colnames(CutTrial)))))
    # Initialize a variable to track duplicates
    Dupe <- 0
    
    # Loop over different cut-off points
    for (cutter in 2:(nrow(PtSurv)-1)){  
      # Repeat the same process with different cut-off points
      PtClinical4 <- merge(PtClinical, PtGeneScore2,  by = "PATIENT_ID", all = TRUE)  
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv[is.na(PtSurv)] <- 0
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      GeneScoreValue <- PtSurv[cutter,3]
      PtSurv[1:cutter,3] <- "High_Gene"
      PtSurv[cutter:(nrow(PtSurv)),3] <- "Low_Gene"
      PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
      
      # Calculate survival difference for the current cut-off point
      CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
      CutTrial <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, GeneScoreValue)
      HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
      CutTrial$HR <- HR
      
      # Add the results to the AggreCut data frame
      AggreCut <- rbind(AggreCut, CutTrial)
    }
    
    # Now we process the AggreCut data frame
    AggreCut <- AggreCut[order(as.numeric(AggreCut[,1]), decreasing = FALSE),]
    AggreCut2 <- AggreCut[AggreCut$groups == "score=High_Gene",]
    AggreCut2 <- AggreCut2[order(as.numeric(AggreCut2[,3]), decreasing = FALSE),]
    AggreCut2$Fraction <- AggreCut2$Freq/nrow(AggreCut2)
    AggreCut2$logp <- -log(AggreCut2$CutTrial.pvalue, base = 10)
    
    #Check for Duplicates
    if (AggreCut2[nrow(AggreCut2), 6] == AggreCut2[nrow(AggreCut2)-1, 6] ){
      Dupe = 1 
    }
    
    
    #Count Dupes
    if (Dupe == 1){
      PreV  <- AggreCut2[1,6]
      AggreCut2[1,10] <- Dupe
      for (n in 2:nrow(AggreCut2)){
        NowV <- AggreCut2[n,6]
        if (PreV==NowV){
          Dupe = Dupe +  1
        } else {
          Dupe = 1
        }
        PreV <- NowV
        AggreCut2[n,10] <- Dupe
      }
    }
    
    DupeStart <- nrow(AggreCut2)-Dupe
    
    #Dupe Penalty
    for (n in DupeStart:nrow(AggreCut2)){
      PreVHR <- AggreCut2[n-1,7] - 0.1
      PreVHR <- max(PreVHR, 1)
      PreVP <- AggreCut2[n-1,9] - 0.1
      PreVP <- max(PreVP, 0)
      AggreCut2[n,7] <- PreVHR
      AggreCut2[n,9] <- PreVP
    }
    
    
    # output data
    write.csv(AggreCut2, file = paste0("AggreCut", Gene_set,cutter, ".csv"))
    
    #Plot LOCC Graph
    {
      df.cut <- subset(AggreCut2, select = c(7,8))
      df.cut <- melt(df.cut, id.vars = "Fraction" , na.rm =  TRUE)
      df.cut2 <- subset(AggreCut2, select = c(8,9))
      df.cut2$logp <- (df.cut2$logp/2)
      df.cut2 <- melt(df.cut2, id.vars = "Fraction" , na.rm =  TRUE)
      df.cut2$variable <- "-Log (p value)"
      df.cut2 <- rbind(df.cut2, df.cut)
      df.cut2 <- df.cut2[df.cut2$Fraction > 0.05,]
      df.cut2 <- df.cut2[df.cut2$Fraction < 0.95,]
      df.cut2$variable <- factor(df.cut2$variable, levels = c("HR", "-Log (p value)"))
      p <- ggplot(df.cut2, aes(x = Fraction, y = value,color = variable)) +geom_point() + geom_line(aes(color = variable)) + scale_color_manual(values=c("black", "orange")) + ggtitle(paste0(Gene_set, " TCGA ", TCGA_Cancer_ID ," LOCC")) + scale_y_continuous(
        # Features of the first axis
        name = "Hazard Ratio (HR)",
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
      ) + scale_x_continuous( name = "Fraction in High Activity Group", limits = c(0, 1), expand = c(0, 0)) + geom_hline(aes(yintercept=1), color = "red") + geom_hline(aes(yintercept=0.65),  color = "green") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30))  + labs(colour="Legend",x="xxx",y="yyy")
      
      
      ggsave(paste0(Gene_set,"_LOCC.pdf"), plot=p)
      # dev.off
    }
    
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    colnames(PtSurv) <- c("time", "status", "score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv$time <- as.numeric(PtSurv$time)
    PtSurv[is.na(PtSurv)] <- 0
    # PtSurv <- LiziSurv2
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    
    {
      #Check if gene is expressed
      if (nrow(PtSurv) > 100 & PtSurv[nrow(PtSurv)/5,3] != PtSurv[nrow(PtSurv)-1,3]){
        
        #Mutation Analysis
        {
          # Subset the dataframe `PtClinical2` to keep only relevant columns
          PtSurvivalCurve <- subset(PtClinical2, select = c(1,30:37,38))
          
          # Further subset `PtSurvivalCurve` to keep only essential columns
          PtSurv <- subset(PtSurvivalCurve, select = c(1,3,2,10))
          
          # Rename the columns of `PtSurv`
          colnames(PtSurv) <- c("ID","time", "status", "mutation")
          
          # Recode the `status` column such that '1:DECEASED' is 1 and '0:LIVING' is 0
          PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
          
          # Remove NA values from `PtSurv`
          PtSurv <- na.omit(PtSurv)
          
          # Replace the mutations with "Mutant"
          PtSurv$mutation[PtSurv$mutation == Mutation1] <- "Mutant"
          
          
          # Order the `PtSurv` dataframe based on `mutation`
          PtSurv <- PtSurv[order(as.character(PtSurv$mutation), decreasing = FALSE), ]
          
          # Keep only distinct rows of `PtSurv` based on `ID`
          PtSurv <- distinct(PtSurv, ID, .keep_all = T)
          
          # Only perform survival analysis if there are more than 5 "Mutant"
          if (length(PtSurv$mutation[PtSurv$mutation == "Mutant"]) > 5){
            # Compute the survival difference
            surv_diff <- survdiff(Surv(time, Status) ~ mutation, data = PtSurv) 
            
            # Compute the survival fit
            surv_fit <- survfit(Surv(time, Status) ~ mutation, data = PtSurv)
            
            # Plot the survival analysis
            r <- survfit2(Surv(time, Status) ~ mutation, data = PtSurv) 
            t <-  ggsurvplot(r, data = PtSurv, risk.table = TRUE, size=1.2, fontsize = 7.5, 
                             risk.table.y.text = FALSE, tables.y.text = FALSE, risk.table.height = 0.35,
                             font.title = c(16, "bold", "black"), font.subtitle = c(16, "bold", "black"),
                             font.caption = c(16, "bold", "black"), font.x = c(16, "bold", "black"),
                             font.y = c(16, "bold", "black"), font.tickslab = c(16, "bold", "black"))
            t$plot <- t$plot + theme(legend.title = element_text(size = 18, color = "black", face = "bold"),
                                     legend.text = element_text(size = 18, color = "black", face = "bold"),
                                     axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                     axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                                     axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                                     axis.title.y = element_text(size = 20, color = "black", face = "bold"))
            t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                                       axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                       axis.title.x = element_text(size = 20, color = "black", face = "bold"))
            
            pdf(paste0(Gene_set,"_mutsurv.pdf"))
            print(t, newpage = FALSE)
            dev.off()
          }
        }
        
        
        
        
        #Cut at best cutpoint
        {
          
          PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
          PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
          colnames(PtSurv) <- c("time", "status", "score")
          PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
          PtSurv <- na.omit(PtSurv)
          PtSurv$score <- as.numeric(PtSurv$score)
          PtSurv$time <- as.numeric(PtSurv$time)
          PtSurv[is.na(PtSurv)] <- 0
          # PtSurv <- LiziSurv2
          PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
          res.cut <- surv_cutpoint(PtSurv, time = "time", event = "Status",
                                   variables = "score", minprop = 0.1)
          summary(res.cut)
          PtSurv$score[as.numeric(PtSurv$score) > res.cut$score[4]] <- as.character("High_Score")
          PtSurv$score[(PtSurv$score) <= res.cut$score[4]] <- as.character("Low_Score")
          
        }
        # Print a summary of `res.cut`
        summary(res.cut)
        
        # Open a new device
        dev.new(width = 200, height = 200, unit = "px")
        
        # Plot the survival curve with risk table
        r <- survfit2(Surv(time, Status) ~ score, data = PtSurv) 
        t <- ggsurvplot(r, data = PtSurv, risk.table = TRUE,
                        size = 1.2,
                        fontsize = 7.5,
                        risk.table.y.text = FALSE,
                        tables.y.text = FALSE,
                        risk.table.height = 0.35,
                        font.title = c(16, "bold", "black"),
                        font.subtitle = c(16, "bold", "black"),
                        font.caption = c(16, "bold", "black"),
                        font.x = c(16, "bold", "black"),
                        font.y = c(16, "bold", "black"),
                        font.tickslab = c(16, "bold", "black")) 
        
        # Add themes to the plot and table
        t$plot <- t$plot + theme(legend.title = element_text(size = 14, color = "black", face = "bold"),
                                 legend.text = element_text(size = 14, color = "black", face = "bold"),
                                 axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                 axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                                 axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                                 axis.title.y = element_text(size = 20, color = "black", face = "bold"))
        t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                                   axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                   axis.title.x = element_text(size = 20, color = "black", face = "bold"),
        )
        
        # Save the plot as a pdf
        pdf(paste0(Gene_set, "survplot.pdf"))
        print(t, newpage = FALSE)
        dev.off()
        
        # Perform Cox proportional hazards model
        fit <- coxph(Surv(time,Status) ~ score, data = PtSurv)
        
        # Compute and save summary of fit, survival difference, and survival fit
        SumFit <- summary(fit)
        Diffsurv <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
        Fitsurv <- survfit(Surv(time, Status) ~ score, data = PtSurv)
        
        capture.output(SumFit, file = paste0(Gene_set,"SumFit.csv"))
        capture.output(Diffsurv, file = paste0(Gene_set,"Diffsurv.csv"))
        capture.output(Fitsurv, file = paste0(Gene_set,"Fitsurv.csv"))
        
        library(ROCR)
        #Calculate and Plot ROC Curve
        {
          #Prepare Dataset
          PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
          PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
          colnames(PtSurv) <- c("time", "status", "score")
          PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
          PtSurv <- na.omit(PtSurv)
          PtSurv$score <- as.numeric(PtSurv$score)
          PtSurv$time <- as.numeric(PtSurv$time)
          PtSurv[is.na(PtSurv)] <- 0
          PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
          PtSurv2 <- PtSurv
          PtSurv2$time <- PtSurv2$time + 0.01
          # Filter by time over 1 month
          # PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
          PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
          PtSurv2 <- rbind(PtSurv2, PtSurv3)
          
          #Time-dependent ROC Curve
          midsurv <- 1000000
          
          df.y <- data.frame(PtSurv$time, PtSurv$Status)
          df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
          colnames(df) <- c("predictions", "time", "status")
          for (o in 1:nrow(df)) {
            df[o,4] <- as.numeric(0)
            if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
              df[o,4] <- 1
            }
          }
          
          
          #Plot ROC Curve
          df <- subset(df, select = c(1,4))
          colnames(df) <- c("predictions", "labels")
          
          
          pred <- prediction(df$predictions, df$labels)
          perf <- performance(pred,"tpr","fpr")
          line = data.frame(1:nrow(df),1:nrow(df))
          line = line/nrow(df)
          colnames(line) <- c("False Positive Rate", "True Positive Rate")
          
          
          pdf(paste0(Gene_set, "AUC.pdf"))
          plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
          abline(a=0, b= 1, col = "red")
          
          
          dev.off()
          auc_ROCR <- performance(pred, measure = "auc")
          auc_ROCR <- auc_ROCR@y.values[[1]]
        }
        
        
        #Select cutoff for activity
        {
          PtScore3 <- merge(PtGeneScore2, PtSpecificMutation, all.x = TRUE)
          bestcutoff <- res.cut$score
          
          write.csv(PtScore3, file = paste0(Gene_set,"vsMut.csv"))
          
          
          bestcutoff <- res.cut$cutpoint$cutpoint
        }
        
        
        #Prepare Activity data
        {
          
          PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
          PtSurv <- subset(PtSurvivalCurve, select = c(1,3,2,10))
          
          colnames(PtSurv) <- c("PATIENT_ID", "time", "status", "score")
          PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
          PtSurv <- na.omit(PtSurv)
          
          df.activity <- data.frame(PtSurv)
          df.activity <- df.activity[order(as.numeric(df.activity[,4]), decreasing = TRUE),]
          df.activity <- df.activity %>% 
            mutate(row_id=row_number())
          df.activity <- merge(df.activity, PtSpecificMutation, by = "PATIENT_ID", keep.all = TRUE, all.x = TRUE)
          df.activity[is.na(df.activity)] <- "WT"
          df.activity <- subset(df.activity, select = c(4,6,7))
          df.activity$N_Score <- as.numeric(df.activity$score)
          df.activity <- df.activity[order(df.activity$score, decreasing = TRUE),]
        }
        
        #Plot Activity Graph
        {
          q <- ggplot(df.activity, aes(x = row_id, y = N_Score,color = Hugo_Symbol)) + geom_point(size = 1.8)  + ggtitle(paste0(Gene_set," TCGA ", TCGA_Cancer_ID,  " Curve")) + scale_y_continuous(
            #+ scale_color_manual(values=c("cornflower blue")) 
            # Features of the first axis
            name = paste0(Gene_set, "Score"), expand = c(0, 0)
            # Add a second axis and specify its features
            #sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
          ) + scale_x_continuous( name = paste0("Ranking by ", Gene_set), expand = c(0, 1))  + geom_hline(aes(yintercept=bestcutoff),  color = "gray50") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30), plot.title = element_text(hjust = 0.5))  + labs(colour="Legend",x="xxx",y="yyy")
          
          # q
          ggsave(paste0(Gene_set," Activityv3.pdf"), plot=q)
          # graphics.off()
        }
        
        
        
        #LOCC Calculations
        {
          
          #Calculate Parameters
          length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])
          max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
          min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
          max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
          #Create Dataframe and Save numbers
          DF.qkscore <- data.frame(Mutation1, max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])/nrow(AggreCut2), max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]),res.cut$cutpoint$cutpoint, Dupe, AggreCut2$HR[AggreCut2$logp == (max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]))], auc_ROCR)
          colnames(DF.qkscore) <- c("Gene", "-log (p value)", "Percentage highly significantly", "Highest HR", "Lowest HR", "Cut","dupe", "significant HR", "auc_ROCR")
          write.csv(DF.qkscore, paste0(Gene_set, "_LOCC_Scorev3.csv"))
        }
      }
    }
  }
  
  graphics.off()
  
}


#Test for cox PH assumption as needed
# fit
# test.ph


# Loop over all genes, will take a long time
for (j in 1:20531){
  Loccfunc(j)
}

Loccfunc(j)
#Create dataframe with all gene LOCC scores
{
  if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set, "_LOCC_Scorev3.csv")))){
    
    df.score <- read.csv(file = paste0(Gene_set, "_LOCC_Scorev3.csv"))
  }
  df.score2 <- df.score
  df.score2 <- df.score2[-1,]
  file.size(paste0(Gene_set, "_LOCC_Scorev3.csv"))
  
  #Combine Dataframes
  x <- 20531
  for (y in 1:x){
    
    Mutation1 <- PtExpressionZ[y,1]
    Gene_set <- paste0(Mutation1,"")
    if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set, "_LOCC_Scorev3.csv")))){
      
      df.score <- read.csv(file = paste0(Gene_set, "_LOCC_Scorev3.csv"))
      df.score2 <- rbind(df.score2, df.score)
    }
  }
  #Calculate LOCC 
  for (y in 1:nrow(df.score2)){
    if (as.numeric(df.score2[y,9]) > 1){
      df.score2[y,11] <- -(df.score2[y,3]*df.score2[y,4]*df.score2[y,9])
    } else {
      df.score2[y,11] <- (df.score2[y,3]*df.score2[y,4]/df.score2[y,9])  
    }
  }
  colnames(df.score2)[11] <- "LOCC"
  write.csv(df.score2, file = "totalLoccScore.csv")
}

# Loop over 10000 random ordered simulation, will take a long time
numRamdom <- 10000
for (r in 1:numRamdom){
  LoccfuncR(r)
}
#Combine Data
df.score <- read.csv(file = paste0(Gene_set, "_LOCC_Scorev3.csv"))
df.score2 <- df.score
df.score2 <- df.score2[-1,]
for (y in 1:numRamdom){
  Gene_set <- paste0('random',y)
  if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set, "_LOCC_Scorev3.csv")))){
    
    df.score <- read.csv(file = paste0(Gene_set, "_LOCC_Scorev3.csv"))
    df.score2 <- rbind(df.score2, df.score)
  }
}
#Calculate LOCC 
{
for (y in 1:nrow(df.score2)){
  if (as.numeric(df.score2[y,9]) > 1){
    df.score2[y,11] <- -(df.score2[y,3]*df.score2[y,4]*df.score2[y,9])
  } else {
    df.score2[y,11] <- (df.score2[y,3]*df.score2[y,4]/df.score2[y,9])  
  }
}
colnames(df.score2)[11] <- "LOCC"
write.csv(df.score2, file = "RandomLoccScores.csv")
}

loccscores <- read.csv(file = 'RandomLoccV1.csv', header = FALSE)
loccscorespvalues <- loccscores[,4]
loccscorespvalues <- loccscorespvalues[2:18789]
loccscorespvalues <- replace(loccscorespvalues, loccscorespvalues >= 1, 1)
loccscorespvalues[18777]
qobj <- qvalue(p = loccscorespvalues, fdr.level = 0.1, pi0.method = "smoother")

qobj$qvalues
options(repr.plot.width=6, repr.plot.height=5)
plot(qobj)

loccscoreqvalues <- qobj$qvalues
write.csv(loccscoreqvalues, file = "LOCCqvalue.csv")

#RISK Score Analysis 
{
  NGS9 <- data.frame(c("ANXA10",
                       "CYP2C9",
                       "HMMR",
                       "KIF20A",
                       "LAPTM4B",
                       "LCAT",
                       "LECT2",
                       "MAGEA6",
                       "RDH16",
                       "SPP1",
                       "TPX2",
                       "TTK"))
  colnames(NGS9) <- "V1"
  Gene_set <- "RISK Original Score"
  {
    PtExpressionN2 <- PtExpression[PtExpression$Hugo_Symbol %in% NGS9$V1,]
    PtExpressionN2 <- na.omit(PtExpressionN2)
    PtExpressionN2[,-c(1)] <- log2(PtExpressionN2[,-c(1)]+1)
    for (p in 3:ncol(PtExpressionN2)){
      PtExpressionN2[13,p]   <- (-0.0627 * PtExpressionN2[1,p]) + ( -0.02365*PtExpressionN2[2,p]) + (0.0565*PtExpressionN2[3,p]) +(0.0386*PtExpressionN2[4,p]) + (0.0015*PtExpressionN2[5,p]) + (-0.02365*PtExpressionN2[6,p]) +(-0.02237*PtExpressionN2[7,p])+(0.05284*PtExpressionN2[8,p])+(-0.01796*PtExpressionN2[9,p])+(0.05124*PtExpressionN2[10,p])+(0.05377*PtExpressionN2[11,p])+(0.09788*PtExpressionN2[12,p])
    }
    PtExpressionN2[nrow(PtExpressionN2),1] <- "Mean"
    
    PtN2Score <- PtExpressionN2
    PtN2Score <- PtExpressionN2[nrow(PtExpressionN2),]
    PtN2Score <- subset(PtN2Score, select = c(-2))
    
    PtN2Score <- rbind(colnames(PtN2Score), PtN2Score)
    PtN2Score2 <- t(PtN2Score)
    PtN2Score2 <- PtN2Score2[2:nrow(PtN2Score2),]
    colnames(PtN2Score2)[1] <- c("PATIENT_ID")
    colnames(PtN2Score2)[1:2] <- c("PATIENT_ID", "N2_Score")
    
    PtN2Score2[,1] <- substr(PtN2Score2[,1],1,nchar(PtN2Score2[,1])-3)
    PtN2Score2[,1] <- str_replace(PtN2Score2[,1], "\\.", "-")
    PtN2Score2[,1] <- str_replace(PtN2Score2[,1], "\\.", "-")
    
  }
  
  {
    cutter <- 22
    PtClinical4 <- merge(PtClinical, PtN2Score2,  by = "PATIENT_ID", .keep_all= TRUE)  
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    
    colnames(PtSurv) <- c("time", "status", "score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv[is.na(PtSurv)] <- 0
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    NScoreValue <- PtSurv[cutter,3]
    PtSurv[1:cutter,3] <- "High_RISK"
    PtSurv[cutter:(nrow(PtSurv)+1),3] <- "Low_RISK"
    PtSurv <- filter(PtSurv, score == "High_RISK" | score == "Low_RISK")
    
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    CutTrial  <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, NScoreValue)
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
  }
  
  # Initialize a data frame to store the results
  AggreCut <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c(colnames(CutTrial)))))
  # Initialize a variable to track duplicates
  Dupe <- 0
  
  # Loop over different cut-off points
  for (cutter in 2:(nrow(PtSurv)-1)){  
    # Repeat the same process with different cut-off points
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    colnames(PtSurv) <- c("time", "status", "score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv[is.na(PtSurv)] <- 0
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    GeneScoreValue <- PtSurv[cutter,3]
    PtSurv[1:cutter,3] <- "High_Gene"
    PtSurv[cutter:(nrow(PtSurv)),3] <- "Low_Gene"
    PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
    
    # Calculate survival difference for the current cut-off point
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    CutTrial <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, GeneScoreValue)
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
    
    # Add the results to the AggreCut data frame
    AggreCut <- rbind(AggreCut, CutTrial)
  }
  
  # Now we process the AggreCut data frame
  AggreCut <- AggreCut[order(as.numeric(AggreCut[,1]), decreasing = FALSE),]
  AggreCut2 <- AggreCut[AggreCut$groups == "score=High_Gene",]
  AggreCut2 <- AggreCut2[order(as.numeric(AggreCut2[,3]), decreasing = FALSE),]
  AggreCut2$Fraction <- AggreCut2$Freq/nrow(AggreCut2)
  AggreCut2$logp <- -log(AggreCut2$CutTrial.pvalue, base = 10)
  
  #Check for Duplicates
  if (AggreCut2[nrow(AggreCut2), 6] == AggreCut2[nrow(AggreCut2)-1, 6] ){
    Dupe = 1 
  }
  
  
  #Count Dupes
  if (Dupe == 1){
    PreV  <- AggreCut2[1,6]
    AggreCut2[1,10] <- Dupe
    for (n in 2:nrow(AggreCut2)){
      NowV <- AggreCut2[n,6]
      if (PreV==NowV){
        Dupe = Dupe +  1
      } else {
        Dupe = 1
      }
      PreV <- NowV
      AggreCut2[n,10] <- Dupe
    }
  }
  
  DupeStart <- nrow(AggreCut2)-Dupe
  
  #Dupe Penalty
  for (n in DupeStart:nrow(AggreCut2)){
    PreVHR <- AggreCut2[n-1,7] - 0.1
    PreVHR <- max(PreVHR, 1)
    PreVP <- AggreCut2[n-1,9] - 0.1
    PreVP <- max(PreVP, 0)
    AggreCut2[n,7] <- PreVHR
    AggreCut2[n,9] <- PreVP
  }
  
  
  # output data
  write.csv(AggreCut2, file = paste0("AggreCut", Gene_set,cutter, ".csv"))
  
  #Plot LOCC Graph
  {
    df.cut <- subset(AggreCut2, select = c(7,8))
    df.cut <- melt(df.cut, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2 <- subset(AggreCut2, select = c(8,9))
    df.cut2$logp <- (df.cut2$logp/2)
    df.cut2 <- melt(df.cut2, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2$variable <- "-Log (p value)"
    df.cut2 <- rbind(df.cut2, df.cut)
    df.cut2 <- df.cut2[df.cut2$Fraction > 0.05,]
    df.cut2 <- df.cut2[df.cut2$Fraction < 0.95,]
    df.cut2$variable <- factor(df.cut2$variable, levels = c("HR", "-Log (p value)"))
    p <- ggplot(df.cut2, aes(x = Fraction, y = value,color = variable)) +geom_point() + geom_line(aes(color = variable)) + scale_color_manual(values=c("black", "orange")) + ggtitle(paste0(Gene_set, " TCGA ", TCGA_Cancer_ID ," LOCC")) + scale_y_continuous(
      # Features of the first axis
      name = "Hazard Ratio (HR)",
      # Add a second axis and specify its features
      sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
    ) + scale_x_continuous( name = "Fraction in High Activity Group", limits = c(0, 1), expand = c(0, 0)) + geom_hline(aes(yintercept=1), color = "red") + geom_hline(aes(yintercept=0.65),  color = "green") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30))  + labs(colour="Legend",x="xxx",y="yyy")
    
    
    ggsave(paste0(Gene_set,"_LOCC.pdf"), plot=p)
    # dev.off
    
    
    
    #Cut at best cutpoint
    {
      
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      # PtSurv <- LiziSurv2
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      res.cut <- surv_cutpoint(PtSurv, time = "time", event = "Status",
                               variables = "score", minprop = 0.1)
      summary(res.cut)
      PtSurv$score[as.numeric(PtSurv$score) > res.cut$score[4]] <- as.character("High_Score")
      PtSurv$score[(PtSurv$score) <= res.cut$score[4]] <- as.character("Low_Score")
      
    }
    # Print a summary of `res.cut`
    summary(res.cut)
    
    # Open a new device
    dev.new(width = 200, height = 200, unit = "px")
    
    # Plot the survival curve with risk table
    r <- survfit2(Surv(time, Status) ~ score, data = PtSurv) 
    t <- ggsurvplot(r, data = PtSurv, risk.table = TRUE,
                    size = 1.2,
                    fontsize = 7.5,
                    risk.table.y.text = FALSE,
                    tables.y.text = FALSE,
                    risk.table.height = 0.35,
                    font.title = c(16, "bold", "black"),
                    font.subtitle = c(16, "bold", "black"),
                    font.caption = c(16, "bold", "black"),
                    font.x = c(16, "bold", "black"),
                    font.y = c(16, "bold", "black"),
                    font.tickslab = c(16, "bold", "black")) 
    
    # Add themes to the plot and table
    t$plot <- t$plot + theme(legend.title = element_text(size = 14, color = "black", face = "bold"),
                             legend.text = element_text(size = 14, color = "black", face = "bold"),
                             axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.y = element_text(size = 20, color = "black", face = "bold"))
    t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                               axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                               axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    )
    
    # Save the plot as a pdf
    pdf(paste0(Gene_set, "survplot.pdf"))
    print(t, newpage = FALSE)
    dev.off()
    
    # Perform Cox proportional hazards model
    fit <- coxph(Surv(time,Status) ~ score, data = PtSurv)
    
    # Compute and save summary of fit, survival difference, and survival fit
    SumFit <- summary(fit)
    Diffsurv <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    Fitsurv <- survfit(Surv(time, Status) ~ score, data = PtSurv)
    
    capture.output(SumFit, file = paste0(Gene_set,"SumFit.csv"))
    capture.output(Diffsurv, file = paste0(Gene_set,"Diffsurv.csv"))
    capture.output(Fitsurv, file = paste0(Gene_set,"Fitsurv.csv"))
    
    library(ROCR)
    #Calculate and Plot ROC Curve (non-time dependent)
    {
      #Prepare Dataset
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      PtSurv2 <- PtSurv
      PtSurv2$time <- PtSurv2$time + 0.01
      # Filter by time over 1 month
      # PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
      PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
      PtSurv2 <- rbind(PtSurv2, PtSurv3)
      
      #Time-dependent ROC Curve
      midsurv <- 1000000
      
      df.y <- data.frame(PtSurv$time, PtSurv$Status)
      df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
      colnames(df) <- c("predictions", "time", "status")
      for (o in 1:nrow(df)) {
        df[o,4] <- as.numeric(0)
        if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
          df[o,4] <- 1
        }
      }
      
      
      #Plot ROC Curve
      df <- subset(df, select = c(1,4))
      colnames(df) <- c("predictions", "labels")
      
      
      pred <- prediction(df$predictions, df$labels)
      perf <- performance(pred,"tpr","fpr")
      line = data.frame(1:nrow(df),1:nrow(df))
      line = line/nrow(df)
      colnames(line) <- c("False Positive Rate", "True Positive Rate")
      
      
      pdf(paste0(Gene_set, "AUC.pdf"))
      plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(a=0, b= 1, col = "red")
      
      
      dev.off()
      auc_ROCR <- performance(pred, measure = "auc")
      auc_ROCR <- auc_ROCR@y.values[[1]]
    }
    
    #Calculate and Plot ROC Curve (1 year)
    {
      #Prepare Dataset
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      PtSurv2 <- PtSurv
      PtSurv2$time <- PtSurv2$time + 0.01
      #Filter by time over 1 month
      PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
      PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
      PtSurv2 <- rbind(PtSurv2, PtSurv3)
      
      #Time-dependent ROC Curve
      midsurv <- 12
      
      df.y <- data.frame(PtSurv$time, PtSurv$Status)
      df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
      colnames(df) <- c("predictions", "time", "status")
      for (o in 1:nrow(df)) {
        df[o,4] <- as.numeric(0)
        if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
          df[o,4] <- 1
        }
      }
      
      
      #Plot ROC Curve
      df <- subset(df, select = c(1,4))
      colnames(df) <- c("predictions", "labels")
      
      
      pred <- prediction(df$predictions, df$labels)
      perf <- performance(pred,"tpr","fpr")
      line = data.frame(1:nrow(df),1:nrow(df))
      line = line/nrow(df)
      colnames(line) <- c("False Positive Rate", "True Positive Rate")
      
      
      pdf(paste0(Gene_set, "1yr_AUC.pdf"))
      plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(a=0, b= 1, col = "red")
      
      
      dev.off()
      auc_ROCR <- performance(pred, measure = "auc")
      auc_ROCR <- auc_ROCR@y.values[[1]]
    }
    
    
    #LOCC Calculations
    {
      
      #Calculate Parameters
      length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])
      max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      #Create Dataframe and Save numbers
      DF.qkscore <- data.frame("Risk Score", max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])/nrow(AggreCut2), max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]),res.cut$cutpoint$cutpoint, AggreCut2$HR[AggreCut2$logp == (max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]))], auc_ROCR)
      colnames(DF.qkscore) <- c("Gene", "-log (p value)", "Percentage highly significantly", "Highest HR", "Lowest HR", "Cut", "significant HR", "auc_ROCR")
      write.csv(DF.qkscore, paste0(Gene_set, "_LOCC_Scorev3.csv"))
      
    }
    
  }
  
}

#AIC Calculations
{
  NGS9 <- data.frame(c("ANXA10",
                       "CYP2C9",
                       "HMMR",
                       "KIF20A",
                       "LAPTM4B",
                       "LCAT",
                       "LECT2",
                       "MAGEA6",
                       "RDH16",
                       "SPP1",
                       "TPX2",
                       "TTK"))
  colnames(NGS9) <- "V1"
  Gene_set <- "RISK Original Score"
  {
    PtExpressionN2 <- PtExpression[PtExpression$Hugo_Symbol %in% NGS9$V1,]
    PtExpressionN2 <- na.omit(PtExpressionN2)
    PtExpressionN2[,-c(1)] <- log2(PtExpressionN2[,-c(1)]+1)
    for (p in 3:ncol(PtExpressionN2)){
      PtExpressionN2[13,p]   <- (-0.0627 * PtExpressionN2[1,p]) + ( -0.02365*PtExpressionN2[2,p]) + (0.0565*PtExpressionN2[3,p]) +(0.0386*PtExpressionN2[4,p]) + (0.0015*PtExpressionN2[5,p]) + (-0.02365*PtExpressionN2[6,p]) +(-0.02237*PtExpressionN2[7,p])+(0.05284*PtExpressionN2[8,p])+(-0.01796*PtExpressionN2[9,p])+(0.05124*PtExpressionN2[10,p])+(0.05377*PtExpressionN2[11,p])+(0.09788*PtExpressionN2[12,p])
    }
    PtExpressionN2[nrow(PtExpressionN2),1] <- "Mean"
    
    PtN2Score <- PtExpressionN2
    # PtN2Score <- PtExpressionN2[nrow(PtExpressionN2),]
    PtN2Score <- subset(PtN2Score, select = c(-2))
    
    PtN2Score <- rbind(colnames(PtN2Score), PtN2Score)
    PtN2Score2 <- t(PtN2Score)
    PtN2Score2 <- PtN2Score2[2:nrow(PtN2Score2),]
    colnames(PtN2Score2)[1] <- c("PATIENT_ID")
    colnames(PtN2Score2)[1:14] <- c("PATIENT_ID",'a','b','c','d','e','f','g','h','i', 'j', 'k','l',"score")
    
    PtN2Score2[,1] <- substr(PtN2Score2[,1],1,nchar(PtN2Score2[,1])-3)
    PtN2Score2[,1] <- str_replace(PtN2Score2[,1], "\\.", "-")
    PtN2Score2[,1] <- str_replace(PtN2Score2[,1], "\\.", "-")
    
  }
{
  cutter <- 22
  PtClinical4 <- merge(PtClinical, PtN2Score2,  by = "PATIENT_ID", .keep_all= TRUE)  
  PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:50))
  PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10:22))
  
  colnames(PtSurv) <- c("time", "status",'a','b','c','d','e','f','g','h','i', 'j', 'k','l',"score")
  PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
  PtSurv <- na.omit(PtSurv)
  PtSurv$score <- as.numeric(PtSurv$score)
  PtSurv[is.na(PtSurv)] <- 0
  PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
  NScoreValue <- PtSurv[cutter,15]
  PtSurv[1:cutter,15] <- "High_RISK"
  PtSurv[cutter:(nrow(PtSurv)+1),15] <- "Low_RISK"
  PtSurv <- filter(PtSurv, score == "High_RISK" | score == "Low_RISK")
  
  CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
  CutTrial  <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, NScoreValue)
  HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
  CutTrial$HR <- HR
}

# Initialize a data frame to store the results
AggreCut <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c(colnames(CutTrial)))))
# Initialize a variable to track duplicates
Dupe <- 0

# Loop over different cut-off points
for (cutter in 2:(nrow(PtSurv)-1)){  
  # Repeat the same process with different cut-off points
  PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
  PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
  colnames(PtSurv) <- c("time", "status", "score")
  PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
  PtSurv <- na.omit(PtSurv)
  PtSurv$score <- as.numeric(PtSurv$score)
  PtSurv[is.na(PtSurv)] <- 0
  PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
  GeneScoreValue <- PtSurv[cutter,3]
  PtSurv[1:cutter,3] <- "High_Gene"
  PtSurv[cutter:(nrow(PtSurv)),3] <- "Low_Gene"
  PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
  
  # Calculate survival difference for the current cut-off point
  CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
  CutTrial <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, GeneScoreValue)
  HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
  CutTrial$HR <- HR
  
  # Add the results to the AggreCut data frame
  AggreCut <- rbind(AggreCut, CutTrial)
}

# Now we process the AggreCut data frame
AggreCut <- AggreCut[order(as.numeric(AggreCut[,1]), decreasing = FALSE),]
AggreCut2 <- AggreCut[AggreCut$groups == "score=High_Gene",]
AggreCut2 <- AggreCut2[order(as.numeric(AggreCut2[,3]), decreasing = FALSE),]
AggreCut2$Fraction <- AggreCut2$Freq/nrow(AggreCut2)
AggreCut2$logp <- -log(AggreCut2$CutTrial.pvalue, base = 10)

#Check for Duplicates
if (AggreCut2[nrow(AggreCut2), 6] == AggreCut2[nrow(AggreCut2)-1, 6] ){
  Dupe = 1 
}


#Count Dupes
if (Dupe == 1){
  PreV  <- AggreCut2[1,6]
  AggreCut2[1,10] <- Dupe
  for (n in 2:nrow(AggreCut2)){
    NowV <- AggreCut2[n,6]
    if (PreV==NowV){
      Dupe = Dupe +  1
    } else {
      Dupe = 1
    }
    PreV <- NowV
    AggreCut2[n,10] <- Dupe
  }
}

DupeStart <- nrow(AggreCut2)-Dupe

#Dupe Penalty
for (n in DupeStart:nrow(AggreCut2)){
  PreVHR <- AggreCut2[n-1,7] - 0.1
  PreVHR <- max(PreVHR, 1)
  PreVP <- AggreCut2[n-1,9] - 0.1
  PreVP <- max(PreVP, 0)
  AggreCut2[n,7] <- PreVHR
  AggreCut2[n,9] <- PreVP
}


# output data
write.csv(AggreCut2, file = paste0("AggreCut", Gene_set,cutter, ".csv"))

#Plot LOCC Graph
{
  df.cut <- subset(AggreCut2, select = c(7,8))
  df.cut <- melt(df.cut, id.vars = "Fraction" , na.rm =  TRUE)
  df.cut2 <- subset(AggreCut2, select = c(8,9))
  df.cut2$logp <- (df.cut2$logp/2)
  df.cut2 <- melt(df.cut2, id.vars = "Fraction" , na.rm =  TRUE)
  df.cut2$variable <- "-Log (p value)"
  df.cut2 <- rbind(df.cut2, df.cut)
  df.cut2 <- df.cut2[df.cut2$Fraction > 0.05,]
  df.cut2 <- df.cut2[df.cut2$Fraction < 0.95,]
  df.cut2$variable <- factor(df.cut2$variable, levels = c("HR", "-Log (p value)"))
  p <- ggplot(df.cut2, aes(x = Fraction, y = value,color = variable)) +geom_point() + geom_line(aes(color = variable)) + scale_color_manual(values=c("black", "orange")) + ggtitle(paste0(Gene_set, " TCGA ", TCGA_Cancer_ID ," LOCC")) + scale_y_continuous(
    # Features of the first axis
    name = "Hazard Ratio (HR)",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
  ) + scale_x_continuous( name = "Fraction in High Activity Group", limits = c(0, 1), expand = c(0, 0)) + geom_hline(aes(yintercept=1), color = "red") + geom_hline(aes(yintercept=0.65),  color = "green") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30))  + labs(colour="Legend",x="xxx",y="yyy")
  
  
  ggsave(paste0(Gene_set,"_LOCC.pdf"), plot=p)
  # dev.off
  
  
  
  #Cut at best cutpoint
  {
    
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:50))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10:22))
    colnames(PtSurv) <- c("time", "status",'a','b','c','d','e','f','g','h','i', 'j', 'k','l',"score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv$time <- as.numeric(PtSurv$time)
    PtSurv[is.na(PtSurv)] <- 0
    # PtSurv <- LiziSurv2
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    res.cut <- surv_cutpoint(PtSurv, time = "time", event = "Status",
                             variables = "score", minprop = 0.1)
    summary(res.cut)
    # PtSurv$score[as.numeric(PtSurv$score) > res.cut$score[4]] <- as.character("High_Score")
    # PtSurv$score[(PtSurv$score) <= res.cut$score[4]] <- as.character("Low_Score")
    
  }
  # Print a summary of `res.cut`
  summary(res.cut)
  
  # Open a new device
  dev.new(width = 200, height = 200, unit = "px")
  
  # Plot the survival curve with risk table
  r <- survfit2(Surv(time, Status) ~ score, data = PtSurv) 
  t <- ggsurvplot(r, data = PtSurv, risk.table = TRUE,
                  size = 1.2,
                  fontsize = 7.5,
                  risk.table.y.text = FALSE,
                  tables.y.text = FALSE,
                  risk.table.height = 0.35,
                  font.title = c(16, "bold", "black"),
                  font.subtitle = c(16, "bold", "black"),
                  font.caption = c(16, "bold", "black"),
                  font.x = c(16, "bold", "black"),
                  font.y = c(16, "bold", "black"),
                  font.tickslab = c(16, "bold", "black")) 
  
  # Add themes to the plot and table
  t$plot <- t$plot + theme(legend.title = element_text(size = 14, color = "black", face = "bold"),
                           legend.text = element_text(size = 14, color = "black", face = "bold"),
                           axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                           axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                           axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                           axis.title.y = element_text(size = 20, color = "black", face = "bold"))
  t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                             axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.x = element_text(size = 20, color = "black", face = "bold"),
  )
  
  # Save the plot as a pdf
  pdf(paste0(Gene_set, "survplot.pdf"))
  print(t, newpage = FALSE)
  dev.off()
  
  #ExtractAIC
  # colnames(PtSurv) <- c("time", "status", 'a','b','c','d','e','f','g','h','i','j', 'k','l',score','Status')
  PtSurv$a <- as.numeric(PtSurv$a)
  PtSurv$b <- as.numeric(PtSurv$b)
  PtSurv$c <- as.numeric(PtSurv$c)
  PtSurv$d <- as.numeric(PtSurv$d)
  PtSurv$e <- as.numeric(PtSurv$e)
  PtSurv$f <- as.numeric(PtSurv$f)
  PtSurv$g <- as.numeric(PtSurv$g)
  PtSurv$h <- as.numeric(PtSurv$h)
  PtSurv$i <- as.numeric(PtSurv$i)
  PtSurv$j <- as.numeric(PtSurv$j)
  PtSurv$k <- as.numeric(PtSurv$k)
  PtSurv$l <- as.numeric(PtSurv$l)
  

  
  fit <- coxph(Surv(time,Status) ~ a + b + c + d  + f  + j + k +l, data = PtSurv)
  
  SumFit <- summary(fit)
  SumFit
  extractAIC(fit)
  newAIC <- as.numeric(extractAIC(fit)) 
  
  fit <- coxph(Surv(time,Status) ~ a + b + c + d + e + f + g + h + i + j + k + l, data = PtSurv)
  
  SumFit <- summary(fit)
  SumFit
  extractAIC(fit)
  oldAIC <- extractAIC(fit) 
}
}



#8-gene Modified RISK Score Analysis 
{
  NGS9 <- data.frame(c("ANXA10",
                       "CYP2C9",
                       "HMMR",
                       "KIF20A",
                       "LAPTM4B",
                       "LCAT",
                       "LECT2",
                       "MAGEA6",
                       "RDH16",
                       "SPP1",
                       "TPX2",
                       "TTK"))
  colnames(NGS9) <- "V1"
  Gene_set <- "RISK New Score"
  {
    PtExpressionN2 <- PtExpression[PtExpression$Hugo_Symbol %in% NGS9$V1,]
    PtExpressionN2 <- na.omit(PtExpressionN2)
    PtExpressionN2[,-c(1)] <- log2(PtExpressionN2[,-c(1)]+1)
    for (p in 3:ncol(PtExpressionN2)){
      PtExpressionN2[13,p]   <- (-0.0627 * PtExpressionN2[1,p]) + ( -0.02365*PtExpressionN2[2,p]) + (0.0565*PtExpressionN2[3,p]) +(0.0386*PtExpressionN2[4,p])  + (-0.02365*PtExpressionN2[6,p])+(0.05124*PtExpressionN2[10,p])+(0.05377*PtExpressionN2[11,p])+(0.09788*PtExpressionN2[12,p])
    }
    PtExpressionN2[nrow(PtExpressionN2),1] <- "Mean"
    
    PtN2Score <- PtExpressionN2
    PtN2Score <- PtExpressionN2[nrow(PtExpressionN2),]
    PtN2Score <- subset(PtN2Score, select = c(-2))
    
    PtN2Score <- rbind(colnames(PtN2Score), PtN2Score)
    PtN2Score2 <- t(PtN2Score)
    PtN2Score2 <- PtN2Score2[2:nrow(PtN2Score2),]
    colnames(PtN2Score2)[1] <- c("PATIENT_ID")
    colnames(PtN2Score2)[1:2] <- c("PATIENT_ID", "N2_Score")
    
    PtN2Score2[,1] <- substr(PtN2Score2[,1],1,nchar(PtN2Score2[,1])-3)
    PtN2Score2[,1] <- str_replace(PtN2Score2[,1], "\\.", "-")
    PtN2Score2[,1] <- str_replace(PtN2Score2[,1], "\\.", "-")
    
  }
  
  {
    cutter <- 22
    PtClinical4 <- merge(PtClinical, PtN2Score2,  by = "PATIENT_ID", .keep_all= TRUE)  
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    
    colnames(PtSurv) <- c("time", "status", "score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv[is.na(PtSurv)] <- 0
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    NScoreValue <- PtSurv[cutter,3]
    PtSurv[1:cutter,3] <- "High_RISK"
    PtSurv[cutter:(nrow(PtSurv)+1),3] <- "Low_RISK"
    PtSurv <- filter(PtSurv, score == "High_RISK" | score == "Low_RISK")
    
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    CutTrial  <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, NScoreValue)
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
  }
  
  # Initialize a data frame to store the results
  AggreCut <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c(colnames(CutTrial)))))
  # Initialize a variable to track duplicates
  Dupe <- 0
  
  # Loop over different cut-off points
  for (cutter in 2:(nrow(PtSurv)-1)){  
    # Repeat the same process with different cut-off points
    PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
    PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
    colnames(PtSurv) <- c("time", "status", "score")
    PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
    PtSurv <- na.omit(PtSurv)
    PtSurv$score <- as.numeric(PtSurv$score)
    PtSurv[is.na(PtSurv)] <- 0
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    GeneScoreValue <- PtSurv[cutter,3]
    PtSurv[1:cutter,3] <- "High_Gene"
    PtSurv[cutter:(nrow(PtSurv)),3] <- "Low_Gene"
    PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
    
    # Calculate survival difference for the current cut-off point
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    CutTrial <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, GeneScoreValue)
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
    
    # Add the results to the AggreCut data frame
    AggreCut <- rbind(AggreCut, CutTrial)
  }
  
  # Now we process the AggreCut data frame
  AggreCut <- AggreCut[order(as.numeric(AggreCut[,1]), decreasing = FALSE),]
  AggreCut2 <- AggreCut[AggreCut$groups == "score=High_Gene",]
  AggreCut2 <- AggreCut2[order(as.numeric(AggreCut2[,3]), decreasing = FALSE),]
  AggreCut2$Fraction <- AggreCut2$Freq/nrow(AggreCut2)
  AggreCut2$logp <- -log(AggreCut2$CutTrial.pvalue, base = 10)
  
  #Check for Duplicates
  if (AggreCut2[nrow(AggreCut2), 6] == AggreCut2[nrow(AggreCut2)-1, 6] ){
    Dupe = 1 
  }
  
  
  #Count Dupes
  if (Dupe == 1){
    PreV  <- AggreCut2[1,6]
    AggreCut2[1,10] <- Dupe
    for (n in 2:nrow(AggreCut2)){
      NowV <- AggreCut2[n,6]
      if (PreV==NowV){
        Dupe = Dupe +  1
      } else {
        Dupe = 1
      }
      PreV <- NowV
      AggreCut2[n,10] <- Dupe
    }
  }
  
  DupeStart <- nrow(AggreCut2)-Dupe
  
  #Dupe Penalty
  for (n in DupeStart:nrow(AggreCut2)){
    PreVHR <- AggreCut2[n-1,7] - 0.1
    PreVHR <- max(PreVHR, 1)
    PreVP <- AggreCut2[n-1,9] - 0.1
    PreVP <- max(PreVP, 0)
    AggreCut2[n,7] <- PreVHR
    AggreCut2[n,9] <- PreVP
  }
  
  
  # output data
  write.csv(AggreCut2, file = paste0("AggreCut", Gene_set,cutter, ".csv"))
  
  #Plot LOCC Graph
  {
    df.cut <- subset(AggreCut2, select = c(7,8))
    df.cut <- melt(df.cut, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2 <- subset(AggreCut2, select = c(8,9))
    df.cut2$logp <- (df.cut2$logp/2)
    df.cut2 <- melt(df.cut2, id.vars = "Fraction" , na.rm =  TRUE)
    df.cut2$variable <- "-Log (p value)"
    df.cut2 <- rbind(df.cut2, df.cut)
    df.cut2 <- df.cut2[df.cut2$Fraction > 0.05,]
    df.cut2 <- df.cut2[df.cut2$Fraction < 0.95,]
    df.cut2$variable <- factor(df.cut2$variable, levels = c("HR", "-Log (p value)"))
    p <- ggplot(df.cut2, aes(x = Fraction, y = value,color = variable)) +geom_point() + geom_line(aes(color = variable)) + scale_color_manual(values=c("black", "orange")) + ggtitle(paste0(Gene_set, " TCGA ", TCGA_Cancer_ID ," LOCC")) + scale_y_continuous(
      # Features of the first axis
      name = "Hazard Ratio (HR)",
      # Add a second axis and specify its features
      sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
    ) + scale_x_continuous( name = "Fraction in High Activity Group", limits = c(0, 1), expand = c(0, 0)) + geom_hline(aes(yintercept=1), color = "red") + geom_hline(aes(yintercept=0.65),  color = "green") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30))  + labs(colour="Legend",x="xxx",y="yyy")
    
    
    ggsave(paste0(Gene_set,"_LOCC.pdf"), plot=p)
    # dev.off
    
    
    
    #Cut at best cutpoint
    {
      
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      # PtSurv <- LiziSurv2
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      res.cut <- surv_cutpoint(PtSurv, time = "time", event = "Status",
                               variables = "score", minprop = 0.1)
      summary(res.cut)
      PtSurv$score[as.numeric(PtSurv$score) > res.cut$score[4]] <- as.character("High_Score")
      PtSurv$score[(PtSurv$score) <= res.cut$score[4]] <- as.character("Low_Score")
      
    }
    # Print a summary of `res.cut`
    summary(res.cut)
    
    # Open a new device
    dev.new(width = 200, height = 200, unit = "px")
    
    # Plot the survival curve with risk table
    r <- survfit2(Surv(time, Status) ~ score, data = PtSurv) 
    t <- ggsurvplot(r, data = PtSurv, risk.table = TRUE,
                    size = 1.2,
                    fontsize = 7.5,
                    risk.table.y.text = FALSE,
                    tables.y.text = FALSE,
                    risk.table.height = 0.35,
                    font.title = c(16, "bold", "black"),
                    font.subtitle = c(16, "bold", "black"),
                    font.caption = c(16, "bold", "black"),
                    font.x = c(16, "bold", "black"),
                    font.y = c(16, "bold", "black"),
                    font.tickslab = c(16, "bold", "black")) 
    
    # Add themes to the plot and table
    t$plot <- t$plot + theme(legend.title = element_text(size = 14, color = "black", face = "bold"),
                             legend.text = element_text(size = 14, color = "black", face = "bold"),
                             axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                             axis.title.y = element_text(size = 20, color = "black", face = "bold"))
    t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                               axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                               axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    )
    
    # Save the plot as a pdf
    pdf(paste0(Gene_set, "survplot.pdf"))
    print(t, newpage = FALSE)
    dev.off()
    
    
    # Perform Cox proportional hazards model
    fit <- coxph(Surv(time,Status) ~ score, data = PtSurv)
    
    # Compute and save summary of fit, survival difference, and survival fit
    SumFit <- summary(fit)
    Diffsurv <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    Fitsurv <- survfit(Surv(time, Status) ~ score, data = PtSurv)
    
    capture.output(SumFit, file = paste0(Gene_set,"SumFit.csv"))
    capture.output(Diffsurv, file = paste0(Gene_set,"Diffsurv.csv"))
    capture.output(Fitsurv, file = paste0(Gene_set,"Fitsurv.csv"))
    
    library(ROCR)
    #Calculate and Plot ROC Curve (non-time dependent)
    {
      #Prepare Dataset
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      PtSurv2 <- PtSurv
      PtSurv2$time <- PtSurv2$time + 0.01
      # Filter by time over 1 month
      # PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
      PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
      PtSurv2 <- rbind(PtSurv2, PtSurv3)
      
      #Time-dependent ROC Curve
      midsurv <- 1000000
      
      df.y <- data.frame(PtSurv$time, PtSurv$Status)
      df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
      colnames(df) <- c("predictions", "time", "status")
      for (o in 1:nrow(df)) {
        df[o,4] <- as.numeric(0)
        if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
          df[o,4] <- 1
        }
      }
      
      
      #Plot ROC Curve
      df <- subset(df, select = c(1,4))
      colnames(df) <- c("predictions", "labels")
      
      
      pred <- prediction(df$predictions, df$labels)
      perf <- performance(pred,"tpr","fpr")
      line = data.frame(1:nrow(df),1:nrow(df))
      line = line/nrow(df)
      colnames(line) <- c("False Positive Rate", "True Positive Rate")
      
      
      pdf(paste0(Gene_set, "AUC.pdf"))
      plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(a=0, b= 1, col = "red")
      
      
      dev.off()
      auc_ROCR <- performance(pred, measure = "auc")
      auc_ROCR <- auc_ROCR@y.values[[1]]
    }
    
    #Calculate and Plot ROC Curve (1 year)
    {
      #Prepare Dataset
      PtSurvivalCurve <- subset(PtClinical4, select = c(1,30:38))
      PtSurv <- subset(PtSurvivalCurve, select = c(3,2,10))
      colnames(PtSurv) <- c("time", "status", "score")
      PtSurv <- mutate(PtSurv, Status = dplyr::recode(PtSurv$status, "1:DECEASED" = 1, "0:LIVING" = 0))
      PtSurv <- na.omit(PtSurv)
      PtSurv$score <- as.numeric(PtSurv$score)
      PtSurv$time <- as.numeric(PtSurv$time)
      PtSurv[is.na(PtSurv)] <- 0
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      PtSurv2 <- PtSurv
      PtSurv2$time <- PtSurv2$time + 0.01
      #Filter by time over 1 month
      PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
      PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
      PtSurv2 <- rbind(PtSurv2, PtSurv3)
      
      #Time-dependent ROC Curve
      midsurv <- 12
      
      df.y <- data.frame(PtSurv$time, PtSurv$Status)
      df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
      colnames(df) <- c("predictions", "time", "status")
      for (o in 1:nrow(df)) {
        df[o,4] <- as.numeric(0)
        if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
          df[o,4] <- 1
        }
      }
      
      
      #Plot ROC Curve
      df <- subset(df, select = c(1,4))
      colnames(df) <- c("predictions", "labels")
      
      
      pred <- prediction(df$predictions, df$labels)
      perf <- performance(pred,"tpr","fpr")
      line = data.frame(1:nrow(df),1:nrow(df))
      line = line/nrow(df)
      colnames(line) <- c("False Positive Rate", "True Positive Rate")
      
      
      pdf(paste0(Gene_set, "1yr_AUC.pdf"))
      plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(a=0, b= 1, col = "red")
      
      
      dev.off()
      auc_ROCR <- performance(pred, measure = "auc")
      auc_ROCR <- auc_ROCR@y.values[[1]]
    }
    
    
    #LOCC Calculations
    {
      
      #Calculate Parameters
      length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])
      max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
      #Create Dataframe and Save numbers
      DF.qkscore <- data.frame("New Risk", max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])/nrow(AggreCut2), max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]),res.cut$cutpoint$cutpoint, AggreCut2$HR[AggreCut2$logp == (max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]))], auc_ROCR)
      colnames(DF.qkscore) <- c("Gene", "-log (p value)", "Percentage highly significantly", "Highest HR", "Lowest HR", "Cut", "significant HR", "auc_ROCR")
      write.csv(DF.qkscore, paste0(Gene_set, "_LOCC_Scorev3.csv"))
      
    }
    
  }
  
}

#Sample Testing
#E2F1
j <- grep("E2F1", PtExpressionZ$Hugo_Symbol)

for (s in 1:100){
  setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018"))
  PtClinical <- read.table(file = "data_clinical_patient.txt", sep = "\t",quote = "", header = TRUE)  
  PtClinical = read.table(file = "data_clinical_patient.txt", sep = "\t",quote = "", header = TRUE)
  PtClinicalv2 <- PtClinical[sample(1:nrow(PtClinical), 186, replace=FALSE),]
  PtClinicalv3 <- PtClinical[!(PtClinical$PATIENT_ID %in% PtClinicalv2$PATIENT_ID),]
  PtClinical <- PtClinicalv3
  setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018/Sample"))
  Loccfunc2(j,s)
}

#Sample Testing
#MTHFR
j <- grep("MTHFR", PtExpressionZ$Hugo_Symbol)

for (s in 1:100){
  setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018"))
  PtClinical <- read.table(file = "data_clinical_patient.txt", sep = "\t",quote = "", header = TRUE)  
  PtClinical = read.table(file = "data_clinical_patient.txt", sep = "\t",quote = "", header = TRUE)
  PtClinicalv2 <- PtClinical[sample(1:nrow(PtClinical), 186, replace=FALSE),]
  PtClinicalv3 <- PtClinical[!(PtClinical$PATIENT_ID %in% PtClinicalv2$PATIENT_ID),]
  PtClinical <- PtClinicalv3
  setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018/Sample"))
  Loccfunc2(j,s)
}

#Sample Testing
#KIF20A
j <- grep("KIF20A", PtExpressionZ$Hugo_Symbol)
for (s in 1:100){
  setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018"))
  PtClinical <- read.table(file = "data_clinical_patient.txt", sep = "\t",quote = "", header = TRUE)  
  PtClinical = read.table(file = "data_clinical_patient.txt", sep = "\t",quote = "", header = TRUE)
  PtClinicalv2 <- PtClinical[sample(1:nrow(PtClinical), 186, replace=FALSE),]
  PtClinicalv3 <- PtClinical[!(PtClinical$PATIENT_ID %in% PtClinicalv2$PATIENT_ID),]
  PtClinical <- PtClinicalv3
  setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018/Sample"))
  Loccfunc2(j,s)
}

#Create E2F1 dataframe with all gene LOCC scores

Gene_set <- "E2F1"
setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018/Sample"))
{
  if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set, "1_LOCC_Scorev3.csv")))){
    
    df.score <- read.csv(file = paste0(Gene_set, "1_LOCC_Scorev3.csv"))
  }
  df.score2 <- df.score
  df.score2 <- df.score2[-1,]
  file.size(paste0(Gene_set, "1_LOCC_Scorev3.csv"))
  
  #Combine Dataframes
  x <- 10
  for (y in 1:x){
    
    if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set,y, "_LOCC_Scorev3.csv")))){
      
      df.score <- read.csv(file = paste0(Gene_set,y, "_LOCC_Scorev3.csv"))
      df.score2 <- rbind(df.score2, df.score)
    }
  }
  #Calculate LOCC 
  for (y in 1:nrow(df.score2)){
    if (as.numeric(df.score2[y,9]) > 1){
      df.score2[y,11] <- -(df.score2[y,3]*df.score2[y,4]*df.score2[y,9])
    } else {
      df.score2[y,11] <- (df.score2[y,3]*df.score2[y,4]/df.score2[y,9])  
    }
  }
  colnames(df.score2)[11] <- "LOCC"
  write.csv(df.score2, file = paste0(Gene_set,"_totalLoccScore.csv"))
}

#Create MTHFR dataframe with all gene LOCC scores

Gene_set <- "MTHFR"
setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018/Sample"))
{
  if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set, "1_LOCC_Scorev3.csv")))){
    
    df.score <- read.csv(file = paste0(Gene_set, "1_LOCC_Scorev3.csv"))
  }
  df.score2 <- df.score
  df.score2 <- df.score2[-1,]
  file.size(paste0(Gene_set, "1_LOCC_Scorev3.csv"))
  
  #Combine Dataframes
  x <- 10
  for (y in 1:x){
    
    if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set,y, "_LOCC_Scorev3.csv")))){
      
      df.score <- read.csv(file = paste0(Gene_set,y, "_LOCC_Scorev3.csv"))
      df.score2 <- rbind(df.score2, df.score)
    }
  }
  #Calculate LOCC 
  for (y in 1:nrow(df.score2)){
    if (as.numeric(df.score2[y,9]) > 1){
      df.score2[y,11] <- -(df.score2[y,3]*df.score2[y,4]*df.score2[y,9])
    } else {
      df.score2[y,11] <- (df.score2[y,3]*df.score2[y,4]/df.score2[y,9])  
    }
  }
  colnames(df.score2)[11] <- "LOCC"
  write.csv(df.score2, file = paste0(Gene_set,"_totalLoccScore.csv"))
}

#Create KIF20A dataframe with all gene LOCC scores

Gene_set <- "KIF20A"
setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018/Sample"))
{
  if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set, "1_LOCC_Scorev3.csv")))){
    
    df.score <- read.csv(file = paste0(Gene_set, "1_LOCC_Scorev3.csv"))
  }
  df.score2 <- df.score
  df.score2 <- df.score2[-1,]
  file.size(paste0(Gene_set, "1_LOCC_Scorev3.csv"))
  
  #Combine Dataframes
  x <- 10
  for (y in 1:x){
    
    if (Gene_set != "_" & !is.na(file.size(paste0(Gene_set,y, "_LOCC_Scorev3.csv")))){
      
      df.score <- read.csv(file = paste0(Gene_set,y, "_LOCC_Scorev3.csv"))
      df.score2 <- rbind(df.score2, df.score)
    }
  }
  #Calculate LOCC 
  for (y in 1:nrow(df.score2)){
    if (as.numeric(df.score2[y,9]) > 1){
      df.score2[y,11] <- -(df.score2[y,3]*df.score2[y,4]*df.score2[y,9])
    } else {
      df.score2[y,11] <- (df.score2[y,3]*df.score2[y,4]/df.score2[y,9])  
    }
  }
  colnames(df.score2)[11] <- "LOCC"
  write.csv(df.score2, file = paste0(Gene_set,"_totalLoccScore.csv"))
}


#LIRI-JP E2F1
{
  
  Mutation1 <- "E2F1"
  Gene_set <- Mutation1
  setwd(paste0("C:/Users/soldi/Documents/R/CTRP2/", TCGA_Cancer_ID,"_tcga_pan_can_atlas_2018"))
  #Pick Gene Mutation Data
  if (all(Mutation1 != "")) {
    
    #TPM of E2F1 from LIRI-JP Dataset
   PtSurv <- read.csv(file = "LIRIJP_E2F1_Ranked.csv", header = T)
    
    # Order by 'score'
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score), decreasing = TRUE), ]
    
    # Save the gene score value at the cut-off point
    cutter <- 100
    GeneScoreValue <- PtSurv[cutter,3]
    
    # Assign high and low gene score groups
    PtSurv[1:cutter,3] <- "High_Gene"
    PtSurv[cutter:(nrow(PtSurv)+1),3] <- "Low_Gene"
    
    # Filter out only the high and low gene score groups
    PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
    
    # Calculate survival difference using the survdiff function
    CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
    
    # Create a dataframe with the results
    CutTrial <- data.frame(CutTrial$pvalue, CutTrial$n, CutTrial$obs, CutTrial$exp, GeneScoreValue)
    
    # Calculate HR
    HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
    CutTrial$HR <- HR
    
    # Initialize AggreCut dataframe with the results from the first cut-off point
    AggreCut <- CutTrial
    
    # Initialize a data frame to store the results
    AggreCut <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c(colnames(CutTrial)))))
    # Initialize a variable to track duplicates
    Dupe <- 0
    
    # Loop over different cut-off points
    for (cutter in 2:(nrow(PtSurv)-1)){  
      # Repeat the same process with different cut-off points
      PtSurv <- read.csv(file = "LIRIJP_E2F1_Ranked.csv", header = T)
      PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
      GeneScoreValue <- PtSurv[cutter,3]
      PtSurv[1:cutter,3] <- "High_Gene"
      PtSurv[cutter:(nrow(PtSurv)),3] <- "Low_Gene"
      PtSurv <- filter(PtSurv, score == "High_Gene" | score == "Low_Gene")
      
      # Calculate survival difference for the current cut-off point
      CutTrial <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
      CutTrial <- data.frame(CutTrial$pvalue,CutTrial$n, CutTrial$obs,CutTrial$exp, GeneScoreValue)
      HR <- as.numeric(CutTrial[1,4])*as.numeric(CutTrial[2,5])/as.numeric(CutTrial[1,5])/as.numeric(CutTrial[2,4])
      CutTrial$HR <- HR
      
      # Add the results to the AggreCut data frame
      AggreCut <- rbind(AggreCut, CutTrial)
    }
    
    # Now we process the AggreCut data frame
    AggreCut <- AggreCut[order(as.numeric(AggreCut[,1]), decreasing = FALSE),]
    AggreCut2 <- AggreCut[AggreCut$groups == "score=High_Gene",]
    AggreCut2 <- AggreCut2[order(as.numeric(AggreCut2[,3]), decreasing = FALSE),]
    AggreCut2$Fraction <- AggreCut2$Freq/nrow(AggreCut2)
    AggreCut2$logp <- -log(AggreCut2$CutTrial.pvalue, base = 10)
    
    #Check for Duplicates
    if (AggreCut2[nrow(AggreCut2), 6] == AggreCut2[nrow(AggreCut2)-1, 6] ){
      Dupe = 1 
    }
    
    
    #Count Dupes
    if (Dupe == 1){
      PreV  <- AggreCut2[1,6]
      AggreCut2[1,10] <- Dupe
      for (n in 2:nrow(AggreCut2)){
        NowV <- AggreCut2[n,6]
        if (PreV==NowV){
          Dupe = Dupe +  1
        } else {
          Dupe = 1
        }
        PreV <- NowV
        AggreCut2[n,10] <- Dupe
      }
    }
    
    DupeStart <- nrow(AggreCut2)-Dupe
    
    #Dupe Penalty
    for (n in DupeStart:nrow(AggreCut2)){
      PreVHR <- AggreCut2[n-1,7] - 0.1
      PreVHR <- max(PreVHR, 1)
      PreVP <- AggreCut2[n-1,9] - 0.1
      PreVP <- max(PreVP, 0)
      AggreCut2[n,7] <- PreVHR
      AggreCut2[n,9] <- PreVP
    }
    
    
    # output data
    write.csv(AggreCut2, file = paste0("LIRI_AggreCut", Gene_set,cutter, ".csv"))
    
    #Plot LOCC Graph
    {
      df.cut <- subset(AggreCut2, select = c(7,8))
      df.cut <- melt(df.cut, id.vars = "Fraction" , na.rm =  TRUE)
      df.cut2 <- subset(AggreCut2, select = c(8,9))
      df.cut2$logp <- (df.cut2$logp/2)
      df.cut2 <- melt(df.cut2, id.vars = "Fraction" , na.rm =  TRUE)
      df.cut2$variable <- "-Log (p value)"
      df.cut2 <- rbind(df.cut2, df.cut)
      df.cut2 <- df.cut2[df.cut2$Fraction > 0.05,]
      df.cut2 <- df.cut2[df.cut2$Fraction < 0.95,]
      df.cut2$variable <- factor(df.cut2$variable, levels = c("HR", "-Log (p value)"))
      p <- ggplot(df.cut2, aes(x = Fraction, y = value,color = variable)) +geom_point() + geom_line(aes(color = variable)) + scale_color_manual(values=c("black", "orange")) + ggtitle(paste0(Gene_set, " LIRI ", TCGA_Cancer_ID ," LOCC")) + scale_y_continuous(
        # Features of the first axis
        name = "Hazard Ratio (HR)",
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
      ) + scale_x_continuous( name = "Fraction in High Activity Group", limits = c(0, 1), expand = c(0, 0)) + geom_hline(aes(yintercept=1), color = "red") + geom_hline(aes(yintercept=0.65),  color = "green") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30))  + labs(colour="Legend",x="xxx",y="yyy")
      
      
      ggsave(paste0(Gene_set,"LIRI_LOCC.pdf"), plot=p)
      # dev.off
    }
    
    PtSurv <- read.csv(file = "LIRIJP_E2F1_Ranked.csv", header = T)
    # PtSurv <- LiziSurv2
    PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
    
    {
      #Check if gene is expressed
      if (nrow(PtSurv) > 100 & PtSurv[nrow(PtSurv)/5,3] != PtSurv[nrow(PtSurv)-1,3]){
        
        
        
        #Cut at best cutpoint
        {
          
          PtSurv <- read.csv(file = "LIRIJP_E2F1_Ranked.csv", header = T)
          PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
          res.cut <- surv_cutpoint(PtSurv, time = "time", event = "Status",
                                   variables = "score", minprop = 0.1)
          summary(res.cut)
          PtSurv$score[as.numeric(PtSurv$score) > res.cut$score[4]] <- as.character("High_Score")
          PtSurv$score[(PtSurv$score) <= res.cut$score[4]] <- as.character("Low_Score")
          
        }
        # Print a summary of `res.cut`
        summary(res.cut)
        
        # Open a new device
        dev.new(width = 200, height = 200, unit = "px")
        
        # Plot the survival curve with risk table
        r <- survfit2(Surv(time, Status) ~ score, data = PtSurv) 
        t <- ggsurvplot(r, data = PtSurv, risk.table = TRUE,
                        size = 1.2,
                        fontsize = 7.5,
                        risk.table.y.text = FALSE,
                        tables.y.text = FALSE,
                        risk.table.height = 0.35,
                        font.title = c(16, "bold", "black"),
                        font.subtitle = c(16, "bold", "black"),
                        font.caption = c(16, "bold", "black"),
                        font.x = c(16, "bold", "black"),
                        font.y = c(16, "bold", "black"),
                        font.tickslab = c(16, "bold", "black")) 
        
        # Add themes to the plot and table
        t$plot <- t$plot + theme(legend.title = element_text(size = 14, color = "black", face = "bold"),
                                 legend.text = element_text(size = 14, color = "black", face = "bold"),
                                 axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                 axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                                 axis.title.x = element_text(size = 20, color = "black", face = "bold"),
                                 axis.title.y = element_text(size = 20, color = "black", face = "bold"))
        t$table <- t$table + theme(plot.title = element_text(size = 16, color = "black", face = "bold"), 
                                   axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                                   axis.title.x = element_text(size = 20, color = "black", face = "bold"),
        )
        
        # Save the plot as a pdf
        pdf(paste0(Gene_set, "LIRI_survplot.pdf"))
        print(t, newpage = FALSE)
        dev.off()
        
        # Perform Cox proportional hazards model
        fit <- coxph(Surv(time,Status) ~ score, data = PtSurv)
        
        # Compute and save summary of fit, survival difference, and survival fit
        SumFit <- summary(fit)
        Diffsurv <- survdiff(Surv(time, Status) ~ score, data = PtSurv) 
        Fitsurv <- survfit(Surv(time, Status) ~ score, data = PtSurv)
        
        capture.output(SumFit, file = paste0(Gene_set,"LIRI_SumFit.csv"))
        capture.output(Diffsurv, file = paste0(Gene_set,"LIRI_Diffsurv.csv"))
        capture.output(Fitsurv, file = paste0(Gene_set,"LIRI_Fitsurv.csv"))
        
        library(ROCR)
        #Calculate and Plot ROC Curve
        {
          #Prepare Dataset
          PtSurv <- read.csv(file = "LIRIJP_E2F1_Ranked.csv", header = T)
          PtSurv <- PtSurv[order(as.numeric(PtSurv$score),decreasing = TRUE), ]
          PtSurv2 <- PtSurv
          PtSurv2$time <- PtSurv2$time + 0.01
          # Filter by time over 1 month
          # PtSurv2 <- PtSurv[PtSurv$time > 1, ] 
          PtSurv3 <- (PtSurv[PtSurv$time < 0.001 & PtSurv$Status == 1, ] )
          PtSurv2 <- rbind(PtSurv2, PtSurv3)
          
          #Time-dependent ROC Curve
          midsurv <- 1000000
          
          df.y <- data.frame(PtSurv$time, PtSurv$Status)
          df <- data.frame(PtSurv2$score, PtSurv2$time, PtSurv2$Status)
          colnames(df) <- c("predictions", "time", "status")
          for (o in 1:nrow(df)) {
            df[o,4] <- as.numeric(0)
            if (as.numeric(df[o,2]) < midsurv & df[o,3] == 1){
              df[o,4] <- 1
            }
          }
          
          
          #Plot ROC Curve
          df <- subset(df, select = c(1,4))
          colnames(df) <- c("predictions", "labels")
          
          
          pred <- prediction(df$predictions, df$labels)
          perf <- performance(pred,"tpr","fpr")
          line = data.frame(1:nrow(df),1:nrow(df))
          line = line/nrow(df)
          colnames(line) <- c("False Positive Rate", "True Positive Rate")
          
          
          pdf(paste0(Gene_set, "LIRI_AUC.pdf"))
          plot(perf,colorize=FALSE, font.size = 32,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
          abline(a=0, b= 1, col = "red")
          
          
          dev.off()
          auc_ROCR <- performance(pred, measure = "auc")
          auc_ROCR <- auc_ROCR@y.values[[1]]
        }
        
        
        #Select cutoff for activity
        {
          PtScore3 <- merge(PtGeneScore2, PtSpecificMutation, all.x = TRUE)
          bestcutoff <- res.cut$score
          
          write.csv(PtScore3, file = paste0(Gene_set,"LIRI_vsMut.csv"))
          
          
          bestcutoff <- res.cut$cutpoint$cutpoint
        }
        
        
        #Prepare Activity data
        {
          PtSurv <- read.csv(file = "LIRIJP_E2F1_Ranked.csv", header = T)
          df.activity <- data.frame(PtSurv)
          df.activity <- df.activity[order(as.numeric(df.activity[,3]), decreasing = TRUE),]
          df.activity <- df.activity %>% 
            mutate(row_id=row_number())
          df.activity$N_Score <- as.numeric(df.activity$score)
          df.activity <- df.activity[order(df.activity$score, decreasing = TRUE),]
        }
        
        #Plot Activity Graph
        {
          q <- ggplot(df.activity, aes(x = row_id, y = N_Score)) + geom_point(size = 1.8)  + ggtitle(paste0(Gene_set," LIRI ", TCGA_Cancer_ID,  " Curve")) + scale_y_continuous(
            #+ scale_color_manual(values=c("cornflower blue")) 
            # Features of the first axis
            name = paste0(Gene_set, "Score"), expand = c(0, 0)
            # Add a second axis and specify its features
            #sec.axis = sec_axis( trans=~.*2, name="-Log (p value)")
          ) + scale_x_continuous( name = paste0("Ranking by ", Gene_set), expand = c(0, 1))  + geom_hline(aes(yintercept=bestcutoff),  color = "gray50") + scale_fill_discrete(labels=c('HR', 'p value')) + theme_classic()+ theme(legend.position="bottom", text=element_text(size=30), axis.text=element_text(size=30), plot.title = element_text(hjust = 0.5))  + labs(colour="Legend",x="xxx",y="yyy")
          
          # q
          ggsave(paste0(Gene_set,"LIRI_Activityv3.pdf"), plot=q)
          # graphics.off()
        }
        
        
        
        #LOCC Calculations
        {
          
          #Calculate Parameters
          length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])
          max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
          min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
          max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9])
          #Create Dataframe and Save numbers
          DF.qkscore <- data.frame(Mutation1, max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), length (AggreCut2$logp[AggreCut2$logp > as.numeric(2)])/nrow(AggreCut2), max(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]), min(AggreCut2$HR[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]),res.cut$cutpoint$cutpoint, Dupe, AggreCut2$HR[AggreCut2$logp == (max(AggreCut2$logp[AggreCut2$Fraction > 0.1 & AggreCut2$Fraction < 0.9]))], auc_ROCR)
          colnames(DF.qkscore) <- c("Gene", "-log (p value)", "Percentage highly significantly", "Highest HR", "Lowest HR", "Cut","dupe", "significant HR", "auc_ROCR")
          write.csv(DF.qkscore, paste0(Gene_set, "_LIRI_LOCC_Scorev3.csv"))
        }
      }
    }
  }
  
  graphics.off()
  
}
