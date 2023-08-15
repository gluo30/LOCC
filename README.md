# README
This program contains code for running LOCC analysis on TCGA LIHC data.  Included in the folder is the TCGA LIHC dataset and the E2F1 expression data from LIRI-JP dataset.  

# Installation
Install R and R studio. Input Code into a new script.  

Install packages, which are located at the top of the code.  

Download the TCGA dataset from [https://www.cbioportal.org/datasets]

Liver Hepatocellular Carcinoma (TCGA, PanCancer Atlas)


Move folder (lihc_tcga_pan_can_atlas_2018) to suitable directory and set directory to that location.  

>#Set Directory
>
>setwd("C:/Users/soldi/Documents/R/CTRP2")

Copy and replace that old directory location to the new directory location for all of the code. 

Add the LIRIJP_E2F1_Ranked.csv to the folder (lihc_tcga_pan_can_atlas_2018).   

# Running the code

The Loccfunc function will calculate LOCC and AUC and save pictures of graphs related to LOCC and AUC for a gene of interest in the TCGA dataset.  The gene of interest will be located at j where j is the row number of the gene in the TCGA expression dataset.  You can query a gene of interest by replacing "E2F1" in grep("E2F1"...) with your gene of interest.

>j <- grep("E2F1", PtExpressionZ$Hugo_Symbol)

The LOCC analysis can will take a long time to run if it attempts analysis LOCC scores for every single gene. 
I would recommend not running this loop until all the other code is run.  If you # out these lines, you can run the rest of the code to make sure it works.  

>#Loop over all genes, will take a long time
>
>for (j in 1:20531){
  Loccfunc(j)
}

The other code include LOCC and AUC analysis of other gene sets (RISK/8-gene RISK) and samples of genes.  
