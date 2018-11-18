#DruGen v1

#Authors:
#Robin Luo – bei.luo@mail.mcgill.ca
#Yuning Bie – yuning.bie@mail.mcgill.ca
#Bessie Luo – bei.x.luo@mail.mcgill.ca

#Made as part of PharmaHackathon – Nov 2018
#Team #4

#MIT-License
#Copyright 2018. Robin Luo, Yuning Bie, Bessie Luo.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED... Always refer to full conditions at https://opensource.org/licenses/MIT

#begin of code

#read all the data files
Expression <- read.delim("./data/sanger1018_brainarray_ensemblgene_rma.txt", row.names=1)
Expression_t <- t(Expression)
library(readxl)
Dose_response <- read_excel("./data/v17.3_fitted_dose_response.xlsx", na = "NA")
Cell_Lines <- read_excel("./data/Cell_Lines_Details.xlsx")
Screened_Compounds <- read_excel("./data/Screened_Compounds.xlsx")

#formatting data
Drug_zscore <- Dose_response[,c(3,5,13)]
for(row in 1:nrow(Drug_zscore))
{
  print("row "+row)
  Drug_zscore[row,1]<-paste("X", Drug_zscore[row,1], sep="")
}

write.csv(Drug_zscore,"./data/Drug_zscore.csv",sep=",")
write.csv(Expression_t,"./data/Expression_t.csv",sep=",")

length(intersect(Drug_zscore[,1],rownames(Expression_t)))

Expression_id <- cbind(rownames(Expression_t),Expression_t)
colnames(Expression_id)[1] <- "COSMIC_ID"

#function to do Drug-gene correlations
DG_corr <- function(DG_table)
{
  results <- c()
  for(c in 4:ncol(DG_table))
  {
    print(c)
    res <- cor(as.numeric(DG_table[,c]),as.numeric(DG_table[,3]))
    results <- c(results,res)
  }

  combined <- rbind(colnames(DG_table),c("NA",DG_table[1,2],"NA",results))
  combined2 <- t(combined)[c(-1,-3),]
  return(combined2[,2])
}


#loop through all the drug-gene data tables
corr_matrix <- c()
for(i in unique(Drug_zscore$DRUG_ID)) #i ended on 163, next time resume from i=163
{
  DrugX <- c()
  print(paste("testing drug ID:",i))
  DrugX <- Drug_zscore[which(Drug_zscore$DRUG_ID==i),]
  temp <- merge(DrugX, Expression_id, by="COSMIC_ID")
  res_temp <- DG_corr(temp)
  corr_matrix <- cbind(corr_matrix,res_temp)
}
header <- colnames(temp)

#write output
write.csv(corr_matrix,"./data/corr_matrix.csv",sep="")
write.csv(header,"./data/corr_matrix_header.csv",sep="")

#took some manual processing here to align headers
corr_matrix_processed <- read_excel("./data/corr_matrix_processed.xlsx", na = "NA")




#order the correlations and list the top5 positive and top5 negative predictors
top_gene_effects <- c()
for(col in 3:ncol(corr_matrix_processed))
{
  print(col)
  temp <- order(corr_matrix_processed[,col])
  res <- corr_matrix_processed[c(head(temp,n=5),tail(temp,n=5)),2]

  top_gene_effects <- cbind(top_gene_effects,res)

}

#output results
colnames(top_gene_effects) <- colnames(corr_matrix_processed)[-c(1,2)]
write.csv(top_gene_effects,"./data/top_gene_effects.csv")



