#########################################
#### Differential abundance analyses #### 
#########################################

library(ggplot2)
library(DESeq2)

## R script to create Figure S7

#input data

exp1_counts_order <- read.csv('micro_ASVs_exp1_counts_order.csv', row.names=1, header=TRUE) #ASV counts grouped at order level, DESeq normalises the data 
exp1_counts_order <- as.data.frame(exp1_counts_order)
metadata_counts <- read.csv('micro_metadata_exp1_all.csv', row.names=1, header=TRUE)
metadata_counts$treatment <- factor(metadata_counts$treatment)
metadata_counts$nutrients <- factor(metadata_counts$nutrients)
metadata_counts$herbivore <- factor(metadata_counts$herbivore)
metadata_counts$dominant <- factor(metadata_counts$dominant)


#list of treatments to compare 

comparison_treatment <- list(
  "AC-AO" = c("AC", "AO"),
  "NO-AO" = c("NO", "AO"),
  "NC-AO" = c("NC", "AO"),
  "AC-NO" = c("AC", "NO"),
  "AC-NC" = c("AC", "NC"),
  "NO-NC" = c("NO", "NC")
)


#differential avundance analysis using DESeq

dsq_exp1_treatment <- DESeqDataSetFromMatrix(countData = exp1_counts_order, colData = metadata_counts, design = ~ treatment)
dsq_exp1_treatment <- estimateSizeFactors(dsq_exp1_treatment)
dsq_exp1_treatment <- DESeq(dsq_exp1_treatment)
for (combo in names(comparison_treatment)) {
  results_treatment.l <- list()
  print(paste(combo))
  combination.v <- comparison_treatment[[combo]]
  results_treatment.r <- results(dsq_exp1_treatment, contrast = c("treatment", combination.v[1], combination.v[2]))
  results_treatment.df <- data.frame(results_treatment.r)
  results_treatment.l[[combo]] <- results_treatment.df
  write.table(x = results_treatment.df, file = paste0(analysis_dir, "/", "deseq_order_exp1_", combo, ".tsv"), append = F, quote = F, sep = "\t", row.names = T, col.names = NA)
}


#create histogram showing log2FC of orders present at 1% relative abundance in at least one of the samples and with DESeq p-value <= 0.05

exp1_dsq <- read.csv('micro_dsq_exp1_treatment.csv', header=TRUE)
exp1_dsq <- as.data.frame(exp1_dsq)
exp1_dsq$order <- factor(exp1_dsq$order)

plot <- ggplot(exp1_dsq, aes(x=log2FC, y=order, fill=treatment)) + 
  geom_bar(stat="identity", color="grey", position=position_dodge()) +
  theme_minimal()

plot + scale_fill_manual(values=c("steelblue1","darkgreen","palegreen3"))


#the histogram was further modified in Illustrator to create Figure S7

