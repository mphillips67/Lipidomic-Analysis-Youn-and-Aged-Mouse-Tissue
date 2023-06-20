setwd("~/Library/CloudStorage/Dropbox/Lipidomics_McReynolds/Kidney_Liver_Data")

#load packages
library(lipidr)   #if you are on mac, you need to download gfortran before installing this package - https://mac.r-project.org/tools/
library(readxl)
library(ggplot2)
library(ggvenn)
library(UpSetR)
library(dbplyr)
library(ggh4x)
library(ggrepel)

#read in data 
liv <- read_excel("Hinton_Kid-Liv_Lipdomics.xlsx", sheet=3) #read in excel sheet with relevant data

#only keep first column and sample columns. will need to be adjusted based on excel sheet
liv <- liv[,c(2,5:ncol(liv))]

#remove "bad" entries (eg. 13-Docosenamide)
liv <- subset(liv, grepl(":", Metabolite.name))

#log transform data
liv[,2:9] <- log2(liv[2:9])

####### MORE PROCESSING AND LIPIDR ANALYSIS ########

#data sets have multiple entires for the same lipids. In these cases, only keep entry with the highest values. This approach does keep duplicates where the rows are just NA. Those are removed at another stp.
liv <- liv %>%
  group_by(Metabolite.name) %>%
  slice_max(across(starts_with("Liver"))) #groups by metabolite name, then only keeps entry with the highest row totals across columns that start with "Liver"

##create tables for  for lipidR
liver_lipids <- na.omit(liv) #remove NAs
colnames(liver_lipids)[1] <- "lipids" #rename first column

nrow(liver_lipids %>%
       group_by(lipids) %>%
       filter(n() > 1)) #if vallue returned is 0, all duplicates were correctly consolidated

#remove rows with same values
same_values_rows <- which(apply(liver_lipids[, 2:8], 1, function(row) all(row == row[1])))
#turns out we have some columns where all entries have the same value. These have to be removed for pca analysis. This line identifies the columns (there is 1 for the liver data).
liver_lipids <- liver_lipids[-same_values_rows, ]


##test for lipid names that don't match lipidr naming system and fix them. to run lipidR, at least 50% of names need to match. 
temp <- liver_lipids$lipids
fix = sub(";(O\\d*)", "(\\1)", temp)
fix = sub(" O-", "O ", fix)
fix = sub(" P-", "P ", fix)
fix = sub("-FA", "/", fix)
fix = sub(";(.*)$", "(\\1)", fix)
fix <- sub("([^\\)])\\|", "\\1(x)|", fix)
fix <- sub("(.*\\(x\\))\\|(.*$)", "\\1|\\2(x2)", fix)
fix <- sub("PE-Cer", "PECer", fix)
fix <- sub("PI-Cer", "PICer", fix)
fix <- sub("LPE-N", "LPEN", fix)

liver_lipids$lipids <- fix #al fixes are done using regex. for a given data set, you will have to look at the names in your orignal table and figure out how to adjust them to the naming system lipidr accepts. 

all_match <- lipidr::annotate_lipids(liver_lipids[[1]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

##create sample annotation table
Sample <- colnames(liver_lipids[,2:ncol(liver_lipids)]) #create a vector with all sample names

SampleType <- c(rep("Young",4), rep("Old",4)) #create vector where each entry describes corresponding sample. 

liver_annotation <- cbind(Sample,SampleType)

write.table(liver_annotation, "Liver_annotation.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart annotation file

###create lipidR objects
liver_lipids_exp <- as_lipidomics_experiment(liver_lipids, normalized = TRUE, logged = TRUE)

#list unmatched lipids  and add sample annotation 
non_parsed_molecules(liver_lipids_exp)
liver_lipids_exp <- add_sample_annotation(liver_lipids_exp,  "Liver_annotation.txt")
liver_lipids_exp <- remove_non_parsed_molecules(liver_lipids_exp)

#pca
mvaresults <- mva(liver_lipids_exp, measure="Area", method="PCA")

plot_mva(mvaresults, color_by="SampleType",hotelling = FALSE, components = c(1,2)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) 

#compare the two groups - univariate analysis
de_old_young <- de_analysis(liver_lipids_exp, Young - Old )
sig <- significant_molecules(de_old_young, p.cutoff = 0.05, logFC.cutoff = 1) #list of sig results
write.table(sig, "Liver_sig_lipids.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart sig lipids
write.table(de_old_young, "Liver_de_results.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart sig lipids

plot_results_volcano(de_old_young, show.labels =TRUE) + xlim(-14,14) + theme(strip.text.x = element_text(size = 20)) + theme(plot.title = element_text(size = 25, face = "bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) + ylim(0,5) # + geom_text_repel(show.legend  = F, data=subset(de_old_young, max.overlaps = 33, adj.P.Val < 0.05 & (logFC > 1 | logFC < -1)))

#enrichment analysis
enrich_results <- lsea(
  de_old_young,
  rank.by = "logFC", min_size = 4, nperm = 100000
)

write.table(enrich_results[,c(1:4,6,8)], "Liver_lipidset.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out enrichment results

sig_lipidsets <- significant_lipidsets(enrich_results)

plot_enrichment(de_old_young, sig_lipidsets, annotation="class")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))

plot_enrichment(de_old_young, sig_lipidsets, annotation="length")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))


#clustering analysis
##all classes
#plot_heatmap(heart_lipids_exp, layout=  list(font = list(size = 8)) ) 

##enriched classes
enriched_classes <- subset(enrich_results, padj <= .05 & grepl("^Class", set))
enriched_classes <- enriched_classes$set
enriched_classes <- sub("Class_","",enriched_classes)
enriched_heatmap <-  rowData(liver_lipids_exp)$Class %in% enriched_classes
liver_lipids_exp_enriched <- liver_lipids_exp[enriched_heatmap,]
liver_lipids_exp_enriched@metadata[["dimnames"]] <- c("Lipids","Sample")
plot_heatmap(liver_lipids_exp_enriched, cluster_col="none", cluster_row = "none", layout = list(font = list(size = 20)))


