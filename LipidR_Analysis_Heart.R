#set working directory to location of files
setwd("~/Dropbox/Lipidomics_McReynolds")

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
heart <- read_excel("Hinton_Muscle-Lipdomics.xlsx", sheet=2) #read in excel sheet with relevant data

#only keep first column and sample columns. will need to be adjusted based on excel sheet
heart <- heart[,c(1,4:ncol(heart))]

#remove "bad" entries (eg. 13-Docosenamide)
heart <- subset(heart, grepl(":", Metabolite.name))

#log transform data
heart[,2:9] <- log2(heart[2:9])

####### MORE PROCESSING AND LIPIDR ANALYSIS ########

#data sets have multiple entires for the same lipids. In these cases, only keep entry with the highest values. This approach does keep duplicates where the rows are just NA. Those are removed at another stp.
heart <- heart %>%
  group_by(Metabolite.name) %>%
  slice_max(across(starts_with("Heart"))) #groups by metabolite name, then only keeps entry with the highest row totals across columns that start with "Heart"


##create tables for  for lipidR
heart_lipids <- na.omit(heart) #remove NAs
colnames(heart_lipids)[1] <- "lipids" #rename first column

nrow(heart_lipids %>%
       group_by(lipids) %>%
       filter(n() > 1)) #if vallue returned is 0, all duplicates were correctly consolidated

#remove rows with same values
same_values_rows <- which(apply(heart_lipids[, 2:8], 1, function(row) all(row == row[1])))
#turns out we have some columns where all samples have the same value. These have to be removed for pca analysis. This line identifies the columns (there are 2 for the heart data).
heart_lipids <- heart_lipids[-same_values_rows, ]

##test for lipid names that don't match lipidr naming system and fix them. to run lipidR, at least 50% of names need to match. 
temp <- heart_lipids$lipids
fix = sub(";(O\\d*)", "(\\1)", temp)
fix = sub(" O-", "O ", fix)
fix = sub(" P-", "P ", fix)
fix = sub("-FA", "/", fix)
fix = sub(";(.*)$", "(\\1)", fix)
fix <- sub("([^\\)])\\|", "\\1(m)|", fix)
fix <- sub("(.*\\(m\\))\\|(.*$)", "\\1|\\2(m2)", fix)
fix <- sub("PE-Cer", "PECer", fix)
fix <- sub("PI-Cer", "PICer", fix)
fix <- sub("LPE-N", "LPEN ", fix)

heart_lipids$lipids <- fix #al fixes are done using regex. for a given data set, you will have to look at the names in your orignal table and figure out how to adjust them to the naming system lipidr accepts. 

all_match <- lipidr::annotate_lipids(heart_lipids[[1]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

##create sample annotation table
Sample <- colnames(heart_lipids[,2:ncol(heart_lipids)]) #create a vector with all sample names

SampleType <- c(rep("Young",4), rep("Old",4)) #create vector where each entry describes corresponding sample. 

heart_annotation <- cbind(Sample,SampleType)

write.table(heart_annotation, "Heart_annotation.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart annotation file

###create lipidR objects
heart_lipids_exp <- as_lipidomics_experiment(heart_lipids, normalized = TRUE, logged = TRUE)

#list unmatched lipids  and add sample annotation 
non_parsed_molecules(heart_lipids_exp)
heart_lipids_exp <- add_sample_annotation(heart_lipids_exp,  "Heart_annotation.txt")
heart_lipids_exp <- remove_non_parsed_molecules(heart_lipids_exp)

#pca
mvaresults <- mva(heart_lipids_exp, measure="Area", method="PCA")

plot_mva(mvaresults, color_by="SampleType",hotelling = FALSE, components = c(1,2)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) 

#compare the two groups - univariate analysis
de_old_young <- de_analysis(heart_lipids_exp, Young - Old )
sig <- significant_molecules(de_old_young, p.cutoff = 0.05, logFC.cutoff = 1) #list of sig results
write.table(sig, "Heart_sig_lipids.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart sig lipids

plot_results_volcano(de_old_young, show.labels =FALSE) + xlim(-8,8) + theme(strip.text.x = element_text(size = 20)) + theme(plot.title = element_text(size = 25, face = "bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) 

#enrichment analysis
enrich_results <- lsea(
  de_old_young,
  rank.by = "logFC", min_size = 4, nperm = 100000
)

write.table(enrich_results[,c(1:4,8)], "Heart_lipidset.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out enrichment results

sig_lipidsets <- significant_lipidsets(enrich_results)

plot_enrichment(de_old_young, sig_lipidsets, annotation="class")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))

plot_enrichment(de_old_young, sig_lipidsets, annotation="length")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))


#clustering analysis
##all classes
plot_heatmap(heart_lipids_exp, layout=  list(font = list(size = 8)) ) 

##enriched classes
enriched_classes <- subset(enrich_results, padj <= .05 & grepl("^Class", set))
enriched_classes <- enriched_classes$set
enriched_classes <- sub("Class_","",enriched_classes)
enriched_heatmap <-  rowData(heart_lipids_exp)$Class %in% enriched_classes
heart_lipids_exp_enriched <- heart_lipids_exp[enriched_heatmap,]
heart_lipids_exp_enriched@metadata[["dimnames"]] <- c("Lipids","Sample")
plot_heatmap(heart_lipids_exp_enriched, cluster_col="none", cluster_row = "none", layout = list(font = list(size = 20)))

