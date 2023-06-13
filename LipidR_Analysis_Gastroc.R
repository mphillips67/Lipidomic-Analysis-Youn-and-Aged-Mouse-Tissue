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
gas <- read_excel("Hinton_Muscle-Lipdomics.xlsx", sheet=3)  #read in excel sheet with relevant data

#only keep first column and sample columns. will need to be adjusted based on excel sheet
gas <- gas[,c(1,3:ncol(gas))]

#remove "bad" entries (eg. 13-Docosenamide)
gas <- subset(gas, grepl(":", Metabolite.name))

#log transform data
gas[,2:9] <- log2(gas[2:9])

####### MORE PROCESSING AND LIPIDR ANALYSIS ########

#data sets have multiple entires for the same lipids. In these cases, only keep entry with the highest values. This approach does keep duplicates where the rows are just NA. Those are removed at another step.
gas <- gas %>%
  group_by(Metabolite.name) %>%
  slice_max(across(starts_with("Gas"))) #groups by metabolite name, then only keeps entry with the highest row totals across columns that start with "Gas"

##create tables for  for lipidR
gas_lipids <- na.omit(gas) #remove NAs
colnames(gas_lipids)[1] <- "lipids" #rename first column

nrow(gas_lipids %>%
       group_by(lipids) %>%
       filter(n() > 1)) #if vallue returned is 0, all duplicates were correctly congasidated

#remove rows with same values
same_values_rows <- which(apply(gas_lipids[, 2:8], 1, function(row) all(row == row[1])))
#turns out we have some columns where all samples have the same value. These have to be removed for pca analysis. This line identifies the columns (there are 2 for the gas data).
gas_lipids <- gas_lipids[-same_values_rows, ]

##test for lipid names that don't match lipidr naming system and fix them. to run lipidR, at least 50% of names need to match. 
temp <- gas_lipids$lipids
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

gas_lipids$lipids <- fix #al fixes are done using regex. for a given data set, you will have to look at the names in your orignal table and figure out how to adjust them to the naming system lipidr accepts. 

all_match <- lipidr::annotate_lipids(gas_lipids[[1]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

##create sample annotation table
Sample <- colnames(gas_lipids[,2:ncol(gas_lipids)]) #create a vector with all sample names

SampleType <- c(rep("Young",4), rep("Old",4)) #create vector where each entry describes corresponding sample. 

gas_annotation <- cbind(Sample,SampleType)

write.table(gas_annotation, "Gastrocnemius_annotation.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out gas annotation file

###create lipidR objects
gas_lipids_exp <- as_lipidomics_experiment(gas_lipids, normalized = TRUE, logged = TRUE)

#list unmatched lipids  and add sample annotation 
non_parsed_molecules(gas_lipids_exp)
gas_lipids_exp <- add_sample_annotation(gas_lipids_exp,  "Gastrocnemius_annotation.txt")
gas_lipids_exp <- remove_non_parsed_molecules(gas_lipids_exp)

#pca
mvaresults <- mva(gas_lipids_exp, measure="Area", method="PCA")

plot_mva(mvaresults, color_by="SampleType",hotelling = FALSE, components = c(1,2)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) 


#compare the two groups - univariate analysis
de_old_young <- de_analysis(gas_lipids_exp, Young - Old )
sig <- significant_molecules(de_old_young, p.cutoff = 0.05, logFC.cutoff = 1) #list of sig results
write.table(sig, "Gastrocnemius_sig_lipids.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out gas sig lipids

plot_results_volcano(de_old_young, show.labels =FALSE) + xlim(-6,6) + theme(strip.text.x = element_text(size = 20)) + theme(plot.title = element_text(size = 25, face = "bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) + ylim(0,2)


#enrichment analysis
enrich_results <- lsea(
  de_old_young,
  rank.by = "logFC", min_size = 4, nperm = 100000
)

write.table(enrich_results[,c(1:4,8)], "Gas_lipidset.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out enrichment results

sig_lipidsets <- significant_lipidsets(enrich_results)

plot_enrichment(de_old_young, sig_lipidsets, annotation="class")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))

plot_enrichment(de_old_young, sig_lipidsets, annotation="length")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))



#clustering analysis
##all classes
plot_heatmap(gas_lipids_exp, layout=  list(font = list(size = 8)) ) 

##enriched classes
enriched_classes <- subset(enrich_results, padj <= .05 & grepl("^Class", set))
enriched_classes <- enriched_classes$set
enriched_classes <- sub("Class_","",enriched_classes)
enriched_heatmap <-  rowData(gas_lipids_exp)$Class %in% enriched_classes
gas_lipids_exp_enriched <- gas_lipids_exp[enriched_heatmap,]
gas_lipids_exp_enriched@metadata[["dimnames"]] <- c("Lipids","Sample")
plot_heatmap(gas_lipids_exp_enriched, cluster_rows="none", cluster_cols = "none", layout = list(font = list(size = 20))) 



