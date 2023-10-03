

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
kid <- read_excel("Hinton_Kid-Liv_Lipdomics.xlsx", sheet=2) #read in excel sheet with relevant data

#only keep first column and sample columns. will need to be adjusted based on excel sheet
kid <- kid[,2:10]

#remove "bad" entries (eg. 13-Docosenamide)
kid <- subset(kid, grepl(":", Metabolite.name))

#log transform data
kid[,2:9] <- log2(kid[2:9])

####### MORE PROCESSING AND LIPIDR ANALYSIS ########

#data sets have multiple entires for the same lipids. In these cases, only keep entry with the highest values. This approach does keep duplicates where the rows are just NA. Those are removed at another stp.
kid <- kid %>%
  group_by(Metabolite.name) %>%
  slice_max(across(starts_with("Kidney"))) #groups by metabolite name, then only keeps entry with the highest row totals across columns that start with "Kidney"

##create tables for  for lipidR
kidney_lipids <- na.omit(kid) #remove NAs
colnames(kidney_lipids)[1] <- "lipids" #rename first column

nrow(kidney_lipids %>%
       group_by(lipids) %>%
       filter(n() > 1)) #if vallue returned is 0, all duplicates were correctly consolidated

#remove rows with same values
same_values_rows <- which(apply(kidney_lipids[, 2:8], 1, function(row) all(row == row[1])))
#turns out we have some columns where all entries have the same value. These have to be removed for pca analysis. This line identifies the columns (there is 3 for the kidney data).
kidney_lipids <- kidney_lipids[-same_values_rows, ]


##test for lipid names that don't match lipidr naming system and fix them. to run lipidR, at least 50% of names need to match. 
temp <- kidney_lipids$lipids
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

kidney_lipids$lipids <- fix #all fixes are done using regex. for a given data set, you will have to look at the names in your orignal table and figure out how to adjust them to the naming system lipidr accepts. 

all_match <- lipidr::annotate_lipids(kidney_lipids[[1]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

##create sample annotation table
Sample <- colnames(kidney_lipids[,2:ncol(kidney_lipids)]) #create a vector with all sample names

SampleType <- c(rep("Young",4), rep("Old",4)) #create vector where each entry describes corresponding sample. 

kidney_annotation <- cbind(Sample,SampleType)

write.table(kidney_annotation, "Kidney_annotation.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart annotation file


###create lipidR objects
kidney_lipids_exp <- as_lipidomics_experiment(kidney_lipids, normalized = TRUE, logged = TRUE)

#list unmatched lipids  and add sample annotation 
non_parsed_molecules(kidney_lipids_exp)
kidney_lipids_exp <- add_sample_annotation(kidney_lipids_exp,  "Kidney_annotation.txt")
kidney_lipids_exp <- remove_non_parsed_molecules(kidney_lipids_exp)

#pca
mvaresults <- mva(kidney_lipids_exp, measure="Area", method="PCA")

plot_mva(mvaresults, color_by="SampleType",hotelling = FALSE, components = c(1,2)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) 

#compare the two groups - univariate analysis
de_old_young <- de_analysis(kidney_lipids_exp, Young - Old )
sig <- significant_molecules(de_old_young, p.cutoff = 0.05, logFC.cutoff = 1) #list of sig results
write.table(sig, "Kidney_sig_lipids.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart sig lipids
write.table(de_old_young, "Kidney_de_results.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart sig lipids

plot_results_volcano(de_old_young, show.labels =TRUE) + xlim(-10,10) + theme(strip.text.x = element_text(size = 20)) + theme(plot.title = element_text(size = 25, face = "bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) + ylim(0,4) # + geom_text_repel(show.legend  = F, data=subset(de_old_young, max.overlaps = 33, adj.P.Val < 0.05 & (logFC > 1 | logFC < -1)))

#enrichment analysis
enrich_results <- lsea(
  de_old_young,
  rank.by = "logFC", min_size = 4, nperm = 100000
)

write.table(enrich_results[,c(1:4,6,8)], "Kidney_lipidset.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out enrichment results

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
enriched_heatmap <-  rowData(kidney_lipids_exp)$Class %in% enriched_classes
kidney_lipids_exp_enriched <- kidney_lipids_exp[enriched_heatmap,]
kidney_lipids_exp_enriched@metadata[["dimnames"]] <- c("Lipids","Sample")
plot_heatmap(kidney_lipids_exp_enriched, cluster_col="none", cluster_row = "none", layout = list(font = list(size = 20)))



