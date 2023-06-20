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

#read in liver data 
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

liver_lipids$converted_lipid_name  <- fix #al fixes are done using regex. for a given data set, you will have to look at the names in your orignal table and figure out how to adjust them to the naming system lipidr accepts. 

all_match <- lipidr::annotate_lipids(liver_lipids[[10]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

write.table(liver_lipids[,c(1,10)], "Lipid_Name_Conversion_Liver_tissue.txt", quote = FALSE, row.names = FALSE, sep="\t")

#read in kidney data 
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

kidney_lipids$converted_lipid_name <- fix #all fixes are done using regex. for a given data set, you will have to look at the names in your orignal table and figure out how to adjust them to the naming system lipidr accepts. 

all_match <- lipidr::annotate_lipids(kidney_lipids[[10]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

write.table(kidney_lipids[,c(1,10)], "Lipid_Name_Conversion_Kidney_tissue.txt", quote = FALSE, row.names = FALSE, sep="\t")


