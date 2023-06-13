#set working directory to location of files
setwd("~/Dropbox/Lipidomics_McReynolds")

#load packages
library(lipidr)   

#read in data 
heart <- read_excel("Hinton_Muscle-Lipdomics.xlsx", sheet=2)
gas <- read_excel("Hinton_Muscle-Lipdomics.xlsx", sheet=3)
sol <- read_excel("Hinton_Muscle-Lipdomics.xlsx", sheet=4)

#only keep first column and sample columns. will need to be adjusted based on excel sheet
heart <- heart[,c(1,4:ncol(heart))]
gas <- gas[,c(1,3:ncol(gas))]
sol <- sol[,c(1,5:ncol(sol))]

#remove "bad" entries (eg. 13-Docosenamide)
heart <- subset(heart, grepl(":", Metabolite.name))
gas <- subset(gas, grepl(":", Metabolite.name))
sol <- subset(sol, grepl(":", Metabolite.name))

#data sets have multiple entires for the same lipids. In these cases, only keep entry with the highest values. This approach does keep duplicates where the rows are just NA. Those are removed at another stp.
heart <- heart %>%
  group_by(Metabolite.name) %>%
  slice_max(across(starts_with("Heart"))) #groups by metabolite name, then only keeps entry with the highest row totals across columns that start with "Heart"

gas <- gas %>%
  group_by(Metabolite.name) %>%
  slice_max(across(starts_with("Gas"))) #groups by metabolite name, then only keeps entry with the highest row totals across columns that start with "Gas"

sol <- sol %>%
  group_by(Metabolite.name) %>%
  slice_max(across(starts_with("Sol"))) #groups by metabolite name, then only keeps entry with the highest row totals across columns that start with "Sol"


#create heart name list 
heart_lipids <- na.omit(heart) #remove NAs
colnames(heart_lipids)[1] <- "lipids" #rename first column

##conversions 
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
fix <- sub("LPE-N ", "LPEN ", fix)

heart_lipids$converted_lipid_name <- fix

all_match <- lipidr::annotate_lipids(heart_lipids[[10]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

write.table(heart_lipids[,c(1,10)], "Lipid_Name_Conversion_Heart_tissue.txt", quote = FALSE, row.names = FALSE, sep="\t")


#create sol list
sol_lipids <- na.omit(sol) #remove NAs
colnames(sol_lipids)[1] <- "lipids" #rename first column

##conversions 
temp <- sol_lipids$lipids
fix = sub(";(O\\d*)", "(\\1)", temp)
fix = sub(" O-", "O ", fix)
fix = sub(" P-", "P ", fix)
fix = sub("-FA", "/", fix)
fix = sub(";(.*)$", "(\\1)", fix)
fix <- sub("([^\\)])\\|", "\\1(m)|", fix)
fix <- sub("(.*\\(m\\))\\|(.*$)", "\\1|\\2(m2)", fix)
fix <- sub("PE-Cer", "PECer", fix)
fix <- sub("PI-Cer", "PICer", fix)
fix <- sub("LPE-N ", "LPEN ", fix)

sol_lipids$converted_lipid_name <- fix

all_match <- lipidr::annotate_lipids(sol_lipids[[10]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

write.table(sol_lipids[,c(1,10)], "Lipid_Name_Conversion_Sol_tissue.txt", quote = FALSE, row.names = FALSE, sep="\t")

#create gas list
gas_lipids <- na.omit(gas) #remove NAs
colnames(gas_lipids)[1] <- "lipids" #rename first column

##conversions 
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
fix <- sub("LPE-N ", "LPEN ", fix)

gas_lipids$converted_lipid_name <- fix

all_match <- lipidr::annotate_lipids(gas_lipids[[10]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

write.table(gas_lipids[,c(1,10)], "Lipid_Name_Conversion_Gas_tissue.txt", quote = FALSE, row.names = FALSE, sep="\t")


