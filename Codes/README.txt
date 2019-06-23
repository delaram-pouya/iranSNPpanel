This is project for finding an optimized snp panel for individual identification based on Iranome database

1. Clean_Iranome.R 
tries to modify the crawled Iranome dataset in the Chr.rds files to gain a clean and easy-to-use data structure

2. MakePanel.R
in this script, we will define and apply some filters to the cleaned iranome dataset, in order to find a set of optimized SNPs for individual identification in iranian population

3. ComparePanel.R
Checks how many of the predicted SNPs are in common with other published SNP panels

4. Functions.R
Includes the functions that are used in the Clean_iranome.R script