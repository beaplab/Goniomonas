library(stringr)
library(openxlsx)
library(data.table)
setwd('/home/dzavadska/Data/Goniomonas/TO_UPLOAD/TREES/tree_viz_rename_pairwise_dist')

headers <- fread("name_dict.csv", header =F)

full_headers <- str_replace_all(headers$V1, ">", "") 
full_headers <- str_replace_all(headers$V1, " 1781 bp", "")
short_headers <- str_replace_all(headers$V2, ">", "") 

tree <- readLines("RAxML_bipartitionsBranchLabels.FIN_23Oct2024_BS_BIPART")
tree <- paste0(tree)
#str(tree)
#full_headers[grep(short_headers[1], full_headers)]
#tree <- str_replace_all(tree, short_headers[1], full_headers[grep(short_headers[1], full_headers)])


for (i in short_headers) {
  if ( length(grep("BEAP", i[1]) ) > 0 ) next
  if ( length(grep("HFCC23", i[1]) ) > 0 ) next
  if ( length(grep("HFCC22", i[1]) ) > 0 ) next
  if ( i == "HFCC" ) next
  tree <- str_replace_all(tree, i[1], full_headers[grep(i[1], full_headers)])
}

tree <- str_replace_all(tree, ">", "")

write(tree, "tree_renamed")


############################################################################################
######################   Custom naming mrbayes
############################################################################################



headers <- fread("name_dict_manual_edit.csv", header =T)

full_headers <- paste(str_replace_all(headers$short, ">", ""), headers$CR_CR_genus , headers$CR_CR_species, headers$special_Separator)
short_headers <- str_replace_all(headers$short, ">", "") 

tree <- readLines("FIN_mafft_gappyout_gt02.nex.con.tre")
tree <- paste0(tree)
#str(tree)
#full_headers[grep(short_headers[1], full_headers)]
#tree <- str_replace_all(tree, short_headers[1], full_headers[grep(short_headers[1], full_headers)])


for (i in short_headers) {
  if ( length(grep("BEAP", i[1]) ) > 0 ) next
  if ( length(grep("HFCC23", i[1]) ) > 0 ) next
  if ( length(grep("HFCC22", i[1]) ) > 0 ) next
  if ( i == "HFCC" ) next
  tree <- str_replace_all(tree, i[1], full_headers[grep(i[1], full_headers)])
  #grep(i[1],tree)
}

tree <- str_replace_all(tree, ">", "")

write(tree, "tree_renamed_custom_mb.nwk")




############################################################################################
######################   Custom naming raxml
############################################################################################


raxml_headers <- data.frame(fread("raxml_IDs.csv", sep = '\n', header = F))
raxml_headers <- raxml_headers$V1
raxml_headers <- str_replace_all(raxml_headers, "[.]", " ")
raxml_headers <- str_replace_all(raxml_headers, "[|]", "_")
raxml_headers <- str_replace_all(raxml_headers, "[+]", "_")
#raxml_headers <- str_replace_all(raxml_headers, "/", "_")
full_headers <- paste(str_replace_all(headers$short, ">", ""), headers$CR_CR_genus , headers$CR_CR_species)

tree <- readLines("RAxML_bipartitionsBranchLabels.FIN_23Oct2024_BS_BIPART")
tree <- paste0(tree)
#str(tree)
#full_headers[grep(short_headers[1], full_headers)]

#tree <- str_replace_all(tree, " ", "_")
tree <- str_replace_all(tree, "[|]", "_")
tree <- str_replace_all(tree, "[+]", "_")

for (i in short_headers) {
  if ( length(raxml_headers[grep(i[1], raxml_headers)]) == 0 ) next
  if ( length(grep("BEAP", i[1]) ) > 0 ) next
  if ( length(grep("HFCC23", i[1]) ) > 0 ) next
  if ( length(grep("HFCC22", i[1]) ) > 0 ) next
  
  tree <- str_replace_all(tree, raxml_headers[grep(i[1], raxml_headers)], full_headers[grep(i[1], full_headers)])
  if ( length(grep("________", tree)) > 0 ) {print(i)}
  #grep(raxml_headers[grep(i[1], raxml_headers)][1],tree)
}
tree


write(tree, "tree_renamed_custom_raxml.nwk")




############################################################################################
######################   CODE DUMP
############################################################################################


headers <- str_replace(headers, "/*", "")
str(headers)

long_h <- fread("FIN_mafft_gappyout.nex.con.tre", header = F)
str(long_h)


tree_ch <- paste0(tree)  
tree_ch <-  toString(tree_ch)

for (i in c(1:length(headers$Gren151F))) {
  print(i)
  full_name <- long_h$V1[grep(headers[i,1], long_h$V1)]

  short_name <- str_remove(toString(headers[i,1]), ">")
  full_name <- str_remove(toString(full_name[1]), ">")
  full_name <- str_remove(toString(full_name[1]), "Eukaryota|Discoba|Euglenozoa|")
  tree_ch <- str_replace(toString(tree_ch),  short_name,  full_name)
}

tree_ch
