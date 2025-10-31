library(ape)

setwd("/home/dzavadska/Data/Goniomonas/Jeuk_micro/TREES_Jeukmicro/tree_viz_rename/")

`%rin%` = function (pattern, list) {
  vapply(pattern, function (p) any(grepl(p, list)), logical(1L), USE.NAMES = FALSE)
}


tree<-read.tree("tree2_0_renamed_custom_raxml.nwk")
PatristicDistMatrix<-cophenetic(tree)
write.csv(PatristicDistMatrix, "RAXML2_0_raw_distance_matrix.csv")


tree<-read.nexus("tree2_0_renamed_custom_mb.nwk")
PatristicDistMatrix<-cophenetic(tree)
write.csv(PatristicDistMatrix, "gtr_mb2_0_raw_distance_matrix.csv")




############################################################################################

headers <- fread("name_dict_manual_edit.csv", header =T)
headers$short <- str_replace_all(headers$short, ">", "") 


mb_matrix <- read.csv("gtr_mb2_0_raw_distance_matrix.csv")

genera <- c("Goniomonas", "Limnogoniomonas", "Aquagoniomonas", "Marigoniomonas", "Neptunogoniomonas", "Thalassogonimonas", "Cosmogoniomonas", "Poseidogoniomonas","Ebisugoniomonas", "Baltigoniomonas")

for(i in genera){
  
ids <- headers$short[grep(i,headers$CR_CR_genus)]

rows <- c()
cols <- c()
for (id in ids) {
  rows <- c(rows, grep(id, mb_matrix$X))
  cols <- c(cols, grep(id, colnames(mb_matrix)))
}

df <- mb_matrix[rows, c(1, cols)]
df$genus <- i
write.csv(df, file = paste0("gtr_mb2_0_raw_distance_matrix_", i))

}




raxml_matrix <- read.csv( "RAXML2_0_raw_distance_matrix.csv")

genera <- c("Goniomonas", "Limnogoniomonas", "Aquagoniomonas", "Marigoniomonas", "Neptunogoniomonas", "Thalassogonimonas", "Cosmogoniomonas", "Poseidogoniomonas","Ebisugoniomonas")

for(i in genera){
  
  ids <- headers$short[grep(i,headers$CR_CR_genus)]
  
  rows <- c()
  cols <- c()
  for (id in ids) {
    rows <- c(rows, grep(id, raxml_matrix$X))
    cols <- c(cols, grep(id, colnames(raxml_matrix)))
  }
  
  df <- raxml_matrix[rows, c(1, cols)]
  df$genus <- i
  write.csv(df, file = paste0("raxml2_0_raw_distance_matrix_", i))
  
}



############################################################################################
######################   Custom naming mrbayes
############################################################################################



headers <- fread("name_dict_manual_edit.csv", header =T)

full_headers <- paste(str_replace_all(headers$short, ">", ""), headers$CR_CR_genus , headers$CR_CR_species)
short_headers <- str_replace_all(headers$short, ">", "") 
short_headers <- str_replace_all(short_headers, "-", ".") 

mb <- read.csv("gtr_mb2_0_raw_distance_matrix.csv")
mb$X <- as.character(mb$X)
#str(tree)
#full_headers[grep(short_headers[1], full_headers)]
#tree <- str_replace_all(tree, short_headers[1], full_headers[grep(short_headers[1], full_headers)])


for (i in short_headers) {
  if ( length(grep("BEAP", i[1]) ) > 0 ) next
  if ( length(grep("HFCC23", i[1]) ) > 0 ) next
  if ( length(grep("HFCC22", i[1]) ) > 0 ) next
  colnames(mb)[grep(i[1],colnames(mb))] <- full_headers[grep(i[1], full_headers)]
  mb$X[grep(i[1],mb$X)][1] <- full_headers[grep(i[1], full_headers)]
  #grep(i[1],tree)
}


write.csv(mb, "gtr_cophenetic_distance_renamed_mb2_0.csv")




############################################################################################
######################   Custom naming raxml
############################################################################################


raxml <- read.csv("RAXML2_0_raw_distance_matrix.csv")
raxml$X <- as.character(mb$X)
#str(tree)
#full_headers[grep(short_headers[1], full_headers)]
#tree <- str_replace_all(tree, short_headers[1], full_headers[grep(short_headers[1], full_headers)])


for (i in short_headers) {
  if ( length(grep("BEAP", i[1]) ) > 0 ) next
  if ( length(grep("HFCC23", i[1]) ) > 0 ) next
  if ( length(grep("HFCC22", i[1]) ) > 0 ) next
  colnames(raxml)[grep(i[1],colnames(raxml))] <- full_headers[grep(i[1], full_headers)]
  raxml$X[grep(i[1],raxml$X)][1] <- full_headers[grep(i[1], full_headers)]
  #grep(i[1],tree)
}


write.csv(raxml, "cophenetic_distance_renamed_raxml2_0.csv")
