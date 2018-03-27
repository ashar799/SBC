### This file loads the Van't Veer 70 gene signature ####
loadVIJVERsignature = function(){
### Load the 70 gene signature
######################## 70 GENE SIGNATURE #################################################
###  Geting Features Which Actually are Part of the 70 Gene Signature 
gns231 <- read.xls(xls ='/home/bit/ashar/ExpressionSets/ONE_VIEW/breastcancer/415530a-s9.xls', skip=0, header=TRUE, stringsAsFactors=FALSE)
###Remove special characters in the colums header, which are due to white spaces present in the Excel files
colnames(gns231) <- gsub("\\.\\.", "", colnames(gns231))
###Remove GO annotation
gns231 <- gns231[, -grep("sp_xref_keyword_list", colnames(gns231))]
###Reorder the genes in decreasing order by absolute correlation
gns231 <- gns231[order(abs(gns231$correlation), decreasing=TRUE),]
###Select the feature identifiers corresponding to the top 231 and 70 genes
gns231$genes231 <- TRUE
gns231$genes70 <- gns231$accession %in% gns231$accession[1:70]
######## Signature Probes ##########################################
vijver.signature <<- gns231$accession[gns231$genes70]
}