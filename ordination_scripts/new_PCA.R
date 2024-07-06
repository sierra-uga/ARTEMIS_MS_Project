## filtering
org_metadata <- read.delim("required_files/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE, row.names="sample.illumina") 

## data culling 
org_metadata <- filter(org_metadata, Sample.Control == "True.Sample") 
org_metadata <- filter(org_metadata, DOC != "NA")
org_metadata <- filter(org_metadata, watertype != "Other")# use tidyr to select "real samples" (non-blanks)
# metadata <- org_metadata[-(which(org_metadata$Station %in% c("STN198", "STN153", "STN056b", "STN012"))),] # removes station 198 and 153 for dotson analysis
distinct <- org_metadata %>% distinct(Depth, Station) # 
fil_metadata <- org_metadata[row.names(org_metadata) %in% row.names(distinct),] # REMOVING DUPLICATES IN THE PCA

metadata <- fil_metadata[ -c( 1, 26:29)] # remove filter-related stuff
metadata <- metadata[-(which(metadata$Station %in% c("STN198", "STN056b"))),] 
metadata <- metadata[ -c( 1, 2:18)] # remove barcode seq/uneeded stuff
metadata <- select(metadata, 1:16, 21, 22, 23, 18) # select numeric + siderophore column + watertype

#metadata[is.na(metadata)] <- 0 # replace NAs with 0's for PCA
#list_true <- replace_na(metadata$True_Flow, "Other") #replace NA with "Other" for coloring
#metadata$True_Flow <- list_true # changing actual column in dataframe

metadata <- metadata %>% mutate_at(1:16, as.numeric) # convert character columns (except sampleID) to numeric for analysis 

# remove unneeded columns (Oxygen, FIECO, Par, Siderophore)
metadata <- select(metadata,-c(Oxygen, FlECO.AFL, CTD_Depth, Par, Chl_a))
metadata <- metadata %>% rename(c(Nitrate = Lab_NO3, 
                                  Nitrite = Lab_NO2, 
                                  Ammonium = Lab_NH4, 
                                  #Chlorophyll = Chl_a, 
                                  Oxygen = Sb_Oxygen,
                                  Phosphate = Lab_PO4))
# remove NA for iron analysis

metadata_above <- filter(metadata, Depth <= 200)
# above 200m 
#metadata_above <- filter(metadata, Depth <= 200)
open <- metadata_above$Location == "Open_polynya"
dotson <- metadata_above$Location == "Dotson"
east <- metadata_above$Location == "Eastern_CC"
#cont <- metadata$Location == "Cont_Shelf"
west <- metadata_above$Location == "Western_CC"
getz <- metadata_above$Location == "Getz"
west_op <- metadata_above$Location == "West_OP"

metadata_above <- metadata_above %>% mutate_at(1:11, as.numeric) 
metadata_above <- dplyr::select(metadata_above, 1:11) # select only numerical columns
# create vectors from orginial metadata file for watertype
metadataPca_above <- prcomp(metadata_above, scale.=TRUE)

sd <- metadataPca_above$sdev
loadings <- metadataPca_above$rotation
rownames(loadings) <- colnames(metadata_above)
scores <- metadataPca_above$x

# use scree plot to reduce dimensionality
var <- sd^2
varPercent <- var/sum(var) * 100
barplot(varPercent, xlab="PC", ylab="Percent Variance", names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col="gray")
abline(h=1/ncol(metadata_above)*100, col="red")

# determine % variance explained
varPercent[1:4]
sum(varPercent[1:4])

# table for the loadings
loadings
sqrt(1/ncol(metadata_above))

scaling <- 4.5
textNudge <- 1.2

metadata <- metadata_above
pdf(file = "final_graphics/PCA_location_with_DOC_upper.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=16, cex=1.25, lwd=2, col="#09A20D", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=16, cex=1.25, lwd=2, col="#5AD0FC", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=16, cex=1.25, col="darkred", bg="darkred")
#points(scores[cont, 1], scores[cont, 2], pch=21, cex=1.25, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=16, cex=1.25, col="red", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=16, cex=1.25, col="#006B93", bg="#006B93")
points(scores[west_op, 1], scores[west_op, 2], pch=16, cex=1.25, col="purple", bg="purple")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
title("Upper 200m")
# add names 
dev.off()

# below 200m
metadata_below <- filter(metadata, Depth >= 200)
open2 <- metadata_below$Location == "Open_polynya"
dotson2 <- metadata_below$Location == "Dotson"
east2 <- metadata_below$Location == "Eastern_CC"
#cont <- metadata$Location == "Cont_Shelf"
west2 <- metadata_below$Location == "Western_CC"
getz2 <- metadata_below$Location == "Getz"
west_op2 <- metadata_below$Location == "West_OP"

metadata_below <- dplyr::select(metadata_below, 1:11) # select only numerical columns
# create vectors from orginial metadata file for watertype
metadataPca_below <- prcomp(metadata_below, scale.=TRUE)

sd2 <- metadataPca_below$sd2ev
loadings2 <- metadataPca_below$rotation
rownames(loadings2) <- colnames(metadata_below)
scores2 <- metadataPca_below$x

# use scree plot to reduce dimensionality
var2 <- sd2^2
var2Percent <- var2/sum(var2) * 100

# determine % var2iance explained
var2Percent[1:4]
sum(var2Percent[1:4])

# table for the loadings2
loadings2
sqrt(1/ncol(metadata_below))

scaling <- 4.5
textNudge <- 1.2

pdf(file = "final_graphics/PCA_location_with_DOC_lower.pdf", width = 6, height = 7) 
plot(scores2[, 1], scores2[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores2[open2, 1], scores2[open2, 2], pch=16, cex=1.25, lwd=2, col="#09A20D", bg="#09A20D")
points(scores2[dotson2, 1], scores2[dotson2, 2], pch=16, cex=1.25, lwd=2, col="#5AD0FC", bg="#5AD0FC")
points(scores2[east2, 1], scores2[east2, 2], pch=16, cex=1.25, col="darkred", bg="darkred")
#points(scores2[cont, 1], scores2[cont, 2], pch=21, cex=1.25, col="#3A3A3A", bg="#A3DCA5")
points(scores2[west2, 1], scores2[west2, 2], pch=16, cex=1.25, col="red", bg="red")
points(scores2[getz2, 1], scores2[getz2, 2], pch=16, cex=1.25, col="#006B93", bg="#006B93")
points(scores2[west_op2, 1], scores2[west_op2, 2], pch=16, cex=1.25, col="purple", bg="purple")
arrows(0, 0, loadings2[, 1]* scaling, loadings2[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings2[, 1]*scaling*textNudge, loadings2[, 2]*scaling*1.3, rownames(loadings2), col="black", cex=0.8)
title("Below 200m")
# add names 
dev.off()

pdf("final_graphics/PCA_location.pdf", width = 7, height = 10)
par(mfrow=c(2,1))
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=16, cex=1.25, lwd=2, col="#09A20D", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=16, cex=1.25, lwd=2, col="#5AD0FC", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=16, cex=1.25, col="darkred", bg="darkred")
#points(scores[cont, 1], scores[cont, 2], pch=21, cex=1.25, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=16, cex=1.25, col="red", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=16, cex=1.25, col="#006B93", bg="#006B93")
points(scores[west_op, 1], scores[west_op, 2], pch=16, cex=1.25, col="purple", bg="purple")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
title("Upper 200m")

plot(scores2[, 1], scores2[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores2[open2, 1], scores2[open2, 2], pch=16, cex=1.25, lwd=2, col="#09A20D", bg="#09A20D")
points(scores2[dotson2, 1], scores2[dotson2, 2], pch=16, cex=1.25, lwd=2, col="#5AD0FC", bg="#5AD0FC")
points(scores2[east2, 1], scores2[east2, 2], pch=16, cex=1.25, col="darkred", bg="darkred")
#points(scores2[cont, 1], scores2[cont, 2], pch=21, cex=1.25, col="#3A3A3A", bg="#A3DCA5")
points(scores2[west2, 1], scores2[west2, 2], pch=16, cex=1.25, col="red", bg="red")
points(scores2[getz2, 1], scores2[getz2, 2], pch=16, cex=1.25, col="#006B93", bg="#006B93")
points(scores2[west_op2, 1], scores2[west_op2, 2], pch=16, cex=1.25, col="purple", bg="purple")
arrows(0, 0, loadings2[, 1]* scaling, loadings2[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings2[, 1]*scaling*textNudge, loadings2[, 2]*scaling*1.3, rownames(loadings2), col="black", cex=0.8)
title("Below 200m")
dev.off()