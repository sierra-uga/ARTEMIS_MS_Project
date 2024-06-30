## filtering
org_metadata <- read.delim("required_files/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE, row.names ="sample_name") 
#org_metadata <- org_metadata_filter
org_metadata <- filter(org_metadata, watertype == "CDW") 
distinct <- org_metadata %>% distinct(Depth, Station) # 
fil_metadata <- org_metadata[row.names(org_metadata) %in% row.names(distinct),] # REMOVING DUPLICATES IN THE PCA


## data culling 
org_metadata <- filter(fil_metadata, Sample.Control == "True.Sample") # use tidyr to select "real samples" (non-blanks)
org_metadata <- filter(org_metadata, Station != "STN056b") 
#org_metadata <- filter(org_metadata, Iron != "NA")
#org_metadata <- filter(org_metadata, DOC != "NA")
org_metadata <- filter(org_metadata, Location != "Cont_Shelf")

# metadata <- org_metadata[-(which(org_metadata$Station %in% c("STN198", "STN153", "STN056b", "STN012"))),] # removes station 198 and 153 for dotson analysis
metadata <- org_metadata[ -c( 1, 2:19)]
# metadata <- org_metadata[-(which(org_metadata$Station %in% c("STN198", "STN153", "STN056b", "STN012"))),] # removes station 198 and 153 for dotson analysis
metadata <- metadata[ -c( 1, 6:10)] # remove filter-related stuff
metadata <- dplyr::select(metadata,-c(Iron, Chl_a, Par, FlECO.AFL, Oxygen))
#metadata <- dplyr::select(org_metadata, 1:12, 17, 19, 21, 23) # select numeric + siderophore column + watertype

metadata <- metadata %>% mutate_at(1:11, as.numeric) 

metadata <- metadata %>% rename(c(Nitrate = Lab_NO3, 
                                  Nitrite = Lab_NO2, 
                                  Ammonium = Lab_NH4, 
                                  Depth = CTD_Depth, 
                                  Oxygen = Sb_Oxygen,
                                  Phosphate = Lab_PO4))
# remove NA for iron analysis


water_type_shapes <- c("AASW"= 8, "AASW-WW"=13, "WW"=23, "WW-CDW"=22, "CDW"=21)
metadata_above <- metadata
# above 200m 
#metadata_above <- filter(metadata, Depth <= 200)
open <- metadata_above$Location == "Open_polynya"
dotson <- metadata_above$Location == "Dotson"
east <- metadata_above$Location == "Eastern_CC"
#cont <- metadata$Location == "Cont_Shelf"
west <- metadata_above$Location == "Western_CC"
getz <- metadata_above$Location == "Getz"
west_op <- metadata$Location == "West_OP"

metadata_above <- metadata_above %>% mutate_at(1:9, as.numeric) 

metadata_above <- dplyr::select(metadata_above, 1:9) # select only numerical columns
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
pdf(file = "ordination_scripts/graphics/PCA_location_with_DOC.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=water_type_shapes[metadata[open, "watertype"]], cex=1, lwd=2, col="#09A20D", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=water_type_shapes[metadata[dotson, "watertype"]], cex=1, lwd=2, col="#5AD0FC", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=water_type_shapes[metadata[east, "watertype"]], cex=1, col="darkred", bg="darkred")
#points(scores[cont, 1], scores[cont, 2], pch=21, cex=1, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=water_type_shapes[metadata[west, "watertype"]], cex=1, col="red", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=water_type_shapes[metadata[getz, "watertype"]], cex=1, col="#006B93", bg="#006B93")
points(scores[west_op, 1], scores[west_op, 2], pch=water_type_shapes[metadata[west_op, "watertype"]], cex=1, col="purple", bg="purple")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
title("Upper 200m")
# add names 
dev.off()

pdf(file = "final_graphics/PCA_CDW.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=16, cex=1, lwd=2, col="#09A20D", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=16, cex=1, lwd=2, col="#5AD0FC", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=16, cex=1, col="darkred", bg="darkred")
#points(scores[cont, 1], scores[cont, 2], pch=21, cex=1, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=16, cex=1, col="red", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=16, cex=1, col="#006B93", bg="#006B93")
points(scores[west_op, 1], scores[west_op, 2], pch=16, cex=1, col="purple", bg="purple")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
title("PCA of Only CDW Samples")
dev.off()

# below 200m
metadata_below <- filter(metadata, Depth >= 200)
open2 <- metadata_below$Location == "Open_polynya"
dotson2 <- metadata_below$Location == "Dotson"
east2 <- metadata_below$Location == "Eastern_CC"
#cont <- metadata$Location == "Cont_Shelf"
west2 <- metadata_below$Location == "Western_CC"
getz2 <- metadata_below$Location == "Getz"

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

pdf(file = "final_graphics/PCA_location_with_DOC.pdf", width = 7, height = 11) 
par(mfrow=c(2,1))
par(mar = c(2, 4, 2, 2)) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=water_type_shapes[metadata[open, "watertype"]], cex=1.3, col="#09A20D", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=water_type_shapes[metadata[dotson, "watertype"]], cex=1.3, col="#5AD0FC", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=water_type_shapes[metadata[east, "watertype"]], cex=1.3, col="darkred", bg="darkred")
#points(scores[cont, 1], scores[cont, 2], pch=21, cex=1.3, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=water_type_shapes[metadata[west, "watertype"]], cex=1.3, col="red", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=water_type_shapes[metadata[getz, "watertype"]], cex=1.3, col="#006B93", bg="#006B93")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
title("Upper 200m")

plot(scores2[, 1], scores2[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores2[open2, 1], scores2[open2, 2], pch=water_type_shapes[metadata[open2, "watertype"]], cex=1.3, col="#09A20D", bg="#09A20D")
points(scores2[dotson2, 1], scores2[dotson2, 2], pch=water_type_shapes[metadata[dotson2, "watertype"]], cex=1.3, col="#5AD0FC", bg="#5AD0FC")
points(scores2[east2, 1], scores2[east2, 2], pch=water_type_shapes[metadata[east2, "watertype"]], cex=1.3, col="darkred", bg="darkred")
#points(scores2[cont, 1], scores2[cont, 2], pch=21, cex=1.3, col="#3A3A3A", bg="#A3DCA5")
points(scores2[west2, 1], scores2[west2, 2], pch=water_type_shapes[metadata[west2, "watertype"]], cex=1.3, col="red", bg="red")
points(scores2[getz2, 1], scores2[getz2, 2], pch=water_type_shapes[metadata[getz2, "watertype"]], cex=1.3, col="#006B93", bg="#006B93")
arrows(0, 0, loadings2[, 1]* scaling, loadings2[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings2[, 1]*scaling*textNudge, loadings2[, 2]*scaling*1.3, rownames(loadings2), col="black", cex=0.8)
title("Lower 200m")
# add names 
dev.off()


plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=21, cex=1, col="#3A3A3A", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=21, cex=1, col="#3A3A3A", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=21, cex=1, col="#3A3A3A", bg="darkred")
#points(scores[cont, 1], scores[cont, 2], pch=21, cex=1, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=21, cex=1, col="#3A3A3A", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=21, cex=1, col="#3A3A3A", bg="#006B93")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="#b54d04")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
title("Upper 200m")