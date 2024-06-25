# load in libraries
library("dplyr")
library("tidyr")
library("phyloseq")

# read in metadata
org_metadata <- read.delim("required_files/temp_metadata.csv", sep=",", header=TRUE, row.names="sample_name") 

#org_metadata <- org_metadata_filter
distinct <- org_metadata %>% distinct(Depth, Station) # 
fil_metadata <- org_metadata[row.names(org_metadata) %in% row.names(distinct),] # REMOVING DUPLICATES IN THE PCA

## data culling 
org_metadata <- filter(fil_metadata, Sample.Control == "True.Sample") # use tidyr to select "real samples" (non-blanks)
org_metadata <- filter(org_metadata, Station != "STN056b") 
# metadata <- org_metadata[-(which(org_metadata$Station %in% c("STN198", "STN153", "STN056b", "STN012"))),] # removes station 198 and 153 for dotson analysis
metadata <- org_metadata[ -c( 1, 26:29)] # remove filter-related stuff
#metadata <- metadata[-(which(metadata$Station %in% c("STN198", "STN056b"))),] 
metadata <- metadata[ -c( 1, 2:18)] # remove barcode seq/uneeded stuff
metadata <- dplyr::select(org_metadata, 1:12, 17, 19, 21, 23) # select numeric + siderophore column + watertype

#metadata[is.na(metadata)] <- 0 # replace NAs with 0's for PCA
#list_true <- replace_na(metadata$True_Flow, "Other") #replace NA with "Other" for coloring
#metadata$True_Flow <- list_true # changing actual column in dataframe
# convert character columns (except sampleID) to numeric for analysis 

# remove unneeded columns (Oxygen, FIECO, Par, Siderophore)
metadata <- dplyr::select(metadata,-c(Iron))
metadata <- metadata %>% mutate_at(1:12, as.numeric) 

metadata <- metadata %>% rename(c(Nitrate = Lab_NO3, 
                                  Nitrite = Lab_NO2, 
                                  Ammonium = Lab_NH4, 
                                  #Chlorophyll = Chl_a, 
                                  Oxygen = Sb_Oxygen,
                                  Phosphate = Lab_PO4))
# remove NA for iron analysis
metadata <- na.omit(metadata)
metadata <- filter(metadata, Iron != "NA")
metadata <- filter(metadata, DOC != "NA")
metadata <- filter(metadata, Location != "Cont_Shelf")


## PCA analysis from stratigrafia.org
# Run the PCA
# set vectors for inflow
Inflow <- metadata$True_Flow == "Inflow"
Outflow <- metadata$True_Flow == "Outflow"
Other <- metadata$True_Flow == "Other"

# set vectors
CDW <- metadata$watertype == "CDW"
WW_CDW <- metadata$watertype == "WW-CDW"
WW <- metadata$watertype == "WW"
AASW_WW <- metadata$watertype == "AASW-WW"
AASW <- metadata$watertype == "AASW"
Other <- metadata$watertype == "Other"

# set vectors for Location
unique(metadata$Location)
open <- metadata$Location == "Open_polynya"
dotson <- metadata$Location == "Dotson"
east <- metadata$Location == "Eastern_CC"
#cont <- metadata$Location == "Cont_Shelf"
west <- metadata$Location == "Western_CC"
getz <- metadata$Location == "Getz"

# set vectors for pos_in_polynya
open <- metadata$Pos_in_polynya == "Open_polynya_surf"
dotson <- metadata$Pos_in_polynya == "Dotson_CC"
east <- metadata$Pos_in_polynya == "Eastern_CC"
#cont <- metadata$Pos_in_polynya == "Cont_Shelf"
west <- metadata$Pos_in_polynya == "Western_CC"
getz <- metadata$Pos_in_polynya == "Getz_CC"
inflow <- metadata$Pos_in_polynya == "Inflow"
outflow <- metadata$Pos_in_polynya == "Outflow"

# high/low iron
high_iron <- metadata$Iron_Level == "High"
low_iron <- metadata$Iron_Level == "Low"

# set vectors for More depth threshold
unique(org_metadata$More_Depth_Threshold)
bottom <- metadata$More_Depth_Threshold == "Bottom"
mid_bottom <- metadata$More_Depth_Threshold == "Mid-Bottom"
mid <- metadata$More_Depth_Threshold == "Mid"
mid_surface <- metadata$More_Depth_Threshold == "Mid-Surface"
surface <- metadata$More_Depth_Threshold == "Surface"

metadata <- dplyr::select(metadata, 1:11) # select only numerical columns
# create vectors from orginial metadata file for watertype
metadataPca <- prcomp(metadata, scale.=TRUE)

# plot sample scores
dev.new(height=7, width=7)
biplot(metadataPca, cex=0.7)
dev.off()

sd <- metadataPca$sdev
loadings <- metadataPca$rotation
rownames(loadings) <- colnames(metadata)
scores <- metadataPca$x

# use scree plot to reduce dimensionality
var <- sd^2
varPercent <- var/sum(var) * 100
barplot(varPercent, xlab="PC", ylab="Percent Variance", names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col="gray")
abline(h=1/ncol(metadata)*100, col="red")

# determine % variance explained
varPercent[1:4]
sum(varPercent[1:4])

# table for the loadings
loadings
sqrt(1/ncol(metadata))

scaling <- 4.5
textNudge <- 1.2

# fix plot
pdf(file = "ordination_scripts/graphics/PCA_watermass_with_DOC.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[CDW, 1], scores[CDW, 2], pch=21, cex=1, col="#3A3A3A", bg="red2")
points(scores[WW_CDW, 1], scores[WW_CDW, 2], pch=21, cex=1, col="#3A3A3A", bg="blueviolet")
points(scores[WW, 1], scores[WW, 2], pch=21, cex=1, col="#3A3A3A", bg="dodgerblue")
points(scores[AASW_WW, 1], scores[AASW_WW, 2], pch=21, cex=1, col="#3A3A3A", bg="aquamarine3")
points(scores[AASW, 1], scores[AASW, 2], pch=21, cex=1, col="#3A3A3A", bg="darkgreen")
points(scores[Other, 1], scores[Other, 2], pch=21, cex=1, col="#3A3A3A", bg="gray")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="black")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*textNudge, rownames(loadings), col="red4", cex=0.7)
# add names 
text(4, 3, "CDW", col="red2")
text(3.4, -2, "WW-CDW", col="blueviolet")
text(-2, 3, "WW", col="dodgerblue")
text(-3, -2, "AASW-WW", col="aquamarine3")
text(-5, 3, "AASW", col="darkgreen")
dev.off()


#Location plot for PCA!
pdf(file = "ordination_scripts/graphics/PCA_location_with_DOC.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=21, cex=1, col="#3A3A3A", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=21, cex=1, col="#3A3A3A", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=21, cex=1, col="#3A3A3A", bg="darkred")
points(scores[cont, 1], scores[cont, 2], pch=21, cex=1, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=21, cex=1, col="#3A3A3A", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=21, cex=1, col="#3A3A3A", bg="#006B93")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="black")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
# add names 
dev.off()
text(4, 3, "CDW", col="red2")
text(3.4, -2, "WW-CDW", col="blueviolet")
text(-2, 3, "WW", col="dodgerblue")
text(-3, -2, "AASW-WW", col="aquamarine3")
text(-5, 3, "AASW", col="darkgreen")


pdf(file = "ordination_scripts/graphics/PCA_pos_in_polynya_with_DOC.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[open, 1], scores[open, 2], pch=21, cex=1, col="#3A3A3A", bg="#09A20D")
points(scores[dotson, 1], scores[dotson, 2], pch=21, cex=1, col="#3A3A3A", bg="#5AD0FC")
points(scores[east, 1], scores[east, 2], pch=21, cex=1, col="#3A3A3A", bg="darkred")
points(scores[cont, 1], scores[cont, 2], pch=21, cex=1, col="#3A3A3A", bg="#A3DCA5")
points(scores[west, 1], scores[west, 2], pch=21, cex=1, col="#3A3A3A", bg="red")
points(scores[getz, 1], scores[getz, 2], pch=21, cex=1, col="#3A3A3A", bg="#006B93")
points(scores[inflow, 1], scores[inflow, 2], pch=21, cex=1, col="#3A3A3A", bg="yellow")
points(scores[outflow, 1], scores[outflow, 2], pch=21, cex=1, col="#3A3A3A", bg="orange")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="black")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
dev.off()

pdf(file = "ordination_scripts/graphics/PCA_high_low_with_DOC.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[high_iron, 1], scores[high_iron, 2], pch=21, cex=1, col="#3A3A3A", bg="brown")
points(scores[low_iron, 1], scores[low_iron, 2], pch=21, cex=1, col="#3A3A3A", bg="yellow")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="black")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*1.3, rownames(loadings), col="black", cex=0.8)
# add names 
dev.off()

pdf(file = "ordination_scripts/graphics/PCA_depth_threshold_below_100m.pdf", width = 6, height = 7) 
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[bottom, 1], scores[bottom, 2], pch=16, cex=0.7, col="darkred")
points(scores[mid_bottom, 1], scores[mid_bottom, 2], pch=16, cex=0.7, col="red")
points(scores[mid, 1], scores[mid, 2], pch=16, cex=0.7, col="dodgerblue")
points(scores[mid_surface, 1], scores[mid_surface, 2], pch=16, cex=0.7, col="aquamarine3")
points(scores[surface, 1], scores[surface, 2], pch=16, cex=0.7, col="darkgreen")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]*scaling, length=0.1, angle=20, col="red4")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*textNudge, rownames(loadings), col="red4", cex=0.7)
# add names 
dev.off()
### Same thing but for > 100m
# select only below 100m
org_metadata <- org_metadata %>% mutate_at("Depth", as.numeric) 
org_metadata_filter <- filter(org_metadata, Depth >= 100)
org_metadata_100 <- filter(org_metadata, CTD_Depth >= 100)
## remove Par/Chl_a (b/c they are all 0)
org_metadata_filter <- select(org_metadata_filter,-Chlorophyll)

metadataPca2 <- prcomp(org_metadata_filter, scale.=TRUE)

sd2 <- metadataPca2$sdev
loadings2 <- metadataPca2$rotation
rownames(loadings) <- colnames(org_metadata_filter)
scores2 <- metadataPca2$x

# create vectors from orginial metadata file for watertype
CDW2 <- org_metadata_100$watertype == "CDW"
WW2 <- org_metadata_100$watertype == "WW"
AASW2 <- org_metadata_100$watertype == "AASW"
Other2 <- org_metadata_100$watertype == "Other"

# fix plot
pdf(file = "graphics/iron_PCA.pdf", width = 6, height = 7) 
plot(scores2[, 1], scores2[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores2[CDW2, 1], scores2[CDW2, 2], pch=16, cex=0.7, col="red2")
points(scores2[WW2, 1], scores2[WW2, 2], pch=16, cex=0.7, col="dodgerblue")
points(scores2[AASW2, 1], scores2[AASW2, 2], pch=16, cex=0.7, col="seagreen")
points(scores2[Other2, 1], scores2[Other2, 2], pch=16, cex=0.7, col="gray")
arrows(0, 0, loadings2[, 1]* scaling, loadings2[, 2]* scaling, length=0.1, angle=20, col="red4")
text(loadings2[, 1]*scaling*textNudge, loadings2[, 2]*scaling*textNudge, rownames(loadings2), col="red4", cex=0.7)
dev.off()

## SPLIT GRAPHIC
pdf("graphics/split_iron_siderophore_PCA.pdf", width = 8, height = 6)
par(mfrow=c(1,2))
## plot 1
plot(scores[, 1], scores[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores[CDW, 1], scores[CDW, 2], pch=16, cex=0.7, col="red2")
points(scores[WW, 1], scores[WW, 2], pch=16, cex=0.7, col="dodgerblue")
points(scores[AASW, 1], scores[AASW, 2], pch=16, cex=0.7, col="seagreen")
points(scores[Other, 1], scores[Other, 2], pch=16, cex=0.7, col="gray")
arrows(0, 0, loadings[, 1]* scaling, loadings[, 2]* scaling, length=0.1, angle=20, col="red4")
text(loadings[, 1]*scaling*textNudge, loadings[, 2]*scaling*textNudge, rownames(loadings), col="red4", cex=0.7)
# add names 
text(4, 3.5, "CDW", col="red2")
text(-2, 3, "WW", col="dodgerblue")
text(-4, 4.5, "AASW", col="seagreen")
text(2.3, -2, "Other", col="darkgray")

## plot 2
plot(scores2[, 1], scores2[, 2], xlab="PC 1", ylab="PC 2", type="n", asp=1, las=1)
points(scores2[CDW2, 1], scores2[CDW2, 2], pch=16, cex=0.7, col="red2")
points(scores2[WW2, 1], scores2[WW2, 2], pch=16, cex=0.7, col="dodgerblue")
points(scores2[AASW2, 1], scores2[AASW2, 2], pch=16, cex=0.7, col="seagreen")
points(scores2[Other2, 1], scores2[Other2, 2], pch=16, cex=0.7, col="gray")
arrows(0, 0, loadings2[, 1]* scaling, loadings2[, 2]* scaling, length=0.1, angle=20, col="red4")
text(loadings2[, 1]*scaling*textNudge, loadings2[, 2]*scaling*textNudge, rownames(loadings2), col="red4", cex=0.7)

dev.off()

pdf("graphics/split_iron_siderophore_PCA.pdf", width = 8, height = 6)
par(mfrow=c(1,2))
siderophore_all
siderophore_dotson 
dev.off()


# for future use 
text(-2.7, .5, "Outflow", col="red2")
text(1, -1.5, "Inflow", col="dodgerblue")
text(2.3, -2, "Other", col="darkgray")