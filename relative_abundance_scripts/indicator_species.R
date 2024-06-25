install.packages("indicspecies")
library(indicspecies)

#### FREE-LIVING #####

data_free <- taxa_names()


data_part <- subset_samples(data_part, Iron_Level != "NA")
data_part <- subset_samples(data_part, Location != "Cont_Shelf")
data_part <- subset_samples(data_part, Location != "Open_polynya")


write.csv(data_part@otu_table,'M3.SP.otus.csv')


taxa_part <- as.data.frame(data_part@tax_table)
taxa_part$taxon <- rownames(taxa_part)
#Import phyloseq OTU table as an OTU table/dataframe
SpOTU<-read.csv('M3.SP.otus.csv')

#do some shuffling of the OTU table
SpOTUFlip <- as.data.frame(t(SpOTU)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip) <- as.matrix(SpOTUFlip[1, ]) # renames columns
SpOTUFlip<- SpOTUFlip[-1, ] #removes first row
SpOTUFlip_num<-as.data.frame(lapply(SpOTUFlip, as.numeric)) #convert from character to number
SpOTUFlip_num$sample_name<-row.names(SpOTUFlip) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata<-ps_part@sam_data #read in metadata
head(metadata) # check

## Join based on SampleID
SpOTU_Final<-left_join(SpOTUFlip_num, metadata, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus = SpOTU_Final[,1:243] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat = SpOTU_Final$Sample.Control #the metadata column group you care about

indisp=multipatt(x=SPotus, cluster=SPwat, func = "r.g", print.perm = TRUE, control = how(nperm=9999))



library(data.table)
#extract table of stats
indisp.sign<-as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
SPind <- indisp.sign[p.value.bh<=0.05, ]

SPind$taxon <- SPind$rn

merged_df <- merge(SPind, taxa_free, by = "taxon")

write.csv(merged_df,'Location_no_cont_open_particle.csv')
summary(SPind)


#### PARTICLE-ASSOCIATED

write.csv(ps_part@otu_table,'M3.SP.otus.csv')

#Import phyloseq OTU table as an OTU table/dataframe
SpOTU<-read.csv('M3.SP.otus.csv')

#do some shuffling of the OTU table
SpOTUFlip <- as.data.frame(t(SpOTU)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip) <- as.matrix(SpOTUFlip[1, ]) # renames columns
SpOTUFlip<- SpOTUFlip[-1, ] #removes first row
SpOTUFlip_num<-as.data.frame(lapply(SpOTUFlip, as.numeric)) #convert from character to number
SpOTUFlip_num$sample_name<-row.names(SpOTUFlip) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata<-ps_part@sam_data #read in metadata
head(metadata) # check

## Join based on SampleID
SpOTU_Final<-left_join(SpOTUFlip_num, metadata, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata
#SpOTU_Final <- subset(SpOTU_Final, Coastal_Current_Name != "NA")

SPotus = SpOTU_Final[,1:1870] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat = SpOTU_Final$Event #the metadata column group you care about

indisp=multipatt(x=SPotus, cluster=SPwat, func = "r.g", print.perm = TRUE, control = how(nperm=9999))


indisp.sign<-as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")] # https://stats.stackexchange.com/questions/370724/indiscpecies-multipatt-and-overcoming-multi-comparrisons/401277#401277
#now can select only the indicators with adjusted significant p-values
indisp.sign <- indisp.sign[p.value.bh<=0.05, ]

write.csv(indisp.sign,'depth_threshold_particle.csv')

summary(SPind)