install.packages("biomformat")
library(biomformat);packageVersion("biomformat")
library(phyloseq)

taxa_names(ps) <- paste0("Seq", seq(ntaxa(ps)))
otu <- as(otu_table(ps), "matrix") # 't' to transform if taxa_are_rows=FALSE
otu_biom <- make_biom(data=otu)
write_biom(otu_biom, "otu.biom")

install.packages("phylotools")
library("phylotools")
fasta.df <- read.fasta("sequences.fna")

fasta.df$seq.name <- row.names(otu) 
dat2fasta(fasta.df, outfile = "rep-seq.fna")



