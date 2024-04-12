library(oce)
#install.packages("ocedata")
library(ocedata)

#importing ctd data from artemis files
base_path <- "~/Documents/Research/Ordination analysis R scripts/ARTEMIS_github" # needed to use full path - so change this or remove from "files" vector to get the actual path
subdirectory <- "ocean_section_plots/required_files/ctd" # if downloaded from github
files <- dir(paste(base_path, subdirectory, sep = "/"), full.names = TRUE, pattern = ".cnv") #works

# loop through list of the files to get the specific cast files
casts <- list()
for (ifile in 1:length(files)) {
  casts[[ifile]] <- read.ctd.sbe(files[ifile])
}
