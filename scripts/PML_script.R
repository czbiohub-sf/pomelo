### PoMeLo pipeline for examining genome streamlining in two bacterial groups

################################################################################################
########################### ONLY RUN THESE LINES THE FIRST TIME ################################
# pkgs = c("igraph","RColorBrewer", "hexbin", "scales","grid", "lattice", "gdata", "gridExtra",
#          "ape","reshape2", "ggplot2", "seqinr", "phangorn", "fs", "hash","ggdendro",
#          "phytools","openxlsx","coop","tidyverse", "rstudioapi", "aplot") # package names
# install.packages(pkgs)
# install.packages("BiocManager", repos = "https://cloud.r-project.org")
# library(BiocManager, warn.conflicts = FALSE)
# BiocManager::install("remotes")
# BiocManager::install("YuLab-SMU/treedataverse")
################################################################################################

library(igraph, warn.conflicts = FALSE)
library(RColorBrewer, warn.conflicts = FALSE)
library(hexbin, warn.conflicts = FALSE)
library(scales, warn.conflicts = FALSE)
library(grid, warn.conflicts = FALSE)
library(lattice, warn.conflicts = FALSE)
suppressPackageStartupMessages(library(gdata, warn.conflicts = FALSE))
library(gridExtra, warn.conflicts = FALSE)
library(ape, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
suppressPackageStartupMessages(library(hash, warn.conflicts = FALSE))
library(ggplot2, warn.conflicts = FALSE)
library(seqinr, warn.conflicts = FALSE)
library(phangorn, warn.conflicts = FALSE)
library(ggdendro, warn.conflicts = FALSE)
suppressPackageStartupMessages(library(phytools, warn.conflicts = FALSE))
library(openxlsx, warn.conflicts = FALSE)
library(coop, warn.conflicts = FALSE)
Sys.sleep(1)
suppressPackageStartupMessages(library(treedataverse, warn.conflicts = FALSE))
library(aplot, warn.conflicts = FALSE)
Sys.sleep(1)
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE))
library(fs, warn.conflicts = FALSE)
library(rstudioapi, warn.conflicts = FALSE)

## better to remove all files from previous runs first
rm(list=ls())

###
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# dir.create("~/pomelo_outputs")
# ## this works ONLY IF YOU HAVE CREATED THE DIRECTORY FIRST...
# setwd("~/pomelo_outputs")
# data_dir0 <- "~/pomelo_outputs"

## allow user to select working directory for outputs (similar code before to avoid tcltk for PCs)
print("Select the directory you wish to use for PoMeLo outputs")
Sys.sleep(1)
data_dir <- rstudioapi::selectDirectory(caption = "Select the directory you wish to use for PoMeLo outputs", label = "Select output directory")
setwd(data_dir)
# request the path to an existing .csv file on disk
# path <- rstudioapi::selectFile(caption = "Select CSV File",
#                                filter = "CSV Files (*.csv)",
#                                existing = TRUE)

## adding options line to remove all of the summarize messages per https://rstats-tips.net/2020/07/31/get-rid-of-info-of-dplyr-when-grouping-summarise-regrouping-output-by-species-override-with-groups-argument/
options(dplyr.summarise.inform = FALSE)

###########################################################################################################################
###########################################################################################################################
######### CODE TO PULL FROM LARGE LIST OF TAB FILES TAKEN FROM BV-BRC VIA FTP 

###  NOW - SPECIFIC DIRECTIONS FOR READING LISTS OF GENOME ID'S FROM BV-BRC DATABASE, THEN DOWNLOADING CSV FILE...
### directions:

# Workflow:
#   Go to BV-BRC website and browse taxonomy for groups you are interested in.
# https://www.bvbrc.org/view/Taxonomy/2
# 
## UPDATE - NO LONGER NEED REFERENCE SET
# This code will compare groups A & B of your choosing, as well as including 'missing' genes by including a
# reference set of genomes from the full genome_summary file 
## NOTE IN CODE target group & non-target group ARE INSTEAD CALLED GROUP A & GROUP B

# For your group A selection, and your group B selection, and your reference selection. Do the following:
#   Select all of your wanted genomes in the page, then download these as a csv file. NOTE
# THAT YOU NEED TO CLICK THE DOWNLOAD BUTTON ON THE RIGHT IN GREEN TO DOWNLOAD ONLY SELECTED GENOMES.
# This will save a file called PATRIC_genome.csv. Change the name of the file to match which dataset it is
# (so each of them to PATRIC_genome_listA.csv, PATRIC_genome_listB.csv, & PATRIC_genome_reference.csv for example).
# For the latter phylogeny steps, we recommend that you also create a Genome Group (or 2 for target & non-target ) in BV-BRC 
# with this full list - this makes it easy to find your list of genomes.
# 
# Next run the R code, and at the prompts select the correct csv file!

# For the last half the pipeline, you will need a corresponding phylogeny of the species in both of your groups.
# The easiest way to do this is to select your Genome Group (or both target & non-target groups if they are separated),
# and from the list of genomes in your group, select a subset - ONLY 1 GENOME PER SPECIES.
# Create a new Genome Group from this selection, and then go to the Phylogenetic Tree tool in BV-BRC, and create a phylogeny of
# just these genomes. You will need to download both this selected set of genomes (as another .csv file), 
# as well as the completed phylogeny (the .nwk file).


## new ftp site with PATRIC name change
# For example, below is the genome directory for Escherichia coli MG1655 genome.
# 
# ftp://ftp.bvbrc.org/genomes/511145.12

## so replace patricbrc.org with bv-brc.org
## ftp://ftp.patricbrc.org/genomes/511145.12
## ftp://ftp.bvbrc.org/genomes/511145.12

#### if re-running, remove files from this run...
# rm(list = ls(pattern='paths_'))
# rm(list = ls(pattern='table_'))
# rm(list = ls(pattern='final_'))
# rm(list = ls(pattern='ranking'))
# 
# rm(list = ls(pattern='combine'))
# rm(list = ls(pattern='diff_'))
# rm(list = ls(pattern='full_dataset'))
# rm(list = ls(pattern='genes_in'))
# rm(list = ls(pattern='mapping_'))
# rm(list = ls(pattern='total_'))
# rm(list = ls(pattern='test_'))
# rm(list = ls(pattern='pathwaysby'))
# 
# ## better to remove all files from previous runs first! see above
# rm(list=ls())

## also removing any pathway files from previous run...change to data_dir via paste0
# unlink("~/pomelo_outputs/listA.*")
# unlink("~/pomelo_outputs/listB.*")
listofAs <- paste0(data_dir,"/listA.*" )
unlink(listofAs)
listofBs <- paste0(data_dir,"/listB.*" )
unlink(listofBs)

### BEGIN CODE TO SELECT GROUPS A & B
## note need to import genome_id as character (pulling from genome summary import below...) & per https://github.com/tidyverse/readr/issues/148
#listAf <- read_csv(tk_choose.files(caption = "Select your List A (target group) .csv file from BV-BRC"), col_types = cols(.default = "c"))

# using rstudioapi - request the path to an existing .csv file on disk
print("Select your List A (target group) .csv file from BV-BRC")
Sys.sleep(1)
listAf <- read_csv(rstudioapi::selectFile(caption = "Select your List A (target group) .csv file from BV-BRC", label = "Select target group .csv", path = data_dir, existing = TRUE, filter = "CSV Files (*.csv)"), col_types = cols(.default = "c"))

## moving select to before the gsub commands...
listA <- listAf %>% select("Genome ID")

##########################################
## idea - to check genome name & species fields, and remove when inconsistent to remove unresolved taxonomic issues...rename
listAf <- listAf %>% rename(genomename = "Genome Name")
listAf$genomename <- gsub("Candidatus ","",listAf$genomename)
listAf$genomename <- gsub("uncultured ","",listAf$genomename)
listAf$genomename <- gsub("sp. ","sp_",listAf$genomename)
## update endosymbiont name changes, OMZ for Treponema & more acronyms
listAf$genomename <- gsub("endosymbiont wPip_Mol of ","endo_",listAf$genomename)
listAf$genomename <- gsub("endosymbiont of ","endo_",listAf$genomename)
listAf$genomename <- gsub("endosymbiont strain TRS of ","endo_",listAf$genomename)
listAf$genomename <- gsub("strain ","strain_",listAf$genomename)
listAf$genomename <- gsub("OMZ ","OMZ_",listAf$genomename)
listAf$genomename <- gsub("ATCC ","ATCC_",listAf$genomename)
listAf$genomename <- gsub("PCC ","PCC_",listAf$genomename)
listAf$genomename <- gsub("FDAARGOS ","FDAARGOS_",listAf$genomename)
listAf$genomename <- gsub("NCTC ","NCTC_",listAf$genomename)
listAf$genomename <- gsub("MAG ","MAG_",listAf$genomename)
listAf$genomename <- gsub("USDA ","USDA_",listAf$genomename)

# ## universal version
# listAf$genomename <- gsub(" ", "_", gsub("endo_(\\S+)", "endo_\\1 ", listAf$genomename))
# listAf$genomename <- gsub("__","_",listAf$genomename)
# # In this code, gsub() function is used twice. The inner gsub() is used to find matches of the pattern "endo_(\S+)" in the 'colla' column of the 'dabba' dataframe. The (\\S+) captures all non-whitespace characters after "endo_" in the match. The outer gsub() is then used to replace the space character with an underscore "" in the matched patterns, using the captured group \\1 to retain the original characters after "endo" and include the space character. The modified 'colla' column is updated in the 'dabba' dataframe. Finally, the modified 'colla' column is printed using print(dabba$colla).
listAf$genomename <- gsub("endo_Brugia ","endo_Brugia_",listAf$genomename)
listAf$genomename <- gsub("endo_Rhopalodia ","endo_Rhopalodia_",listAf$genomename)
listAf$genomename <- gsub("endo_Epithemia ","endo_Epithemia_",listAf$genomename)
listAf$genomename <- gsub("endo_Onchocerca ","endo_Onchocerca_",listAf$genomename)
listAf$genomename <- gsub("endo_Drosophila ","endo_Drosophila_",listAf$genomename)
listAf$genomename <- gsub("endo_Amblyomma ","endo_Amblyomma_",listAf$genomename)
listAf$genomename <- gsub("endo_Culex ","endo_Culex_",listAf$genomename)
listAf$genomename <- gsub("endo_Aedes ","endo_Aedes_",listAf$genomename)
listAf$genomename <- gsub("endo_Donacia ","endo_Donacia_",listAf$genomename)
listAf$genomename <- gsub("endo_Macroplea ","endo_Macroplea_",listAf$genomename)
listAf$genomename <- gsub("endo_Plateumaris ","endo_Plateumaris_",listAf$genomename)

listAf$genomename <- gsub("[","",listAf$genomename, fixed = TRUE)
listAf$genomename <- gsub("]","",listAf$genomename, fixed = TRUE)
listAf$genomename <- gsub('"',"",listAf$genomename, fixed = TRUE)
listAf$genomename <- gsub("'","",listAf$genomename, fixed = TRUE)

unique(listAf$genomename)

## then pull just genus & species
## create new column with shorter description
listAf <- listAf %>% separate(genomename, into = c("genus", "species"), sep = " ", remove = FALSE, convert = TRUE, extra = "drop")
listAf <- listAf %>% unite("genusspecies", genus, species, sep = "_", remove = TRUE)

## now repeat for Species field...
listAf$Species <- gsub("Candidatus ","",listAf$Species)
listAf$Species <- gsub("uncultured ","",listAf$Species)
listAf$Species <- gsub("sp. ","sp_",listAf$Species)
## update endosymbiont name changes, OMZ for Treponema & more acronyms
listAf$Species <- gsub("endosymbiont wPip_Mol of ","endo_",listAf$Species)
listAf$Species <- gsub("endosymbiont of ","endo_",listAf$Species)
listAf$Species <- gsub("endosymbiont strain TRS of ","endo_",listAf$Species)
listAf$Species <- gsub("strain ","strain_",listAf$Species)
listAf$Species <- gsub("OMZ ","OMZ_",listAf$Species)
listAf$Species <- gsub("ATCC ","ATCC_",listAf$Species)
listAf$Species <- gsub("PCC ","PCC_",listAf$Species)
listAf$Species <- gsub("FDAARGOS ","FDAARGOS_",listAf$Species)
listAf$Species <- gsub("NCTC ","NCTC_",listAf$Species)
listAf$Species <- gsub("MAG ","MAG_",listAf$Species)
listAf$Species <- gsub("USDA ","USDA_",listAf$Species)

# ## universal version
# listAf$Species <- gsub(" ", "_", gsub("endo_(\\S+)", "endo_\\1 ", listAf$Species))
# listAf$Species <- gsub("__","_",listAf$Species)
# # In this code, gsub() function is used twice. The inner gsub() is used to find matches of the pattern "endo_(\S+)" in the 'colla' column of the 'dabba' dataframe. The (\\S+) captures all non-whitespace characters after "endo_" in the match. The outer gsub() is then used to replace the space character with an underscore "" in the matched patterns, using the captured group \\1 to retain the original characters after "endo" and include the space character. The modified 'colla' column is updated in the 'dabba' dataframe. Finally, the modified 'colla' column is printed using print(dabba$colla).

listAf$Species <- gsub("endo_Brugia ","endo_Brugia_",listAf$Species)
listAf$Species <- gsub("endo_Rhopalodia ","endo_Rhopalodia_",listAf$Species)
listAf$Species <- gsub("endo_Epithemia ","endo_Epithemia_",listAf$Species)
listAf$Species <- gsub("endo_Onchocerca ","endo_Onchocerca_",listAf$Species)
listAf$Species <- gsub("endo_Drosophila ","endo_Drosophila_",listAf$Species)
listAf$Species <- gsub("endo_Culex ","endo_Culex_",listAf$Species)
listAf$Species <- gsub("endo_Aedes ","endo_Aedes_",listAf$Species)
listAf$Species <- gsub("endo_Amblyomma ","endo_Amblyomma_",listAf$Species)
listAf$Species <- gsub("endo_Donacia ","endo_Donacia_",listAf$Species)
listAf$Species <- gsub("endo_Macroplea ","endo_Macroplea_",listAf$Species)
listAf$Species <- gsub("endo_Plateumaris ","endo_Plateumaris_",listAf$Species)
## also just in Species field change XXX sp. at end of field to null, to force use of the genome name - specify the metacharacter $ after the sp_ to signify the end of the string
listAf$Species <- gsub(".*sp.$","null",listAf$Species)

listAf$Species <- gsub("[","",listAf$Species, fixed = TRUE)
listAf$Species <- gsub("]","",listAf$Species, fixed = TRUE)
listAf$Species <- gsub('"',"",listAf$Species, fixed = TRUE)
listAf$Species <- gsub("'","",listAf$Species, fixed = TRUE)

unique(listAf$Species)

listAf <- listAf %>% mutate(Species = na_if(Species, "null"))

#listAf$Species <- gsub("null ","",listAf$Species)
## instead of gsub, use coalesce??
# dfABy %>% mutate(A = coalesce(A,B)) this is replacing all NAs in cregion with the value from locus
listAf <- listAf %>% mutate(Species = coalesce(Species,genomename))


## then pull just genus & species
## create new column with shorter description
listAf <- listAf %>% separate(Species, into = c("genus", "species"), sep = " ", remove = FALSE, convert = TRUE, extra = "drop")
listAf <- listAf %>% unite("genusspecies2", genus, species, sep = "_", remove = TRUE) %>% relocate(genusspecies2, .after = genusspecies)

#listAf <- listAf %>% mutate(genusspecies2 = na_if(genusspecies2, "NA_NA"))

## now remove if genusspecies2 doesn't = genusspecies update do not remove, but use the latter genusspecies2 in every case!
## further down we will replace the collated paths_groupAandB$genusspecies with genusspecies2
listAfmismatch <- subset(listAf , genusspecies != genusspecies2)
#listAf2 <- subset(listAf , genusspecies == genusspecies2)

## select only Genome ID, then convert to list of characters for downloading:
## moving select to before the gsub commands...
#listA <- listAf %>% select("Genome ID")
# listA <- listAf2 %>% select("Genome ID")
listA <- as.character(unlist(listA))


listAurl <- paste0("ftp://ftp.bvbrc.org/genomes/" ,listA, "/",listA,".PATRIC.pathway.tab" )
destinations <- paste0("listA.",listA,".PATRIC.pathway.tab" )

## time estimate
print("Downloading these files may take as long as")
print(sum(table(destinations)) / 30)
print("minutes")


## modified loop will continue past errors that can't find a particular file for downloading...
for(i in seq(listAurl)){
  skip_to_next <- FALSE
  tryCatch(download.file(listAurl[i], destinations[i], mode="wb", quiet = TRUE, cacheOK = FALSE), error = function(e) { skip_to_next <<- TRUE})
}

#### list B

# listBf <- read_csv(tk_choose.files(caption = "Select your List B (non-target group) .csv file from BV-BRC"), col_types = cols(.default = "c"))

# using rstudioapi - request the path to an existing .csv file on disk
print("Select your List B (non-target group) .csv file from BV-BRC")
Sys.sleep(1)
listBf <- read_csv(rstudioapi::selectFile(caption = "Select your List B (non-target group) .csv file from BV-BRC", label = "Select non-target group .csv", path = data_dir, existing = TRUE, filter = "CSV Files (*.csv)"), col_types = cols(.default = "c"))

## moving select to before the gsub commands...
listB <- listBf %>% select("Genome ID")

##########################################
## idea - to check genome name & species fields, and remove when inconsistent to remove unresolved taxonomic issues...rename
listBf <- listBf %>% rename(genomename = "Genome Name")
listBf$genomename <- gsub("Candidatus ","",listBf$genomename)
listBf$genomename <- gsub("uncultured ","",listBf$genomename)
listBf$genomename <- gsub("sp. ","sp_",listBf$genomename)
## update endosymbiont name changes, OMZ for Treponema & more acronyms
listBf$genomename <- gsub("endosymbiont wPip_Mol of ","endo_",listBf$genomename)
listBf$genomename <- gsub("endosymbiont of ","endo_",listBf$genomename)
listBf$genomename <- gsub("endosymbiont strain TRS of ","endo_",listBf$genomename)
listBf$genomename <- gsub("strain ","strain_",listBf$genomename)
listBf$genomename <- gsub("OMZ ","OMZ_",listBf$genomename)
listBf$genomename <- gsub("ATCC ","ATCC_",listBf$genomename)
listBf$genomename <- gsub("PCC ","PCC_",listBf$genomename)
listBf$genomename <- gsub("FDAARGOS ","FDAARGOS_",listBf$genomename)
listBf$genomename <- gsub("NCTC ","NCTC_",listBf$genomename)
listBf$genomename <- gsub("MAG ","MAG_",listBf$genomename)
listBf$genomename <- gsub("USDA ","USDA_",listBf$genomename)

# ## universal version
# listBf$genomename <- gsub(" ", "_", gsub("endo_(\\S+)", "endo_\\1 ", listBf$genomename))
# listBf$genomename <- gsub("__","_",listBf$genomename)
# # In this code, gsub() function is used twice. The inner gsub() is used to find matches of the pattern "endo_(\S+)" in the 'colla' column of the 'dabba' dataframe. The (\\S+) captures all non-whitespace characters after "endo_" in the match. The outer gsub() is then used to replace the space character with an underscore "" in the matched patterns, using the captured group \\1 to retain the original characters after "endo" and include the space character. The modified 'colla' column is updated in the 'dabba' dataframe. Finally, the modified 'colla' column is printed using print(dabba$colla).

listBf$genomename <- gsub("endo_Brugia ","endo_Brugia_",listBf$genomename)
listBf$genomename <- gsub("endo_Rhopalodia ","endo_Rhopalodia_",listBf$genomename)
listBf$genomename <- gsub("endo_Epithemia ","endo_Epithemia_",listBf$genomename)
listBf$genomename <- gsub("endo_Onchocerca ","endo_Onchocerca_",listBf$genomename)
listBf$genomename <- gsub("endo_Drosophila ","endo_Drosophila_",listBf$genomename)
listBf$genomename <- gsub("endo_Culex ","endo_Culex_",listBf$genomename)
listBf$genomename <- gsub("endo_Aedes ","endo_Aedes_",listBf$genomename)
listBf$genomename <- gsub("endo_Amblyomma ","endo_Amblyomma_",listBf$genomename)
listBf$genomename <- gsub("endo_Donacia ","endo_Donacia_",listBf$genomename)
listBf$genomename <- gsub("endo_Macroplea ","endo_Macroplea_",listBf$genomename)
listBf$genomename <- gsub("endo_Plateumaris ","endo_Plateumaris_",listBf$genomename)

listBf$genomename <- gsub("[","",listBf$genomename, fixed = TRUE)
listBf$genomename <- gsub("]","",listBf$genomename, fixed = TRUE)
listBf$genomename <- gsub('"',"",listBf$genomename, fixed = TRUE)
listBf$genomename <- gsub("'","",listBf$genomename, fixed = TRUE)

unique(listBf$genomename)

## then pull just genus & species
## create new column with shorter description
listBf <- listBf %>% separate(genomename, into = c("genus", "species"), sep = " ", remove = FALSE, convert = TRUE, extra = "drop")
listBf <- listBf %>% unite("genusspecies", genus, species, sep = "_", remove = TRUE)

## now repeat for Species field...
listBf$Species <- gsub("Candidatus ","",listBf$Species)
listBf$Species <- gsub("uncultured ","",listBf$Species)
listBf$Species <- gsub("sp. ","sp_",listBf$Species)
## update endosymbiont name changes, OMZ for Treponema & more acronyms
listBf$Species <- gsub("endosymbiont wPip_Mol of ","endo_",listBf$Species)
listBf$Species <- gsub("endosymbiont of ","endo_",listBf$Species)
listBf$Species <- gsub("endosymbiont strain TRS of ","endo_",listBf$Species)
listBf$Species <- gsub("strain ","strain_",listBf$Species)
listBf$Species <- gsub("OMZ ","OMZ_",listBf$Species)
listBf$Species <- gsub("ATCC ","ATCC_",listBf$Species)
listBf$Species <- gsub("PCC ","PCC_",listBf$Species)
listBf$Species <- gsub("FDAARGOS ","FDAARGOS_",listBf$Species)
listBf$Species <- gsub("NCTC ","NCTC_",listBf$Species)
listBf$Species <- gsub("MAG ","MAG_",listBf$Species)
listBf$Species <- gsub("USDA ","USDA_",listBf$Species)
## also just in Species field change XXX sp. at end of field to null, to force use of the genome name - specify the metacharacter $ after the sp_ to signify the end of the string
listBf$Species <- gsub(".*sp.$","null",listBf$Species)

# ## universal version
# listBf$Species <- gsub(" ", "_", gsub("endo_(\\S+)", "endo_\\1 ", listBf$Species))
# listBf$Species <- gsub("__","_",listBf$Species)
# # In this code, gsub() function is used twice. The inner gsub() is used to find matches of the pattern "endo_(\S+)" in the 'colla' column of the 'dabba' dataframe. The (\\S+) captures all non-whitespace characters after "endo_" in the match. The outer gsub() is then used to replace the space character with an underscore "" in the matched patterns, using the captured group \\1 to retain the original characters after "endo" and include the space character. The modified 'colla' column is updated in the 'dabba' dataframe. Finally, the modified 'colla' column is printed using print(dabba$colla).

listBf$Species <- gsub("endo_Brugia ","endo_Brugia_",listBf$Species)
listBf$Species <- gsub("endo_Rhopalodia ","endo_Rhopalodia_",listBf$Species)
listBf$Species <- gsub("endo_Epithemia ","endo_Epithemia_",listBf$Species)
listBf$Species <- gsub("endo_Onchocerca ","endo_Onchocerca_",listBf$Species)
listBf$Species <- gsub("endo_Drosophila ","endo_Drosophila_",listBf$Species)
listBf$Species <- gsub("endo_Culex ","endo_Culex_",listBf$Species)
listBf$Species <- gsub("endo_Aedes ","endo_Aedes_",listBf$Species)
listBf$Species <- gsub("endo_Amblyomma ","endo_Amblyomma_",listBf$Species)
listBf$Species <- gsub("endo_Donacia ","endo_Donacia_",listBf$Species)
listBf$Species <- gsub("endo_Macroplea ","endo_Macroplea_",listBf$Species)
listBf$Species <- gsub("endo_Plateumaris ","endo_Plateumaris_",listBf$Species)

listBf$Species <- gsub("[","",listBf$Species, fixed = TRUE)
listBf$Species <- gsub("]","",listBf$Species, fixed = TRUE)
listBf$Species <- gsub('"',"",listBf$Species, fixed = TRUE)
listBf$Species <- gsub("'","",listBf$Species, fixed = TRUE)
unique(listBf$Species)

listBf <- listBf %>% mutate(Species = na_if(Species, "null"))

# listBf <- listBf %>%
#   mutate(Species,
#          ~ replace(.,   str_detect(., "null$"), NA))

#listBf$Species <- gsub("null ","",listBf$Species)
## instead of gsub, use coalesce??
# dfABy %>% mutate(A = coalesce(A,B)) this is replacing all NAs in cregion with the value from locus
listBf <- listBf %>% mutate(Species = coalesce(Species,genomename))

## then pull just genus & species
## create new column with shorter description
listBf <- listBf %>% separate(Species, into = c("genus", "species"), sep = " ", remove = FALSE, convert = TRUE, extra = "drop")
listBf <- listBf %>% unite("genusspecies2", genus, species, sep = "_", remove = TRUE) %>% relocate(genusspecies2, .after = genusspecies)


## now remove if genusspecies2 doesn't = genusspecies update do not remove, but use the latter genusspecies2 in every case!
## further down we will replace the collated paths_groupAandB$genusspecies with genusspecies2
#listBf2 <- subset(listBf , genusspecies == genusspecies2)
listBfmismatch <- subset(listBf , genusspecies != genusspecies2)

##########################################


## select only Genome ID, then convert to list of characters for downloading:
# listB <- listB %>% select("Genome.ID") ## first was called "Genome.ID" then "Genome ID"
# listB <- as.character(unlist(listB$Genome.ID))

## moving select to before the gsub commands...
#listB <- listBf %>% select("Genome ID")
#listB <- listBf2 %>% select("Genome ID")
listB <- as.character(unlist(listB))

listBurl <- paste0("ftp://ftp.bvbrc.org/genomes/" ,listB, "/",listB,".PATRIC.pathway.tab" )
destinationsB <- paste0("listB.",listB,".PATRIC.pathway.tab" )
## sum(table(destinations)) + sum(table(destinationsB))

## time estimate
print("Downloading these files may take as long as")
print(sum(table(destinationsB)) / 30)
print("minutes")

for(i in seq(listBurl)){
  skip_to_next <- FALSE
  tryCatch(download.file(listBurl[i], destinationsB[i], mode="wb", quiet = TRUE, cacheOK = FALSE), error = function(e) { skip_to_next <<- TRUE})
}

## new after directly downloading via R
tab_filesA <- fs::dir_ls(data_dir, regexp = "listA.")
tab_filesB <- fs::dir_ls(data_dir, regexp = "listB.")

## this code takes the large list of tab files I pulled and imports them into R, then combines them into 2 large dataframes
## taken from https://www.gerkelab.com/blog/2018/09/import-directory-csv-purrr-readr/

pathways_groupA <- tab_filesA %>% 
  map_dfr(read_tsv, show_col_types = FALSE, col_types = cols("pathway_id" = "i", "refseq_locus_tag" = "c", "genome_id" = "c"))

head(pathways_groupA)

pathways_groupB <- tab_filesB %>% 
  map_dfr(read_tsv, show_col_types = FALSE, col_types = cols("pathway_id" = "i", "refseq_locus_tag" = "c", "genome_id" = "c"))

pathways_groupA$alt_locus_tag <- NULL
pathways_groupB$alt_locus_tag <- NULL

paths_groupAandB <- bind_rows(pathways_groupA, pathways_groupB, .id = "group")

## now very early remove ECs with dashes
paths_groupAandB <- paths_groupAandB[ grep(".-", paths_groupAandB$ec_number, invert = TRUE) , ]


paths_groupAandB$group <- gsub("1","group A",paths_groupAandB$group)
paths_groupAandB$group <- gsub("2","group B",paths_groupAandB$group)
paths_groupAandB$genome_id <- as.character(paths_groupAandB$genome_id)

###################################################################################################################
## one-off code to get 4 missing pathways added back - after line 390...
# paths_groupAandB_topull <- paths_groupAandB %>% group_by(genome_name,ec_number) %>%
#   summarize_all(first) ## update to use across instead of _all
paths_groupAandB_topull <- paths_groupAandB %>% group_by(genome_name,ec_number) %>%
  summarize(across(everything(), first))
## need to pull every instance, so need a big filter for each of the four
paths_groupAandB_topull220 <- paths_groupAandB_topull %>%
  dplyr::filter((ec_number == "1.2.1.38")|
                  (ec_number == "1.4.1.2")|
                  (ec_number == "1.4.1.4")|
                  (ec_number == "2.1.3.3")|
                  (ec_number == "2.3.1.1")|
                  (ec_number == "2.3.1.35")|
                  (ec_number == "2.6.1.1")|
                  (ec_number == "2.6.1.11")|
                  (ec_number == "2.6.1.2")|
                  (ec_number == "2.7.2.2")|
                  (ec_number == "2.7.2.8")|
                  (ec_number == "3.5.1.14")|
                  (ec_number == "3.5.1.16")|
                  (ec_number == "3.5.1.2")|
                  (ec_number == "3.5.1.5")|
                  (ec_number == "3.5.1.54")|
                  (ec_number == "3.5.3.1")|
                  (ec_number == "3.5.3.6")|
                  (ec_number == "4.3.2.1")|
                  (ec_number == "6.3.1.2")|
                  (ec_number == "6.3.4.5")|
                  (ec_number == "6.3.4.6"))
paths_groupAandB_topull220$pathway_id <- 220
paths_groupAandB_topull220$pathway_name <- "Arginine biosynthesis"

paths_groupAandB_topull270 <- paths_groupAandB_topull %>%
  dplyr::filter((ec_number == "1.1.1.37")|
                  (ec_number == "1.1.1.95")|
                  (ec_number == "2.6.1.42")|
                  (ec_number == "2.6.1.52")|
                  (ec_number == "3.5.99.7")|
                  (ec_number == "6.3.2.2")|
                  (ec_number == "6.3.2.3"))
paths_groupAandB_topull270$pathway_id <- 270
paths_groupAandB_topull270$pathway_name <- "Cysteine and methionine metabolism"

paths_groupAandB_topull470 <- paths_groupAandB_topull %>%
  dplyr::filter((ec_number == "1.2.1.26")|
                  (ec_number == "1.4.1.12")|
                  (ec_number == "1.4.3.3")|
                  (ec_number == "2.6.1.21")|
                  (ec_number == "3.5.1.2")|
                  (ec_number == "3.5.4.22")|
                  (ec_number == "4.1.1.20")|
                  (ec_number == "4.3.1.18")|
                  (ec_number == "4.4.1.15")|
                  (ec_number == "5.1.1.13")|
                  (ec_number == "5.1.1.7")|
                  (ec_number == "5.1.1.8")|
                  (ec_number == "6.3.2.4")|
                  (ec_number == "6.3.2.8")|
                  (ec_number == "6.3.2.9"))
paths_groupAandB_topull470$pathway_id <- 470
paths_groupAandB_topull470$pathway_name <- "D-Amino acid metabolism"

paths_groupAandB_topull541 <- paths_groupAandB_topull %>%
  dplyr::filter((ec_number == "1.1.1.133")|
                  (ec_number == "1.1.1.136")|
                  (ec_number == "1.1.1.22")|
                  (ec_number == "1.1.1.271")|
                  (ec_number == "1.1.1.336")|
                  (ec_number == "2.3.1.201")|
                  (ec_number == "2.5.1.56")|
                  (ec_number == "2.5.1.97")|
                  (ec_number == "2.7.7.13")|
                  (ec_number == "2.7.7.23")|
                  (ec_number == "2.7.7.24")|
                  (ec_number == "2.7.7.33")|
                  (ec_number == "2.7.7.43")|
                  (ec_number == "2.7.7.9")|
                  (ec_number == "4.2.1.115")|
                  (ec_number == "4.2.1.45")|
                  (ec_number == "4.2.1.46")|
                  (ec_number == "4.2.1.47")|
                  (ec_number == "5.1.3.13")|
                  (ec_number == "5.1.3.14")|
                  (ec_number == "5.1.3.2")|
                  (ec_number == "5.1.3.23")|
                  (ec_number == "5.3.1.8")|
                  (ec_number == "5.4.2.8"))
paths_groupAandB_topull541$pathway_id <- 541
paths_groupAandB_topull541$pathway_name <- "O-Antigen nucleotide sugar biosynthesis"

##  note we're simply replacing the pathway_id & pathway_name for each of these four
## FINALLY APPEND THESE TO paths_groupAandB
paths_groupAandBn <- bind_rows(paths_groupAandB, paths_groupAandB_topull220, paths_groupAandB_topull270, paths_groupAandB_topull470, paths_groupAandB_topull541)
rm(paths_groupAandB)
paths_groupAandB <- paths_groupAandBn
Sys.sleep(2)
## end one-off code to get 4 missing pathways added back - after line 390...
###################################################################################################################

### this makes new column pathwayid not number but character p1-p1000
paths_groupAandB$pathwayid0 <- "p"
paths_groupAandB <- paths_groupAandB %>% unite("pathwayid", pathwayid0, pathway_id, sep = "", remove = FALSE)
paths_groupAandB$pathwayid0 <- NULL

## also removing some pathways that have been deleted from KEGG: https://www.genome.jp/kegg/docs/upd_map.html
# 1058, 471, 472, 473, 72, 231 - note to avoid removing 720, using name of p72
## note there are still a few pathways >1000 - 1040,1051,1053,1055,1056,1057,1059
paths_groupAandB <- paths_groupAandB[ grep("p1058", paths_groupAandB$pathwayid, invert = TRUE) , ]
paths_groupAandB <- paths_groupAandB[ grep("p471", paths_groupAandB$pathwayid, invert = TRUE) , ]
paths_groupAandB <- paths_groupAandB[ grep("p472", paths_groupAandB$pathwayid, invert = TRUE) , ]
paths_groupAandB <- paths_groupAandB[ grep("p473", paths_groupAandB$pathwayid, invert = TRUE) , ]
paths_groupAandB <- paths_groupAandB[ grep("Synthesis and degradation of ketone bodies", paths_groupAandB$pathway_name, invert = TRUE) , ]
paths_groupAandB <- paths_groupAandB[ grep("p231", paths_groupAandB$pathwayid, invert = TRUE) , ]
paths_groupAandB <- paths_groupAandB[ grep("p4070", paths_groupAandB$pathwayid, invert = TRUE) , ]
paths_groupAandB <- paths_groupAandB[ grep("p4150", paths_groupAandB$pathwayid, invert = TRUE) , ]

###########################################################################################################################
### REFERENCE GENOMES - UPDATE, PULLING FROM MAPPING FILE BELOW
### now remove all of these files...
# unlink("~/pomelo_outputs/listA.*")
# unlink("~/pomelo_outputs/listB.*")

###########################################################################################################################
###########################################################################################################################
## instead of splitting (see below), start analyzing each gene in terms of counts in A vs. B - simply plotting the heat maps
# idea is to recapitulate the PATRIC pathway heatmaps....

## this creates new column with combined pathwayid to the gene ec_number
paths_groupAandB <- paths_groupAandB %>% unite("full_ec_number", pathwayid, ec_number, sep = ".", remove = FALSE)

## ALSO GET GENUS + SPECIES FOR ALL (AND TRANSFORM ENDOSYMBIONT OF & CANDIDATUS)
## ADDING MODIFICATION TO NOT LUMP ALL SPP...
paths_groupAandB$genome_name <- gsub("Candidatus ","",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("uncultured ","",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("sp. ","sp_",paths_groupAandB$genome_name)
## update endosymbiont name changes, OMZ for Treponema & more acronyms
paths_groupAandB$genome_name <- gsub("endosymbiont wPip_Mol of ","endo_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endosymbiont of ","endo_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endosymbiont strain TRS of ","endo_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("strain ","strain_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("OMZ ","OMZ_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("ATCC ","ATCC_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("PCC ","PCC_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("FDAARGOS ","FDAARGOS_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("NCTC ","NCTC_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("MAG ","MAG_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("USDA ","USDA_",paths_groupAandB$genome_name)

unique(paths_groupAandB$genome_name)
# ## universal version
# paths_groupAandB$genome_name <- gsub(" ", "_", gsub("endo_(\\S+)", "endo_\\1 ", paths_groupAandB$genome_name))
# paths_groupAandB$genome_name <- gsub("__","_",paths_groupAandB$genome_name)
# # In this code, gsub() function is used twice. The inner gsub() is used to find matches of the pattern "endo_(\S+)" in the 'colla' column of the 'dabba' dataframe. The (\\S+) captures all non-whitespace characters after "endo_" in the match. The outer gsub() is then used to replace the space character with an underscore "" in the matched patterns, using the captured group \\1 to retain the original characters after "endo" and include the space character. The modified 'colla' column is updated in the 'dabba' dataframe. Finally, the modified 'colla' column is printed using print(dabba$colla).

paths_groupAandB$genome_name <- gsub("endo_Brugia ","endo_Brugia_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Rhopalodia ","endo_Rhopalodia_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Epithemia ","endo_Epithemia_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Onchocerca ","endo_Onchocerca_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Drosophila ","endo_Drosophila_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Culex ","endo_Culex_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Aedes ","endo_Aedes_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Amblyomma ","endo_Amblyomma_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Donacia ","endo_Donacia_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Macroplea ","endo_Macroplea_",paths_groupAandB$genome_name)
paths_groupAandB$genome_name <- gsub("endo_Plateumaris ","endo_Plateumaris_",paths_groupAandB$genome_name)

paths_groupAandB$genome_name <- gsub("[","",paths_groupAandB$genome_name, fixed = TRUE)
paths_groupAandB$genome_name <- gsub("]","",paths_groupAandB$genome_name, fixed = TRUE)
paths_groupAandB$genome_name <- gsub('"',"",paths_groupAandB$genome_name, fixed = TRUE)
paths_groupAandB$genome_name <- gsub("'","",paths_groupAandB$genome_name, fixed = TRUE)

unique(paths_groupAandB$genome_name)

## then pull just genus & species
## create new column with shorter description
paths_groupAandB <- paths_groupAandB %>% separate(genome_name, into = c("genus", "species"), sep = " ", remove = FALSE, convert = TRUE, extra = "drop")
paths_groupAandB <- paths_groupAandB %>% unite("genusspecies", genus, species, sep = "_", remove = FALSE)
paths_groupAandB$species <- NULL

unique(paths_groupAandB$genusspecies)
n_distinct(paths_groupAandB$genusspecies) #119 here

## 2023 if we are going to pull in names from species column...do it here
listAtojoin <- listAf
listBtojoin <- listBf
listAtojoin <- listAtojoin %>% rename(genome_id = "Genome ID")
listBtojoin <- listBtojoin %>% rename(genome_id = "Genome ID")

listAtojoin1 <- listAtojoin %>% select(genome_id, genusspecies2)
listBtojoin1 <- listBtojoin %>% select(genome_id, genusspecies2)
listAandBtojoin <- bind_rows(listAtojoin1, listBtojoin1)
## include genome size to use instead of genome_summary below
listAtojoin2 <- listAtojoin %>% select(genome_id, Size)
listBtojoin2 <- listBtojoin %>% select(genome_id, Size)
listAandBtojoin2 <- bind_rows(listAtojoin2, listBtojoin2) %>% 
  rename(genome_length = Size)
listAandBtojoin2$genome_length <- as.integer(listAandBtojoin2$genome_length)


## now join
paths_groupAandB <- left_join(paths_groupAandB, listAandBtojoin)
paths_groupAandB <- paths_groupAandB %>% 
  rename(genusspecies1 = genusspecies) %>%
  rename(genusspecies = genusspecies2) %>%
  relocate(genusspecies, .before = genusspecies1)
## replace NAs...
paths_groupAandB <- paths_groupAandB %>% mutate(genusspecies = coalesce(genusspecies,genusspecies1))


## now new columns based on updated species
paths_groupAandB <- paths_groupAandB %>% unite("genusspeciespathway", genusspecies, pathwayid, sep = "_", remove = FALSE) %>% relocate(genusspeciespathway, .after = genusspecies)
paths_groupAandB <- paths_groupAandB %>% unite("genuspathway", genus, pathwayid, sep = "_", remove = FALSE) %>% relocate(genuspathway, .after = genus)


## some stats
# paths_groupAandB %>% group_by(group) %>% summarise_at(vars(genusspecies), n_distinct) ## update to use across instead of _at
paths_groupAandB %>% group_by(group) %>% summarise(across(genusspecies, n_distinct))
#btw summarise(across(genusspecies), n_distinct) & summarise(across(everything), first) are wrong, need to move the parenthses

#or for summary of all
#paths_groupAandB %>% group_by(group) %>% summarise_all(n_distinct)
paths_groupAandB %>% group_by(group) %>% summarize(across(everything(), n_distinct))
paths_groupAandB_counts <- paths_groupAandB %>% group_by(group) %>% summarize(across(everything(), n_distinct))

paths_groupAandB <- paths_groupAandB %>% ungroup()
paths_groupAandB_counts$genusspecies[1]
paths_groupAandB_counts$genusspecies[2]
unique(paths_groupAandB$genusspecies)
n_distinct(paths_groupAandB$genusspecies)


################################################################################################################
### this is summarizing above table into fewer rows (by group + genusspecies + ec_number) - update to only group_by species & ec number...
## update to summarize_if: summarize_if(is.numeric, diff) becomes summarize(across(where(is.numeric), diff), .groups = "keep")

paths_groupAandB_stats0 <- paths_groupAandB %>%
  group_by(genusspecies, full_ec_number) %>%
  summarize(across(where(is.numeric), Mode))
# summarize_if(is.numeric, Mode)

paths_groupAandB_stats1 <- paths_groupAandB %>%
  group_by(genusspecies, full_ec_number) %>%
  summarize(across(where(is.character), Mode))
paths_groupAandB_stats <- inner_join(paths_groupAandB_stats1, paths_groupAandB_stats0)

rm(paths_groupAandB_stats0)
rm(paths_groupAandB_stats1)

## create new column with shorter description
paths_groupAandB_stats <- paths_groupAandB_stats %>% separate(product, into = c("productshort"), sep = "\\(EC", remove = FALSE, convert = TRUE, extra = "drop")

## creating new column with number and short description
paths_groupAandB_prestats <- paths_groupAandB_stats %>% unite("ec_numberanddescription", ec_number, ec_description, sep = "_", remove = FALSE)
rm(paths_groupAandB_stats)

#########################################################
### ADDING GENOME LENGTHS HERE BEFORE MAKING DOWNSTREAM FILES...
## ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary

## new tack - instead just pull sizes from our listAandBtojoin2

# genome_summaryurl <- paste0("ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary" )
# genome_summarydest <- paste0("genome_summary.tab" )
# 
# download.file(genome_summaryurl, genome_summarydest, mode="wb", quiet = TRUE, cacheOK = FALSE)
# 
# genome_summary0 <- read_tsv("genome_summary.tab", col_types = list(
#   genome_id = col_character(), genome_name = col_character(), taxon_id = col_character()))
# 
# genome_summary <- genome_summary0 %>% select(genome_id, genome_length, taxon_id)

## needs to not be an inner_join, but a left_join to be sure not missing values
## replacing genome_summary with listAandBtojoin2
paths_groupAandB_stats <- left_join(paths_groupAandB_prestats, listAandBtojoin2)

### MAKING A Mbp genome length column
paths_groupAandB_stats$genome_size_mbp <- paths_groupAandB_stats$genome_length / 1000000
paths_groupAandB_stats$genome_size_mbp <- round(paths_groupAandB_stats$genome_size_mbp, 1)

## CHECK HERE TO SEE IF THERE ARE MISSING GENOME LENGTHS..
## check to make sure this doesn't have NAs - if it does, need to examine genome_summary file to see if there are missing genomes...

if (nrow(paths_groupAandB_stats %>% dplyr::filter(genome_size_mbp == "NA")) == 0) {
  print("You have no missing genome values!")
} else {
  print("You have some missing genome values:")
  print(paths_groupAandB_stats %>% dplyr::filter(genome_size_mbp == "NA"))
}

Sys.sleep(2)

# rm(genome_summary)
# rm(genome_summary0)
#unlink("~/pomelo_outputs/genome_summary.tab")

#########################################################
## TESTS TO MAKE ADDITIONAL COLUMNS FOR GENE COUNTS...

paths_groupAandB_statsb1 <- paths_groupAandB %>%
  ungroup() %>%
  group_by(genusspecies, ec_number) %>%
  summarize(across(where(is.character), Mode)) %>% select(genusspecies, group, ec_number) # %>% add_tally() %>% rename(gene_countbygenusspecies = n)
paths_groupAandB_statsb1$seen <- 1
# summarize_if(is.character, Mode) becomes summarize(across(where(is.character), Mode))

## first change genusspecies order of statsb1 (moved) by group A alphabetical, then B alphabetical
paths_groupAandB_speciesorder <- paths_groupAandB_statsb1 %>% arrange(group,genusspecies)
orderingspecies <- paths_groupAandB_speciesorder[["genusspecies"]]
orderingspecies <- unique(orderingspecies)
paths_groupAandB_statsb1$genusspecies <- factor(paths_groupAandB_statsb1$genusspecies, levels = orderingspecies, ordered = TRUE)
## also make a separate file without group...
paths_groupAandB_statsa1 <- paths_groupAandB_statsb1 %>%
  select(-group)

unique(paths_groupAandB$genusspecies)
n_distinct(paths_groupAandB$genusspecies)
unique(paths_groupAandB_statsb1$genusspecies)
n_distinct(paths_groupAandB_statsb1$genusspecies)

## THIS IS THE PIVOT_WIDER COMMAND - MAKES A FILE TO GET THE ACTUAL STATS FROM...
paths_groupAandB_statswide <- paths_groupAandB_statsb1 %>%
  pivot_wider(names_from = genusspecies, values_from = seen, values_fill = 0, names_sort = TRUE)

paths_groupAandB_statswide2 <- paths_groupAandB_statswide %>%
  mutate(genecount_group_by_ecnumber = select(., 3:ncol(.)) %>% rowSums(na.rm = TRUE)) %>% select(group, ec_number, genecount_group_by_ecnumber)

paths_groupAandB_statswide2a <- paths_groupAandB_statswide2 %>% dplyr::filter(group == "group A")
paths_groupAandB_statswide2b <- paths_groupAandB_statswide2 %>% dplyr::filter(group == "group B")

## NOW MORE GENERALLY GETTING SPECIES COUNTS IN GROUP A AND GROUP B - FROM ABOVE...
paths_groupAandB_statswide2a$genepercentage_group_by_ecnumber <- (paths_groupAandB_statswide2a$genecount_group_by_ecnumber / paths_groupAandB_counts$genusspecies[1]) * 100
paths_groupAandB_statswide2b$genepercentage_group_by_ecnumber <- (paths_groupAandB_statswide2b$genecount_group_by_ecnumber / paths_groupAandB_counts$genusspecies[2]) * 100
paths_groupAandB_statswide2a$genepercentage_group_by_ecnumber <- round(paths_groupAandB_statswide2a$genepercentage_group_by_ecnumber, 0)
paths_groupAandB_statswide2b$genepercentage_group_by_ecnumber <- round(paths_groupAandB_statswide2b$genepercentage_group_by_ecnumber, 0)
paths_groupAandB_statswide2a <- paths_groupAandB_statswide2a %>% select(-group)
paths_groupAandB_statswide2b <- paths_groupAandB_statswide2b %>% select(-group)

## 2023 to get stats per species and per genus, need to repeat above steps but change from by group...or do it in section of line 575?
## what we need are in PML_bypathway to have similar columns but instead of by group, by each unique genusspecies, and also by genus

######################
### USING MAPPING FILE FOR SET OF REFERENCE PATHWAYS IN PATRIC

## ALLOW USER TO FIND THE FILE - these first steps are now added earlier
# paths_ref2b <- read_tsv(tk_choose.files(caption = "Find the 'mapping_GO_to_ecgene_and_ecpathway_toPATRIC.tab' file"), show_col_types = FALSE)

# using rstudioapi - request the path to an existing .csv file on disk
print("Select the 'mapping_BVBRC_allECs.tab' file, found in the pomelo-main/scripts folder")
Sys.sleep(1)
paths_ref2b <- read_tsv(rstudioapi::selectFile(caption = "Find the 'mapping_BVBRC_allECs.tab' file", label = "Select mapping_BVBRC_allECs.tab", path = data_dir, existing = TRUE, filter = "Tab Files (*.tab)"), show_col_types = FALSE)


## immediately remove .- as well..
paths_ref2b <- paths_ref2b[ grep(".-", paths_ref2b$ec_number, invert = TRUE) , ]

## removing some pathways that have been deleted from KEGG: https://www.genome.jp/kegg/docs/upd_map.html
# 1058, 471, 472, 473, 72, 231
## note there are still a few pathways >1000 - 1040,1051,1053,1055,1056,1057,1059
paths_ref2b <- paths_ref2b[ grep(".-", paths_ref2b$pathwayid, invert = TRUE) , ]
paths_ref2b <- paths_ref2b[ grep("p471", paths_ref2b$pathwayid, invert = TRUE) , ]
paths_ref2b <- paths_ref2b[ grep("p472", paths_ref2b$pathwayid, invert = TRUE) , ]
paths_ref2b <- paths_ref2b[ grep("p473", paths_ref2b$pathwayid, invert = TRUE) , ]
paths_ref2b <- paths_ref2b[ grep("Synthesis and degradation of ketone bodies", paths_ref2b$pathway_name, invert = TRUE) , ]
paths_ref2b <- paths_ref2b[ grep("p231", paths_ref2b$pathwayid, invert = TRUE) , ]
paths_ref2b <- paths_ref2b[ grep("p4070", paths_ref2b$pathwayid, invert = TRUE) , ]
paths_ref2b <- paths_ref2b[ grep("p4150", paths_ref2b$pathwayid, invert = TRUE) , ]


## now steps to pull back in any missing EC stats very early...check 1.1.1.133 & 5.4.3.5 !! also 1.6.99.5
paths_ref2b_a <- paths_ref2b %>% dplyr::filter(group == "group A")
paths_ref2b_b <- paths_ref2b %>% dplyr::filter(group == "group B")

paths_ref2b_a2 <- full_join(paths_ref2b_a, paths_groupAandB_statswide2a) %>% ungroup()
paths_ref2b_b2 <- full_join(paths_ref2b_b, paths_groupAandB_statswide2b) %>% ungroup()

## was just
#paths_ref2b2 <- bind_rows(paths_ref2b_a2, paths_ref2b_b2)

# fix to add ec_number and group, then but then also use a coalesce 
paths_ref2b2 <- bind_rows(paths_ref2b_a2, paths_ref2b_b2) %>%
  unite("ec_numberandgroup", ec_number, group, sep = "_", remove = FALSE) %>%
  group_by(ec_numberandgroup)
paths_groupAandB_stats <- paths_groupAandB_stats %>%
  unite("ec_numberandgroup", ec_number, group, sep = "_", remove = FALSE)

# ## corrected to only this, but removes too much...maybe not just go back to full join?? seems to help but not really
## instead of removing pathwayids and collapsing, we need to use a coalesce with this separate file...
paths_ref2b2_replacemissing <- bind_rows(paths_ref2b_a2, paths_ref2b_b2) %>%
  select(group,ec_number,genecount_group_by_ecnumber,genepercentage_group_by_ecnumber) %>%
  unite("ec_numberandgroup", ec_number, group, sep = "_", remove = TRUE) %>%
  group_by(ec_numberandgroup) %>%
  summarize(across(everything(), first)) %>%
  rename(genecount_group_by_ecnumber2 = genecount_group_by_ecnumber) %>%
  rename(genepercentage_group_by_ecnumber2 = genepercentage_group_by_ecnumber)

## split into A then B, then left join stats2wide 2a & 2b separately, then bind_rows
##stats2 had 27467, stats3 had 37524...
## NOW COMBINE PREVIOUS STATS with genecount_group_by_ecnumber AND genepercentage_group_by_ecnumber
#paths_groupAandB_stats2 <- left_join(paths_groupAandB_stats, paths_groupAandB_statswide2) %>% ungroup()
## now create stats2 & stats3...but stats3 first, then remove NA to get stats2...

## to avoid missed pathways adding NAs, drop pathways from paths_ref2b2
## was full join - keeping but adding next command!
paths_groupAandB_stats3 <- full_join(paths_groupAandB_stats, paths_ref2b2) %>% ungroup()

## tried to change to left join
# paths_groupAandB_stats3 <- left_join(paths_groupAandB_stats, paths_ref2b2) %>% ungroup() %>% 
#   select(-ec_numberandgroup)
## then full join after fix above - still removes NAs
# paths_groupAandB_stats3 <- full_join(paths_groupAandB_stats, paths_ref2b2) %>% ungroup() %>% 
#   select(-ec_numberandgroup)

## first a left join of new columns, then coalesce with genecount_group_by_ecnumber2 AND genepercentage_group_by_ecnumber2
paths_groupAandB_stats3t <- left_join(paths_groupAandB_stats3, paths_ref2b2_replacemissing) %>% 
  mutate(genecount_group_by_ecnumber = coalesce(genecount_group_by_ecnumber,genecount_group_by_ecnumber2)) %>% 
  mutate(genepercentage_group_by_ecnumber = coalesce(genepercentage_group_by_ecnumber,genepercentage_group_by_ecnumber2))
  
## note using stats3t just briefly.  ## update to use across instead of _at
paths_groupAandB_stats3 <- paths_groupAandB_stats3t %>%
  mutate(across(genusspecies, ~ replace_na(., "reference_set")))
#  mutate_at(vars(genusspecies), ~replace_na(., "reference_set"))

paths_groupAandB_stats2 <- paths_groupAandB_stats3 %>%
  dplyr::filter(genusspecies != "reference_set")

rm(paths_groupAandB_stats3t)
# write.table(paths_groupAandB_statswide, file = paste("paths_groupAandB_statswide",Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)
## counts of genes
n_distinct(paths_groupAandB_statswide2$ec_number)
n_distinct(paths_groupAandB_stats2$ec_number)
n_distinct(paths_ref2b$ec_number)

rm(paths_groupAandB_statswide2a)
rm(paths_groupAandB_statswide2b)
rm(paths_ref2b_a)
rm(paths_ref2b_b)
rm(paths_ref2b_a2)
rm(paths_ref2b_b2)

unique(paths_groupAandB_stats$genusspecies)
n_distinct(paths_groupAandB_stats$genusspecies)

unique(paths_groupAandB_stats3$genusspecies)
n_distinct(paths_groupAandB_stats3$genusspecies)

unique(paths_groupAandB_stats2$genusspecies)
n_distinct(paths_groupAandB_stats2$genusspecies)
## lost 3 species in Treponema run...

# want to re-order by genome size - maybe import the genome data then sort?
## then will need to order with a factor like so:

## removing this section for generalized code just by group A & group B
### INSTEAD ADDING GENOME LENGTHS EARLIER BEFORE MAKING DOWNSTREAM FILES...

### THIS WILL AUTOMATICALLY ORDER ALL SUBSEQUENT PLOTS WITH GENUSSPECIES TO BE IN ORDER OF GROUP, THEN GENOME SIZE
## update to use tree above
paths_groupAandB_stats3 <- paths_groupAandB_stats3 %>% arrange(group,genome_length)
ordering <- paths_groupAandB_stats3[["genusspecies"]]
ordering <- unique(ordering)
paths_groupAandB_stats3$genusspecies <- factor(paths_groupAandB_stats3$genusspecies, levels = ordering, ordered = TRUE)

##############################################################################
## end reordering
##############################################################################

###################################################################################################################################################################################
################################################# CODE FOR HIGH PRIORITY PATHWAYS  ################################################################################################
###################################################################################################################################################################################

### Formatting ###
# Change all NA's to 0 & also removing EC numbers with '-' - moving last part to very early now...
paths_groupAandB_stats3nonas <- paths_groupAandB_stats3
paths_groupAandB_stats3nonas <- paths_groupAandB_stats3nonas[ grep(".-", paths_groupAandB_stats3nonas$ec_number, invert = TRUE) , ]
paths_groupAandB_stats3nonas <- replace_na(paths_groupAandB_stats3nonas, list(genepercentage_group_by_ecnumber = 0))
## THIS ONE IS OKAY, BUT CHANGE STATS2_NONAS - fixing early now

paths_groupAandB_stats2_nonas <- paths_groupAandB_stats2
paths_groupAandB_stats2_nonas <- replace_na(paths_groupAandB_stats2_nonas, list(genepercentage_group_by_ecnumber = 0))
## INSTEAD OF REPLACING ALL NAs WITH 0, USE COALESCE (moving below to do this...) - fixing early now


##### Pathway Ranking Code #####
# List all pathways to rank #
paths_to_rank <- paths_groupAandB_stats3nonas %>%
  ungroup() %>%
  distinct(pathwayid)

# Calculate total number of genes in each pathway with and without ref set #
## number with reference set
total_genes_list <- paths_groupAandB_stats3nonas %>%
  ungroup() %>%
  group_by(pathwayid) %>%
  summarize(total_genes_inpathway = n_distinct(ec_number))
table_with_total <- dplyr::left_join(paths_groupAandB_stats3nonas, total_genes_list, by="pathwayid")

## number without reference set - MOVE FURTHER DOWN okay fixed above
## but still minor change - take larger set stats3 but remove only genecounts are NAs....
## in treponema seeing NAs in genecount_group_by_ecnumber that shouldn't be...
total_genes_list2 <- paths_groupAandB_stats3nonas %>%
  dplyr::filter(genecount_group_by_ecnumber != "NA") %>%
  ungroup() %>%
  group_by(pathwayid) %>%
  summarize(bvbrc_genes_inpathway = n_distinct(ec_number))
table_with_total_2 <- dplyr::left_join(table_with_total, total_genes_list2, by="pathwayid")

# Calculate promiscuity index for each gene (EC number) & avg promiscuity index for each pathway #
prom_index <- paths_groupAandB_stats3nonas %>%
  ungroup() %>%
  group_by(ec_number) %>%
  summarize(ec_pindex = n_distinct(pathwayid))
table_withprom1 <- dplyr::left_join(table_with_total_2, prom_index, by="ec_number")

avg_prom_index <- table_withprom1 %>%
  ungroup()%>%
  select(pathwayid, pathway_name, ec_number, total_genes_inpathway, bvbrc_genes_inpathway, ec_pindex) %>% 
  distinct()%>%
  group_by(pathwayid) %>%
  summarize(avg_p_index_bypathway = mean(ec_pindex, na.rm=TRUE))
table_withprom2 <- dplyr::left_join(table_withprom1, avg_prom_index, by="pathwayid")

ecs_in_more_than_one <- table_withprom2 %>%
  ungroup %>%
  select(pathwayid, pathway_name, ec_number, total_genes_inpathway, bvbrc_genes_inpathway, ec_pindex, avg_p_index_bypathway) %>%
  distinct()%>%
  group_by(pathwayid) %>%
  summarize(ecs_in_greater_than_one_pathway= sum(ec_pindex>1))
table_withprom3 <- dplyr::left_join(table_withprom2, ecs_in_more_than_one, by="pathwayid")

table_withprom4 <- table_withprom3 %>%
  mutate(perc_ecs_in_greater_than_one_pathway = ((ecs_in_greater_than_one_pathway/total_genes_inpathway)*100))

# Calculate number of genes in each pathway found in A/B #
## MOVE bvbrc_genes_inpathway CALCULATION HERE - UPDATE, FIXED ABOVE!!!
## ALSO BETTER WAY TO GRAB THESE IS TO USE COALESCE TO REPLACE genecount_group_by_ecnumber & genepercentage_group_by_ecnumber VALUES WHEN NAS...
## first have to split into A & B first to correctly coalesce the right value...

## ADDING EXTRA dplyr::filter TO REMOVE REFERENCE_SET - change instead remove only when genecounts are NAs....
genes_in_A <- paths_groupAandB_stats3nonas %>%
  ungroup() %>%
  dplyr::filter(group == "group A") %>%
  # dplyr::filter(genusspecies != "reference_set") %>%
  dplyr::filter(genecount_group_by_ecnumber != "NA") %>%
  group_by(pathwayid) %>%
  summarize(genes_found_A_bypathway = n_distinct(ec_number))
table_withA <- dplyr::left_join(table_withprom4, genes_in_A, by="pathwayid")

genes_in_B <- paths_groupAandB_stats3nonas %>%
  ungroup() %>%
  dplyr::filter(group == "group B") %>%
  # dplyr::filter(genusspecies != "reference_set") %>%
  dplyr::filter(genecount_group_by_ecnumber != "NA") %>%
  group_by(pathwayid) %>%
  summarize(genes_found_B_bypathway = n_distinct(ec_number))
table_withB <- dplyr::left_join(table_withA, genes_in_B, by="pathwayid")

## 2023 tests for counts by speciesname & bygenus (alternately - use genusspecies2?)
test_genes_by_species <- paths_groupAandB_stats3nonas %>% unite("genusspeciespathway", genusspecies, pathwayid, sep = "_", remove = FALSE) %>% relocate(genusspeciespathway, .after = genusspecies)
test_genes_by_genus <- paths_groupAandB_stats3nonas %>% unite("genuspathway", genus, pathwayid, sep = "_", remove = FALSE) %>% relocate(genuspathway, .after = genus)

test_genes_by_species <- test_genes_by_species %>%
  ungroup() %>%
  group_by(genusspeciespathway) %>%
  summarize(genes_found_byspecies_bypathway = n_distinct(ec_number))

test_genes_by_genus <- test_genes_by_genus %>%
  ungroup() %>%
  group_by(genuspathway) %>%
  summarize(genes_found_bygenus_bypathway = n_distinct(ec_number))

# create new column and use that to join...do earlier
# table_withB <- table_withB %>% unite("genusspeciespathway", genusspecies, pathwayid, sep = "_", remove = FALSE) %>% relocate(genusspeciespathway, .after = genusspecies)
# table_withB <- table_withB %>% unite("genuspathway", genus, pathwayid, sep = "_", remove = FALSE) %>% relocate(genuspathway, .after = genus)

# test_genes_by_species$genusspecies <- NULL
# test_genes_by_genus$genus <- NULL
# test_genes_by_species$pathwayid <- NULL
# test_genes_by_genus$pathwayid <- NULL

## can't join these because there are too many combinations. Instead add these at the end...plan B adding new combined column to join

## adding this as well
table_withgenes <- dplyr::left_join(table_withB, test_genes_by_species)
table_withgenes2 <- dplyr::left_join(table_withgenes, test_genes_by_genus)

table_withgenes3 <- replace_na(table_withgenes2, list(genes_found_A_bypathway = 0))
table_withgenes4 <- replace_na(table_withgenes3, list(genes_found_B_bypathway = 0))
table_withgenes5 <- replace_na(table_withgenes4, list(genes_found_byspecies_bypathway = 0))
table_withgenes6 <- replace_na(table_withgenes5, list(genes_found_bygenus_bypathway = 0))

# Calculate percent of genes missing in each pathway in A/B #
table_with_perc0 <- table_withgenes6 %>%
  mutate(perc_missing_in_A_bypathway = (100-((genes_found_A_bypathway/bvbrc_genes_inpathway)*100))) %>%
  mutate(perc_missing_in_B_bypathway = (100-((genes_found_B_bypathway/bvbrc_genes_inpathway)*100)))

## now same calculation for genes_found_bygenus_bypathway & genes_found_byspecies_bypathway
table_with_perc <- table_with_perc0 %>%
  mutate(perc_in_species_bypathway = (genes_found_byspecies_bypathway/bvbrc_genes_inpathway)*100) %>%
  mutate(perc_in_genus_bypathway = (genes_found_bygenus_bypathway/bvbrc_genes_inpathway)*100)


table_with_perc_na <- replace_na(table_with_perc, list(perc_missing_in_A_bypathway = 100))
table_with_perc_na2 <- replace_na(table_with_perc_na, list(perc_missing_in_B_bypathway = 100))
table_with_perc_na3 <- replace_na(table_with_perc_na2, list(perc_in_species_bypathway = 0))
table_with_perc_na4 <- replace_na(table_with_perc_na3, list(perc_in_genus_bypathway = 0))

# Calculate Percentage Differential between groups A/B #
table_with_perc_2 <- table_with_perc_na4 %>%
  mutate(perc_differential_A_bypathway = ((perc_missing_in_A_bypathway)-(perc_missing_in_B_bypathway))) %>%
  mutate(perc_differential_B_bypathway = ((perc_missing_in_B_bypathway)-(perc_missing_in_A_bypathway))) %>%
  dplyr::filter(genusspecies != "reference_set")
## INSTEAD  remove only when genecounts are NAs?? but see above, should be okay - 
  # dplyr::filter(genusspecies != "reference_set")
  # dplyr::filter(genecount_group_by_ecnumber != "NA") maybe changing to genepercentage?? NO
  # dplyr::filter(genepercentage_group_by_ecnumber != "NA")


# Determine Differences at gene level using gene_percentaage_group # 
test_A <- table_with_perc_2 %>% dplyr::filter(group == "group A")
test_B <- table_with_perc_2 %>% dplyr::filter(group == "group B")

#total_genomes <- paths_groupAandB %>% ungroup() %>% group_by(group) %>% summarise_at(vars(genusspecies), n_distinct) 
total_genomes <- paths_groupAandB %>% ungroup() %>% group_by(group) %>% summarise(across(genusspecies, n_distinct))
total_genomes_A <- as.integer(total_genomes[1,2])
total_genomes_B <- as.integer(total_genomes[2,2])

## changing from summarize to mutate which keeps other columns (but then have dplyr::filter) df %>% mutate(z = x + y, .keep = "all")
##  also pull genepercentage_group_by_ecnumber & use this updated value instead of gene_percentage_
## another update, we do not need to keep gene_count at all
combine_A <- test_A %>%
  ungroup() %>%
  group_by(ec_number) %>%
  # summarize(gene_count_A = n_distinct(genusspecies)) %>%
  # mutate(gene_count_A = n_distinct(genusspecies)) %>%
  # mutate(gene_percentage_A = (gene_count_A/total_genomes_A)) %>%
  # select(ec_number, gene_count_A, gene_percentage_A) %>%
  select(ec_number, genepercentage_group_by_ecnumber) %>%
  rename(gene_percentage_A = genepercentage_group_by_ecnumber) %>%
  distinct()

combine_B <- test_B %>%
  ungroup() %>%
  group_by(ec_number) %>%
  # summarize(gene_count_B = n_distinct(genusspecies)) %>%
  # mutate(gene_count_B = n_distinct(genusspecies)) %>%
  # mutate(gene_percentage_B=(gene_count_B/total_genomes_B)) %>%
  # select(ec_number, gene_count_B, gene_percentage_B) %>%
  select(ec_number, genepercentage_group_by_ecnumber) %>%
  rename(gene_percentage_B = genepercentage_group_by_ecnumber) %>%
  distinct()

## change here, also pull genepercentage_group_by_ecnumber & use this updated value instead of gene_percentage_
table_with_diff <- dplyr::left_join(table_with_perc_2, combine_A, by="ec_number")
table_with_diff_2 <- dplyr::left_join(table_with_diff, combine_B, by="ec_number")

table_with_diff_3 <- table_with_diff_2 %>%
  mutate(differential_bygene=0)

table_with_diff_4 <- replace_na(table_with_diff_3, list(gene_percentage_A = 0))
table_with_diff_5 <- replace_na(table_with_diff_4, list(gene_percentage_B = 0))

num_rows <- nrow(table_with_diff_5)

# calculating raw differential_bygene & keeping separate for _A & _B..., also round to 1 decimal point
## SOME GENES ARE BEING MISSED...fix is above in combine_A and _B

# table_with_diff_5A <- table_with_diff_5 %>% dplyr::filter(group == "group A")
# table_with_diff_5B <- table_with_diff_5 %>% dplyr::filter(group == "group B")
# table_with_diff_5A$genepercentage_groupA_by_ecnumber <- table_with_diff_5A$genepercentage_group_by_ecnumber
# table_with_diff_5B$genepercentage_groupB_by_ecnumber <- table_with_diff_5B$genepercentage_group_by_ecnumber
# table_with_diff_5 <- full_join(table_with_diff_5A,table_with_diff_5B)
# table_with_diff_5 <- replace_na(table_with_diff_5, list(genepercentage_groupA_by_ecnumber = 0))
# table_with_diff_5 <- replace_na(table_with_diff_5, list(genepercentage_groupB_by_ecnumber = 0))
# table_with_diff_5$differential_bygene_A <- (table_with_diff_5$genepercentage_groupA_by_ecnumber - table_with_diff_5$genepercentage_groupB_by_ecnumber)
# table_with_diff_5$differential_bygene_B <- (table_with_diff_5$genepercentage_groupB_by_ecnumber - table_with_diff_5$genepercentage_groupA_by_ecnumber) ## like perc_differential will just be -ve of _A

table_with_diff_5$differential_bygene_A <- (table_with_diff_5$gene_percentage_A - table_with_diff_5$gene_percentage_B)
table_with_diff_5$differential_bygene_B <- (table_with_diff_5$gene_percentage_B - table_with_diff_5$gene_percentage_A) ## like perc_differential will just be -ve of _A

# table_with_diff_5$differential_bygene_A <- (table_with_diff_5$gene_percentage_A - table_with_diff_5$gene_percentage_B) * 100
# table_with_diff_5$differential_bygene_B <- (table_with_diff_5$gene_percentage_B - table_with_diff_5$gene_percentage_A) * 100 ## like perc_differential will just be -ve of _A
# table_with_diff_5$differential_bygene_A <- round(table_with_diff_5$differential_bygene_A, 1)
# table_with_diff_5$differential_bygene_B <- round(table_with_diff_5$differential_bygene_B, 1)

## convert negative values to zero
table_with_diff_5 <- table_with_diff_5 %>% mutate(differential_bygene_A = if_else(differential_bygene_A < 0, 0, differential_bygene_A)) %>% 
  mutate(differential_bygene_B = if_else(differential_bygene_B < 0, 0, differential_bygene_B))
## will keep separate for downstream PML_A & PML_B score, but can also combine to a single value
table_with_diff_5 <- table_with_diff_5  %>% mutate(differential_bygene=pmax(differential_bygene_A,differential_bygene_B))

table_with_e <- table_with_diff_5 %>%
  ungroup() %>%
  group_by(pathwayid)

## DO WE NEED TO PULL OUT ONE GROUP HERE?? or a coalesce with 5A & 5B...
# Grab factors from larger table to print to new data frame #
ranking_table_1 <- dplyr::right_join(table_with_e, paths_to_rank, by="pathwayid") %>%
  ungroup() %>%
  distinct(group,pathwayid,pathway_id,pathway_name,ec_number,total_genes_inpathway,bvbrc_genes_inpathway, genes_found_A_bypathway,genes_found_B_bypathway, 
           perc_missing_in_A_bypathway,perc_missing_in_B_bypathway,avg_p_index_bypathway, differential_bygene_A, differential_bygene_B, differential_bygene, perc_differential_A_bypathway, 
           perc_differential_B_bypathway, gene_percentage_A, gene_percentage_B)

# Calculate number of differential genes for each pathway #
## updating to divide into _A and _B
## need to not count genes twice!! - note new grouping & first, then just a mutate
diff_counts <- ranking_table_1 %>%
  # dplyr::filter(group == "group A") %>% 
  group_by(pathwayid,ec_number) %>%
  summarize(across(everything(), first)) %>%
  # summarize_all(first) %>%
  ungroup() %>%
  group_by(pathwayid) %>% 
  mutate(total_differences_A_bygene = sum(differential_bygene_A)) %>% 
  mutate(total_differences_B_bygene = sum(differential_bygene_B))

# diff_countsB <- ranking_table_1 %>%
#   # dplyr::filter(group == "group A") %>% 
#   group_by(pathwayid,ec_number) %>%
#   summarize_all(first) %>%
#   ungroup() %>%
#   group_by(pathwayid) %>% 
#   mutate(total_differences_B_bygene = sum(differential_bygene_B))
# diff_counts <- full_join(diff_countsA,diff_countsB)

# ranking_table_2 <- dplyr::left_join(ranking_table_1, diff_counts, by="pathwayid")
ranking_table_2 <- dplyr::left_join(ranking_table_1, diff_counts)

# Make smaller table for PML calculation #
## need to grab max values here

ranking_table_3x <- ranking_table_2 %>%
  select(pathwayid,pathway_id,pathway_name,total_genes_inpathway,bvbrc_genes_inpathway, genes_found_A_bypathway, genes_found_B_bypathway, perc_missing_in_A_bypathway, 
         perc_missing_in_B_bypathway, avg_p_index_bypathway, total_differences_A_bygene, total_differences_B_bygene) %>%
  distinct() %>% 
  dplyr::filter(genes_found_B_bypathway != "NA" & total_differences_B_bygene != "NA")
  # group_by(pathwayid) %>%
  # summarize_all(max)
  # dplyr::filter(genecount_group_by_ecnumber != "NA") !is.na(mean(x$v_identity))

ranking_table_3y <- ranking_table_2 %>%
  select(pathwayid,pathway_id,pathway_name,total_genes_inpathway,bvbrc_genes_inpathway, genes_found_A_bypathway, genes_found_B_bypathway, perc_missing_in_A_bypathway, 
         perc_missing_in_B_bypathway, avg_p_index_bypathway, total_differences_A_bygene, total_differences_B_bygene) %>%
  distinct() %>% 
  dplyr::filter(is.na(pathway_name))
  # dplyr::filter(pathway_name == "NA")
ranking_table_3 <- rbind(ranking_table_3x,ranking_table_3y)

##### Pathway Metabolic Loss (PML) score calculation (was Metabolic Loss Score MLS) #####

# Turn all NAs to 0 #
ranking_table_4 <- replace_na(ranking_table_3, list(total_genes_inpathway=0, bvbrc_genes_inpathway=0, genes_found_A_bypathway=0, genes_found_B_bypathway=0,
                                                    perc_missing_in_A_bypathway=0, perc_missing_in_B_bypathway=0, perc_in_genus_bypathway=0, perc_in_genus_bypathway=0, avg_p_index_bypathway=0, total_differences_A_bygene=0, total_differences_B_bygene=0))

# Assign points based on number of differential genes #
ranking_table_PML <- ranking_table_4 %>%
  mutate(PML_score_A=0) %>%
  mutate(PML_score_Bo=0) %>%
  arrange(total_differences_A_bygene)

num_rows <- nrow(ranking_table_PML)

## update
ranking_table_PML$PML_score_A <- (ranking_table_PML$total_differences_A_bygene) * 10
ranking_table_PML$PML_score_Bo <- (ranking_table_PML$total_differences_B_bygene) * 10

# Assign points based on number percentage differential #
## first calculating new measures
ranking_table_PML_2 <- ranking_table_PML
# ranking_table_PML_2$PML_score <- NULL

ranking_table_PML_2$genes_absent_A_bypathway <- ranking_table_PML_2$bvbrc_genes_inpathway - ranking_table_PML_2$genes_found_A_bypathway
ranking_table_PML_2$genes_absent_B_bypathway <- ranking_table_PML_2$bvbrc_genes_inpathway - ranking_table_PML_2$genes_found_B_bypathway
ranking_table_PML_2$genes_absent_differential <- ranking_table_PML_2$genes_absent_A_bypathway - ranking_table_PML_2$genes_absent_B_bypathway

# ranking_table_PML_2$PML_score_Bp <- ranking_table_PML_2$PML_score_Bo
# ranking_table_PML_2$PML_score_B <- NULL
ranking_table_PML_2$PML_score_B <- ranking_table_PML_2$PML_score_Bo + (ranking_table_PML_2$genes_absent_differential) * 50

# Assign points based on total number of genes #
ranking_table_PML_4 <- ranking_table_PML_2 %>%
  arrange(bvbrc_genes_inpathway)

ranking_table_PML_4$PML_score_A <- ranking_table_PML_4$PML_score_A + (ranking_table_PML_4$bvbrc_genes_inpathway * 5)
ranking_table_PML_4$PML_score_B <- ranking_table_PML_4$PML_score_B + (ranking_table_PML_4$bvbrc_genes_inpathway * 5) 

# remove negatives from PML_score_B
ranking_table_PML_5 <- ranking_table_PML_4 %>%
  mutate(PML_score_B = ifelse(PML_score_B < 0, 0, PML_score_B))

# drop PML_score_Bo,PML_score_Bp from final tables & change names - PML_score_B to PML& PML_score_A to loss_innontarget_group_score
ranking_table_PML_6 <- ranking_table_PML_5 %>% select (-PML_score_Bo)

## RENAME scores here, rest later...
## update - PML_score_B is now PML_composite, PML_score_A is loss_innontarget_group_composite
## new PML & loss_innontarget_group_score now defined below
ranking_table_PML_6 <- ranking_table_PML_6 %>%
  # rename(PML = PML_score_B) %>%
  # rename(loss_innontarget_group_score = PML_score_A) %>%
  # relocate(loss_innontarget_group_score, .before = PML)
  rename(PML_composite = PML_score_B) %>%
  rename(loss_innontarget_group_composite = PML_score_A) %>%
  relocate(loss_innontarget_group_composite, .before = PML_composite)

################################################################################################
################################################################################################

# Add PML score to larger data frame # but first have to remove pathway_id now...
ranking_table_PML_6b <- ranking_table_PML_6 %>% select (-pathway_id)

paths_groupAandB_stats3nonas_PML <- dplyr::left_join(paths_groupAandB_stats3nonas, ranking_table_PML_6b, by="pathwayid")

## change columns here!! check for any downstream use of stats2 THERE ARE SOME - ALSO WATCH FOR ANY JOINS WITH UNCHANGED FILES...
## renamed
paths_groupAandB_stats3nonas_PML <- paths_groupAandB_stats3nonas_PML %>% rename(genes_found_target_group_bypathway = genes_found_A_bypathway) %>%
  rename(genes_found_nontarget_group_bypathway = genes_found_B_bypathway) %>%
  rename(perc_missing_in_target_group_bypathway = perc_missing_in_A_bypathway) %>%
  rename(perc_missing_in_nontarget_group_bypathway = perc_missing_in_B_bypathway) %>%
  # rename(sum_differences_bygene_target_group_perpathway = total_differences_A_bygene) %>% ## keep this & next column downstream
  # rename(sum_differences_bygene_nontarget_group_perpathway = total_differences_B_bygene) %>%
  rename(loss_innontarget_group_score = total_differences_A_bygene) %>% ## new PML scores are just these
  rename(PML = total_differences_B_bygene) %>%
  rename(genes_absent_target_group_bypathway = genes_absent_A_bypathway) %>%
  rename(genes_absent_nontarget_group_bypathway = genes_absent_B_bypathway)

## and row change...
paths_groupAandB_stats3nonas_PML$group <- gsub('group A','target',paths_groupAandB_stats3nonas_PML$group)
paths_groupAandB_stats3nonas_PML$group <- gsub('group B','non-target',paths_groupAandB_stats3nonas_PML$group)
# 
# x$group <- gsub('group A','target',x$group)
# x$group <- gsub('group B','non-target',x$group)

# name changes
# genes_found_A_bypathway
# genes_found_B_bypathway
# perc_missing_in_A_bypathway
# perc_missing_in_B_bypathway
# perc_differential_A_bypathway
# perc_differential_B_bypathway
# total_differences_A_bygene
# total_differences_B_bygene
# genes_absent_A_bypathway
# genes_absent_B_bypathway
# ## to
# genes_found_target_group_bypathway
# genes_found_nontarget_group_bypathway
# perc_missing_in_target_group_bypathway
# perc_missing_in_nontarget_group_bypathway
# perc_differential_target_group_bypathway
# perc_differential_nontarget_group_bypathway
# total_differences_target_group_bygene           now sum_differences_bygene_target_group_perpathway 
# total_differences_nontarget_group_bygene        now sum_differences_bygene_nontarget_group_perpathway 
# genes_absent_target_group_bypathway
# genes_absent_nontarget_group_bypathway

## also change any select columns & relocate

## remove duplicate pathway name
paths_groupAandB_stats3nonas_PML$pathway_name.y <- NULL
paths_groupAandB_stats3nonas_PML$pathway_name <- paths_groupAandB_stats3nonas_PML$pathway_name.x
paths_groupAandB_stats3nonas_PML$pathway_name.x <- NULL

## round some values...
paths_groupAandB_stats3nonas_PML$avg_p_index_bypathway <- round(paths_groupAandB_stats3nonas_PML$avg_p_index_bypathway, 2)
paths_groupAandB_stats3nonas_PML$perc_missing_in_nontarget_group_bypathway <- round(paths_groupAandB_stats3nonas_PML$perc_missing_in_nontarget_group_bypathway, 1)
paths_groupAandB_stats3nonas_PML$perc_missing_in_target_group_bypathway <- round(paths_groupAandB_stats3nonas_PML$perc_missing_in_target_group_bypathway, 1)

# paths_groupAandB_stats3nonas_PML$perc_in_species_bypathway <- round(paths_groupAandB_stats3nonas_PML$perc_in_species_bypathway, 1)
# paths_groupAandB_stats3nonas_PML$perc_in_genus_bypathway <- round(paths_groupAandB_stats3nonas_PML$perc_in_genus_bypathway, 1)
## 2023 only now bring in perc_in_species_bypathway & perc_in_genus_bypathway from table_with_e
pull_perc_in_species <- table_with_e %>%
  ungroup() %>%
  select(genusspeciespathway,perc_in_species_bypathway) %>% 
  distinct()
pull_perc_in_genus <- table_with_e %>%
  ungroup() %>%
  select(genuspathway,perc_in_genus_bypathway) %>% 
  distinct()
## now join just these two columns to paths_groupAandB_stats3nonas_PML
paths_groupAandB_stats3nonas_PML01 <- dplyr::left_join(paths_groupAandB_stats3nonas_PML, pull_perc_in_species)
paths_groupAandB_stats3nonas_PML02 <- dplyr::left_join(paths_groupAandB_stats3nonas_PML01, pull_perc_in_genus)
paths_groupAandB_stats3nonas_PML <- paths_groupAandB_stats3nonas_PML02
Sys.sleep(1)
rm(paths_groupAandB_stats3nonas_PML02)
rm(paths_groupAandB_stats3nonas_PML01)

## examine final scores by pathway
PML_bypathway <- paths_groupAandB_stats3nonas_PML %>%
  ungroup() %>%
  group_by(pathwayid) %>%
  # summarize_all(first) %>%
  summarize(across(everything(), first)) %>%
  arrange(desc(PML),desc(loss_innontarget_group_score), pathway_id) %>%
  ## adding simpler PML scores, keeping composite values here but add to very end, PML right after pathway...also removing the genes_absent_target_group_bypathway & genes_absent_nontarget_group_bypathway columns...
  select(pathwayid,pathway_id,pathway_name,PML,total_genes_inpathway,bvbrc_genes_inpathway,genes_found_target_group_bypathway,genes_found_nontarget_group_bypathway,perc_missing_in_target_group_bypathway,perc_missing_in_nontarget_group_bypathway,avg_p_index_bypathway,genes_absent_differential,loss_innontarget_group_score,PML_composite,loss_innontarget_group_composite)

## ADDING FACTOR OF LATEST PML_SCORE FILE WITH GENOME SIZE ORDER - this is changed, use desc
paths_groupAandB_stats3nonas_PMLorder <- paths_groupAandB_stats3nonas_PML %>% arrange(desc(group),genome_length)
ordering2 <- paths_groupAandB_stats3nonas_PMLorder[["genusspecies"]]
ordering2 <- unique(ordering2)

paths_groupAandB_stats3nonas_PML$genusspecies <- factor(paths_groupAandB_stats3nonas_PML$genusspecies, levels = ordering2, ordered = TRUE)

##############################
## IF INSTEAD ORDERING BY PHYLO...see further below for separate dataframe with phylo ordering
# paths_groupAandB_stats3nonas_PML$genusspecies <- factor(paths_groupAandB_stats3nonas_PML$genusspecies, levels = ordering_byphylo, ordered = TRUE)

## then make new plots with x & y reversed from above...
##############################

## moving plot_title naming here...
unlist(strsplit(paths_groupAandB_stats3nonas_PMLorder$genome_name[5], " "))[1]
plot_title <- unlist(strsplit(paths_groupAandB_stats3nonas_PMLorder$genome_name[5], " "))[1]

print("We've guessed that your taxon plot title should be")
print(plot_title)
print("Please rename if this is a bad guess!")

write.table(PML_bypathway, file = paste("Summary_of_ranked_pathways",plot_title,Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(final_ranked_AnoNAs, file = paste("Ranked_pathways_A",plot_title,Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(final_ranked_BnoNAs, file = paste("Ranked_pathways_B",plot_title,Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################################

### why 170.....removes pathways with nearly all missing data (>80% missing in both A & B)
## adding second dplyr::filter to remove pathways with bvbrc_genes_inpathway = 0
paths_groupAandB_stats3nonasnofullymissing_PML <- paths_groupAandB_stats3nonas_PML %>% dplyr::filter((perc_missing_in_target_group_bypathway + perc_missing_in_nontarget_group_bypathway < 170)) %>% 
  dplyr::filter(bvbrc_genes_inpathway != "0")


## note difference when you look at full dataset vs one row per pathway...
median(paths_groupAandB_stats3nonas_PML$loss_innontarget_group_score)
median(paths_groupAandB_stats3nonas_PML$PML)
median(PML_bypathway$loss_innontarget_group_score)
median(PML_bypathway$PML)

median(paths_groupAandB_stats3nonasnofullymissing_PML$loss_innontarget_group_score)
median(paths_groupAandB_stats3nonasnofullymissing_PML$PML)

## add PML score to pathway name for plots below...
paths_groupAandB_stats3nonasnofullymissing_PML$xtext <- "PML score"
paths_groupAandB_stats3nonasnofullymissing_PML$xtext2 <- "loss in non-target group score"
paths_groupAandB_stats3nonasnofullymissing_PML <- paths_groupAandB_stats3nonasnofullymissing_PML %>% unite("pathway_name_with_loss_innontarget_group_score", pathway_name, xtext2, loss_innontarget_group_score, sep = " ", remove = FALSE)
paths_groupAandB_stats3nonasnofullymissing_PML <- paths_groupAandB_stats3nonasnofullymissing_PML %>% unite("pathway_name_with_PML_score", pathway_name, xtext, PML, sep = " ", remove = FALSE)
paths_groupAandB_stats3nonasnofullymissing_PML$xtext <- NULL
paths_groupAandB_stats3nonasnofullymissing_PML$xtext2 <- NULL

################
## late in calculations, combining perc_missing_target and _nontarget into a single column for by_pathway plots
paths_groupAandB_stats3nonasnofullymissing_PML0 <- paths_groupAandB_stats3nonasnofullymissing_PML %>%
  ungroup() %>% 
  group_by(group,pathwayid) %>% 
  # summarize_all(first) %>%
  summarize(across(everything(), first)) %>%
  select(group,pathwayid,perc_missing_in_target_group_bypathway,perc_missing_in_nontarget_group_bypathway)

paths_groupAandB_stats3nonasnofullymissing_PMLa <- paths_groupAandB_stats3nonasnofullymissing_PML0 %>% dplyr::filter(group == "target")
paths_groupAandB_stats3nonasnofullymissing_PMLb <- paths_groupAandB_stats3nonasnofullymissing_PML0 %>% dplyr::filter(group == "non-target")

paths_groupAandB_stats3nonasnofullymissing_PMLa <- paths_groupAandB_stats3nonasnofullymissing_PMLa %>%
  ungroup() %>% 
  rename(perc_missing_in_group_bypathway = perc_missing_in_target_group_bypathway) %>%
  select(pathwayid,perc_missing_in_group_bypathway)
paths_groupAandB_stats3nonasnofullymissing_PMLb <- paths_groupAandB_stats3nonasnofullymissing_PMLb %>%
  ungroup() %>% 
  rename(perc_missing_in_group_bypathway = perc_missing_in_nontarget_group_bypathway) %>%
  select(pathwayid,perc_missing_in_group_bypathway)

paths_groupAandB_stats3nonasnofullymissing_PMLaf <- paths_groupAandB_stats3nonasnofullymissing_PML %>% dplyr::filter(group == "target")
paths_groupAandB_stats3nonasnofullymissing_PMLbf <- paths_groupAandB_stats3nonasnofullymissing_PML %>% dplyr::filter(group == "non-target")

paths_groupAandB_stats3nonasnofullymissing_PMLa2 <- left_join(paths_groupAandB_stats3nonasnofullymissing_PMLaf, paths_groupAandB_stats3nonasnofullymissing_PMLa)
paths_groupAandB_stats3nonasnofullymissing_PMLb2 <- left_join(paths_groupAandB_stats3nonasnofullymissing_PMLbf, paths_groupAandB_stats3nonasnofullymissing_PMLb)

paths_groupAandB_stats3nonasnofullymissing_PMLx <- bind_rows(paths_groupAandB_stats3nonasnofullymissing_PMLa2, paths_groupAandB_stats3nonasnofullymissing_PMLb2) %>% 
  replace_na(list(perc_missing_in_group_bypathway = 100))
paths_groupAandB_stats3nonasnofullymissing_PMLx$perc_in_group_bypathway <- (100 - paths_groupAandB_stats3nonasnofullymissing_PMLx$perc_missing_in_group_bypathway)

## now replace PML with PMLx, then remove PMLx
paths_groupAandB_stats3nonasnofullymissing_PML <- paths_groupAandB_stats3nonasnofullymissing_PMLx
Sys.sleep(1)
rm(paths_groupAandB_stats3nonasnofullymissing_PMLa)
rm(paths_groupAandB_stats3nonasnofullymissing_PMLb)
rm(paths_groupAandB_stats3nonasnofullymissing_PML0)
rm(paths_groupAandB_stats3nonasnofullymissing_PMLx)
rm(paths_groupAandB_stats3nonasnofullymissing_PMLaf)
rm(paths_groupAandB_stats3nonasnofullymissing_PMLbf)
rm(paths_groupAandB_stats3nonasnofullymissing_PMLa2)
rm(paths_groupAandB_stats3nonasnofullymissing_PMLb2)
paths_groupAandB_stats3nonasnofullymissing_PML$perc_missing_in_group_bypathway <- NULL
######################## end combining perc_missing_target and _nontarget into a single colum

## focus plots on pathways with PML score > median - CORRECTED
paths_enrichedinA <- paths_groupAandB_stats3nonasnofullymissing_PML %>% 
  ungroup() %>% 
  dplyr::filter(loss_innontarget_group_score > median(loss_innontarget_group_score, na.rm = TRUE))

median(paths_groupAandB_stats3nonasnofullymissing_PML$PML)
median(paths_groupAandB_stats3nonasnofullymissing_PML$loss_innontarget_group_score)

paths_enrichedinB <- paths_groupAandB_stats3nonasnofullymissing_PML %>% 
  ungroup() %>% 
  dplyr::filter(PML > median(PML, na.rm = TRUE))

## do we need these??
# paths_enrichedinA_stats <- paths_enrichedinA %>%
#   ungroup() %>%
#   distinct(pathwayid, .keep_all = TRUE) %>%
#   select(pathwayid,pathway_name,total_genes_inpathway,bvbrc_genes_inpathway,perc_missing_in_A_bypathway,perc_missing_in_B_bypathway,perc_in_species_bypathway,perc_in_genus_bypathway,genes_absent_A_bypathway,genes_absent_B_bypathway,genes_absent_differential,avg_p_index_bypathway,total_differences_A_bygene,total_differences_B_bygene,loss_innontarget_group_score,PML,pathway_name_with_loss_innontarget_group_score)
# 
# paths_enrichedinB_stats <- paths_enrichedinB %>%
#   ungroup() %>%
#   distinct(pathwayid, .keep_all = TRUE) %>%
#   select(pathwayid,pathway_name,total_genes_inpathway,bvbrc_genes_inpathway,perc_missing_in_A_bypathway,perc_missing_in_B_bypathway,perc_in_species_bypathway,perc_in_genus_bypathway,genes_absent_A_bypathway,genes_absent_B_bypathway,genes_absent_differential,avg_p_index_bypathway,total_differences_A_bygene,total_differences_B_bygene,loss_innontarget_group_score,PML,pathway_name_with_PML)
# write.table(paths_enrichedinA_stats, file = paste("pathways_enriched_in_groupA_",plot_title,Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(paths_enrichedinB_stats, file = paste("pathways_enriched_in_groupB_",plot_title,Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)

## extracting the ec pathways...
paths_enrichedinA_unique <- as.data.frame(unique(paths_enrichedinA$pathwayid))
colnames(paths_enrichedinA_unique)[1] <- "pathwayid"

paths_enrichedinB_unique <- as.data.frame(unique(paths_enrichedinB$pathwayid))
colnames(paths_enrichedinB_unique)[1] <- "pathwayid"

write.table(paths_groupAandB_stats3nonasnofullymissing_PML, file = paste("pathway_stats_",plot_title,Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)


## - CHECK IN BY SPECIES PLOTS - LOOK CURRENTLY BY ALPHABET...WANT BY GROUP, THEN BY GEN SIZE OR PHYLO

## ADD SUBSET DATA FOR GEOM_TEXT data=subset(mtcars, wt > 4 | mpg > 25)  data = subset(paths_enrichedinA, genepercentage_group_by_ecnumber > 0.2) or after aes portion: , subset = .(wt > 4 | mpg > 25)
## test with paths_rickettsiales_stats.d900 <- paths_groupAandB_stats4 %>% dplyr::filter(pathway_id == "900")
## also adding ordering facets by PML score instead of alphabetical
paths_enrichedinA <- paths_enrichedinA %>% arrange(desc(loss_innontarget_group_score), pathway_id)
ordering_enrichedinA <- paths_enrichedinA[["pathway_name_with_loss_innontarget_group_score"]]
ordering_enrichedinA <- unique(ordering_enrichedinA)
paths_enrichedinA$pathway_name_with_loss_innontarget_group_score <- factor(paths_enrichedinA$pathway_name_with_loss_innontarget_group_score, levels = ordering_enrichedinA, ordered = TRUE)


paths_enrichedinB <- paths_enrichedinB %>% arrange(desc(PML), desc(loss_innontarget_group_score), pathway_id)
ordering_enrichedinB <- paths_enrichedinB[["pathway_name_with_PML_score"]]
ordering_enrichedinB <- unique(ordering_enrichedinB)
paths_enrichedinB$pathway_name_with_PML_score <- factor(paths_enrichedinB$pathway_name_with_PML_score, levels = ordering_enrichedinB, ordered = TRUE)


### if we want to have smaller plots with pathways only found in reference set, remove fullec is NA...
### removing rows with EC number is missing...in other words, removing data only in the reference set!!
## only use in phylo commands further below...
paths_groupAandB_stats3nonasnofullymissing_PMLnoref <- paths_groupAandB_stats3nonasnofullymissing_PML %>%
  dplyr::filter(full_ec_number != "NA")
paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genusspecies <- factor(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genusspecies, levels = ordering2, ordered = TRUE)

paths_enrichedinAnonas <- paths_enrichedinA %>%
  dplyr::filter(full_ec_number != "NA")
paths_enrichedinBnonas <- paths_enrichedinB %>%
  dplyr::filter(full_ec_number != "NA")

##### to export table of full dataset
## making table of wide stats by pathway...then EXPORT TABLES...  https://stackoverflow.com/questions/25378184/extract-data-from-a-ggplot
paths_ref2b_colstoadd <- paths_groupAandB_stats3nonasnofullymissing_PMLnoref %>%
  ungroup() %>% 
  group_by(ec_number) %>% 
  # summarize_all(first) %>%
  summarize(across(everything(), first)) %>%
  select(ec_number,ec_description)

## now join paths_ref2b to get ec_description - but first change group names
paths_ref2b_renamed <- paths_ref2b
paths_ref2b_renamed$group <- gsub('group A','target',paths_ref2b_renamed$group)
paths_ref2b_renamed$group <- gsub('group B','non-target',paths_ref2b_renamed$group)

paths_ref2c <- left_join(paths_ref2b_renamed, paths_ref2b_colstoadd) %>%
  dplyr::filter(ec_description != "NA") %>% 
  dplyr::filter(group == "target") %>% 
  select(-group)

## save mapping stats
mapping_stats0 <- left_join(paths_ref2c, prom_index)
mapping_stats <- left_join(mapping_stats0, avg_prom_index) %>%
  arrange(pathway_id, ec_number) %>%
  relocate(ec_pindex, .after = ec_number) %>%
  relocate(ec_description, .after = ec_number) %>%
  relocate(avg_p_index_bypathway, .after = pathway_name)
mapping_stats$avg_p_index_bypathway <- round(mapping_stats$avg_p_index_bypathway, 2)

write.table(mapping_stats, file = paste("ec_number_pathway_stats",plot_title,Sys.Date(),".tab"), sep = "\t", row.names = FALSE, quote = FALSE)

## now add genecount_group_by_ecnumber & genepercentage_group_by_ecnumber from stats2
# unique(paths_groupAandB_stats3$ec_number)
# unique(paths_groupAandB_stats2$ec_number)
# n_distinct(paths_ref2c$ec_number)

## renaming stats2 here - only need to change group names (but also do this for every file that is joined!)
##  row change...
paths_groupAandB_stats2$group <- gsub('group A','target',paths_groupAandB_stats2$group)
paths_groupAandB_stats2$group <- gsub('group B','non-target',paths_groupAandB_stats2$group)

paths_groupAandB_stats2_ecstatsonly0 <- paths_groupAandB_stats2 %>%
  ungroup() %>% 
  group_by(group,ec_number) %>% 
  # summarize_all(first) %>%
  summarize(across(everything(), first)) %>%
  select(group,ec_number,genecount_group_by_ecnumber,genepercentage_group_by_ecnumber)

n_distinct(paths_groupAandB_stats2_ecstatsonly0$ec_number)

paths_groupAandB_stats2_ecstatsonlya <- paths_groupAandB_stats2_ecstatsonly0 %>% dplyr::filter(group == "target")
paths_groupAandB_stats2_ecstatsonlyb <- paths_groupAandB_stats2_ecstatsonly0 %>% dplyr::filter(group == "non-target")

paths_groupAandB_stats2_ecstatsonlya <- paths_groupAandB_stats2_ecstatsonlya %>%
  ungroup() %>% 
  rename(genecount_target_group_by_ecnumber = genecount_group_by_ecnumber) %>%
  rename(genepercentage_target_group_by_ecnumber = genepercentage_group_by_ecnumber) %>%
  select(-group)
paths_groupAandB_stats2_ecstatsonlyb <- paths_groupAandB_stats2_ecstatsonlyb %>%
  ungroup() %>% 
  rename(genecount_nontarget_group_by_ecnumber = genecount_group_by_ecnumber) %>%
  rename(genepercentage_nontarget_group_by_ecnumber = genepercentage_group_by_ecnumber) %>%
  select(-group)

paths_groupAandB_stats2_ecstatsonly <- full_join(paths_groupAandB_stats2_ecstatsonlya, paths_groupAandB_stats2_ecstatsonlyb) %>% 
  replace_na(list(genecount_target_group_by_ecnumber = 0)) %>% 
  replace_na(list(genepercentage_target_group_by_ecnumber = 0)) %>% 
  replace_na(list(genecount_nontarget_group_by_ecnumber = 0)) %>% 
  replace_na(list(genepercentage_nontarget_group_by_ecnumber = 0))

n_distinct(paths_groupAandB_stats2_ecstatsonly$ec_number)
rm(paths_groupAandB_stats2_ecstatsonlya)
rm(paths_groupAandB_stats2_ecstatsonlyb)
rm(paths_groupAandB_stats2_ecstatsonly0)

## instead of ref2c adding prom index from mappingstats0
paths_ref2d <- left_join(mapping_stats0, paths_groupAandB_stats2_ecstatsonly)
## now add new shortened pathway name column for excel...
paths_ref2d$pathway_namesh <- substring(paths_ref2d$pathway_name, 1,25)
paths_ref2d <- paths_ref2d %>% unite("pathway_nameshort", pathway_namesh, pathwayid, sep = "_", remove = FALSE)
## need to remove / from pathway names...
paths_ref2d$pathway_nameshort <- gsub("\\/","\\-",paths_ref2d$pathway_nameshort)
paths_ref2d$pathway_namesh <- NULL

## next a pivot, without pulling into groups first
paths_groupAandB_statswidea <- paths_groupAandB_statsa1 %>%
  pivot_wider(names_from = genusspecies, values_from = seen, values_fill = 0, names_sort = TRUE)


## try full join of paths_ref2c & paths_groupAandB_statswide ??
paths_groupAandB_statswiderstill <- full_join(paths_ref2d,paths_groupAandB_statswidea)

## then another full join of PML_bypathway!!! first remove all but one matching column (pathwayid)
PML_bypathway2 <- PML_bypathway %>%
  ungroup() %>% 
  select(-pathway_id,-pathway_name)


## REMOVING SOME NAs & 0s also updating some missing EC stats very early on
paths_groupAandB_statswidest <- full_join(PML_bypathway2,paths_groupAandB_statswiderstill, multiple = "all")
paths_groupAandB_statswidest <- paths_groupAandB_statswidest %>%
  dplyr::filter(ec_description != "NA") %>% 
  dplyr::filter(bvbrc_genes_inpathway > 0)

## now move colummns...
paths_groupAandB_statswidest <- paths_groupAandB_statswidest %>%
  relocate(ec_number, .before = pathwayid) %>%
  relocate(pathway_nameshort, .before = pathwayid) %>%
  relocate(PML, .after = pathwayid) %>%
  relocate(genepercentage_nontarget_group_by_ecnumber, .after = pathwayid) %>%
  relocate(genepercentage_target_group_by_ecnumber, .after = pathwayid) %>%
  relocate(genecount_nontarget_group_by_ecnumber, .after = pathwayid) %>%
  relocate(genecount_target_group_by_ecnumber, .after = pathwayid) %>%
  relocate(ec_description, .after = pathwayid) %>%
  relocate(ec_pindex, .after = avg_p_index_bypathway) %>% 
  relocate(loss_innontarget_group_score, .after = genes_absent_differential) %>%
  relocate(PML_composite, .after = loss_innontarget_group_score) %>%
  relocate(loss_innontarget_group_composite, .after = PML_composite)
#paths_groupAandB_statswidest$avg_p_index_bypathway <- round(paths_groupAandB_statswidest$avg_p_index_bypathway, 2)

## to split into multiple Excel sheets....https://stackoverflow.com/questions/60834065/split-large-dataframe-in-r-and-output-into-separate-sheets-in-a-single-excel-wor

full_dataset_part2 <- paths_groupAandB_statswidest
full_dataset_part2 <- full_dataset_part2 %>% arrange(desc(PML), desc(loss_innontarget_group_score), pathway_id)
ordering_full_dataset_part2b <- full_dataset_part2[["pathway_nameshort"]]
ordering_full_dataset_part2b <- unique(ordering_full_dataset_part2b)
full_dataset_part2$pathway_nameshort <- factor(full_dataset_part2$pathway_nameshort, levels = ordering_full_dataset_part2b, ordered = TRUE)
## remove the perc_by_
output <- split(full_dataset_part2, full_dataset_part2$pathway_nameshort)

wb <- createWorkbook()
for (i in 1:length(output)) {
  addWorksheet(wb, sheetName=names(output[i]))
  writeData(wb, sheet=names(output[i]), x=output[[i]]) # Note [[]]
}

saveWorkbook(wb, file = paste("PML_fulldata_bypathway_",plot_title,Sys.Date(),".xlsx"), overwrite = TRUE)


##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

######################################################
## PLOTS OF PATHWAY PRESENCE (PATHWAYS ON X-AXIS)
######################################################

### FIRST FOCUS PATHWAYS CLUSTERED BY PML SCORE

## want to add perc_A and perc_B 100- perc_missing?
# paths_groupAandB_stats3nonasnofullymissing_PML
# paths_enrichedinBnonas
#theme(panel.background = element_rect(fill = 'black'))
#paths_groupAandB_stats3nonas_PMLorder <- paths_groupAandB_stats3nonas_PML %>% arrange(desc(group),genome_length)
## reminder this will order by group, then by genome size...
orderinggenus <- paths_groupAandB_stats3nonas_PMLorder[["genus"]]
orderinggenus <- unique(orderinggenus)

paths_groupAandB_stats3nonasnofullymissing_PML$genus <- factor(paths_groupAandB_stats3nonasnofullymissing_PML$genus, levels = orderinggenus, ordered = TRUE)
paths_enrichedinBnonas$genus <- factor(paths_enrichedinBnonas$genus, levels = orderinggenus, ordered = TRUE)

## manual 10-color GnBu facsimile: #fafcf6 #ecf4e9 #deecdb #cbe2cf #b1d4c9 #9bc8cc #86b7c9 #6f9fba #5784a5 #446181 #446181
mutedGnBu <- c("#fafcf6", "#ecf4e9", "#deecdb","#cbe2cf", "#b1d4c9","#9bc8cc", "#86b7c9", "#6f9fba", "#5784a5", "#446181", "#446181") 
## finally for plot aesthetics replacing all zeros in perc_in_species_bypathway, perc_in_genus_bypathway, perc_in_group_bypathway with 0.001
## also changing theme from theme_bw to theme_classic & adding background fill #fafcf6
paths_enrichedinBnonas$perc_in_species_bypathway[paths_enrichedinBnonas$perc_in_species_bypathway<0.01] <- 0.01
paths_enrichedinBnonas$perc_in_genus_bypathway[paths_enrichedinBnonas$perc_in_genus_bypathway<0.01] <- 0.01
paths_enrichedinBnonas$perc_in_group_bypathway[paths_enrichedinBnonas$perc_in_group_bypathway<0.01] <- 0.01
paths_groupAandB_stats3nonasnofullymissing_PMLnoref$perc_in_species_bypathway[paths_groupAandB_stats3nonasnofullymissing_PMLnoref$perc_in_species_bypathway<0.01] <- 0.01
paths_groupAandB_stats3nonasnofullymissing_PMLnoref$perc_in_genus_bypathway[paths_groupAandB_stats3nonasnofullymissing_PMLnoref$perc_in_genus_bypathway<0.01] <- 0.01
paths_groupAandB_stats3nonasnofullymissing_PMLnoref$perc_in_group_bypathway[paths_groupAandB_stats3nonasnofullymissing_PMLnoref$perc_in_group_bypathway<0.01] <- 0.01

# ## using the following ggplot settings:
# ## manual 10-color GnBu facsimile: #fafcf6 #ecf4e9 #deecdb #cbe2cf #b1d4c9 #9bc8cc #86b7c9 #6f9fba #5784a5 #446181 #446181
# pathwaysbygenus <- ggplot(paths_enrichedinBnonas, aes(x=pathway_name_with_PML_score, y=genus)) + scale_y_discrete(limits=rev)
# pathwaysbygenus <- pathwaysbygenus + geom_tile(aes(fill = perc_in_genus_bypathway), color="white", linewidth=0.1)
# pathwaysbygenus <- pathwaysbygenus + scale_fill_gradientn(name = "Pathway \npresence", colors = mutedGnBu)
# pathwaysbygenus <- pathwaysbygenus + theme(axis.ticks=element_blank())
# pathwaysbygenus <- pathwaysbygenus + theme_classic(base_family="Helvetica")
# ## background not white but lightest shade of GnBu
# pathwaysbygenus <- pathwaysbygenus + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#fafcf6'))
# pathwaysbygenus

## change labels + ylab("Species") + xlab("Pathway Name with PML")
pathwaysbyspecies <- ggplot(paths_enrichedinBnonas, aes(x=pathway_name_with_PML_score, y=genusspecies)) + scale_y_discrete(limits=rev)
#pathwaysbyspecies <- pathwaysbyspecies + geom_tile(aes(fill = perc_in_species_bypathway), color="thistle1", linewidth=0.1)
pathwaysbyspecies <- pathwaysbyspecies + geom_tile(aes(fill = perc_in_species_bypathway), color="white", linewidth=0.1)
#pathwaysbyspecies <- pathwaysbyspecies + geom_tile(aes(fill = perc_in_species_bypathway), color="white", linewidth=0.1) + geom_text(data = subset(paths_enrichedinBnonas, perc_in_species_bypathway > 0.2), aes(label = perc_in_species_bypathway), color="blue", size = 2)
pathwaysbyspecies <- pathwaysbyspecies + scale_fill_gradientn(name = "Pathway \npresence", colors = mutedGnBu)
##other palette options
## note we can reverse the direction
#pathwaysbyspecies <- pathwaysbyspecies + scale_fill_viridis_c(name = "Pathway \npresence", na.value = "transparent", direction = -1)
## if using Reds color pallete
#pathwaysbyspecies <- pathwaysbyspecies + scale_fill_distiller(name = "Pathway \npresence", palette = "Reds", direction = 1, na.value = "#fff5f0")
## YlGnBu palette ffffd9
#pathwaysbyspecies <- pathwaysbyspecies + scale_fill_distiller(name = "Pathway \npresence", palette = "YlGnBu", direction = 1, na.value = "#ffffd9")
## best GnBu palette f7fcf0
#pathwaysbyspecies <- pathwaysbyspecies + scale_fill_distiller(name = "Pathway \npresence", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
## manual version? not really...
#pathwaysbyspecies <- pathwaysbyspecies + scale_fill_gradient2(name = "Pathway \npresence", low = "#f0f9e8", mid = "#7bccc4", high = "#08589e", midpoint = 60, na.value = "#f7fcf0")
pathwaysbyspecies <- pathwaysbyspecies + theme(axis.ticks=element_blank()) + ylab("Species") + xlab("Pathway Name with PML")
#pathwaysbyspecies <- pathwaysbyspecies + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for focus", plot_title, "pathways missing in target group"))
#pathwaysbyspecies <- pathwaysbyspecies + theme_dark(base_family="Helvetica")

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  pathwaysbyspecies <- pathwaysbyspecies + theme_classic(base_family="Helvetica")
} else {
  pathwaysbyspecies <- pathwaysbyspecies + theme_classic(base_family="sans")
}
## background not white but lightest shade of mutedGnBu
pathwaysbyspecies <- pathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
## other backgrounds
#pathwaysbyspecies <- pathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), axis.text.y = element_text(size = 6))
## dark background (if using theme_dark)
#pathwaysbyspecies <- pathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 45, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#440154', color = '#440154'))
## background not white but lightest shade of GnBu
#pathwaysbyspecies <- pathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 45, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#f7fcf0', color = '#f7fcf0'))
### for lightest shade of YlGnBu #ffffd9
#pathwaysbyspecies <- pathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 45, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#ffffd9', color = '#ffffd9'))
#pathwaysbyspecies <- pathwaysbyspecies + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
pathwaysbyspecies

pathwaysbygenus <- ggplot(paths_enrichedinBnonas, aes(x=pathway_name_with_PML_score, y=genus)) + scale_y_discrete(limits=rev)
pathwaysbygenus <- pathwaysbygenus + geom_tile(aes(fill = perc_in_genus_bypathway), color="white", linewidth=0.1)
#pathwaysbygenus <- pathwaysbygenus + geom_tile(aes(fill = perc_in_species_bypathway), color="white", linewidth=0.1) + geom_text(data = subset(paths_enrichedinBnonas, perc_in_species_bypathway > 0.2), aes(label = perc_in_species_bypathway), color="blue", size = 2)
pathwaysbygenus <- pathwaysbygenus + scale_fill_gradientn(name = "Pathway \npresence", colors = mutedGnBu)
#pathwaysbygenus <- pathwaysbygenus + scale_fill_viridis_c(name = "Pathway \npresence", na.value = "transparent", direction = -1)
## red palettes (second removes strongest colors)
#pathwaysbygenus <- pathwaysbygenus + scale_fill_distiller(name = "Pathway \npresence", palette = "Reds", direction = 1, na.value = "#fff5f0")
#pathwaysbygenus <- pathwaysbygenus + scale_fill_gradient(name = "Pathway \npresence", low="#fff5f0", high="#ef3b2c", na.value = "#fff5f0")
## GnBu palette f7fcf0
#pathwaysbygenus <- pathwaysbygenus + scale_fill_distiller(name = "Pathway \npresence", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
## YlGnBu palette ffffd9
#pathwaysbygenus <- pathwaysbygenus + scale_fill_distiller(name = "Pathway \npresence", palette = "YlGnBu", direction = 1, na.value = "#ffffd9")
pathwaysbygenus <- pathwaysbygenus + theme(axis.ticks=element_blank()) + ylab("Genus") + xlab("Pathway Name with PML")
#pathwaysbygenus <- pathwaysbygenus + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for focus", plot_title, "pathways missing in target group"))

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  pathwaysbygenus <- pathwaysbygenus + theme_classic(base_family="Helvetica")
} else {
  pathwaysbygenus <- pathwaysbygenus + theme_classic(base_family="sans")
}

pathwaysbygenus <- pathwaysbygenus + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
#pathwaysbygenus <- pathwaysbygenus + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), axis.text.y = element_text(size = 6), panel.background = element_rect(fill = '#ffffd9', color = '#ffffd9'))
#pathwaysbygenus <- pathwaysbygenus + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
pathwaysbygenus

## also by group
pathwaysbygroup <- ggplot(paths_enrichedinBnonas, aes(x=pathway_name_with_PML_score, y=group)) + scale_y_discrete()
pathwaysbygroup <- pathwaysbygroup + geom_tile(aes(fill = perc_in_group_bypathway), color="white", linewidth=0.1)
pathwaysbygroup <- pathwaysbygroup + scale_fill_gradientn(name = "Pathway \npresence", colors = mutedGnBu)
## GnBu palette f7fcf0
#pathwaysbygroup <- pathwaysbygroup + scale_fill_distiller(name = "Pathway \npresence", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
pathwaysbygroup <- pathwaysbygroup + theme(axis.ticks=element_blank()) + ylab("Group") + xlab("Pathway Name with PML")

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  pathwaysbygroup <- pathwaysbygroup + theme_classic(base_family="Helvetica")
} else {
  pathwaysbygroup <- pathwaysbygroup + theme_classic(base_family="sans")
}
pathwaysbygroup <- pathwaysbygroup + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
pathwaysbygroup

## to manually make a less 'saturated' version of Reds colorBrewer palette: https://ggplot2.tidyverse.org/reference/scale_hue.html
## or scale_fill_gradient scale_fill_gradient(low="fff5f0", high="#ef3b2c") (tried alpha = 0.8, doesn't work)

### PML PLOTS
## all pathways
#PML_bypathway$pathway_name <- factor(PML_bypathway$pathway_name, levels = ordering_enrichedinB, ordered = TRUE)
## for later aplot combinations, try pulling pathway_name_with_PML
PML_bypathway_all <- PML_bypathway
PML_bypathway0 <- distinct(paths_groupAandB_stats3nonasnofullymissing_PMLnoref, pathway_name)
PML_bypathway <- inner_join(PML_bypathway, PML_bypathway0)

PML_bypathway$xtext <- "PML score"
PML_bypathway <- PML_bypathway %>% unite("pathway_name_with_PML_score", pathway_name, xtext, PML, sep = " ", remove = FALSE)
PML_bypathway$xtext <- NULL
unique(PML_bypathway$pathway_name_with_PML_score)

d <- ggplot(PML_bypathway, aes(x=reorder(pathway_name_with_PML_score, -PML), y=PML))
d <- d + theme_bw() + theme(axis.ticks=element_blank()) + ylab("PML score")
d <- d + geom_col() + theme(axis.text.x = element_blank(), panel.border = element_blank()) + labs(x=NULL)
d
## focus pathways...
PML_bypathway_focus0 <- distinct(paths_enrichedinBnonas, pathway_name)
PML_bypathway_focus <- inner_join(PML_bypathway, PML_bypathway_focus0)
unique(PML_bypathway_focus$pathway_name_with_PML_score)

df <- ggplot(PML_bypathway_focus, aes(x=reorder(pathway_name_with_PML_score, -PML), y=PML))
df <- df + theme_bw() + theme(axis.ticks=element_blank()) + ylab("PML score")
df <- df + geom_col() + theme(axis.text.x = element_blank(), panel.border = element_blank()) + labs(x=NULL)
df

## max widths for alignment in grid.arrange
pathwaysbyspecies_grob <- ggplotGrob(pathwaysbyspecies)
pathwaysbygenus_grob <- ggplotGrob(pathwaysbygenus)
pathwaysbygroup_grob <- ggplotGrob(pathwaysbygroup)
df_grob <- ggplotGrob(df)
d_grob <- ggplotGrob(d)
maxWidth = grid::unit.pmax(pathwaysbyspecies_grob$widths, pathwaysbygenus_grob$widths, pathwaysbygroup_grob$widths)
pathwaysbyspecies_grob$widths <- maxWidth
pathwaysbygenus_grob$widths <- maxWidth
pathwaysbygroup_grob$widths <- maxWidth
df_grob$widths <- maxWidth
d_grob$widths <- maxWidth

## use aplot here?? NO THAT IS BEST FOR ADDING A PHYLO... USE GRID.ARRANGE
## NEED TO SET MAXWIDTH, THEN IS ALIGNED...
#pathwaysbyspecies %>% insert_top(df, height = 0.2)
PMLaboveheatmap <- rbind(c(1,1,1,1,1,1,1),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2),
                        c(2,2,2,2,2,2,2))
pathwaysbyspeciesPML <- grid.arrange(df_grob,pathwaysbyspecies_grob, layout_matrix = PMLaboveheatmap, top = textGrob(paste("Pathway presences by species for focus", plot_title, "pathways missing in target group"),gp=gpar(fontsize=14,font=1)))
pathwaysbygenusPML <- grid.arrange(df_grob,pathwaysbygenus_grob, layout_matrix = PMLaboveheatmap, top = textGrob(paste("Pathway presences by genera for focus", plot_title, "pathways missing in target group"),gp=gpar(fontsize=14,font=1)))
pathwaysbygroupPML <- grid.arrange(df_grob,pathwaysbygroup_grob, layout_matrix = PMLaboveheatmap, top = textGrob(paste("Pathway presences in target group vs. non-target group for focus", plot_title, "pathways missing in target group"),gp=gpar(fontsize=14,font=1)))

## to add titles - can do easily in grid.arrange (for aplot, need to add to topmost plot)
## https://stackoverflow.com/questions/14726078/changing-title-in-multiplot-ggplot2-using-grid-arrange

## start saving all .png to supplemental as well
## save allpathways into subfolder
# data_dir <- "~/pomelo_outputs"
## use two separate plot folders _taxon_by_pathway & _ec_by_taxon_per_pathway
# dir.create("~/pomelo_outputs/supplemental_plots_taxon_by_pathway")
# dir.create("~/pomelo_outputs/supplemental_plots_ec_by_taxon_per_pathway")

## note new way to do it more generically 
suppdir1tocreate <- paste0(data_dir,"/supplemental_plots_taxon_by_pathway" )
suppdir2tocreate <- paste0(data_dir,"/supplemental_plots_ec_by_taxon_per_pathway" )
suppressWarnings(dir.create(suppdir1tocreate))
suppressWarnings(dir.create(suppdir2tocreate))

ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_byspecies_and_pathway_",plot_title,Sys.Date(),".png", sep=""), pathwaysbyspeciesPML, width = 32, height = 16, units = "in", limitsize = FALSE)
ggsave(filename = paste("heatmaps_focuspathways_byspecies_and_pathway_",plot_title,Sys.Date(),".pdf", sep=""), pathwaysbyspeciesPML, width = 32, height = 16, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygenus_and_pathway_",plot_title,Sys.Date(),".png", sep=""), pathwaysbygenusPML, width = 16, height = 8, units = "in", limitsize = FALSE)
## moving by genus & by group to supplemental folder
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygenus_and_pathway_",plot_title,Sys.Date(),".pdf", sep=""), pathwaysbygenusPML, width = 16, height = 8, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygroup_and_pathway_",plot_title,Sys.Date(),".png", sep=""), pathwaysbygroupPML, width = 16, height = 8, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygroup_and_pathway_",plot_title,Sys.Date(),".pdf", sep=""), pathwaysbygroupPML, width = 16, height = 8, units = "in", limitsize = FALSE)


####################################################
### FOCUS PATHWAYS ONLY CLUSTERED BY PATHWAY PRESENCE OVERALL SIMILARITY

## instead try clustering?? https://uc-r.github.io/hc_clustering
# # Hierarchical clustering using Complete Linkage
# hc1 <- hclust(d, method = "complete" )
# ## from https://rstudio-pubs-static.s3.amazonaws.com/707309_cec3abaa920b4b42afedcb002bbc36a6.html
## might want to try pheatmap package...https://btep.ccr.cancer.gov/docs/data-visualization-with-r/Lesson5_intro_to_ggplot/

## better example here https://stackoverflow.com/questions/25528059/cluster-data-in-heat-map-in-r-ggplot
## first make tibble of just x & y we want from above...
pathwaysbyspecies_tibble <- paths_enrichedinBnonas %>% 
  select(pathway_name_with_PML_score, genusspecies, perc_in_species_bypathway) %>% 
  distinct()
pathwaysbygenus_tibble <- paths_enrichedinBnonas %>% 
  select(pathway_name_with_PML_score, genus, perc_in_genus_bypathway) %>% 
  distinct()

## now pivot wider, and make rownames for the hclust command, then hclust
pathwaysbyspecies_wider <- pathwaysbyspecies_tibble %>% 
  pivot_wider(
    names_from = pathway_name_with_PML_score, 
    values_from = perc_in_species_bypathway,
    values_fill = 0)
pathwaysbyspecies_wider <- as.data.frame(pathwaysbyspecies_wider)
rownames(pathwaysbyspecies_wider) <- pathwaysbyspecies_wider$genusspecies
pathwaysbyspecies_wider$genusspecies <- NULL
## default hclust settings - but using ward.D instead of complete
hc.cols <- hclust(dist(t(pathwaysbyspecies_wider)), method = "ward.D")
#plot(hc.cols)

## to use cosine distances
# # https://community.rstudio.com/t/r-code-calculating-cosine-similarity-pairwise/125015/2 & https://stackoverflow.com/questions/52391558/hierarchical-clustering-using-cosine-distance-in-r

# don't transform, but divide by 100 - (pathwaysbyspecies_wider) / 100
#cosine_dist <- log2(pathwaysbyspecies_wider)
# cosine_dist <- cosine((pathwaysbyspecies_wider/100))
# cosine_dist <- as.dist(cosine_dist)
# #sq_cosine_dist <- sqrt(cosine_dist)
# hc.cols <- hclust(cosine_dist, method = "ward.D")

## can skip this: just use species level info for clustering...there is more information there anyway!
# pathwaysbygenus_wider <- pathwaysbygenus_tibble %>% 
#   pivot_wider(
#     names_from = pathway_name_with_PML_score, 
#     values_from = perc_in_genus_bypathway,
#     values_fill = 0)
# pathwaysbygenus_wider <- as.data.frame(pathwaysbygenus_wider)
# rownames(pathwaysbygenus_wider) <- pathwaysbygenus_wider$genus
# pathwaysbygenus_wider$genus <- NULL
# hc.cols2 <- hclust(dist(t(pathwaysbygenus_wider)), method = "complete")

## now make new levels for plots
orderingclust <- hc.cols$order
orderingclust
# levels = rownames(pathwaysbygenus_wider)[orderingclust]

## re ordering pathways...just for this plot!
paths_enrichedinBnonas2 <- paths_enrichedinBnonas
paths_enrichedinBnonas2$pathway_name_with_PML_score <- factor(paths_enrichedinBnonas2$pathway_name_with_PML_score, levels = colnames(pathwaysbyspecies_wider)[orderingclust], ordered = TRUE)
## if fails need to repeat
# paths_enrichedinBnonas <- paths_enrichedinB %>%
#   dplyr::filter(full_ec_number != "NA")
# paths_enrichedinBnonas$genus <- factor(paths_enrichedinBnonas$genus, levels = orderinggenus, ordered = TRUE)

## repeating but now sorted by clusters
pathwaysbyspecies <- ggplot(paths_enrichedinBnonas2, aes(x=pathway_name_with_PML_score, y=genusspecies)) + scale_y_discrete(limits=rev)
pathwaysbyspecies <- pathwaysbyspecies + geom_tile(aes(fill = perc_in_species_bypathway), color="white", linewidth=0.1)
pathwaysbyspecies <- pathwaysbyspecies + scale_fill_gradientn(name = paste("Pathway presences\nby species for focus\n", plot_title, "pathways\nmissing in target group\n\n"), colors = mutedGnBu)
pathwaysbyspecies <- pathwaysbyspecies + theme(axis.ticks=element_blank()) + ylab("Species") + xlab("Pathway Name with PML")
## trying to add plot title to bottom as second lined of xlabel...not ideal - instead adding to legend???
#pathwaysbyspecies <- pathwaysbyspecies + theme(axis.ticks=element_blank()) + ylab("Species") + xlab(paste("Pathway Name with PML\n\nPathway presences by species for focus", plot_title, "pathways missing in target group"))
#pathwaysbyspecies <- pathwaysbyspecies + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for focus", plot_title, "pathways missing in target group"))

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  pathwaysbyspecies <- pathwaysbyspecies + theme_classic(base_family="Helvetica")
} else {
  pathwaysbyspecies <- pathwaysbyspecies + theme_classic(base_family="sans")
}

pathwaysbyspecies <- pathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
#pathwaysbyspecies <- pathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 45, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#440154', color = '#440154'))
#pathwaysbyspecies <- pathwaysbyspecies + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
pathwaysbyspecies

pathwaysbygenus <- ggplot(paths_enrichedinBnonas2, aes(x=pathway_name_with_PML_score, y=genus)) + scale_y_discrete(limits=rev)
pathwaysbygenus <- pathwaysbygenus + geom_tile(aes(fill = perc_in_genus_bypathway), color="white", linewidth=0.1)
pathwaysbygenus <- pathwaysbygenus + scale_fill_gradientn(name = paste("Pathway presences\nby genus for focus\n", plot_title, "pathways\nmissing in target group\n\n"), colors = mutedGnBu)
pathwaysbygenus <- pathwaysbygenus + theme(axis.ticks=element_blank()) + ylab("Genus") + xlab("Pathway Name with PML")
#pathwaysbygenus <- pathwaysbygenus + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for focus", plot_title, "pathways missing in target group"))

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  pathwaysbygenus <- pathwaysbygenus + theme_classic(base_family="Helvetica")
} else {
  pathwaysbygenus <- pathwaysbygenus + theme_classic(base_family="sans")
}
pathwaysbygenus <- pathwaysbygenus + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
#pathwaysbygenus <- pathwaysbygenus + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
pathwaysbygenus

## also by group
pathwaysbygroup <- ggplot(paths_enrichedinBnonas2, aes(x=pathway_name_with_PML_score, y=group)) + scale_y_discrete()
pathwaysbygroup <- pathwaysbygroup + geom_tile(aes(fill = perc_in_group_bypathway), color="white", linewidth=0.1)
## GnBu palette f7fcf0
pathwaysbygroup <- pathwaysbygroup + scale_fill_gradientn(name = paste("Pathway presences\nin target vs. non-target\ngroup for focus\n", plot_title, "pathways\nmissing in target group\n\n"), colors = mutedGnBu)
pathwaysbygroup <- pathwaysbygroup + theme(axis.ticks=element_blank()) + ylab("Group") + xlab("Pathway Name with PML")

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  pathwaysbygroup <- pathwaysbygroup + theme_classic(base_family="Helvetica")
} else {
  pathwaysbygroup <- pathwaysbygroup + theme_classic(base_family="sans")
}
pathwaysbygroup <- pathwaysbygroup + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
pathwaysbygroup

## PML PLOTS SORTED BY CLUSTERS...
### just PML plots but after clustering....
# ## all pathways
# PML_bypathwayc <- PML_bypathway
# PML_bypathwayc$pathway_name <- factor(PML_bypathwayc$pathway_name, levels = unique(PML_bypathwayc$pathway_name)[orderingclust2], ordered = TRUE)
# d2 <- ggplot(PML_bypathwayc, aes(x=pathway_name, y=PML))
# d2 <- d2 + theme_bw() + theme(axis.ticks=element_blank()) + ylab("PML score")
# d2 <- d2 + geom_col() + theme(axis.text.x = element_blank(), panel.border = element_blank()) + labs(x=NULL)
# d2

## focus pathways. for all colnames(pathwaysbyspecies_wider2)
PML_bypathwayc_focus <- PML_bypathway_focus
PML_bypathwayc_focus$pathway_name_with_PML_score <- factor(PML_bypathwayc_focus$pathway_name_with_PML_score, levels = colnames(pathwaysbyspecies_wider)[orderingclust], ordered = TRUE)
d2f <- ggplot(PML_bypathwayc_focus, aes(x=pathway_name_with_PML_score, y=PML))
d2f <- d2f + theme_bw() + theme(axis.ticks=element_blank()) + ylab("PML score")
d2f <- d2f + geom_col() + theme(axis.text.x = element_blank(), panel.border = element_blank()) + labs(x=NULL)
d2f

## only really need clustered_withdendro!
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_byspecies_and_pathway_clustered_",plot_title,Sys.Date(),".png", sep=""), pathwaysbyspecies, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_byspecies_and_pathway_clustered_",plot_title,Sys.Date(),".pdf", sep=""), pathwaysbyspecies, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygenus_and_pathway_clustered_",plot_title,Sys.Date(),".png", sep=""), pathwaysbygenus, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygenus_and_pathway_clustered_",plot_title,Sys.Date(),".pdf", sep=""), pathwaysbygenus, width = 16, height = 8, units = "in", limitsize = FALSE)

####################################################
### ALL PATHWAYS CLUSTERED BY PML SCORE

## now all pathways using paths_groupAandB_stats3nonasnofullymissing_PML or paths_groupAandB_stats3nonasnofullymissing_PMLnoref
paths_groupAandB_stats3nonasnofullymissing_PMLnoref <- paths_groupAandB_stats3nonasnofullymissing_PMLnoref %>% arrange(desc(PML), desc(loss_innontarget_group_score), pathway_id)
ordering_fullpathways <- paths_groupAandB_stats3nonasnofullymissing_PMLnoref[["pathway_name_with_PML_score"]]
ordering_fullpathways <- unique(ordering_fullpathways)
paths_groupAandB_stats3nonasnofullymissing_PMLnoref$pathway_name_with_PML_score <- factor(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$pathway_name_with_PML_score, levels = ordering_fullpathways, ordered = TRUE)

allpathwaysbyspecies <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PMLnoref, aes(x=pathway_name_with_PML_score, y=genusspecies)) + scale_y_discrete(limits=rev)
allpathwaysbyspecies <- allpathwaysbyspecies + geom_tile(aes(fill = perc_in_species_bypathway), color="white", linewidth=0.1)
allpathwaysbyspecies <- allpathwaysbyspecies + scale_fill_gradientn(name = "Pathway \npresence", colors = mutedGnBu)
allpathwaysbyspecies <- allpathwaysbyspecies + theme(axis.ticks=element_blank()) + ylab("Species") + xlab("Pathway Name with PML")
#allpathwaysbyspecies <- allpathwaysbyspecies + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for all", plot_title, "pathways missing in target group"))

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  allpathwaysbyspecies <- allpathwaysbyspecies + theme_classic(base_family="Helvetica")
} else {
  allpathwaysbyspecies <- allpathwaysbyspecies + theme_classic(base_family="sans")
}
allpathwaysbyspecies <- allpathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
#allpathwaysbyspecies <- allpathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 45, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#440154', color = '#440154'))
#allpathwaysbyspecies <- allpathwaysbyspecies + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
allpathwaysbyspecies

paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genus <- factor(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genus, levels = orderinggenus, ordered = TRUE)

allpathwaysbygenus <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PMLnoref, aes(x=pathway_name_with_PML_score, y=genus)) + scale_y_discrete(limits=rev)
allpathwaysbygenus <- allpathwaysbygenus + geom_tile(aes(fill = perc_in_genus_bypathway), color="white", linewidth=0.1)
allpathwaysbygenus <- allpathwaysbygenus + scale_fill_gradientn(name = "Pathway \npresence", colors = mutedGnBu)
allpathwaysbygenus <- allpathwaysbygenus + theme(axis.ticks=element_blank()) + ylab("Genus") + xlab("Pathway Name with PML")
#allpathwaysbygenus <- allpathwaysbygenus + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for all", plot_title, "pathways missing in target group"))

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  allpathwaysbygenus <- allpathwaysbygenus + theme_classic(base_family="Helvetica")
} else {
  allpathwaysbygenus <- allpathwaysbygenus + theme_classic(base_family="sans")
}
allpathwaysbygenus <- allpathwaysbygenus + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
#allpathwaysbygenus <- allpathwaysbygenus + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
allpathwaysbygenus

## also by group
allpathwaysbygroup <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PMLnoref, aes(x=pathway_name_with_PML_score, y=group)) + scale_y_discrete()
allpathwaysbygroup <- allpathwaysbygroup + geom_tile(aes(fill = perc_in_group_bypathway), color="white", linewidth=0.1)
## GnBu palette f7fcf0
allpathwaysbygroup <- allpathwaysbygroup + scale_fill_gradientn(name = "Pathway \npresence", colors = mutedGnBu)
allpathwaysbygroup <- allpathwaysbygroup + theme(axis.ticks=element_blank()) + ylab("Group") + xlab("Pathway Name with PML")

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  allpathwaysbygroup <- allpathwaysbygroup + theme_classic(base_family="Helvetica")
} else {
  allpathwaysbygroup <- allpathwaysbygroup + theme_classic(base_family="sans")
}
allpathwaysbygroup <- allpathwaysbygroup + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
allpathwaysbygroup

### grobs & maxwidths
allpathwaysbyspecies_grob <- ggplotGrob(allpathwaysbyspecies)
allpathwaysbygenus_grob <- ggplotGrob(allpathwaysbygenus)
allpathwaysbygroup_grob <- ggplotGrob(allpathwaysbygroup)
d_grob <- ggplotGrob(d)
maxWidth = grid::unit.pmax(allpathwaysbyspecies_grob$widths, allpathwaysbygenus_grob$widths, allpathwaysbygroup_grob$widths)
allpathwaysbyspecies_grob$widths <- maxWidth
allpathwaysbygenus_grob$widths <- maxWidth
allpathwaysbygroup_grob$widths <- maxWidth
d_grob$widths <- maxWidth
allpathwaysbyspeciesPML <- grid.arrange(d_grob,allpathwaysbyspecies_grob, layout_matrix = PMLaboveheatmap, top = textGrob(paste("Pathway presences by species for all", plot_title, "pathways"),gp=gpar(fontsize=14,font=1)))
allpathwaysbygenusPML <- grid.arrange(d_grob,allpathwaysbygenus_grob, layout_matrix = PMLaboveheatmap, top = textGrob(paste("Pathway presences by genera for all", plot_title, "pathways"),gp=gpar(fontsize=14,font=1)))
allpathwaysbygroupPML <- grid.arrange(d_grob,allpathwaysbygroup_grob, layout_matrix = PMLaboveheatmap, top = textGrob(paste("Pathway presences in target group vs. non-target group for all", plot_title, "pathways"),gp=gpar(fontsize=14,font=1)))


ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_byspecies_and_pathway_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbyspeciesPML, width = 32, height = 16, units = "in", limitsize = FALSE)
ggsave(filename = paste("heatmaps_allpathways_byspecies_and_pathway_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbyspeciesPML, width = 32, height = 16, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygenus_and_pathway_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbygenusPML, width = 16, height = 8, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygenus_and_pathway_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbygenusPML, width = 16, height = 8, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygroup_and_pathway_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbygroupPML, width = 16, height = 8, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygroup_and_pathway_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbygroupPML, width = 16, height = 8, units = "in", limitsize = FALSE)

####################################################
### ALL PATHWAYS CLUSTERED BY PATHWAY PRESENCE OVERALL SIMILARITY

## need to re-cluster all of the pathways...
## first unfactor...
paths_groupAandB_stats3nonasnofullymissing_PMLnoref <- paths_groupAandB_stats3nonasnofullymissing_PML %>%
  dplyr::filter(full_ec_number != "NA")

pathwaysbyspecies_tibble2 <- paths_groupAandB_stats3nonasnofullymissing_PMLnoref %>% 
  select(pathway_name_with_PML_score, genusspecies, perc_in_species_bypathway) %>% 
  distinct()
pathwaysbygenus_tibble2 <- paths_groupAandB_stats3nonasnofullymissing_PMLnoref %>% 
  select(pathway_name_with_PML_score, genus, perc_in_genus_bypathway) %>% 
  distinct()

## now pivot wider, and make rownames for the hclust command, then hclust
pathwaysbyspecies_wider2 <- pathwaysbyspecies_tibble2 %>% 
  pivot_wider(
    names_from = pathway_name_with_PML_score, 
    values_from = perc_in_species_bypathway,
    values_fill = 0)
pathwaysbyspecies_wider2 <- as.data.frame(pathwaysbyspecies_wider2)
rownames(pathwaysbyspecies_wider2) <- pathwaysbyspecies_wider2$genusspecies
pathwaysbyspecies_wider2$genusspecies <- NULL
hc.cols22 <- hclust(dist(t(pathwaysbyspecies_wider2)), method = "ward.D")

## if using cosine distance...
# cosine_dist22 <- cosine((pathwaysbyspecies_wider2/100))
# cosine_dist22 <- as.dist(cosine_dist22)
# hc.cols22 <- hclust(cosine_dist22, method = "ward.D")
# orderingclust2 <- hc.cols22$order
# orderingclust2

#plot(hc.cols22)
orderingclust2 <- hclust(dist(t(pathwaysbyspecies_wider2)), method = "ward.D")$order
orderingclust2
#tmatrix <- t(pathwaysbyspecies_wider2)

## PML PLOTS OF ALL PATHWAYS SORTED BY CLUSTERS...
### just PML plots but after clustering....
## all pathways
PML_bypathwayc <- PML_bypathway
PML_bypathwayc$pathway_name_with_PML_score <- factor(PML_bypathwayc$pathway_name_with_PML_score, levels = unique(PML_bypathwayc$pathway_name_with_PML_score)[orderingclust2], ordered = TRUE)
d2 <- ggplot(PML_bypathwayc, aes(x=pathway_name_with_PML_score, y=PML))
d2 <- d2 + theme_bw() + theme(axis.ticks=element_blank()) + ylab("PML score")
d2 <- d2 + geom_col() + theme(axis.text.x = element_blank(), panel.border = element_blank()) + labs(x=NULL)
d2
## BE CAREFUL - LOOKS INCORRECT ON ITS OWN BUT OKAY IN APLOT BELOW...

# ## focus pathways
# PML_bypathwayc_focus <- PML_bypathway_focus
# PML_bypathwayc_focus$pathway_name_with_PML_score <- factor(PML_bypathwayc_focus$pathway_name_with_PML_score, levels = unique(PML_bypathwayc$pathway_name_with_PML_score)[orderingclust], ordered = TRUE)
# d2f <- ggplot(PML_bypathwayc_focus, aes(x=pathway_name_with_PML_score, y=PML))
# d2f <- d2f + theme_bw() + theme(axis.ticks=element_blank()) + ylab("PML score")
# d2f <- d2f + geom_col() + theme(axis.text.x = element_blank()) + labs(x=NULL)
# d2f
# rm(PML_bypathway_focus0)


## replotting, sorted by clusters
paths_groupAandB_stats3nonasnofullymissing_PMLnoref$pathway_name_with_PML_score <- factor(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$pathway_name_with_PML_score, levels = colnames(pathwaysbyspecies_wider2)[orderingclust2], ordered = TRUE)
#paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genus <- factor(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genus, levels = orderinggenus, ordered = TRUE)

allpathwaysbyspecies <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PMLnoref, aes(x=pathway_name_with_PML_score, y=genusspecies)) + scale_y_discrete(limits=rev)
allpathwaysbyspecies <- allpathwaysbyspecies + geom_tile(aes(fill = perc_in_species_bypathway), color="white", linewidth=0.1)
allpathwaysbyspecies <- allpathwaysbyspecies + scale_fill_gradientn(name = paste("Pathway presences\nby species for all\n", plot_title, "\npathways"), colors = mutedGnBu)
allpathwaysbyspecies <- allpathwaysbyspecies + theme(axis.ticks=element_blank()) + ylab("Species") + xlab("Pathway Name with PML")
#allpathwaysbyspecies <- allpathwaysbyspecies + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for all", plot_title, "pathways missing in target group"))

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  allpathwaysbyspecies <- allpathwaysbyspecies + theme_classic(base_family="Helvetica")
} else {
  allpathwaysbyspecies <- allpathwaysbyspecies + theme_classic(base_family="sans")
}
allpathwaysbyspecies <- allpathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
#allpathwaysbyspecies <- allpathwaysbyspecies + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 45, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#440154', color = '#440154'))
#allpathwaysbyspecies <- allpathwaysbyspecies + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
allpathwaysbyspecies


paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genus <- factor(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genus, levels = orderinggenus, ordered = TRUE)

allpathwaysbygenus <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PMLnoref, aes(x=pathway_name_with_PML_score, y=genus)) + scale_y_discrete(limits=rev)
allpathwaysbygenus <- allpathwaysbygenus + geom_tile(aes(fill = perc_in_genus_bypathway), color="white", linewidth=0.1)
allpathwaysbygenus <- allpathwaysbygenus + scale_fill_gradientn(name = paste("Pathway presences\nby genera for all\n", plot_title, "\npathways"), colors = mutedGnBu)
allpathwaysbygenus <- allpathwaysbygenus + theme(axis.ticks=element_blank()) + ylab("Genus") + xlab("Pathway Name with PML")
#allpathwaysbygenus <- allpathwaysbygenus + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for all", plot_title, "pathways missing in target group"))

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  allpathwaysbygenus <- allpathwaysbygenus + theme_classic(base_family="Helvetica")
} else {
  allpathwaysbygenus <- allpathwaysbygenus + theme_classic(base_family="sans")
}
allpathwaysbygenus <- allpathwaysbygenus + theme(plot.title=element_text(hjust=0), 
                                                 axis.text.x = element_text(angle = 80, hjust=1), 
                                                 panel.background = element_rect(fill = '#fafcf6'))
#allpathwaysbygenus <- allpathwaysbygenus + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
allpathwaysbygenus

## also by group
allpathwaysbygroup <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PMLnoref, aes(x=pathway_name_with_PML_score, y=group)) + scale_y_discrete()
allpathwaysbygroup <- allpathwaysbygroup + geom_tile(aes(fill = perc_in_group_bypathway), color="white", linewidth=0.1)
## GnBu palette f7fcf0
allpathwaysbygroup <- allpathwaysbygroup + scale_fill_gradientn(name = paste("Pathway presences\nin target vs. non-target\ngroup for all", plot_title, "\npathways"), colors = mutedGnBu)
allpathwaysbygroup <- allpathwaysbygroup + theme(axis.ticks=element_blank()) + ylab("Group") + xlab("Pathway Name with PML")

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  allpathwaysbygroup <- allpathwaysbygroup + theme_classic(base_family="Helvetica")
} else {
  allpathwaysbygroup <- allpathwaysbygroup + theme_classic(base_family="sans")
}
allpathwaysbygroup <- allpathwaysbygroup + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
allpathwaysbygroup

#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_byspecies_and_pathway_clustered_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbyspecies, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_byspecies_and_pathway_clustered_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbyspecies, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygenus_and_pathway_clustered_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbygenus, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygenus_and_pathway_clustered_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbygenus, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygroup_and_pathway_clustered_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbygroup, width = 16, height = 8, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygroup_and_pathway_clustered_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbygroup, width = 16, height = 8, units = "in", limitsize = FALSE)

## see about adding the tree to the plot?? can use ggdendro 
ggdendrogram(hc.cols22, rotate = FALSE, size = 2)
#allpathways_dendro <- ggdendrogram(data = hc.cols22, rotate = FALSE)
allpathways_dendronolabels <- ggdendrogram(data = hc.cols22, rotate = FALSE, theme_dendro = TRUE, leaf_labels = FALSE, labels = FALSE) ## dendrogram faces downward
allpathways_dendronolabels1 <- allpathways_dendronolabels + theme(axis.title = element_text(size = 0), axis.text = element_text(size = 0))
allpathways_dendronolabels2 <- allpathways_dendronolabels + theme(plot.margin = margin(0, 1.9, 0, 2.8, "cm"), axis.title = element_text(size = 0), axis.text = element_text(size = 0))
## continue and try more ggplot bits..

## ALTERNATE TACK - SAVE HCLUST AS A TREE, THEN USE GGTREE TYPE COMMANDS...https://stackoverflow.com/questions/21727820/how-to-create-a-newick-file-from-a-cluster-in-r
allpathways_tree <- as.phylo(hc.cols22) 

#allpathways_tree$edge.length <- (log(allpathways_tree$edge.length)+(min(log(allpathways_tree$edge.length))*-1))

# allpathways_ggtree <- ggtree(allpathways_tree) + geom_tiplab()
# p <- ggtree(allpathways_dendronolabels, branch.length = "none") + 
#   geom_tiplab() + theme(legend.position='none')

## now try aplot ...further below should be able to use this with heatmaps + phylogeny as well!
## or aplot? this combines 1 heatmap & one tree...want to add heatmap & PML plot & tree though!!!
#p2 %>% insert_left(g) %>% insert_right(p1, width=.5)
## note if using cosine distances, modify to  ggtree(allpathways_tree, branch.length='none')
allpathways_ggtree2 <- ggtree(allpathways_tree) + layout_dendrogram()

## tree & heatmap - works!
allpathwaysbyspecies %>% insert_top(allpathways_ggtree2, height = 0.2)
## PML with tree above & heatmap below
## allpathwaysbyspeciesPML %>% insert_top(allpathways_ggtree2, height = 0.2)

### tree + PML + heatmap - THIS WORKS, USE
allpathwaysbyspecies_withdendroontop <- allpathwaysbyspecies %>% insert_top(d2, height = 0.2) %>% insert_top(allpathways_ggtree2, height = 0.2)
allpathwaysbygenus_withdendroontop <- allpathwaysbygenus %>% insert_top(d2, height = 0.2) %>% insert_top(allpathways_ggtree2, height = 0.2)
allpathwaysbygroup_withdendroontop <- allpathwaysbygroup %>% insert_top(d2, height = 0.4) %>% insert_top(allpathways_ggtree2, height = 0.2)

#allpathwaysbygenusandspecies_withdendroontop <- grid.arrange(allpathways_dendronolabels2,allpathwaysbygenus,allpathwaysbyspecies, layout_matrix = treeaboveheatmap3)
# twoplots_onlyphylo <- grid.arrange(gsizeforphylo,dendro_plotref_forphyloheatmap, ncol=1, nrow=2, widths = unit(4.5, "in"))
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_byspecies_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbyspecies_withdendroontop, width = 32, height = 16, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_byspecies_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbyspecies_withdendroontop, width = 32, height = 16, units = "in", limitsize = FALSE)
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygenus_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), allpathwaysbygenus_withdendroontop, width = 16, height = 8, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygenus_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbygenus_withdendroontop, width = 16, height = 8, units = "in", limitsize = FALSE)

ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygroup_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".pdf", sep=""), allpathwaysbygroup_withdendroontop, width = 16, height = 8, units = "in", limitsize = FALSE)

## trying alternate way to save as png this works
png(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_byspecies_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), width = 32, height = 16, units = "in", res = 300)
Sys.sleep(2)
allpathwaysbyspecies %>% insert_top(d2, height = 0.2) %>% insert_top(allpathways_ggtree2, height = 0.2)
Sys.sleep(2)
dev.off()
Sys.sleep(2)

png(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygenus_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), width = 16, height = 8, units = "in", res = 300)
Sys.sleep(2)
allpathwaysbygenus %>% insert_top(d2, height = 0.2) %>% insert_top(allpathways_ggtree2, height = 0.2)
Sys.sleep(2)
dev.off()
Sys.sleep(2)

png(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_allpathways_bygroup_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), width = 16, height = 8, units = "in", res = 300)
Sys.sleep(2)
allpathwaysbygroup %>% insert_top(d2, height = 0.4) %>% insert_top(allpathways_ggtree2, height = 0.2)
Sys.sleep(2)
dev.off()
Sys.sleep(2)

## other alternatives to grid.arrange - removed

#################
## focus pathways...USE hc.cols
## not using this code except rename
# focuspathways_dendronolabels <- ggdendrogram(data = hc.cols, rotate = FALSE, theme_dendro = TRUE, leaf_labels = FALSE, labels = FALSE) ## dendrogram faces downward
# focuspathways_dendronolabels1 <- focuspathways_dendronolabels + theme(axis.title = element_text(size = 0), axis.text = element_text(size = 0))
# focuspathways_dendronolabels2 <- focuspathways_dendronolabels + theme(plot.margin = margin(0, 1.9, 0, 2.8, "cm"), axis.title = element_text(size = 0), axis.text = element_text(size = 0))
# ## continue and try more ggplot bits..
# 
# #focuspathways_dendro_forphyloheatmap <- print(focuspathways_dendro) + scale_y_reverse() + theme(legend.position="left", legend.key.width=unit(0.2,"cm"))
# # using this one!
# #focuspathways_dendro_forphyloheatmap <- print(focuspathways_dendronolabels) + scale_y_reverse() ## dendrogram faces upward
focuspathwaysbyspecies <- pathwaysbyspecies
focuspathwaysbygenus <- pathwaysbygenus
focuspathwaysbygroup <- pathwaysbygroup

# focuspathwaysbyspecies_grob <- ggplotGrob(focuspathwaysbyspecies)
# focuspathwaysbygenus_grob <- ggplotGrob(pathwaysbygenus)
# focuspathwaysbygroup_grob <- ggplotGrob(pathwaysbygroup)
# # grid::grid.newpage()
# # grid::grid.draw(rbind(gA, gB))
# maxWidth = grid::unit.pmax(focuspathwaysbyspecies_grob$widths, focuspathwaysbygenus_grob$widths, focuspathwaysbygroup_grob$widths)
# # gA$widths[2:5] <- as.list(maxWidth)
# # gB$widths[2:5] <- as.list(maxWidth)
# # gC$widths[2:5] <- as.list(maxWidth)
# focuspathwaysbyspecies_grob$widths <- maxWidth
# focuspathwaysbygenus_grob$widths <- maxWidth
# focuspathwaysbygroup_grob$widths <- maxWidth
# ## adding margin instead?? plot + theme(plot.margin = margin(2, 2, 2, 2, "cm") see above for where it is added
# focuspathwaysbyspecies_withdendroontop <- grid.arrange(focuspathways_dendronolabels2,focuspathwaysbyspecies, layout_matrix = treeaboveheatmap2)
# focuspathwaysbygenus_withdendroontop <- grid.arrange(focuspathways_dendronolabels2,focuspathwaysbygenus_grob, layout_matrix = treeaboveheatmap2)
# focuspathwaysbygroup_withdendroontop <- grid.arrange(focuspathways_dendronolabels2,focuspathwaysbygroup_grob, layout_matrix = treeaboveheatmap2)

## see about adding the tree to the plot?? can use ggdendro 
ggdendrogram(hc.cols, rotate = FALSE, size = 2)
#focuspathways_dendro <- ggdendrogram(data = hc.cols22, rotate = FALSE)
focuspathways_dendronolabels <- ggdendrogram(data = hc.cols, rotate = FALSE, theme_dendro = TRUE, leaf_labels = FALSE, labels = FALSE) ## dendrogram faces downward
focuspathways_dendronolabels1 <- focuspathways_dendronolabels + theme(axis.title = element_text(size = 0), axis.text = element_text(size = 0))
focuspathways_dendronolabels2 <- focuspathways_dendronolabels + theme(plot.margin = margin(0, 1.9, 0, 2.8, "cm"), axis.title = element_text(size = 0), axis.text = element_text(size = 0))
## continue and try more ggplot bits..

## ALTERNATE TACK - SAVE HCLUST AS A TREE, THEN USE GGTREE TYPE COMMANDS...https://stackoverflow.com/questions/21727820/how-to-create-a-newick-file-from-a-cluster-in-r
focuspathways_tree <- as.phylo(hc.cols) 
#focuspathways_tree$edge.length <- (log(focuspathways_tree$edge.length)+(min(log(focuspathways_tree$edge.length))*-1))

# focuspathways_ggtree <- ggtree(focuspathways_tree, branch.length='none') + geom_tiplab()
# p <- ggtree(focuspathways_dendronolabels, branch.length = "none") + 
#   geom_tiplab() + theme(legend.position='none')

## now try aplot ...further below should be able to use this with heatmaps + phylogeny as well!
## or aplot? this combines 1 heatmap & one tree...want to add heatmap & PML plot & tree though!!!

## note if using cosine distances, modify to  ggtree(focuspathways_tree, branch.length='none')
focuspathways_ggtree2 <- ggtree(focuspathways_tree) + layout_dendrogram()

## tree & heatmap - works!
focuspathwaysbyspecies %>% insert_top(focuspathways_ggtree2, height = 0.2)
## PML with tree above & heatmap below. d2   focuspathwaysbyspeciesPML
## focuspathwaysbyspeciesPML %>% insert_top(focuspathways_ggtree2, height = 0.2)

### tree + PML + heatmap - THIS WORKS, USE
focuspathwaysbyspecies_withdendroontop <- focuspathwaysbyspecies %>% insert_top(d2f, height = 0.2) %>% insert_top(focuspathways_ggtree2, height = 0.2)
focuspathwaysbygenus_withdendroontop <- focuspathwaysbygenus %>% insert_top(d2f, height = 0.2) %>% insert_top(focuspathways_ggtree2, height = 0.2)
focuspathwaysbygroup_withdendroontop <- focuspathwaysbygroup %>% insert_top(d2f, height = 0.4) %>% insert_top(focuspathways_ggtree2, height = 0.2)


#focuspathwaysbygenusandspecies_withdendroontop <- grid.arrange(focuspathways_dendronolabels2,focuspathwaysbygenus,focuspathwaysbyspecies, layout_matrix = treeaboveheatmap3)
# twoplots_onlyphylo <- grid.arrange(gsizeforphylo,dendro_plotref_forphyloheatmap, ncol=1, nrow=2, widths = unit(4.5, "in"))
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_byspecies_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), focuspathwaysbyspecies_withdendroontop, width = 32, height = 16, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_byspecies_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".pdf", sep=""), focuspathwaysbyspecies_withdendroontop, width = 32, height = 16, units = "in", limitsize = FALSE)

## moving by genus and by group to supplemental_plots_taxon_by_pathway/ folder
#ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygenus_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), focuspathwaysbygenus_withdendroontop, width = 16, height = 8, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygenus_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".pdf", sep=""), focuspathwaysbygenus_withdendroontop, width = 16, height = 8, units = "in", limitsize = FALSE)

ggsave(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygroup_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".pdf", sep=""), focuspathwaysbygroup_withdendroontop, width = 16, height = 8, units = "in", limitsize = FALSE)

png(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_byspecies_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), width = 32, height = 16, units = "in", res = 300)
Sys.sleep(2)
focuspathwaysbyspecies %>% insert_top(d2f, height = 0.2) %>% insert_top(focuspathways_ggtree2, height = 0.2)
Sys.sleep(2)
dev.off()
Sys.sleep(2)

png(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygenus_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), width = 16, height = 8, units = "in", res = 300)
Sys.sleep(2)
focuspathwaysbygenus %>% insert_top(d2f, height = 0.2) %>% insert_top(focuspathways_ggtree2, height = 0.2)
Sys.sleep(2)
dev.off()
Sys.sleep(2)

png(filename = paste("supplemental_plots_taxon_by_pathway/heatmaps_focuspathways_bygroup_and_pathway_clusteredwithdendro_",plot_title,Sys.Date(),".png", sep=""), width = 16, height = 8, units = "in", res = 300)
Sys.sleep(2)
focuspathwaysbygroup %>% insert_top(d2f, height = 0.4) %>% insert_top(focuspathways_ggtree2, height = 0.2)
Sys.sleep(2)
dev.off()
Sys.sleep(2)

####################################################################################################################################
####################################################################################################################################
######################################################################   MAIN PLOTS   ##############################################
####################################################################################################################################
####################################################################################################################################

###############################  plots not including pathways only found in reference set  ##################### 
############################################ FOR PLOTS BY GROUP WITH TARGET & NON-TARGET, NEED TO REVERSE X_SCALE + scale_x_reverse() OR -GROUP?, scale_x_discrete(limits=rev).
## for main plots - switching to mutedGnBu & theme_classic, but not trying to fill in missing...SCRATCH THAT, FOR FACETS KEEPING GnBu & theme_bw
## ENRICHED IN GROUP B PLOTS
ggtestallb <- ggplot(paths_enrichedinBnonas, aes(x=group, y=ec_number)) + scale_x_discrete(limits=rev) + scale_y_discrete(limits=rev)
ggtestallb <- ggtestallb + geom_tile(aes(fill = genepercentage_group_by_ecnumber), color="white", linewidth=0.1) + geom_text(data = subset(paths_enrichedinBnonas, genepercentage_group_by_ecnumber > 0.2), aes(label = genepercentage_group_by_ecnumber), color="red", size = 2)
#ggtestallb <- ggtestallb + scale_fill_viridis_c(name = "Species+Gene \n Percentage of Group", na.value = "transparent")
ggtestallb <- ggtestallb + scale_fill_distiller(name = "Species+Gene \n Percentage of Group", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
ggtestallb <- ggtestallb + theme(axis.ticks=element_blank())
ggtestallb <- ggtestallb + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for", plot_title, "focus pathways missing in target group"))
#ggtestallb <- ggtestallb + theme_bw(base_family="Helvetica")

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  ggtestallb <- ggtestallb + theme_bw(base_family="Helvetica")
} else {
  ggtestallb <- ggtestallb + theme_bw(base_family="sans")
}
#ggtestallb <- ggtestallb + theme_dark(base_family="Helvetica")
ggtestallb <- ggtestallb + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(angle = 45, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#fafcf6'))
ggtestallb <- ggtestallb + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
#ggtestallb

#ggsave("heatmaps_highpriorityrefgenes_enrichedinB.png", ggtestallb, width = 28, height = 28, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_focuspathways_bygroup_without_refECs_",plot_title,Sys.Date(),".png", sep=""), ggtestallb, width = 28, height = 28, units = "in", limitsize = FALSE)
ggsave(filename = paste("heatmaps_focuspathways_bygroup_without_refECs_",plot_title,Sys.Date(),".pdf", sep=""), ggtestallb, width = 28, height = 28, units = "in", limitsize = FALSE)

###########################
### with all species not just group a vs b - note removing labels inside boxes...also making plot a bit wider
## THIS PLOT MAY BE VERY BIG...ONLY SAVE
### NOTE WE ALSO DON'T WANT FREE_X IN FACET
## note this order is by genome size...for phylo ordering see near bottom...
ggtestallb2 <- ggplot(paths_enrichedinBnonas, aes(x=genusspecies, y=ec_number)) + scale_y_discrete(limits=rev)
ggtestallb2 <- ggtestallb2 + geom_tile(aes(fill = genepercentage_group_by_ecnumber), color="white", linewidth=0.1) # + geom_text(aes(label = genepercentage_group_by_ecnumber), color="red", size=2)
ggtestallb2 <- ggtestallb2 + scale_fill_distiller(name = "Species+Gene \n Percentage of Group", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
ggtestallb2 <- ggtestallb2 + theme(axis.ticks=element_blank())
ggtestallb2 <- ggtestallb2 + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for", plot_title, "focus pathways missing in target group"))
#ggtestallb2 <- ggtestallb2 + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  ggtestallb2 <- ggtestallb2 + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)
} else {
  ggtestallb2 <- ggtestallb2 + theme_bw(base_family="sans") + geom_vline(xintercept = total_genomes_A + 0.5)
}
ggtestallb2 <- ggtestallb2 + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(size = 4, angle = 80, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#fafcf6'))
ggtestallb2 <- ggtestallb2 + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
#ggtestallb2

ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_focuspathways_byspecies_without_refECs_",plot_title,Sys.Date(),".png", sep=""), ggtestallb2, width = 32, height = 28, units = "in", limitsize = FALSE)
ggsave(filename = paste("heatmaps_focuspathways_byspecies_without_refECs_",plot_title,Sys.Date(),".pdf", sep=""), ggtestallb2, width = 32, height = 28, units = "in", limitsize = FALSE)

###########################
## everything plots - all pathways (no including PML score in the name)
ggtestallc <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PML, aes(x=group, y=ec_number)) + scale_x_discrete(limits=rev) + scale_y_discrete(limits=rev)
ggtestallc <- ggtestallc + geom_tile(aes(fill = genepercentage_group_by_ecnumber), color="white", linewidth=0.1) + geom_text(data = subset(paths_groupAandB_stats3nonasnofullymissing_PML, genepercentage_group_by_ecnumber > 0.2), aes(label = genepercentage_group_by_ecnumber), color="red", size = 1.5)
ggtestallc <- ggtestallc + scale_fill_distiller(name = "Species+Gene \n Percentage of Group", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
ggtestallc <- ggtestallc + theme(axis.ticks=element_blank())
ggtestallc <- ggtestallc + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for all", plot_title, "pathways"))
#ggtestallc <- ggtestallc + theme_bw(base_family="Helvetica")

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  ggtestallc <- ggtestallc + theme_bw(base_family="Helvetica")
} else {
  ggtestallc <- ggtestallc + theme_bw(base_family="sans")
}
ggtestallc <- ggtestallc + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(size = 6, angle = 45, hjust=1), axis.text.y = element_text(size = 2), panel.background = element_rect(fill = '#fafcf6'))
ggtestallc <- ggtestallc + facet_wrap(~ pathway_name, scales = "free_y", labeller = labeller(pathway_name = label_wrap_gen(40))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
#ggtestallc
ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_allpathways_bygroup_with_refECs_",plot_title,Sys.Date(),".png", sep=""), ggtestallc, width = 28, height = 28, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_allpathways_bygroup_with_refECs_",plot_title,Sys.Date(),".pdf", sep=""), ggtestallc, width = 28, height = 28, units = "in", limitsize = FALSE)

ggtestallc2 <- ggplot(paths_groupAandB_stats3nonasnofullymissing_PML, aes(x=genusspecies, y=ec_number)) + scale_y_discrete(limits=rev)
ggtestallc2 <- ggtestallc2 + geom_tile(aes(fill = genepercentage_group_by_ecnumber), color="white", linewidth=0.1) # + geom_text(aes(label = genepercentage_group_by_ecnumber), color="blue", size=2)
ggtestallc2 <- ggtestallc2 + scale_fill_distiller(name = "Species+Gene \n Percentage of Group", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
ggtestallc2 <- ggtestallc2 + theme(axis.ticks=element_blank())
ggtestallc2 <- ggtestallc2 + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for all", plot_title, "pathways"))
#ggtestallc2 <- ggtestallc2 + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  ggtestallc2 <- ggtestallc2 + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)
} else {
  ggtestallc2 <- ggtestallc2 +theme_bw(base_family="sans") + geom_vline(xintercept = total_genomes_A + 0.5)
}
ggtestallc2 <- ggtestallc2 + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(size = 2, angle = 80, hjust=1), axis.text.y = element_text(size = 2), panel.background = element_rect(fill = '#fafcf6'))
ggtestallc2 <- ggtestallc2 + facet_wrap(~ pathway_name, scales = "free_y", labeller = labeller(pathway_name = label_wrap_gen(40))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
#ggtestallc2
ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_allpathways_byspecies_with_refECs_",plot_title,Sys.Date(),".png", sep=""), ggtestallc2, width = 28, height = 28, units = "in", limitsize = FALSE)
ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_allpathways_byspecies_with_refECs_",plot_title,Sys.Date(),".pdf", sep=""), ggtestallc2, width = 28, height = 28, units = "in", limitsize = FALSE)

## end everything plots

## removing "reversed" plots
######################################################
## adding one genome size plot to each enriched...
## NEED TO CHANGE TO ALWAYS PICK THE FIRST 
#paths_enrichedinB_gsize <- paths_enrichedinB %>% dplyr::filter(pathway_id == paths_enrichedinB$pathway_id[1])
paths_groupAandB_stats3_gsize <- paths_groupAandB_stats3nonasnofullymissing_PML %>% dplyr::filter(pathway_id == paths_groupAandB_stats3nonasnofullymissing_PML$pathway_id[1])
paths_groupAandB_stats3_gsize$pathway_name <- "genome_size"
paths_groupAandB_stats3_gsize$ec_number <- "Mbp"

## also a thinning step to remove multiple genome sizes...
paths_groupAandB_stats3_gsize0 <- paths_groupAandB_stats3_gsize %>%
  group_by(genusspecies, group) %>%
  summarize(across(where(is.numeric), Mode))
# summarize_if(is.numeric, Mode)

paths_groupAandB_stats3_gsize1 <- paths_groupAandB_stats3_gsize %>%
  group_by(genusspecies, group) %>%
  summarize(across(where(is.character), Mode))
  # summarize_if(is.character, Mode)
paths_groupAandB_stats3_gsize <- inner_join(paths_groupAandB_stats3_gsize1, paths_groupAandB_stats3_gsize0)

## adding a dplyr::filter to remove reference set & change ec_number to "Genome Size"
paths_groupAandB_stats3_gsize <- paths_groupAandB_stats3_gsize %>% arrange(desc(group),genome_length) %>% 
  dplyr::filter(genusspecies != "reference_set")
paths_groupAandB_stats3_gsize$pathway_name <- gsub("genome_size","Genome Size",paths_groupAandB_stats3_gsize$pathway_name)


orderinggs <- paths_groupAandB_stats3_gsize[["genusspecies"]]
orderinggs <- unique(orderinggs)
paths_groupAandB_stats3_gsize$genusspecies <- factor(paths_groupAandB_stats3_gsize$genusspecies, levels = orderinggs, ordered = TRUE)

### NOTE NEW IDEA - CHANGE ORDER FROM GENOME SIZE TO PHYLO ORDER, THEN BELOW YOU WILL PLOT PHYLO...

#write.table(paths_groupAandB_stats, "paths_groupAandB_stats.tab", sep = "\t", row.names = FALSE, quote = FALSE)
rm(paths_groupAandB_stats3_gsize0)
rm(paths_groupAandB_stats3_gsize1)


## then plot all as before, but add new plot to end??
## use ggtestall from before!! see new code with grid.arrange
gsize <- ggplot(paths_groupAandB_stats3_gsize, aes(x=genusspecies, y=ec_number)) + scale_y_discrete(limits=rev)
gsize <- gsize + geom_tile(aes(fill = genome_size_mbp), color="white", linewidth=0.1) + geom_text(aes(label = genome_size_mbp), color="red", size=2.5, position=position_jitter(width=0,height=0.2)) #height was 1 but multiple sizes
gsize <- gsize + scale_fill_distiller(name = "Mbp", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
#gsize <- gsize + scale_fill_distiller(name = "Genome Size, Mbp", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
#gsize <- gsize + scale_fill_viridis_c(name = "Genome Size, \nMbp", na.value = "transparent")
gsize <- gsize + theme(axis.ticks=element_blank())
gsize <- gsize + labs(x=NULL, y=NULL, title="")
#gsize <- gsize + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  gsize <- gsize + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)
} else {
  gsize <- gsize + theme_bw(base_family="sans") + geom_vline(xintercept = total_genomes_A + 0.5)
}
size <- gsize + theme(plot.title=element_text(hjust=0), axis.ticks=element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
gsize <- gsize + facet_wrap(~ pathway_name, scales = "free_y") + theme(aspect.ratio = 1)
gsize

layouthkl3 <- rbind(c(1,1,1,NA),
                    c(1,1,1,2),
                    c(1,1,1,NA))


twoplots_enrichedinB <- grid.arrange(ggtestallb2,gsize, layout_matrix = layouthkl3)

print("Main heatmap plots are complete!")
print("The rest of the pipeline is to combine heatmaps with a phylogeny")
Sys.sleep(2)

# ggsave(filename = paste("heatmaps_focuspathways_byspecies_withgenomesize_",plot_title,Sys.Date(),".png", sep=""), twoplots_enrichedinB, width = 32, height = 28, units = "in", limitsize = FALSE)
# ggsave(filename = paste("heatmaps_focuspathways_byspecies_withgenomesize_",plot_title,Sys.Date(),".pdf", sep=""), twoplots_enrichedinB, width = 32, height = 28, units = "in", limitsize = FALSE)
# no longer saving this file, instead saving version with phylo tree - see below

##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

#####################################################################################################################
#####################################################################################################################
### BELOW COMMANDS ONLY IF USER IS MAKING A PHYLOGENY
#####################################################################################################################
#####################################################################################################################


################################################################################################################################
################################################################################################################################
# ### directions:
## First set of code works to print and examine any phylogeny made using the BV-BRC Phylogenetic Tree Tool
## Second set of code will prune the phylogeny to the target and non-target species and allow the user to show the phylogeny along with the heatmaps...

# ##   - First create a genome group that includes all of the species you are interested in - ONLY ONE GENOME PER SPECIES THOUGH!
# ##   - Next go to the Phylogenetic Tree Building (under Services): https://bvbrc.org/app/PhylogeneticTree
# ##   - Add the genome group, select an output name, then submit the job - NOTE THE NUMBER OF SPECIES MUST BE BETWEEN 4-100
# ##   - When run is complete, find the output name and search through the outputs to find the Newick tree file (.nwk)
# 
# # Next run the following R code, at the prompt select the correct nwk file

#### FOR PIC ANALYSES - START HERE RUNNING MANUALLY

## IF YOU ARE IMPORTING A CUSTOM-BUILD TREE USING THE PATRIC TREE TOOL...
# phylo <- read.tree(tk_choose.files(caption = "Select your phylotree (.nwk) made using PATRIC/BV-BRC"))
## using rstudioapi
print("Select your phylotree (.nwk) made using BV-BRC")
Sys.sleep(1)
phylo <- read.tree(rstudioapi::selectFile(caption = "Select your phylotree (.nwk) made using BV-BRC", label = "Select phylotree .nwk", path = data_dir, existing = TRUE, filter = "Newick Files (*.nwk)"))

phylo$tip.label
plot.phylo(phylo)

Sys.sleep(2)
##  also import the same custom datasets's accompanying CSV file...

print("Select your downloaded table (.csv) of the genome group matching the phylotree .nwk file")
Sys.sleep(1)
phylo.data <- read_csv(rstudioapi::selectFile(caption = "Select your downloaded table (.csv) of the genome group matching the phylotree .nwk file", label = "Select matching .csv table", path = data_dir, existing = TRUE, filter = "CSV Files (*.csv)")
                       , col_types = list(
                         "Genome ID" = col_character(), "Genome Name" = col_character(), "NCBI Taxon ID" = col_character()))


phylo.data <- phylo.data %>%
  rename(genome_id = "Genome ID") %>%
  rename(genome_name = "Genome Name") %>%
  rename(genome_length = Size) %>%
  rename(taxon_id = "NCBI Taxon ID") %>%
  select(genome_id, genome_name, Species, taxon_id, genome_length)

## some correcting (not thinning)
phylo.data$genome_name <- gsub("Candidatus ","",phylo.data$genome_name)
phylo.data$genome_name <- gsub("uncultured ","",phylo.data$genome_name)
phylo.data$genome_name <- gsub("sp. ","sp_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("strain ","strain_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("OMZ ","OMZ_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("ATCC ","ATCC_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("PCC ","PCC_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("FDAARGOS ","FDAARGOS_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("NCTC ","NCTC_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("MAG ","MAG_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("USDA ","USDA_",phylo.data$genome_name)

phylo.data$genome_name <- gsub("[","",phylo.data$genome_name, fixed = TRUE)
phylo.data$genome_name <- gsub("]","",phylo.data$genome_name, fixed = TRUE)
phylo.data$genome_name <- gsub('"',"",phylo.data$genome_name, fixed = TRUE)
phylo.data$genome_name <- gsub("'","",phylo.data$genome_name, fixed = TRUE)

phylo.data$genome_name <- gsub("endosymbiont wPip_Mol of ","endo_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endosymbiont of ","endo_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endosymbiont strain_TRS of ","endo_",phylo.data$genome_name)

# ## universal version
# phylo.data$genome_name <- gsub(" ", "_", gsub("endo_(\\S+)", "endo_\\1 ", phylo.data$genome_name))
# phylo.data$genome_name <- gsub("__","_",phylo.data$genome_name)
# # In this code, gsub() function is used twice. The inner gsub() is used to find matches of the pattern "endo_(\S+)" in the 'colla' column of the 'dabba' dataframe. The (\\S+) captures all non-whitespace characters after "endo_" in the match. The outer gsub() is then used to replace the space character with an underscore "" in the matched patterns, using the captured group \\1 to retain the original characters after "endo" and include the space character. The modified 'colla' column is updated in the 'dabba' dataframe. Finally, the modified 'colla' column is printed using print(dabba$colla).

phylo.data$genome_name <- gsub("endo_Brugia ","endo_Brugia_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Rhopalodia ","endo_Rhopalodia_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Epithemia ","endo_Epithemia_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Onchocerca ","endo_Onchocerca_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Drosophila ","endo_Drosophila_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Culex ","endo_Culex_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Aedes ","endo_Aedes_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Amblyomma ","endo_Amblyomma_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Donacia ","endo_Donacia_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Macroplea ","endo_Macroplea_",phylo.data$genome_name)
phylo.data$genome_name <- gsub("endo_Plateumaris ","endo_Plateumaris_",phylo.data$genome_name)

## repeat for Species field
phylo.data$Species <- gsub("Candidatus ","",phylo.data$Species)
phylo.data$Species <- gsub("uncultured ","",phylo.data$Species)
phylo.data$Species <- gsub("sp. ","sp_",phylo.data$Species)
phylo.data$Species <- gsub("strain ","strain_",phylo.data$Species)
phylo.data$Species <- gsub("OMZ ","OMZ_",phylo.data$Species)
phylo.data$Species <- gsub("ATCC ","ATCC_",phylo.data$Species)
phylo.data$Species <- gsub("PCC ","PCC_",phylo.data$Species)
phylo.data$Species <- gsub("FDAARGOS ","FDAARGOS_",phylo.data$Species)
phylo.data$Species <- gsub("NCTC ","NCTC_",phylo.data$Species)
phylo.data$Species <- gsub("MAG ","MAG_",phylo.data$Species)
phylo.data$Species <- gsub("USDA ","USDA_",phylo.data$Species)
## also just in Species field change XXX sp. at end of field to null, to force use of the genome name - specify the metacharacter $ after the sp_ to signify the end of the string
phylo.data$Species <- gsub(".*sp.$","null",phylo.data$Species)

phylo.data$Species <- gsub("[","",phylo.data$Species, fixed = TRUE)
phylo.data$Species <- gsub("]","",phylo.data$Species, fixed = TRUE)
phylo.data$Species <- gsub('"',"",phylo.data$Species, fixed = TRUE)
phylo.data$Species <- gsub("'","",phylo.data$Species, fixed = TRUE)


phylo.data$Species <- gsub("endosymbiont wPip_Mol of ","endo_",phylo.data$Species)
phylo.data$Species <- gsub("endosymbiont of ","endo_",phylo.data$Species)
phylo.data$Species <- gsub("endosymbiont strain_TRS of ","endo_",phylo.data$Species)

# ## universal version
# phylo.data$Species <- gsub(" ", "_", gsub("endo_(\\S+)", "endo_\\1 ", phylo.data$Species))
# phylo.data$Species <- gsub("__","_",phylo.data$Species)
# # In this code, gsub() function is used twice. The inner gsub() is used to find matches of the pattern "endo_(\S+)" in the 'colla' column of the 'dabba' dataframe. The (\\S+) captures all non-whitespace characters after "endo_" in the match. The outer gsub() is then used to replace the space character with an underscore "" in the matched patterns, using the captured group \\1 to retain the original characters after "endo" and include the space character. The modified 'colla' column is updated in the 'dabba' dataframe. Finally, the modified 'colla' column is printed using print(dabba$colla).

phylo.data$Species <- gsub("endo_Brugia ","endo_Brugia_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Rhopalodia ","endo_Rhopalodia_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Epithemia ","endo_Epithemia_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Onchocerca ","endo_Onchocerca_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Drosophila ","endo_Drosophila_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Culex ","endo_Culex_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Aedes ","endo_Aedes_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Amblyomma ","endo_Amblyomma_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Donacia ","endo_Donacia_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Macroplea ","endo_Macroplea_",phylo.data$Species)
phylo.data$Species <- gsub("endo_Plateumaris ","endo_Plateumaris_",phylo.data$Species)

phylo.data <- phylo.data %>% mutate(Species = na_if(Species, "null"))
#phylo.data$Species <- gsub("null ","",phylo.data$Species)
## instead of gsub, use coalesce??
# dfABy %>% mutate(A = coalesce(A,B)) this is replacing all NAs in cregion with the value from locus
phylo.data <- phylo.data %>% mutate(Species = coalesce(Species,genome_name))

## then pull just genus & species
## create new column with shorter description
phylo.data <- phylo.data %>% separate(genome_name, into = c("genus", "species"), sep = " ", remove = FALSE, convert = TRUE, extra = "drop")
phylo.data <- phylo.data %>% unite("genusspecies", genus, species, sep = "_", remove = TRUE)

## then pull just genus & species
## create new column with shorter description
phylo.data <- phylo.data %>% separate(Species, into = c("genus", "species"), sep = " ", remove = FALSE, convert = TRUE, extra = "drop")
phylo.data <- phylo.data %>% unite("genusspecies2", genus, species, sep = "_", remove = TRUE) %>% relocate(genusspecies2, .after = genusspecies)

## now remove if genusspecies2 doesn't = genusspecies update do not remove, but use the latter genusspecies2 in every case!
## further down we will replace the collated paths_groupAandB$genusspecies with genusspecies2
phylo.data_mismatch <- subset(phylo.data , genusspecies != genusspecies2)

phylo.data <- phylo.data %>% 
  rename(genusspecies1 = genusspecies) %>%
  rename(genusspecies = genusspecies2) %>%
  relocate(genusspecies, .before = genusspecies1)
## replace NAs...
phylo.data <- phylo.data %>% mutate(genusspecies = coalesce(genusspecies,genusspecies1))


phylotipslist <- phylo$tip.label
phylotipsdf <- as.data.frame(phylotipslist, col.names = "genome_id")
phylotipsdf <- rename(phylotipsdf, genome_id = phylotipslist)
phylotips_withname_guidetree <- left_join(phylotipsdf, phylo.data) ## need full list, not curated either, but genome_summary_phylo (but still some missing, see plan B)

phylotips_withname_guidetree$genome_size_mbp <- phylotips_withname_guidetree$genome_length / 1000000
phylotips_withname_guidetree$genome_size_mbp <- round(phylotips_withname_guidetree$genome_size_mbp, 2)

### ADD A CHECK HERE TO MAKE SURE NO NAS - IF SO NEED TO GIVE USER AN ERROR...
if (nrow(phylotips_withname_guidetree %>% dplyr::filter(is.na(.[["taxon_id"]])|is.na(.[["genome_size_mbp"]]))) == 0) {
  print("Your BV-BRC phylogeny has no missing genome values!")
} else {
  print("Your BV-BRC phylogeny has some missing genome values:")
  print(phylotips_withname_guidetree %>% dplyr::filter(is.na(.[["taxon_id"]])|is.na(.[["genome_size_mbp"]])))
  stop()
}
Sys.sleep(1)

phylotips_withname_guidetree2 <- phylotips_withname_guidetree %>% unite("genome_name_length", genusspecies, genome_size_mbp, sep = " ", remove = FALSE)
## adding gsub to remove ':', also check to make sure no duplicates
phylotips_withname_guidetree2$genusspecies <- gsub(":","-",phylotips_withname_guidetree2$genusspecies)
phylotips_withname_guidetree2$genome_name_length <- gsub(":","-",phylotips_withname_guidetree2$genome_name_length)
phylotips_withname_guidetree2$genome_name <- gsub(":","-",phylotips_withname_guidetree2$genome_name)

# phylotips_withname_guidetree2 %>% dplyr::filter(genusspecies %in% unique(.[["genusspecies"]][duplicated(.[["genusspecies"]])]))

if (nrow(phylotips_withname_guidetree2 %>% dplyr::filter(genusspecies %in% unique(.[["genusspecies"]][duplicated(.[["genusspecies"]])]))) == 0) {
  print("Your BV-BRC phylogeny has no duplicate taxon names!")
} else {
  print("Your BV-BRC phylogeny may still some duplicate taxon names:")
  print(phylotips_withname_guidetree2 %>% dplyr::filter(genusspecies %in% unique(.[["genusspecies"]][duplicated(.[["genusspecies"]])])))
}
Sys.sleep(1)

## IF YOU HAVE DUPLICATES, REMOVE HERE...
phylotips_withname_guidetree2 <- phylotips_withname_guidetree2 %>%
  dplyr::filter(genusspecies1 != "Rickettsia_heilongjiangensis") %>%
  dplyr::filter(genusspecies1 != "Rickettsia_raoultii") %>%
  dplyr::filter(genusspecies1 != "Rickettsia_amblyommii")
phylotips_withname_guidetree2 <- phylotips_withname_guidetree2 %>%
  dplyr::filter(genusspecies1 != "Strawberry_lethal") %>%
  dplyr::filter(genusspecies1 != "Mycoplasma_ovipneumoniae")
phylotips_withname_guidetree2 <- phylotips_withname_guidetree2 %>%
  dplyr::filter(genome_id != "1660071.3") %>%
  dplyr::filter(genome_id != "1660070.3") %>%
  dplyr::filter(genome_id != "1244528.3") %>%
  dplyr::filter(genome_id != "1244531.3") %>%
  dplyr::filter(genome_id != "1660063.4")
phylotips_withname_guidetree2 <- phylotips_withname_guidetree2 %>%
  dplyr::filter(genome_id != "1458985.3")

#and the tree!
phylo <- ape::drop.tip(phylo, c("226665.5", "369822.3", "1105111.3"))
phylo <- ape::drop.tip(phylo, c("1795874.3", "1439853.3"))
phylo <- ape::drop.tip(phylo, c("980422.3", "29562.18"))
phylo <- ape::drop.tip(phylo, c("1660071.3", "1660070.3", "1244528.3","1660063.4","1244531.3"))
phylo <- ape::drop.tip(phylo, c("1458985.3"))

## new tack - if not true above, take phylotips_withname_guidetree2 and group by genusspecies, then take only first...
## issue is we lose the order in guidetree22...need to save original order as new column, then reorder by this
phylotips_withname_guidetree2 <- phylotips_withname_guidetree2 %>% dplyr::mutate(row = row_number())
phylotips_withname_guidetree22 <- phylotips_withname_guidetree2 %>% 
  group_by(genusspecies) %>% 
  summarize(across(everything(), first)) %>% ungroup() %>% relocate(genusspecies, .before = genusspecies1) %>% arrange(row)
## BUT NEED TO NOTE THE GENOME_ID OF THOSE DROPPED AND ALSO DROP THEM FROM PHYLO!!
phylotips_withname_guidetree_todrop <- anti_join(phylotips_withname_guidetree2,phylotips_withname_guidetree22)
phylo_todrop <- phylotips_withname_guidetree_todrop$genome_id
###
## only drop tips if there are tips to drop... check phylo$tip.label[1], maybe root by this?
## issue is we lose the order in guidetree22...need to save original order as new column, then reorder by this
if (n_distinct(phylotips_withname_guidetree_todrop$genome_id) > 0) {
  phylo <- ape::drop.tip(phylo, phylo_todrop)
  rm(phylotips_withname_guidetree2)
  phylotips_withname_guidetree2 <- phylotips_withname_guidetree22
  rm(phylotips_withname_guidetree22)
  phylo <- ape::root(phylo, 1)
  print("Your BV-BRC phylogeny had some duplicate taxon names that have been removed.")
  print("Carefully check your resulting phylogeny to make sure the topology is correct.")
  print("If there are changes, we recommend removing the duplicates in BV-BRC and re-calculating your phylogeny.")
  print("Then re-run the pipeline")
} else {
}
Sys.sleep(3)

## check again for duplicates
if (nrow(phylotips_withname_guidetree2 %>% dplyr::filter(genusspecies %in% unique(.[["genusspecies"]][duplicated(.[["genusspecies"]])]))) == 0) {
  print("Your BV-BRC phylogeny now has no duplicate taxon names!")
} else {
  print("Your BV-BRC phylogeny still has some duplicate taxon names:")
  print(phylotips_withname_guidetree2 %>% dplyr::filter(genusspecies %in% unique(.[["genusspecies"]][duplicated(.[["genusspecies"]])])))
  print("You will have to manually add two commands (one with dplyr::filter & one with ape::drop.tip) to filter out the genome_id numbers of all duplicates!")
  print("Look for the code snippets after 'IF YOU HAVE DUPLICATES, REMOVE HERE' and replicate")
  print("Then re-run the pipeline")
  stop()
}
Sys.sleep(1)

### printing...
species <-  phylotips_withname_guidetree2$genome_name_length
# plot_title <- unlist(strsplit(species[14], "_"))[1]
# print("We've guessed that your taxon plot title should be")
# print(plot_title)
# print("Please rename if this is a bad guess!")

write.table(phylotips_withname_guidetree2, file = paste("supplemental_plots_ec_by_taxon_per_pathway/phylogeny_focus_",plot_title,Sys.Date(),".tab", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

## also save version with just genus+species for pulling A+B species names below...
phylo0 <- phylo

### check to stop if tip labels are not equal (will crash R if not)
#n_distinct(phylo$tip.label) == n_distinct(phylotips_withname_guidetree2$genusspecies)
## check if true!
if ((n_distinct(phylo$tip.label) == n_distinct(phylotips_withname_guidetree2$genusspecies)) == TRUE) {
  phylo0$tip.label <- phylotips_withname_guidetree2$genusspecies
  phylo$tip.label <- phylotips_withname_guidetree2$genome_name_length
  phylo$tip.label[1]
  plot.phylo(phylo)
  phylo0$tip.label[1]
  Sys.sleep(2)
} else {
  print("Your phylotree (.nwk) file and your downloaded table (.csv) of the genome group do not have the same number of taxa.")
  print("Please check these, re-run in BV-BRC if necessary, then re-run the pipeline")
  stop()
}

# Margins area
par(oma=c(1,2,1,2)) # left & right have 3 lines of space
# par(mar=c(2,4,2,4) + 0.1)
# 
# par(oma=c(1,5,1,5)) # left & right have 3 lines of space
# par(mar=c(2,4,2,6)) midpoint(tree)
# plot.phylo(phylo, edge.width=2, cex = 0.5, align.tip.label = TRUE, adj = 1, label.offset = -0.8)
# plot.phylo(phylo, edge.width=2, cex = 0.5, align.tip.label = TRUE, adj = 0, label.offset = 0)
plot.phylo(midpoint(phylo), edge.width=1.4, cex = 0.6, align.tip.label = TRUE, adj = 0, label.offset = 0.01, no.margin = FALSE)

## outer margins seem to help more!!
# ## bottom, left, top, right
par(oma=c(1,4,1,5))

# par(“mai”)/par(“mar”)
# [1] 0.2 0.2 0.2 0.2

### saving a plot...
png(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/phylogeny_focus_",plot_title,Sys.Date(),".png", sep=""), width = 16, height = 24, units = "in", res = 300)
# ggsave(filename = paste("phylogeny_thintips_withlabels8x8_",plot_title,Sys.Date(),".png", sep=""), dendro_plot2c, width = 8, height = 8, units = "in", limitsize = FALSE)
par(oma=c(1,4,1,5))

if (Ntip(phylo) < 500) {
  # plot.phylo(phylo, edge.width=1, cex = 2.2, align.tip.label = TRUE, adj = 0, label.offset = 0, no.margin = FALSE)
  plot.phylo(midpoint(phylo), edge.width=1.4, cex = 0.6, align.tip.label = TRUE, adj = 0, label.offset = 0.01, no.margin = FALSE)
} else {
  print("Your phylogeny has more than 500 tips, we recommend thinning!")
  plot.phylo(midpoint(phylo), edge.width=1.4, cex = 0.6, align.tip.label = TRUE, adj = 0, label.offset = 0.01, no.margin = FALSE)
}
dev.off()

Sys.sleep(2)
##########
## FOR PIC ANALYSES - use midpoint rooted phylo!
write.tree(midpoint(phylo), file = paste("supplemental_plots_ec_by_taxon_per_pathway/PATRICphylotree_midpointrooted_withlabels_",plot_title,Sys.Date(),".nwk", sep=""))

#### END PIC ANALYSIS CODE
################################################################################################################################
################################################################################################################################

################################################################################################################################
################################################################################################################################
### making phylogeny for ggdendrogram

tree.unmatched <- multi2di(phylo, random=FALSE) 
#Let's plot our new tree.  Huzzah!
#plotTree(tree.unmatched,fsize=0.6) 
plot.phylo(tree.unmatched)

tree.unmatched$tip.label[1]

is.binary(tree.unmatched)

# Mode(phylotips_withname$genome_name)
# unlist(strsplit(phylotips_withname$genome_name[200], " "))[1]
#First, we’ll make the tree ultrametric using the chronos function in APE.
phyloguide.pruned_um <- chronos(tree.unmatched, lambda=0)

phyloguide.pruned_um
phyloguide.pruned_d <- as.dendrogram(as.hclust.phylo(phyloguide.pruned_um))
is.ultrametric(phyloguide.pruned_um)

dendro_plot <- ggdendrogram(data = phyloguide.pruned_d, rotate = FALSE)
#dendro_plot

dendro_plot2 <- ggdendrogram(data = phyloguide.pruned_d, rotate = TRUE)

################################################################################################################################
################################################################################################################################
### to pull just species in your groups A & B from previously loaded tree (including check to make sure no tips are missing!)

## to make phylogeny of your species in groups A & B...
n_distinct(paths_enrichedinAnonas$genusspecies)
n_distinct(paths_enrichedinBnonas$genusspecies)

unique(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genusspecies)
n_distinct(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genusspecies)

n_distinct(paths_enrichedinAnonas$genusspecies) == n_distinct(paths_enrichedinBnonas$genusspecies)
#n_distinct(paths_enrichedinAnonas$genusspecies) == n_distinct(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genusspecies)
## check if true!
if ((n_distinct(paths_enrichedinAnonas$genusspecies) == n_distinct(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genusspecies)) == TRUE) {
  print("You are ready to make more trees!")
  Sys.sleep(2)
} else {
  print("You are missing taxa!")
  stop()
}

## if true, then
species.aplusb <- unique(paths_groupAandB_stats3nonasnofullymissing_PMLnoref$genusspecies)

## new tack - iteratively try rooting first with midpoint rooting, then by non-target, then by target group
## this way if at all possible the tree will have a root that splits target & non-target species
## first try midpoint rooting 
phylo0_mp <- phytools::midpoint.root(phylo0)
plot.phylo(phylo0_mp)
Sys.sleep(1)
## note that midpoint rooting usually but not always results in target & non-target groups being monophyletic.

## secondly try rooting by non-target group
species.bonly <- paths_groupAandB_stats %>% dplyr::filter(group == "group B") %>% 
  select(genusspecies) %>% unique()
species.bonly <- species.bonly[["genusspecies"]]


#mrca.b <- ape::getMRCA(phylo0, species.bonly)
# testrootb <- try(ape::root(phylo0, node = mrca.b, resolve.root = TRUE), silent = TRUE)
# phylo0_mp <- ape::root(phylo0, node = mrca.b, resolve.root = TRUE)
# plot.phylo(phylo0_mp)
# Sys.sleep(1)

## thirdly try rooting by target group
species.aonly <- paths_groupAandB_stats %>% dplyr::filter(group == "group A") %>% 
  select(genusspecies) %>% unique()
species.aonly <- species.aonly[["genusspecies"]]

## only make mrca.b if 3 or more species (otherwise will error out)
if (n_distinct(species.aonly) < 3)  {
  Sys.sleep(2)
} else {
  mrca.b <- ape::getMRCA(phylo0, species.bonly)
}

## only make mrca.a if 3 or more species (otherwise will error out)
if (n_distinct(species.bonly) < 3)  {
  Sys.sleep(2)
} else {
  mrca.a <- ape::getMRCA(phylo0, species.aonly)
}

#mrca.a <- ape::getMRCA(phylo0, species.aonly)
## try to root with MRCA of group A...
# testroota <- try(ape::root(phylo0, node = mrca.a, resolve.root = TRUE), silent = TRUE)
# phylo0_mp <- ape::root(phylo0, node = mrca.a, resolve.root = TRUE)
# plot.phylo(phylo0_mp)
# Sys.sleep(1)

## now if then statement depending on whether these rootings result in errors or not
testrootb <- try(ape::root(phylo0, node = mrca.b, resolve.root = TRUE), silent = TRUE)
testroota <- try(ape::root(phylo0, node = mrca.a, resolve.root = TRUE), silent = TRUE)

if ((n_distinct(species.aonly) < 2) & (n_distinct(species.bonly) < 2)) {
  print("Your target and non-target group are both a single species, you cannot make a phylogeny!")
  stop()
} else {
  Sys.sleep(2)
}


if (class(testrootb) == "phylo") {
  print("Phylogeny is rooted by the non-target group.")
  phylo0_mp <- ape::root(phylo0, node = mrca.b, resolve.root = TRUE)
} else if (class(testroota) == "phylo") {
  print("Phylogeny is rooted by the target group.")
  phylo0_mp <- ape::root(phylo0, node = mrca.a, resolve.root = TRUE)
} else if (n_distinct(species.aonly) < 2) {
  print("Phylogeny is rooted by the target group, as it contains a single species.")
  phylo0_mp <- ape::root(phylo0, outgroup = species.aonly, resolve.root = TRUE)
} else if (n_distinct(species.bonly) < 2) {
  print("Phylogeny is rooted by the non-target group, as it contains a single species.")
  phylo0_mp <- ape::root(phylo0, outgroup = species.bonly, resolve.root = TRUE)
} else {
  print("We were unable to root your phylogeny by either the target group or the non-target group.")
  print("As a result, we are using midpoint rooting.")
  print("Check to see if your target & non-target groups are monophyletic!")
  phylo0_mp <- phytools::midpoint.root(phylo0)
}

## now prune to just the species in the heatmaps
phylo.pruned <- ape::drop.tip(phylo0_mp, setdiff(phylo0$tip.label, species.aplusb));
plot.phylo(phylo.pruned)
## phylo.pruned is needed to make the ordering for final phylo plots!


### CODE TO AUTOMATICALLY RE-ALIGN TREE SO THAT CORRECT TAXON IS FIRST (OR LAST)
## first try by manually identifying nodes...fails
# print("Click on the phylogeny the ancestral node of this taxon:")
# print(ordering[1])
# pickedclade1 <- identify(phylo.pruned)
# 
# print("Now click on the phylogeny the ancestral node of the sister clade that needs to be swapped:")
# print(ordering[1])
# pickedclade2 <- identify(phylo.pruned)
# 
# print("Finally click on the phylogeny the ancestor of both of the previous two nodes:")
# print(ordering[1])
# pickedroot <- identify(phylo.pruned)
# 
# exchangedclades <- c(pickedclade1,pickedclade2)
# # ape::identify.phylo(phylo.pruned, tips = TRUE, labels = TRUE)
# # identify(phylo.pruned, tips = TRUE, labels = TRUE)
# 
# phylo.prtest <- ape::rotate(phylo.pruned, pickedroot, exchangedclades)
# ## then rotate nodes based on this node?? phytools::rotateNodes or better ape::rotate - look into this...using constraint option!
# plot.phylo(phylo.prtest)

## this didn't work...new idea is to take the list of phylo.pruned$tip.label & make a new list that starts with ordering[1], then wraps to get everything
## like this https://stackoverflow.com/questions/47292561/reorder-a-vector-with-wrap-around-in-r
#### that would be the correct order for the constraint like in the example: co <- tr2$tip.label[order(match(tr2$tip.label, tr1$tip.label))]. with:
# rotateConstr(phy, constraint)

clade_name0 <- phylo.pruned$tip.label
shiftname <- ordering[1]
#shiftnum <- which.edge(phylo.pruned, shiftname) doesn't work
## this works: shiftnum <- 3 # this works, which.edge is wrong

# https://stackoverflow.com/questions/57090663/how-to-find-word-index-or-position-in-a-given-string-using-r-programming
# Output_text <- c("applicable to any future potential contract termination disputes as the tepco dispute was somewhat unique")
# words <- unlist(str_split(Output_text, " "))
# which(words == "termination")
shiftnum <- which(clade_name0 == shiftname)
clade_name1 <-c(clade_name0[shiftnum:length(clade_name0)],clade_name0[1:shiftnum-1])

#c(clade_name0[shiftnum:length(clade_name0)],clade_name0[1:shiftnum-1])
#or relevel(x, ref, ...) - needs to be a factor, adds strangeness clade_name1 <- relevel(clade_name0, shiftnum)
# clade_name1 <-c(clade_name0[shiftnum:length(clade_name0)],clade_name0[1:shiftnum-1])

phylo.prunedt <- rotateConstr(phylo.pruned, clade_name1)
plot.phylo(phylo.prunedt)

## not sure why, but new phylo plots okay but only is in correct order if we export then re-import?
write.tree(phylo.prunedt, file = "phylogeny_realigned.nwk")
phylo.pruned <- read.tree("phylogeny_realigned.nwk")
unlink("phylogeny_realigned.nwk")
rm(clade_name0)
rm(clade_name1)
rm(shiftnum)
rm(shiftname)
rm(phylo.prunedt)


####################################################################################################################
## IF YOU NEED TO RE-ALIGN THE TIPS SO A TARGET GROUP TAXON IS FIRST, DO THAT IN GENEIOUS, THEN RE-IMPORT - note automatic fix just above!
## for last - sum(table(phylo.pruned$tip.label))     phylo.pruned$tip.label[sum(table(phylo.pruned$tip.label))]
##. ordering[sum(table(ordering)) ]
if ((phylo.pruned$tip.label[1] == ordering[1]) == TRUE) {
  print("Your pruned phylogeny has a target group taxon first, great!!")
} else if ((phylo.pruned$tip.label[sum(table(phylo.pruned$tip.label))] == ordering[1]) == TRUE) {
  print("Your pruned phylogeny has a target group taxon last, great!!")
} else {
  print("Your pruned phylogeny does not have a target group taxon either first or last.")
  print("As a result, heatmaps + phylo will not parse well.")
  print("Consider re-aligning the 'PATRICphylotree_pruned_withoutlabels_' .nwk file so that")
  print(ordering[1])
  print(" is either the first, or the last taxon!")
  write.tree(phylo.pruned, file = paste("PATRICphylotree_pruned_withoutlabels_",plot_title,Sys.Date(),".nwk", sep=""))
  print("After re-aligning that phylogeny, you will be prompted to re-input it.")
  readline(prompt="Press [enter] to continue")
  # phylo.pruned <- read.tree(tk_choose.files(caption = "Select your re-rooted phylotree made from PATRIC/BV-BRC"))
  phylo.pruned <- read.tree(rstudioapi::selectFile(caption = "Select your re-rooted phylotree made using BV-BRC", label = "Select re-rooted phylotree .nwk", path = data_dir, existing = TRUE, filter = "Newick Files (*.nwk)"))
  phylo.pruned$tip.label <- gsub("'","",phylo.pruned$tip.label, fixed = TRUE)
  Sys.sleep(3)
#  stop()
}

####################################################################################################################

## moving this to below to change tree order depending on whether phylogeny has target species first or last...
# if ((phylo.pruned$tip.label[1] == ordering[1]) == TRUE) {
#   print("Your phylogeny has a target group taxon first, great!!")
# } else {
#   print("Your phylogeny still does not have a target group taxon first, your heatmaps + phylo will not parse well.")
#   stop()
# }

## phylo.pruned is needed to make the ordering for final phylo plots!

if ((phylo.pruned$tip.label[1] == ordering[1]) == TRUE) {
  print("Your pruned phylogeny has a target group taxon first, great!!")
  tree.unmatched <- multi2di(phylo.pruned, random=FALSE) 
  phylo.pruned_um <- chronos(tree.unmatched, lambda=0)
  phylo.pruned_d <- as.dendrogram(as.hclust.phylo(phylo.pruned_um))
  clade_order <- order.dendrogram(phylo.pruned_d)
  clade_name <- labels(phylo.pruned_d)
  dendro_plotref1 <- ggdendrogram(data = phylo.pruned_d, rotate = FALSE)  
} else if ((phylo.pruned$tip.label[sum(table(phylo.pruned$tip.label))] == ordering[1]) == TRUE) {
  print("Your pruned phylogeny has a target group taxon last, great!!")
  tree.unmatched <- multi2di(phylo.pruned, random=FALSE) 
  phylo.pruned_um <- chronos(tree.unmatched, lambda=0)
  phylo.pruned_d <- as.dendrogram(as.hclust.phylo(phylo.pruned_um))
  clade_order <- order.dendrogram(phylo.pruned_d)
  ## if your target taxon is last...
  clade_name <-  rev(labels(phylo.pruned_d))
  dendro_plotref1 <- ggdendrogram(data = rev(phylo.pruned_d), rotate = FALSE)
  } else {
  print("Your pruned phylogeny still does not have a target group taxon first & heatmaps + phylo will not parse well.")
  stop()
}

## below commands moved into if else depending on order...
### back to regular commands

clade_position <- data.frame(clade_name,
                             clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
ordering_byphylo <- clade_name

dendro_plotref_flipped1 <- print(dendro_plotref1) + scale_y_reverse()
#plot.phylo(rev(phylo.pruned_um), direction = "upwards")

## for heatmaps with phylo...
#print(phylo.pruned_d)
#dendro_plotref_flipped <- print(dendro_plotref) + scale_y_reverse()
#plot.phylo(phylo.pruned_um, direction = "upwards")

#dendro_plotref <- ggdendrogram(data = rev(phylo.pruned_gs_d), rotate = FALSE)
## no dendro_plotref_forphyloheatmap <- print(dendro_plotref) + scale_x_reverse()

## try adding margins - here is without for pdf version
dendro_plotref_forphyloheatmap <- print(dendro_plotref1) + scale_y_reverse() + theme(legend.position="left", legend.key.width=unit(0.2,"cm"))

## too much margin, trying less
#dendro_plotref_forphyloheatmap <- print(dendro_plotref1) + scale_y_reverse() + theme(legend.position="left", legend.key.width=unit(0.2,"cm"), plot.margin = margin(0, 3.4, 0, 4.3, "cm"))
dendro_plotref_forphyloheatmapwithmargin <- print(dendro_plotref1) + scale_y_reverse() + theme(legend.position="left", legend.key.width=unit(0.2,"cm"), plot.margin = margin(0, 2.9, 0, 3.8, "cm"))
#dendro_plotref_forphyloheatmapwithbigmargin <- print(dendro_plotref1) + scale_y_reverse() + theme(legend.position="left", legend.key.width=unit(0.2,"cm"), plot.margin = margin(0, 3.4, 0, 4.8, "cm"))


## this instead? or as.phylo?? yikes changes a lot, do not use
# dendroplotref_ggtree <- ggtree(phylo.pruned_d) + geom_tiplab()
# dendroplotref_ggtree2 <- ggtree(phylo.pruned_d) + coord_flip()

# allpathways_dendronolabels <- ggdendrogram(data = hc.cols22, rotate = FALSE, theme_dendro = TRUE, leaf_labels = FALSE, labels = FALSE) ## dendrogram faces downward
# allpathways_dendronolabels1 <- allpathways_dendronolabels + theme(axis.title = element_text(size = 0), axis.text = element_text(size = 0))
# allpathways_dendronolabels2 <- allpathways_dendronolabels + theme(plot.margin = margin(0, 1.9, 0, 2.8, "cm"), axis.title = element_text(size = 0), axis.text = element_text(size = 0))


##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

#########################
## plots with phylogenies included

## NOW REPEATING LATEST PLOT COMMANDS BUT ORDERED FOR PHYLOGENY...
## new files...
############ ADDING FACTOR OF LATEST PML FILE WITH GENOME SIZE ORDER ##########

paths_enrichedinBnonas_plotswithphylo <- paths_enrichedinBnonas %>% ungroup()
unique(paths_enrichedinBnonas_plotswithphylo$genusspecies)

################################################
## CHECK TO MAKE SURE ORDERING & ORDERING_BYPHYLO ARE THE SAME
## note when removing reference set from "ordering" do not subtract 1
# sum(table(ordering)) - 1
# sum(table(ordering_byphylo))

if (((sum(table(ordering)) - 1) == sum(table(ordering_byphylo))) == TRUE) {
# if ((sum(table(ordering)) == sum(table(ordering_byphylo))) == TRUE) {
  print("Your pruned phylogeny has the same number of taxa as your heatmaps, great!!")
} else {
  print("Your pruned phylogeny is missing some taxa in your heatmaps. You will need to re-run your phylogeny in BV-BRC with the missing taxa")
  # print("Your heatmaps contain")
  # print(sum(table(ordering)) - 1)
  # print("taxa")
  print(paste0("Your heatmaps contain ", sum(table(ordering)) - 1, " taxa"))
  # print(paste0("Your heatmaps contain ", sum(table(ordering)), " taxa"))
  print(paste0("While your pruned phylogeny only contains ", sum(table(ordering_byphylo)), " taxa"))
  # print("While your phylogeny only contains")
  # print(sum(table(ordering_byphylo)))
  # print("taxa")
  print("You will need to re-run your pruned phylogeny in BV-BRC with the missing taxa, then re-run the pipeline")
  stop()
}
################################################

paths_enrichedinBnonas_plotswithphylo$genusspecies <- factor(paths_enrichedinBnonas_plotswithphylo$genusspecies, levels = ordering_byphylo, ordered = TRUE)
unique(paths_enrichedinBnonas_plotswithphylo$genusspecies)

## do not want gsize, that was a small dataset, use paths_enrichedinBnonas
paths_groupAandB_stats3_gsize_plotswithphylo <- paths_groupAandB_stats3_gsize %>% ungroup()
paths_groupAandB_stats3_gsize_plotswithphylo$genusspecies <- factor(paths_groupAandB_stats3_gsize_plotswithphylo$genusspecies, levels = ordering_byphylo, ordered = TRUE)
unique(paths_groupAandB_stats3_gsize_plotswithphylo$genusspecies)

## now change order to ordering_by_phylo
ggtestallb2withphylo <- ggplot(paths_enrichedinBnonas_plotswithphylo, aes(x=genusspecies, y=ec_number)) + scale_y_discrete(limits=rev)
ggtestallb2withphylo <- ggtestallb2withphylo + geom_tile(aes(fill = genepercentage_group_by_ecnumber), color="white", linewidth=0.1) # + geom_text(aes(label = genepercentage_group_by_ecnumber), color="red", size=2)
ggtestallb2withphylo <- ggtestallb2withphylo + scale_fill_distiller(name = "Species+Gene \n Percentage of Group", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
ggtestallb2withphylo <- ggtestallb2withphylo + theme(axis.ticks=element_blank())
ggtestallb2withphylo <- ggtestallb2withphylo + labs(x=NULL, y=NULL, title=paste("Species+Gene percentages in target group vs. non-target group \n     for", plot_title, "focus pathways missing in target group"))
#ggtestallb2withphylo <- ggtestallb2withphylo + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5) # ...lastly below changing from angle 80 & size 4 to angle 90 & size 6 (no, only for gsizeforphylo plot!)

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  ggtestallb2withphylo <- ggtestallb2withphylo + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)
} else {
  ggtestallb2withphylo <- ggtestallb2withphylo + theme_bw(base_family="sans") + geom_vline(xintercept = total_genomes_A + 0.5)
}
ggtestallb2withphylo <- ggtestallb2withphylo + theme(plot.title=element_text(hjust=0), axis.text.x = element_text(size = 4, angle = 80, hjust=1), axis.text.y = element_text(size = 4), panel.background = element_rect(fill = '#fafcf6'))
ggtestallb2withphylo <- ggtestallb2withphylo + facet_wrap(~ pathway_name_with_PML_score, scales = "free_y", labeller = labeller(pathway_name_with_PML_score = label_wrap_gen(56))) + theme(aspect.ratio = 1) + theme(strip.text.x = element_text(size = 6))
#ggtestallb2withphylo

#ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_focuspathways_byspecies_without_refECs_forphylo_",plot_title,Sys.Date(),".png", sep=""), ggtestallb2withphylo, width = 32, height = 28, units = "in", limitsize = FALSE)
## do not even need, just make combined plot with phylo!
#ggsave(filename = paste("heatmaps_focuspathways_byspecies_without_refECs_forphylo_",plot_title,Sys.Date(),".pdf", sep=""), ggtestallb2withphylo, width = 32, height = 28, units = "in", limitsize = FALSE)

######################################################
## adding one genome size plot to each enriched...
### NOTE NEW IDEA - CHANGE ORDER FROM GENOME SIZE TO PHYLO ORDER, THEN BELOW YOU WILL PLOT PHYLO...
gsizeforphylo <- ggplot(paths_groupAandB_stats3_gsize_plotswithphylo, aes(x=genusspecies, y=ec_number)) + scale_y_discrete(limits=rev)
## also lately removing labels in next line
gsizeforphylo <- gsizeforphylo + geom_tile(aes(fill = genome_size_mbp), color="white", linewidth=0.1) # + geom_text(aes(label = genome_size_mbp), color="red", size=2.5, position=position_jitter(width=0,height=0.2)) #height was 1 but multiple sizes
#gsizeforphylo <- gsizeforphylo + scale_fill_gradientn(name = "Mbp", colors = mutedGnBu)
gsizeforphylo <- gsizeforphylo + scale_fill_distiller(name = "Mbp", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
#gsizeforphylo <- gsizeforphylo + scale_fill_viridis_c(name = "", na.value = "transparent")
# gsizeforphylo <- gsizeforphylo + scale_fill_viridis_c(name = "Genome \nSize, \nMbp", na.value = "transparent")
#gsizeforphylo <- gsizeforphylo + scale_fill_viridis_c(name = "Genome Size, Mbp", na.value = "transparent") \nSHM
#gsizeforphylo <- gsizeforphylo + theme(legend.position="none") 
# gsizeforphylo <- gsizeforphylo + theme(axis.ticks=element_blank())  - theme(legend.key.width=unit(4,"cm"))
gsizeforphylo <- gsizeforphylo + labs(x=NULL, y=NULL, title="")
gsizeforphylo <- gsizeforphylo + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)

## adding an if/then for Macs vs PCs - former use Helvetica font, latter just specify 'sans' to avoid font warnings
if (grepl("apple", sessionInfo()[2]) == TRUE) {
  gsizeforphylo <- gsizeforphylo + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)
} else {
  gsizeforphylo <- gsizeforphylo + theme_bw(base_family="sans") + geom_vline(xintercept = total_genomes_A + 0.5)
}
#gsizeforphylo <- gsizeforphylo + theme(plot.title=element_text(hjust=0), axis.ticks=element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(angle = 80, hjust=1), panel.background = element_rect(fill = '#fafcf6'))
## need this command to line up phylogeny, moving legend to the left...lastly changing from angle 80 & size 4 to angle 90 & size 6, 6 too small trying 8
gsizeforphylo <- gsizeforphylo + theme(legend.position="left", legend.key.width=unit(0.1,"cm"), axis.ticks=element_blank(), plot.title=element_text(hjust=0), axis.text.x = element_text(size = 8, angle = 90, hjust=1), axis.text.y = element_text(size = 0), panel.background = element_rect(fill = '#fafcf6'))
# gsizeforphylo <- gsizeforphylo + theme(legend.position="left", legend.key.width=unit(0.1,"cm"), axis.ticks=element_blank(), plot.title=element_text(hjust=0), axis.text.x = element_text(size = 4, angle = 45, hjust=1), axis.text.y = element_text(size = 4))
# gsizeforphylo <- gsizeforphylo + theme(legend.position="none", axis.ticks=element_blank(), plot.title=element_text(hjust=0), axis.text.x = element_text(size = 4, angle = 45, hjust=1), axis.text.y = element_text(size = 4))
gsizeforphylo <- gsizeforphylo + facet_wrap(~ pathway_name, scales = "free_y") + theme(aspect.ratio = 1)
gsizeforphylo

## try gsizeforphylo without a legend
# paths_groupAandB_stats3_gsize_plotswithphylo$pathway_name <- gsub("Genome Size","Genome Size, Mbp",paths_groupAandB_stats3_gsize$pathway_name)
# gsizeforphylonolegend <- ggplot(paths_groupAandB_stats3_gsize_plotswithphylo, aes(x=genusspecies, y=ec_number)) + scale_y_discrete(limits=rev)
# gsizeforphylonolegend <- gsizeforphylonolegend + geom_tile(aes(fill = genome_size_mbp), color="white", linewidth=0.1) + geom_text(aes(label = genome_size_mbp), color="red", size=2.5, position=position_jitter(width=0,height=0.2)) #height was 1 but multiple sizes
# gsizeforphylonolegend <- gsizeforphylonolegend + scale_fill_distiller(name = "Mbp", palette = "GnBu", direction = 1, na.value = "#f7fcf0")
# gsizeforphylonolegend <- gsizeforphylonolegend + labs(x=NULL, y=NULL, title="")
# gsizeforphylonolegend <- gsizeforphylonolegend + theme_bw(base_family="Helvetica") + geom_vline(xintercept = total_genomes_A + 0.5)
# ## no legend
# gsizeforphylonolegend <- gsizeforphylonolegend + theme(legend.position="none", axis.ticks=element_blank(), plot.title=element_text(hjust=0), axis.text.x = element_text(size = 4, angle = 80, hjust=1), axis.text.y = element_text(size = 0), panel.background = element_rect(fill = '#fafcf6'))
# gsizeforphylonolegend <- gsizeforphylonolegend + facet_wrap(~ pathway_name, scales = "free_y") + theme(aspect.ratio = 1)
# gsizeforphylonolegend

## new layout to fit phylogeny as well...
## use dendro_plotref_flipped BELOW last plot
layouthkl4 <- rbind(c(1,1,1,2),
                    c(1,1,1,3),
                    c(1,1,1,NA))

layouthkl22 <- rbind(c(1,2),
                    c(1,3))

## grid arrange
twoplots_enrichedinB_withphylo <- grid.arrange(ggtestallb2withphylo,gsizeforphylo,dendro_plotref_forphyloheatmap, layout_matrix = layouthkl4)

#ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_focuspathways_byspecies_withphylogeny_",plot_title,Sys.Date(),".", sep=""), twoplots_enrichedinB_withphylo, width = 32, height = 28, units = "in", limitsize = FALSE)
## best dimensions here are 48x42 - for 2 plots with grobs 17.5x12 or 40x24
ggsave(filename = paste("heatmaps_focuspathways_byspecies_withphylogeny_",plot_title,Sys.Date(),".pdf", sep=""), twoplots_enrichedinB_withphylo, width = 48, height = 42, units = "in", limitsize = FALSE)

png(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/heatmaps_focuspathways_byspecies_withphylogeny_",plot_title,Sys.Date(),".png", sep=""), width = 48, height = 42, units = "in", res = 300)
twoplots_enrichedinB_withphylo <- grid.arrange(ggtestallb2withphylo,gsizeforphylo,dendro_plotref_forphyloheatmap, layout_matrix = layouthkl4)
Sys.sleep(2)
dev.off()
Sys.sleep(2)

## ALSO OUTPUT just gsizeforphylo,dendro_plotref_forphyloheatmap
# layouthkl2 <- rbind(c(1),
#                     c(2))
layouthkl2b <- rbind(c(1,1,1),
                    c(NA,2,NA))

#twoplots_onlyphylo <- grid.arrange(gsizeforphylo,dendro_plotref_forphyloheatmap, layout_matrix = layouthkl2b)
twoplots_onlyphylo2 <- grid.arrange(gsizeforphylo,dendro_plotref_forphyloheatmapwithmargin, layout_matrix = layouthkl2b)

#ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/suppfig_plot_onlywithgenomesizeandphylo_",plot_title,Sys.Date(),".png", sep=""), twoplots_onlyphylo, width = 24, height = 21, units = "in", limitsize = FALSE)
## dimensions here were 48x42 - for 2 plots with grobs WITH dendro_plotref_forphyloheatmapwithmargin best is 17.5x12 or 40x24
ggsave(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/suppfig_plot_onlywithgenomesizeandphylo_",plot_title,Sys.Date(),".pdf", sep=""), twoplots_onlyphylo2, width = 40, height = 24, units = "in", limitsize = FALSE)

png(filename = paste("supplemental_plots_ec_by_taxon_per_pathway/suppfig_plot_onlywithgenomesizeandphylo_",plot_title,Sys.Date(),".png", sep=""), width = 40, height = 24, units = "in", res = 300)
twoplots_onlyphylo2 <- grid.arrange(gsizeforphylo,dendro_plotref_forphyloheatmapwithmargin, layout_matrix = layouthkl2b)
Sys.sleep(2)
dev.off()
Sys.sleep(2)
print("Pipeline is complete!")
print("Output plots will be in your specified directory, with additional plots in the subfolders /supplemental_plots_ec_by_taxon_per_pathway & /supplemental_plots_taxon_by_pathway")

print("Pathways not annotated by BV-BRC (if any) will be listed below.")
paths_groupAandB_statsmissing <- paths_groupAandB_stats3 %>%
  dplyr::filter(genusspecies == "reference_set" & genecount_group_by_ecnumber != "NA") %>% 
  add_count(pathwayid) %>% 
  dplyr::filter(n > 6)
#print(unique(paths_groupAandB_statsmissing$pathway_name))

# paths_groupAandB_statsmissingtopull <- paths_groupAandB_stats3 %>%
#   dplyr::filter(genusspecies == "reference_set" & genecount_group_by_ecnumber != "NA") %>% 
#   expand(nesting(ec_number,pathway_name)) %>% 
#   add_count(pathway_name) %>% 
#   dplyr::filter(n > 4)
# print(unique(paths_groupAandB_statsmissingtopull$pathway_name))
# print("Thus they will have PML scores of 0 but may still be important.")

paths_groupAandB_statsmissingtopull <- paths_groupAandB_stats3 %>% 
  dplyr::filter(genusspecies == "reference_set" & genecount_group_by_ecnumber != "NA") %>% 
  group_by(ec_number,pathway_name) %>% 
  # summarize_all(first) %>%  # replacing all summarize_all commands
  summarize(across(everything(), first)) %>% 
  select(ec_number, pathwayid, pathway_id, pathway_name) %>% 
  add_count(pathway_name) %>% 
  dplyr::filter(n > 4)
print(unique(paths_groupAandB_statsmissingtopull$pathway_name))
print("These pathways will have PML scores of 0 but may still be important.")

## also removing any pathway files at end of run...change to data_dir via paste0
# unlink("~/pomelo_outputs/listA.*")
# unlink("~/pomelo_outputs/listB.*")
listofAs <- paste0(data_dir,"/listA.*" )
unlink(listofAs)
listofBs <- paste0(data_dir,"/listB.*" )
unlink(listofBs)

####################################################################################################################
####################################################################################################################
