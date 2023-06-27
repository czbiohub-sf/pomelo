![](pomelo_logo2.png)

# PoMeLo

PoMeLo is a novel R-based bioinformatics approach to identifying metabolic vulnerabilities of bacterial pathogens to host-directed therapeutics.

## Dependencies & R Session Info

To run PoMeLo, R & RStudio must be installed. PoMeLo was developed using R v4.2.2 with RStudio 2022.12.0.

Attached base R packages:	_stats, graphics, grDevices, utils, datasets, methods, base_

Other attached packages:
	_fs_1.6.1, lubridate_1.9.2, forcats_1.0.0, stringr_1.5.0, purrr_1.0.1, readr_2.1.4, tidyr_1.3.0, tibble_3.1.8, tidyverse_2.0.0, aplot_0.1.10, ggtreeExtra_1.8.1, ggtree_3.6.2, treeio_1.22.0, tidytree_0.4.2, dplyr_1.1.0, treedataverse_0.0.1, coop_0.6-3, openxlsx_4.2.5.2, phytools_1.5-1, maps_3.4.1, ggdendro_0.1.23, phangorn_2.11.1, seqinr_4.2-23, ggplot2_3.4.1, hash_2.2.6.2, reshape2_1.4.4, ape_5.7, gridExtra_2.3, gdata_2.18.0.1, lattice_0.20-45, scales_1.2.1, hexbin_1.28.2, RColorBrewer_1.1-3, igraph_1.4.0_

## Running PoMeLo

To run PoMeLo, clone this repo, open the ```PML_script.R``` file in RStudio, and then click "Source".

PoMeLo is closely linked to analyses using the BV-BRC website (https://www.bv-brc.org/), and requires files generated via BV-BRC. Here is a typical workflow incorporating both BV-BRC & PoMeLo:

* Users will first need to use BV-BRC to select genomes of two related, preferably monophyletic, groups of bacterial species (_target & non-target_). If the user wishes to determine their groupings via analysis of genome loss via a phylogenetically independent contrast approach, see the PIC section below.
    + When examining genomes in BV-BRC (under the 'genomes' tab) of interested taxa, we recommend including only genomes which are Complete (no WGS), and of Good Quality.
    + When examining BV-BRC genomes we also recommend to add the "Species" column as one of the viewed columns when viewing and selecting genomes. This column commonly has updated taxonomic information and is superior to the name used in the initial "Genome Name" column.
    + The genomes should be placed into separate BV-BRC Genome Groups. From the page listing the users' Genome Groups (https://www.bv-brc.org/workspace/yourusername@patricbrc.org/home/Genome%20Groups), next download both the target group & non-target group. Be sure to download each table as .csv and also be sure to not modify these .csv files after downloading, as Excel will drop zeros in genome IDs.

* It is not required, but the last portion of the PoMeLo pipeline incorporates a BV-BRC phylogenetic tree. If the user would like to include a phylogeny of the two groups, first create a separate Genome Group in BV-BRC incorporating the species in your target and non-target groups - but for this Genome Group, be sure to only include **one genome per "Species" column**. The user will also need to save this Genome Group as a .csv file.
    + Within BV-BRC, use the Bacterial Genome Tree tool (https://www.bv-brc.org/app/PhylogeneticTree) under the _Tools & Services_ menu to create the phylogeny. Select your Genome Group with **one genome per species**, but otherwise use default settings (e.g. 100 genes) and have the output saved to your Genome Groups folder. After ~2-20 hours depending on your group size, the completed phylogeny will be saved to a subfolder.
    + Find the folder with the completed phylogeny, and download the tree file ending in "_tree.nwk" - this file along with the .csv of the Genome Group will be used by PoMeLo.

* With these BV-BRC files in hand, start the PoMeLo pipeline. Users will be prompted to select the required files. Note the file ```mapping_GO_to_ecgene_and_ecpathway_toPATRIC.tab``` can be found in the scripts subfolder.

* All of the output files will be saved in your ```~/code``` directory, with additional plots in the subfolders ```/supplemental_plots_ec_by_taxon_per_pathway``` & ```/supplemental_plots_taxon_by_pathway```. There are prompts throughout the PoMeLo pipeline that will indicate progress, and if there are any intermediate problems.

## PIC analysis

To generate a phylogeny of bacterial species calculating genome size change via a PIC approach, we include a subpipeline: the ```PML_PICanalysis.Rmd``` file. This script incorporates both Python and R code, and will output a phylogeny as a .pdf file.  This script requires the same BV-BRC phylogeny (.nwk) and associated table of genomes (.csv) as desribed earlier.

## Tips

Users need only click on "Source", and pop-up windows will appear when required input files need to be selected. The script will also provide notes as it progreses, indicating progress or warnings if there are errors. For example based on the input species names, the pipeline will select a taxonomic name that is included in all output filenames - the script informs the user of this name and suggests chaning the filenames if it is incorrect. If the user does not wish to include a phylogeny in the final steps, they can simply select the 'cancel' button in the last two prompts for .nwk and .csv files, and the script will end there. Also note that if there are only 2 species in the analysis, the last steps incorporating a phylogeny will fail.

## Citation
If you use PoMeLo, please cite our publication:

**PoMeLo: A systematic computational approach to predicting metabolic loss in pathogen genomes**

Abigail Leigh Glascock, Eric Waltari, Gytis Dudas, Joan Wong, Vida Ahyong

tk 2023; doi: https://doi.org/tk

The two folders called "supplemental_materials..." contain both the inputs and output plots for the Treponema example in Glascock et al & the Rickettsia example in the companion paper by Medicielo et al.:

**Evolutionary genomics identifies host-directed therapeutics to treat intracellular bacterial infections**

Josette Medicielo, Eric Waltari, Abigail Leigh Glascock, Gytis Dudas, Brian DeFelice, Ira Gray, Cristina M. Tato, Joan Wong, Vida Ahyong
tk 2023
