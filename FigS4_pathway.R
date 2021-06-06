###########################################################################################
###Project: Dietary fiber intake, the gut microbiome, and chronic systemic inflammation
###Lead author: Wenjie Ma
###Study population: Men's Lifestyle Validation Study (MLVS)
###Data: 
#Recent (7DDRs) total fiber and pectin 
#Cumulative average of total dietary fiber and fiber from cereals, fruits, and vegetables (HPFS FFQs from 1986 to 2010)
#Circulating C-reactive protein
#Stool microbiome 
###Last update date: 5/28/2021 for Genome Medicine
#1) species-level taxonomic composition
#2) MetaCyc pathways (DNA and RNA/DNA ratio)
#3) EC analysis (DNA and RNA/DNA ratio)
#4) CAZy analysis (DNA and RNA/DNA ratio)

###########################################################################################
# Figure S4. multivariable associations between fiber and metagenomic pathways from MaAsLin2 #
###########################################################################################
###########################################################################################
                                      # read in data #
###########################################################################################
rm(list=ls())
getwd()
setwd("/Users/wenjiema/Dropbox (Partners HealthCare)/MGH/AC/AC/microbiome/mlvs_fiber_crp/")
source("./code/final/Rstart.R")

###########################################################################################
                                      # MaAsLin2 #
###########################################################################################
gpd_file<-c('gpd.qc', 'gpd.all')
metadata_file <- c( 'aofib_ffqcum', 'frtaf_ffqcum', 'vegaf_ffqcum', 'ceraf_ffqcum','aofib_ffq', 'aofib_ddr', 'aofib_ddr_food', 'pect_ddr_food', 'biomarker_logCRP','insol_ddr','sol_ddr')



# arcsin sqrt transformation
for (a in 1:length(metadata_file))
{
  # (a) with abundance/prevalence filter
  for (b in 1:length(gpd_file))
  {
    Maaslin2(input_data     = paste('./maaslin2/pcl/', gpd_file[b], '.pcl', sep=''), 
             input_metadata   = paste('./maaslin2/pcl/', metadata_file[a], '.pcl', sep=''),
             output           = paste('./maaslin2/output/ast_tss/gpd/', metadata_file[a], '_', gpd_file[b], '/', sep=''),
             normalization    = 'NONE', 
             standardize      = 'FALSE',
             transform        = 'AST', 
             analysis_method  = 'LM', 
             random_effects   = 'participant',
             min_abundance    = 0, 
             min_prevalence   = 0, 
             cores            = 4) 
  }
}  
