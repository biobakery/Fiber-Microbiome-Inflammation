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
      # Figure 3. multivariable associations between fiber and species from MaAsLin2 #
###########################################################################################
###########################################################################################
                                      # read in data #
###########################################################################################

rm(list=ls())
getwd()
setwd("/Users/wenjiema/Dropbox (Partners HealthCare)/MGH/AC/AC/microbiome/mlvs_fiber_crp/")
source("./code/final/Rstart.R")

# read in metadata (925)
mlvs_metadata_925 <- read.table('./maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# read in species
species_all <- read.table('./maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# merge metadata and species
species_metadata <- merge(x = mlvs_metadata_925, y = species_all, by = "row.names", all = TRUE)
rownames(species_metadata) <- species_metadata[ ,1]
species_metadata[ ,1] <- NULL

ddr_plot <- species_metadata
checkit <- colnames(ddr_plot)
View(checkit)

for(i in 1:ncol(ddr_plot)) 
{ 
  ddr_plot[ , i] <- as.numeric(as.character(ddr_plot[, i])) 
}


###########################################################################################
                                     # Fig 3A. MaAsLin2 #
###########################################################################################

bug_file<-c('species_all')
metadata_file <- c('aofib_ffqcum', 'frtaf_ffqcum', 'vegaf_ffqcum', 'ceraf_ffqcum','aofib_ffq', 'aofib_ddr', 'aofib_ddr_food', 'pect_ddr_food','insol_ddr','sol_ddr'
                   ,'biomarker_logcrp', 'biomarker_logcrp_complete', 'aofib_ddr_logcrp')


# arcsin sqrt transformation  
for (a in 1:length(metadata_file))
{
  # (1) with no abundance/prevalence filter since i already did this myself
  
  for (b in 1:length(bug_file))
  {
    Maaslin2(input_data     = paste('./maaslin2/pcl/', bug_file[b], '.pcl', sep=''), 
             input_metadata   = paste('./maaslin2/pcl/', metadata_file[a], '.pcl', sep=''),
             output           = paste('./maaslin2/output/ast_tss/bugs/', metadata_file[a], '_', bug_file[b], '/', sep=''),
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


###########################################################################################
                                    # Figure 3B/C. scatter plots #
###########################################################################################
# Function of positive associations
scatter_plot_positive <- function(x,y,xlab,ylab) {
  plot <- ggplot(
    data=ddr_plot, aes_string(x=x, y=y)) +
    geom_point( aes(), fill = '#D9EF8B', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
    stat_smooth(method = "glm", color ='black')+ 
    guides(alpha='none')+labs("")+
    xlab(xlab) +  ylab(ylab)  +
    nature_theme +
    guides(legend.position=NULL)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.ticks =element_blank(),
          plot.title = element_blank())
  filepath <- paste('./result/figures/scatter_species/', y, '_', x, '.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 5, height = 5, dpi = 600) 
}# end of function

scatter_plot_positive('a_aofib_fs_dr','eubacterium_eligens',"Short-term dietary fiber\n(g/d)","Eubacterium eligens\n(relative abundance")
scatter_plot_positive('a_aofib_fs_dr','faecalibacterium_prausnitzii',"Short-term dietary fiber\n(g/d)","Faecalibacterium prausnitzii\n(relative abundance")
scatter_plot_positive('logcrp','bacteroides_uniformis',"Log CRP","Bacteroides uniformis\n(relative abundance)")


# Function of negative associations
scatter_plot_negative <- function(x,y,xlab,ylab) {
  plot <- ggplot(
    data=ddr_plot, aes_string(x=x, y=y)) +
    geom_point( aes(), fill = '#C92D39', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
    stat_smooth(method = "glm", color ='black')+ 
    guides(alpha='none')+labs("")+
    xlab(xlab) +  ylab(ylab)  +
    nature_theme +
    guides(legend.position=NULL)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.ticks =element_blank(),
          plot.title = element_blank())
  filepath <- paste('./result/figures/scatter_species/', y, '_', x, '.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 5, height = 5, dpi = 600) 
}# end of function

scatter_plot_negative('a_aofib_fs_dr','lachnospiraceae_bacterium_1_4_56faa',"Short-term dietary fiber\n(g/d)","Lachnospiraceae bacterium 1 4 56FAA\n(relative abundance")
scatter_plot_negative('a_aofib_fs_dr','ruminococcus_torques',"Short-term dietary fiber\n(g/d)","Ruminococcus torques\n(relative abundance")
scatter_plot_negative('logcrp','eubacterium_eligens',"Short-term dietary fiber\n(g/d)","Eubacterium eligens\n(relative abundance")


