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
# Figure S6. multivariable associations between fiber and CAZy RNA/DNA ratios from MaAsLin2 #
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

# read in cazyrnadna.qc
cazyrnadnaec.qc <- read.table('./maaslin2/pcl/cazyrnadnaec.qc.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(cazyrnadnaec.qc) <- cazyrnadnaec.qc$sample
cazyrnadnaec.qc <- cazyrnadnaec.qc[,-1]
cazyrnadnaec.qc <- as.data.frame(t(cazyrnadnaec.qc))

#Manually Log2 transformation
#merge with ddr fiber for scatter plot
#Manually Log2 transformation
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log2(y))
}
cazyrnadnaec.qc_log <- as.data.frame(apply(cazyrnadnaec.qc, 2, LOG))
cazyrnadna.qc_ddr <- merge (x=cazyrnadnaec.qc_log, y=mlvs_metadata_925, by = "row.names", all.x = TRUE)
rownames(cazyrnadna.qc_ddr) <- cazyrnadna.qc_ddr [ ,1]
cazyrnadna.qc_ddr [ ,1] <- NULL


ddr_plot <- cazyrnadna.qc_ddr
checkit <- colnames(cazyrnadna.qc_ddr)
View(checkit)

for(i in 1:ncol(ddr_plot)) 
{ 
  ddr_plot[ , i] <- as.numeric(as.character(ddr_plot[, i])) 
}


###########################################################################################
                                      # MaAsLin2 #
###########################################################################################

############### log transformation 
bug_file<-c('cazyrnadnaec', 'cazyrnadnaec.qc')
metadata_file <- c( 'aofib_ffqcum', 'frtaf_ffqcum', 'vegaf_ffqcum', 'ceraf_ffqcum', 'aofib_ffq', 'aofib_ddr', 'pect_ddr_food', 'biomarker_logCRP')

# log-transform (no such thing as AST for RNA/DNA ratio)  
for (a in 1:length(metadata_file))
  
{
  for (c in 1:length(bug_file))
  {
    Maaslin2(input_data     = paste('./maaslin2/pcl/', bug_file[c], '.pcl', sep=''), 
             input_metadata   = paste('./maaslin2/pcl/', metadata_file[a], '.pcl', sep=''),
             output           = paste('./maaslin2/output/log_tss/cazy/cazy_rna_dna_ratio/', metadata_file[a], '_', bug_file[c], '/', sep=''),
             normalization    = 'NONE', 
             standardize      = 'FALSE',
             transform        = 'LOG', 
             analysis_method  = 'LM', 
             random_effects   = 'participant',
             min_abundance    = 0, 
             min_prevalence   = 0, 
             cores            = 4)
  }
}

###########################################################################################
                                      # scatter plots #
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
  filepath <- paste('./result/figures/cazy/scatter/', y, '_', x, '.cazyrnadna.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 6, height = 6, dpi = 600) 
}# end of function

scatter_plot_positive('a_aofib_fs_dr','PL9',"Dietary fiber intake (g/day)","PL9 (log(RNA/DNA ratio))")
scatter_plot_positive('a_pect_fo_dr','PL9',"Pectin intake (g/day)","PL9 (log(RNA/DNA ratio))")

