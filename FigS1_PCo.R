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
        # Figure S1. PcoA on species-level taxonomy#
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


###########################################################################################
          # PCoA plots of samples ordinated on all bugs decorated with FFQ cum fiber #
###########################################################################################
# metadata
meta_pcoa <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('aofib10v', 'frtaf10v', 'ceraf10v', 'vegaf10v', 'a_aofib_fs_dr')] 

bugs_pcoa <- species_all
bugs_pcoa <- capscale(t(bugs_pcoa) ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )

# append metadata to scores
ord.bug.scores.meta <- merge(ord.bug.scores, meta_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.bug.scores.meta) <- ord.bug.scores.meta[ ,1]
ord.bug.scores.meta[ ,1] <- NULL

# number of samples
n_samples <- nrow(ord.bug.scores.meta)

# Generate ordination plot objects and store all in list object
ord.plot.list <- list()
library(viridis)

ordcolors <-  scale_colour_gradientn(colours = viridis(10))

pco_titles <- c('long-term total fiber', 'long-term fruit fiber', 'long-term cereal fiber', 'long-term vegetable fiber', 'short-term total fiber')

for( i in 1 : length( colnames( meta_pcoa ) ) ) {
  
  # plot ordination
  print(colnames(meta_pcoa)[i])
  pcoa.plot <- 
    ggplot( ord.bug.scores.meta, aes(MDS1, MDS2) ) + 
    geom_point(size = 2, alpha = 0.75, aes_string( color = colnames(meta_pcoa)[i] ) ) + 
    #geom_text(aes(label=rownames(ord.scores.meta))) +
    theme_classic() + 
    theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
          axis.ticks = element_blank(), axis.text = element_blank(),
          axis.title = element_text( size = 16, face = 'bold'), 
          legend.text = element_text(size=12), legend.title = element_blank(),
          legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
          plot.title = element_text( size = 20, face = 'bold')
    ) +
    guides( color = guide_legend(nrow=2, byrow = TRUE) ) +
    xlab(paste0("PCo1 (",x_variance.bugs,"%)")) + ylab(paste0("PCo2 (",y_variance.bugs,"%)")) + 
    ggtitle(pco_titles[i])
  
  if( is.numeric(meta_pcoa[,i]) ) {
    pcoa.plot <- pcoa.plot +   ordcolors + guides(colour="colourbar")
  }
  ord.plot.list[[i]] <- pcoa.plot
  rm(pcoa.plot, i)
}

# save ordinations to pdf, four ordinations per page
pdf( "./result/figures/pcoa/pcoa_ordination.pdf", onefile = TRUE)
marrangeGrob( ord.plot.list, nrow=1, ncol=1, top = NULL )
dev.off()
