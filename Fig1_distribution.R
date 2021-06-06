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
         # Figure 1. Distribution of data and PcoA on species-level taxonomy#
###########################################################################################


###########################################################################################
             # read in data#
###########################################################################################

rm(list=ls())
getwd()
setwd("/Users/wenjiema/Dropbox (Partners HealthCare)/MGH/AC/AC/microbiome/mlvs_fiber_crp/")
source("./code/final/Rstart.R")

# read in metadata (307)
z_metadata_all <- read.table('./maaslin2/pcl/z_metadata_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(z_metadata_all) <- z_metadata_all$sample
z_metadata_all <- z_metadata_all[,-1]
meta_dis <- as.data.frame(t(z_metadata_all))

# read in species
species_all <- read.table('./maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# read in metadata (925)
mlvs_metadata_925 <- read.table('./maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

###########################################################################################
              # Fig 1B. Distribution of DDR fiber and CRP#
###########################################################################################
graphics.off()

####distribution of ddr fiber, count##########
aofib_ddr_distribution<-ggplot(meta_dis, aes(x=a_aofib_fs_dr_w1))+
  geom_histogram(color="black", fill="grey", binwidth=5)+
  labs(x="Dietary fiber intake (g/day)", y = "Count")+
  xlim(5,55)+
  theme(axis.title=element_text(size=30,face = 'bold'),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(aofib_ddr_distribution)
ggsave(filename='./result/figures/distribution/aofib_ddr_histogram_count.png', plot=aofib_ddr_distribution, width = 12, height = 5, dpi = 600)

####distribution of logcrp, count##########
summary(meta_dis$logcrp_w1)
logcrp_distribution<-ggplot(meta_dis, aes(x=logcrp_w1))+
  geom_histogram(color="black", fill="grey", binwidth=0.5)+
  labs(x="log C-reactive protein", y = "Count")+
  xlim(-2,3)+
  theme(axis.title=element_text(size=30, face = 'bold'),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(logcrp_distribution)
ggsave(filename='./result/figures/distribution/logcrp_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)


###########################################################################################
            # Fig 1C. Correlation between fiber and age/BMI#
###########################################################################################
####scatter plots of ddr fiber and age####
ddrfiber_age <- ggplot(
  data=meta_dis, aes_string(x='agemlvs', y='a_aofib_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab('Age (years)') +  ylab("Dietary fiber intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=30, face = 'bold'),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )+
  annotate("text", x = 78, y=50, label = "spearman r=-0.02", size=10)
plot(ddrfiber_age)
ggsave(filename='./result/figures/distribution/ddrfiber_age.png', plot=ddrfiber_age, width = 10, height = 10, dpi = 600) 

cor.test(z_metadata_all$a_aofib_fs_dr_w1, z_metadata_all$agemlvs,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)



####scatter plots of ddr fiber and BMI####
ddrfiber_bmi <- ggplot(
  data=meta_dis, aes_string(x='bmi12', y='a_aofib_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  #xlab(expression(paste("Body mass index (kg/m"^"2"*")"))) +  #ylab("Dietary fiber intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        #axis.title=element_text(size=30),
        #axis.title.x=element_text(size=30, face = 'bold'),
        #axis.title.y = element_blank(),
        axis.title = element_blank(),
        #axis.text.x = element_text(size=30),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=30),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )+
  annotate("text", x = 35, y=50, label = "spearman r=-0.24", size=10)
plot(ddrfiber_bmi)
ggsave(filename='./result/figures/distribution/ddrfiber_bmi.png', plot=ddrfiber_bmi, width = 10, height = 10, dpi = 600) 

cor.test(z_metadata_all$a_aofib_fs_dr_w1, z_metadata_all$bmi12,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)


###########################################################################################
            # Fig 1D. stacked barplots of fiber sources#
##########################################################################################
barplot_cum <- meta_dis[,c('aofib10v','frtaf10v','ceraf10v','vegaf10v')]
barplot_cum <- barplot_cum %>% 
  mutate(other = aofib10v - frtaf10v - ceraf10v - vegaf10v)
rownames(barplot_cum) <- rownames(z_metadata_cum)
colnames(barplot_cum) <- c('total', 'fruit','cereal','vegetable','other')

# sort by aofib10v
barplot_cum_sort <- barplot_cum[order(-barplot_cum$total),]

# keep only fiber subsets
barplot_cum_sort_stratified <- barplot_cum_sort[,2:5]

# transform to relative abundances
#barplot_cum_sort_stratified_a <- sweep(barplot_cum_sort_stratified, 1, rowSums(barplot_cum_sort_stratified), `/`)

# transform to samples in columns and metadata in rows
t_barplot_cum_sort_stratified <- as.data.frame(t(barplot_cum_sort_stratified))


t_barplot_cum_sort_stratified$sources <- row.names(t_barplot_cum_sort_stratified)
dim(t_barplot_cum_sort_stratified) # 4 308
# melt
t_barplot_cum_sort_stratified_melt <- melt(t_barplot_cum_sort_stratified, id.vars = 'sources')
dim(t_barplot_cum_sort_stratified_melt) # 1228 3


top30_s_colors <- scale_fill_manual(values = c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D53E4F", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D", "#CCCCCC"))

bar_cum<-ggplot(data=t_barplot_cum_sort_stratified_melt, aes(x=variable, y=value, fill=sources)) +
  geom_bar(stat="identity")+
  theme_bw(base_size = 30) + 
  theme(plot.title = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title=element_text(face = 'bold'),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=30),
        axis.ticks.x=element_blank()) + 
  top30_s_colors + 
  xlab("Participant") + 
  ylab("Dietary fiber intake (g/day)") 
#scale_y_continuous(labels = percent_format(), expand = c(0, 0)) 
print(bar_cum)

ggsave(filename='./result/figures/distribution/fiber_subsets.png', plot=bar_cum, width = 15, height = 6, dpi = 600) 

###########################################################################################
          # Fig 1E. PCoA of species-level taxonomy (decorated with fruit fiber and CRP) #
##########################################################################################
meta_pcoa <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('frtaf10v','logcrp')] 

# create quartiles for frtaf10v
quantile(meta_pcoa$frtaf10v)
# 0%         25%         50%         75%        100% 
# 0.04428571  3.19000000  4.43000000  6.05571429 15.22000000 

meta_pcoa <- meta_pcoa %>%
  mutate(frtaf10vq = case_when(
    0 < frtaf10v & frtaf10v <= 3.19 ~ 'Q1',
    3.19 < frtaf10v & frtaf10v<= 4.43 ~ 'Q2',
    4.43 < frtaf10v & frtaf10v<= 6.05571429 ~ 'Q3',
    6.05571429 < frtaf10v ~ 'Q4'))
table(meta_pcoa$frtaf10vq)
meta_pcoa$frtaf10vq<- factor(meta_pcoa$frtaf10vq, 
                             levels = c('Q1', 'Q2', 'Q3', 'Q4'))

# create quartiles for logcrp
quantile(meta_pcoa$logcrp, na.rm=TRUE)
# 0%        25%        50%        75%       100% 
# -1.6094379 -0.9162907 -0.1053605  0.6418539  2.9704145 

meta_pcoa <- meta_pcoa %>%
  mutate(logcrpq = case_when(
    -2.0 <= logcrp & logcrp <= -0.9162907 ~ 'Q1',
    -0.9162907 < logcrp & logcrp<= -0.1053605 ~ 'Q2',
    -0.1053605 < logcrp & logcrp<= 0.6418539 ~ 'Q3',
    0.6418539 < logcrp ~ 'Q4'))
table(meta_pcoa$logcrpq)

library(tidyverse)
meta_pcoa$logcrpq <- fct_explicit_na(meta_pcoa$logcrpq, na_level = "Missing")

rownames(meta_pcoa) <- rownames(mlvs_metadata_925)  



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

# remove logCRP = NA
ord.bug.scores.meta <- ord.bug.scores.meta[complete.cases(ord.bug.scores.meta), ]

# number of samples
n_samples <- nrow(ord.bug.scores.meta)


library(viridis)
ordcolors <- scale_colour_gradientn(colours = viridis(10), limits=c(0,16))

pcoa.plot <- 
  ggplot( ord.bug.scores.meta, aes(MDS1, MDS2) ) + 
  geom_point(size=3, alpha = 0.75, aes( color = frtaf10v, shape = factor(logcrpq)))+
  scale_shape_manual(values=c(3, 16, 17, 15), name = "CRP quartiles")+
  labs(color = "Fruit fiber (g/day)")+
  theme_classic() + 
  #theme_bw(base_size = 30) +
  theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 30, face = 'bold'), 
        legend.text = element_text(size=20), legend.title = element_text(size=30),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'vertical',
        plot.title = element_text( size = 30, face = 'bold')
  ) +
  guides( color = guide_legend(nrow=2, byrow = TRUE) ) +
  xlab(paste0("PCo1 (",x_variance.bugs,"%)")) + ylab(paste0("PCo2 (",y_variance.bugs,"%)"))

pcoa.plot <- pcoa.plot +   ordcolors + guides(colour="colourbar")

print(pcoa.plot)

ggsave(filename='./result/figures/pcoa/pcoa_ordination_frtaf_crp.png', plot=pcoa.plot, width = 12, height = 10, dpi = 600) 

