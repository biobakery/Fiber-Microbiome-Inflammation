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
           # Figure 2. interaction between DDR fiber and P. copri on logCRP#
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

# read in MDS1 and MDS2
bugs_pcoa <- species_all
bugs_pcoa <- capscale(t(bugs_pcoa) ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )

# merge the data
association <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('participant', 'agemlvs','abx', 'bmi12', 'logcrp', 'logtchdl', 'logtc', 'loghdl', 'logtg', 'aofib10v', 'frtaf10v', 'ceraf10v', 'vegaf10v', 'a_aofib_fs_dr', 'calor_fs_dr', 'bristolcat', "act10", "a_pect_fo_dr","a_fat_fs_dr")] 

#merge association data with MDS bug scores
association_mds <- merge(x = association, y = ord.bug.scores, by = "row.names", all = TRUE)

# set first column as rownames and remove it
rownames(association_mds) <- association_mds[ ,1]
association_mds[ ,1] <- NULL

# merge with the bug data
mns <- colMeans(species_all, na.rm=TRUE)
order(-mns)
species_all_sort <- species_all[,order(-mns)]
colMeans(species_all_sort)

association_bugs <- merge(x = association_mds, y = species_all_sort, by = 'row.names', all = TRUE)
rownames(association_bugs) <- association_bugs[ ,1]
association_bugs[ ,1] <- NULL

head(association_bugs)
#sort by MDS1
association_bugs_sort <- association_bugs[order(association_bugs$MDS1),]

association_bugs$prevotella_copri2c <- with(association_bugs, ifelse(prevotella_copri>0, 1, 0))
table(association_bugs$prevotella_copri2c)
#0   1 
#702 223


library(lmerTest)
###########################################################################################
                                    # interaction model with MDS1 #
###########################################################################################
#crude interaction model: n=788, p-int=0.02
lm.crp <- lmer(formula = logcrp ~ a_aofib_fs_dr + MDS1 + a_aofib_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#further adjusting for age: n=788, p-int=0.025
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_aofib_fs_dr + MDS1 + a_aofib_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#main model, try further adjusting for some covariates: age, antibiotics, calorie; n=788, p-int=0.024
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_aofib_fs_dr + MDS1 + a_aofib_fs_dr*MDS1 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#additionally adjusting for bristolcat or physical activity did not change the results; n=767, p-int=0.0157
association_bugs$bristolcat <- as.factor(association_bugs$bristolcat)
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_aofib_fs_dr + MDS1 + a_aofib_fs_dr*MDS1 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)


###########################################################################################
                                    # interaction model with P copri #
###########################################################################################
#model adjusting for age only, p-int=0.01 ########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_aofib_fs_dr + prevotella_copri2c + a_aofib_fs_dr*prevotella_copri2c + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)


#model adjusting for other covariates, p-int=0.01 ########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_aofib_fs_dr + prevotella_copri2c + a_aofib_fs_dr*prevotella_copri2c + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)


#further adjusting for bristol score and physical activity, p-int=0.02 #########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_aofib_fs_dr + prevotella_copri2c + a_aofib_fs_dr*prevotella_copri2c + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)


###########################################################################################
                                   # Scatter plots by P copri #
###########################################################################################
#scatter plot with smooth line 
association_bugs$prevotella_copri2c <- as.factor(association_bugs$prevotella_copri2c)
head(association_bugs)

ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_aofib_fs_dr, y=logcrp, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color=prevotella_copri2c), size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary fiber intake (g/day)", y = "Log-transformed CRP")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title=element_text(size=20,face='bold'),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank())+
  scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=608)", "Yes (n=180)"))
#annotate("text", x = 50, y = 3, label = "P-interaction=0.01"))

plot(ggp_by_prevotella_copri2c)

ggsave(filename='./result/figures/interaction/ddr.crp.by.prevotella_copri2c.png', plot=ggp_by_prevotella_copri2c, width = 10, height = 5, dpi = 600) 



