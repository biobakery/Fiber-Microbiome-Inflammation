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
        # Figure 4. multivariable associations between fiber and CAZy DNA from MaAsLin2 #
###########################################################################################
###########################################################################################
                                      # Fig 4A. MaAsLin2 #
###########################################################################################
dna_file<-c('cazydna925.qc')
metadata_file <- c( 'aofib_ffqcum', 'frtaf_ffqcum', 'vegaf_ffqcum', 'ceraf_ffqcum','aofib_ffq', 'aofib_ddr', 'aofib_ddr_food', 'pect_ddr_food', 'biomarker_logCRP')


# arcsin sqrt transformation
for (a in 1:length(metadata_file))
{
  # (a) with abundance/prevalence filter
  for (b in 1:length(dna_file))
  {
    Maaslin2(input_data     = paste('./maaslin2/pcl/', dna_file[b], '.pcl', sep=''), 
             input_metadata   = paste('./maaslin2/pcl/', metadata_file[a], '.pcl', sep=''),
             output           = paste('./maaslin2/output/ast_tss/cazy/cazydna/', metadata_file[a], '_', dna_file[b], '/', sep=''),
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




bug_colors <- list(
  # Bacteroidetes (red/yellow)
  "Bacteroides_sp_1_1_6"         = "gold",
  "Bacteroides_thetaiotaomicron" = "coral",
  "Bacteroides_vulgatus"         = "pink",
  "Bacteroides_dorei"            = "red1",
  "Bacteroides_fragilis"         = "darkgoldenrod3",
  "Bacteroides_caccae"           = "peru",
  "Bacteroides_ovatus"           = "chocolate",
  "Bacteroides_massiliensis"     = "violetred1",
  
  
  # Firmicutes (blue/purple)
  "Coprococcus_catus"            = "aquamarine3",
  "Coprococcus_sp_ART55_1"       = "cadetblue2",
  "Eubacterium_eligens"          = "dodgerblue1",
  "Eubacterium_siraeum"          = "deepskyblue2",
  "Eubacterium_rectale"          = "cyan2",
  "Eubacterium_limosum"          = "royalblue3",
  "Faecalibacterium_prausnitzii" = "darkblue",
  "Roseburia_intestinalis"       = "slateblue3",
  "Roseburia_inulinivorans"      = "skyblue4",
  "Ruminococcus_champanellensis" = "darkorchid2",
  "Anaerotruncus_colihominis"    = "cornflowerblue",
  "Clostridium_sp_L2_50"         = "mediumorchid3",
  "Subdoligranulum_variabile"    = "paleturquoise3",
  
  # Actinobacteria (green)
  "Akkermansia_muciniphila"      = "green3",
  "Bifidobacterium_longum"       = "yellowgreen",
  "Bifidobacterium_bifidum"      = "chartreuse2",
  
  "Other"                        = "#555555",
  "unclassified"                 = "#aaaaaa")

bugcolor <- scale_fill_manual(values = unlist(bug_colors))


###########################################################################################
                                    # Fig 4B. Stacked barplot (DNA) #
###########################################################################################
# beginning of function plot_stratified_distribution
plot_stratified_distribution <- function(pwy, name, pcl, label) {
  pwy <- grepl(name, rownames(pcl), fixed=T)
  sum_stratified <- pcl[pwy,]
  t_sum_stratified <- as.data.frame(t(sum_stratified))
  #dim(t_sum_stratified)
  
  # split into two data
  t_sum <- as.data.frame(t_sum_stratified[,1])
  rownames(t_sum) <- rownames(t_sum_stratified)
  colnames(t_sum) <- c(name)
  
  t_stratified <- t_sum_stratified[,2:nrow(sum_stratified)]
  # shorten species names
  checkit <- colnames(t_stratified)
  View(checkit)
  colnames(t_stratified) <- gsub(".*s__","", colnames(t_stratified))
  colnames(t_stratified) <- gsub(paste(name, ".", sep = ""), "", colnames(t_stratified))
  
  # step 5: transform to relative abundance by sum normalization and filtering
  t_stratified[ is.na(t_stratified) ] <- NA
  # make numeric
  for(i in 1:ncol(t_stratified)) 
  { 
    t_stratified[ ,i] <- as.numeric(as.character(t_stratified[,i])) 
  }
  dim(t_stratified) #372 5
  
  # sum normalize - transform to relative abundance and filtering
  rowSums(t_stratified)
  t_stratified <- sweep(t_stratified, 1, rowSums(t_stratified), `/`)
  t_stratified[is.na(t_stratified)] <- 0
  stratified <- as.data.frame(t(t_stratified)) 
  sum <- as.data.frame(t(t_sum)) 
  
  #sort by mean relative abundance
  mns <- rowMeans(stratified, na.rm=TRUE)
  order(-mns)
  stratified <- stratified[order(-mns),]
  rowMeans(stratified)
  
  # step 6: keep 10 species, and group species beyond 10 to "other"
  if (nrow(stratified)>10)
  {
    other <- as.data.frame(t(colSums(stratified[11:nrow(stratified),])))
    rownames(other) <- "Other"
    stratified <- rbind(stratified[1:10,], other)
  }
  
  ###########################
  ####stacked barplot########
  ###########################
  stratified$species <- row.names(stratified)
  stratified_names <- c(rownames(stratified))
  stratified$species <- factor(stratified$species, ordered = TRUE, levels = stratified_names)
  
  # melt
  stratified_melt <- melt(stratified, id.vars = 'species')
  dim(stratified_melt) # 1705 3
  
  #top30_s_colors <- scale_fill_manual(values = c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D53E4F", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D", "#CCCCCC"))
  
  bar<-ggplot(data=stratified_melt, aes(x=variable, y=value, fill=species)) +
    geom_bar(stat="identity")+
    theme_bw(base_size = 20) + 
    theme(plot.title = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=20),
          axis.ticks.x=element_blank()) + 
    guides(fill = guide_legend(ncol = 1,keyheight = 0.55)) + 
    bugcolor + 
    xlab("Sample") + 
    ylab("DNA") +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) 
  print(bar)
  
  filepath <- paste('./result/figures/cazy/barplot/', name, '_', label, '_barplot.png', sep ="" )
  ggsave(filename=filepath, plot=bar, width = 12, height = 5, dpi = 600) 
}# end of stacked barplot function

########################################################################
############################DNA###################################
# read in stratified cazydna file
cazydna372_sortpect <- read.table('./maaslin2/pcl/cazydna372_sortpect_stratified.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(cazydna372_sortpect) <- cazydna372_sortpect$sample
cazydna372_sortpect <- cazydna372_sortpect[,-1]


plot_stratified_distribution(pl9, "PL9", cazydna372_sortpect, 'dna')
plot_stratified_distribution(gh25, "GH25", cazydna372_sortpect, 'dna')
plot_stratified_distribution(gh29, "GH29", cazydna372_sortpect, 'dna')
plot_stratified_distribution(cbm13, "CBM13", cazydna372_sortpect, 'dna')


########################################################################
############################RNA###################################
# read in stratified cazyrna file
cazyrna372_sortpect <- read.table('./maaslin2/pcl/cazyrna372_sortpect_stratified.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(cazyrna372_sortpect) <- cazyrna372_sortpect$sample
cazyrna372_sortpect <- cazyrna372_sortpect[,-1]

plot_stratified_distribution(pl9, "PL9", cazyrna372_sortpect, 'rna')
plot_stratified_distribution(gh25, "GH25", cazyrna372_sortpect, 'rna')
plot_stratified_distribution(gh29, "GH29", cazyrna372_sortpect, 'rna')
plot_stratified_distribution(cbm13, "CBM13", cazyrna372_sortpect, 'rna')


