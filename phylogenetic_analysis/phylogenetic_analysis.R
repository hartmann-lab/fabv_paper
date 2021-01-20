#
'phylogenetic_analysis.R contains all the code used to generate figures and data that generated phylogenetic tree or phylogenetics-based analyses'
#
#
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
#plotting
library(cowplot)
library(ggsci)
library(scales)
library(RColorBrewer)
#Tree
library(ape)
library(ggrepel)
library(ggtree)
library(ggstance)
library(phytools)
library(phylogram)
library(dendextend)
library(treespace)
#ASR
library(ggimage)

main_path <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/manuscript/ismej_submission/github/fabv_paper/phylogenetic_analysis'

#### Section 0 ####
#### 34 gene MLSA phylogenetic tree construction
section_0_path <- paste(main_path,'section_0',sep='/')
setwd(section_0_path)
#
#
## taxonomic assignments based off hesse et al
df_taxassignments <- read.csv('pseudomonas_groups1.csv',colClasses=c("character"))%>%select(species,clade)
## metdata for 160 pseudomonas type strains, 30 study isolates, and c japonicus outgroup
df_meta <- read.csv('processed_genomes_metadata.csv',colClasses=c("character"))%>%select(speciesID,refseqSpecies,genID,source)
## fabv blast search results. used to map to pseudomonas with fabv
df_fabv <- read.csv('pao1_fabv_filteredhomologs.csv',colClasses=c("character"))
## merge all datasets
df_meta <- merge(df_meta,df_taxassignments,all.x=TRUE,by.x='speciesID',by.y='species')
df_meta <- df_meta%>%mutate(fabv_pos = case_when(speciesID%in%df_fabv$speciesID==TRUE~1,TRUE~0))
## read in 34 mlsa concatenated phyogeny
tree_species <- read.tree('RAxML_bipartitions.conctat34mlsa')
tree_species <- ape::ladderize(tree_species,F)
tree_species <- ape::root(tree_species,outgroup='cjaponicus')
## make dataframe from newick tree and add meta data from df_meta
df_treespecies <- data.frame('speciesID'=tree_species$tip.label,'order'=seq(1:length(tree_species$tip.label)))
df_treespecies <- merge(df_meta%>%select(speciesID,clade,source,fabv_pos),df_treespecies,by='speciesID')
df_treespecies <- df_treespecies%>%arrange(order)%>%mutate(source2=na_if(source,'REFSEQ'))
## remove all bootstraps that are better than 99%. remove 'Root' label
tree_species$node.label[as.numeric(tree_species$node.label)>=99] <- ''
tree_species$node.label[as.character(tree_species$node.label)=='Root'] <- ''
#
## artificially lower the branch length of the C japonicus outgroup to emphasize associations between Pseudomonas
xtree <- tree_species
xtree[["edge.length"]][380] <- mean(xtree[["edge.length"]])
p <- ggtree(xtree,branch.length = 'rate')+geom_nodelab(aes(label=label),nudge_x=-0.02,nudge_y=.4,size=2.5)
## plot phylogeny. distinguish between study isolates and type strains
p1 <- p%<+%df_treespecies+geom_tiplab(size=3.5,offset=0,align=FALSE)+scale_color_manual(values=c('maroon','white'))+scale_size_manual(values=3)+
  geom_treescale(y=-10,fontsize=5,linesize=1,offset=1)+geom_tippoint(aes(color=as.factor(source2)))
#
## Add in colors for clades
colors_ <- colorRampPalette(brewer.pal(12,'Paired'))
df_cladecolors <- df_treespecies%>%select(speciesID,clade)
df_cladecolors <- df_treespecies%>%select(clade)
rownames(df_cladecolors) <- df_treespecies$speciesID
## Finish plotting and save
plot_out <- gheatmap(p1, df_cladecolors, width = 0.1, color = "black",offset=.1)+
  scale_fill_manual(values=colors_(length(unique(df_treespecies$clade))))
#
save_plot('sec0_34mlsa_phylogeny.png',plot_out, base_height = 22, base_aspect_ratio = .4,bg='white')



#### Section 1 ####
#### Ancestral state reconstruction
section_1_path <- paste(main_path,'section_1',sep='/')
setwd(section_1_path)
#
## fabv blast search results. used to map to pseudomonas with fabv
df_fabv <- read.csv('pao1_fabv_filteredhomologs.csv',colClasses=c("character"))
## read in 34 mlsa concatenated phyogeny
tree_species <- read.tree('RAxML_bestTree.conctat34mlsa')
tree_species <- ape::root(tree_species,outgroup='cjaponicus')
tree_species <- ape::ladderize(tree_species,F)
## map fabv metadata to dataframe generated from newick
df_treespecies <- data.frame('speciesID'=tree_species$tip.label,'order'=seq(1:length(tree_species$tip.label)))
df_treespecies <- df_treespecies%>%mutate(fabv_pos=case_when(speciesID%in%df_fabv$speciesID~'fabvPos',TRUE~'fabvNeg'))%>%arrange(order)
#
## Use ACE for ancestral state reconstruction
char1 <- df_treespecies$fabv_pos 
names(char1) <- tree_species$tip.label
tree_species$edge.length[tree_species$edge.length == 0] <- 1e-7
#
fitER <- ace(char1, tree_species, type="discrete",
             method="ML",CI=TRUE,model="ARD",
             scaled=TRUE,kappa=1,
             corStruct = NULL,
             ip = 0.1,
             use.expm = FALSE,
             use.eigen = TRUE,
             marginal = FALSE)
#
ancstats <- as.data.frame(fitER$lik.anc)
ancstats$node <- 1:Nnode(tree_species)+Ntip(tree_species)
## convert pies with values <.1 fabV pos to NA
ancstats <- ancstats%>%mutate(fabvPos = case_when(fabvPos>.9~NA_real_,
                                                  fabvPos<.1~NA_real_,
                                                  TRUE~fabvPos),
                              fabvNeg = case_when(fabvNeg>.9~NA_real_,
                                                  fabvNeg<.1~NA_real_,
                                                  TRUE~fabvNeg))


#
#
## drmatically shortening C japonicus outgroup branch length to focus on Pseudomonas relationships
tree_species[["edge.length"]][1] <- 0.01
pies <- nodepie(ancstats,cols=c('fabvNeg','fabvPos'),alpha=0.8)#, cols=1:3)
pies <- lapply(pies, function(g) g + scale_fill_manual(values = c('blue','maroon')))
#
pg <- ggtree(tree_species)+xlim(0,5)+geom_tiplab(size=2.6,align = TRUE,offset=.05)
#
pg2 <- pg %<+% df_treespecies+
  geom_tippoint(aes(fill = fabv_pos,x=x+0.025), shape=22, size=1,color='white')+
  scale_fill_manual(values = c('blue','maroon'))+
  geom_inset(pies,width=.15,height=.15,hjust=.025)+
  theme(legend.position='none')+
  theme_transparent()
##
save_plot('sec1_ace_reconstruction.png',pg2,base_height = 16, base_aspect_ratio = .8,bg='white')
## Ancestral state reconstruction stats
fitER$se
fitER$loglik 
fitER$rates
fitER$lik.anc # these are the values that are plotted in figure


#### Section 2 ####
#### Make species-gene tree using the 34 mlsa concat newick and fabv newick for the type strain & study isolate genomes
section_2_path <- paste(main_path,'section_2',sep='/')
setwd(section_2_path)
#
## Load in fabv hits and speciesmetadata
df_fabv <- read.csv('pao1_fabv_filteredhomologs.csv',colClasses=c("character"))
df_species <- read.csv('processed_genomes_metadata.csv')%>%select(cleaned_filename,source,species,genID,speciesID)%>%unique()
## read in species tree
tree_species <- read.tree('RAxML_bipartitions.conctat34mlsa')
## root and prune the species tree to keep only genomes with fabV
df_drop <- df_species%>%filter(!speciesID%in%df_fabv$speciesID)
tree_species <- drop.tip(tree_species,as.character(df_drop$speciesID))
tree_species <- ape::root(tree_species,outgroup='cjaponicus')
## save the pruned species tree
write.tree(tree_species,'pruned_RAxML_bipartitions.conctat34mlsa')
#
#
## re-order and root the fabV gene tree so that it matches the species tree.
df_gene <- df_fabv
tree_gene <- read.tree('RAxML_bipartitions.pao1_fabv_filteredhomologs')
df_dropgene <- df_gene%>%filter(!speciesID%in%df_fabv$speciesID)
tree_gene <- drop.tip(tree_gene,as.character(df_dropgene$inputtedDelimID))
#
df_geneorder <- data.frame('inputtedDelimID'=tree_gene$tip.label,'order'=seq(1:length(tree_gene$tip.label)))
df_geneorder <- merge(df_geneorder,df_gene,by='inputtedDelimID')
df_geneorder <- df_geneorder%>%arrange(order)
tree_gene$tip.label <- as.character(df_geneorder$speciesID)
tree_gene <- ape::root(tree_gene,outgroup='cjaponicus')
## save the modified gene tree
write.tree(tree_gene,'pruned_RAxML_bipartitions.pao1_fabv_filteredhomologs')
#
#
## Make the tree difference plot
## read in pruned species and gene trees
tree_species <- read.tree('pruned_RAxML_bipartitions.conctat34mlsa')
tree_gene <- read.tree('pruned_RAxML_bipartitions.pao1_fabv_filteredhomologs')
## modfy into dendograms so that its compatabile with plotTreeDiff
dend_species <-  as.dendrogram.phylo(tree_species)
dend_gene <- as.dendrogram.phylo(tree_gene)
combined_dend <- dendlist(dend_species,dend_gene)
#
## Code for tanglegram is hashed out
#combined_dend <- untangle(combined_dend,method='step1side') # value used by tanglegram
#entanglement(combined_dend)  # metric used to optimize tanglegram
#png(file=paste(geneoutfile,'_tangle.png',sep=''),width=1400,height=900)
#tanglegram(combined_dend,margin_inner=5,margin_outer=5,axes=FALSE,dLeaf=0,sort=FALSE)
#dev.off()
#
## Tree Difference 
png(file='fabv_concat34mlsa_diff.png',width=1400,height=900,bg='white')
plotTreeDiff(tr1=tree_species,tr2=tree_gene,col1='#4575b4',col2='#d73027',baseCol = 'black',
             font=1,edge.width=1.5,cex=1.2)
dev.off()
#
#
## Pair-wise topography and distance metric comparisons between fabV and species tree
## Read in pruned gene and species trees
tree_gene <- read.tree('pruned_RAxML_bipartitions.conctat34mlsa')
tree_species <- read.tree('pruned_RAxML_bipartitions.pao1_fabv_filteredhomologs')
## Calculate robinson-foulds, kuhner-felsenstein, kendall-colijn 
a <- phangorn::RF.dist(tree_species,tree_gene,rooted=TRUE,normalize=TRUE) ## This is normalized robinson-foulds distance reported in the manuscript
b <- phangorn::RF.dist(tree_species,tree_gene,rooted=FALSE,normalize=FALSE)
c <- phangorn::wRF.dist(tree_species,tree_gene,rooted=FALSE,normalize=FALSE)
d <- phangorn::KF.dist(tree_species,tree_gene,rooted=FALSE)
e <- treeDist(tree_species,tree_gene) ## This is the Kendall-Colijn metric reported in the manuscript
f <- phangorn::treedist(tree_species,tree_gene)
df_stats <- data.frame('gene'='fabv','robinson_foulds_normal'=a,'robinson_foulds'=b,'robinson_foulds_weighted'=c,'kuhner_felsenstein'=d,'kendall_colijn'=e)
## write all distances to a single file
write.csv(df_stats,'fabv_concat34mlsa_tree_distances.csv',row.names=FALSE)


#### Section 3 ####
#### Graph the number of genomes that have have fabV, separated by the species group assigned by ANI
section_3_path <- paste(main_path,'section_3',sep='/')
setwd(section_3_path)
#
#
## read in metadata with all genome info
df_meta <- read.csv('processed_genomes_metadata.csv')%>%select(filename,cleaned_filename)
## read in species group assignments generated by ANI
df_tax <- read.csv('species_assignments.csv')
## read in fabv blast results
df_fabv <- read.csv('pao1_fabv_homologs.csv')
## merge all datasets
df_fabv <- merge(df_fabv,df_meta,by='cleaned_filename')
df_fabv <- merge(df_fabv,df_tax,by.x='filename',by.y='subject',all.y=TRUE)
## generate presence absence dataset. remove instances where the coverage is less than 80
df_fabv <- df_fabv%>%mutate(fabv_pa=case_when(is.na(pident)==TRUE~0,TRUE~1))
df_fabv <- df_fabv%>%filter(qcovs>80|fabv_pa==0)
#
write.csv(df_fabv,'fabv_homologs_processed.csv',row.names=FALSE)
#
#
## Create a per-species group percentage dataframe of fabv count categories
#
df_fabv <- read.csv('fabv_homologs_processed.csv')
#
df_fabvx <- df_fabv%>%group_by(filename)%>%mutate(genome_fabvcount=sum(fabv_pa))%>%ungroup()%>%
  select(filename,mod_species,ani,clade,fabv_pa,genome_fabvcount)%>%unique()%>%ungroup()%>%
  group_by(mod_species)%>%mutate(species_count=length(filename))%>%ungroup()%>%
  group_by(mod_species,genome_fabvcount)%>%mutate(category_count=length(filename))%>%ungroup()%>%
  select(mod_species,clade,fabv_pa,genome_fabvcount,species_count,category_count)%>%unique()%>%ungroup()%>%
  group_by(mod_species,genome_fabvcount)%>%
  mutate(category_perc=100*(category_count/species_count))%>%ungroup()%>%
  arrange(clade,species_count)%>%mutate(plot_order=seq(1:length(mod_species)))
## write this to file
write.csv(df_fabvx,paste(path_largeanalysis,'fabv_pa.csv',sep='/'),row.names=FALSE)
#
#
## integrating fabv p/a (fabv_pa.csv) plot with species counts plot (reorderd_speciesassignments.csv generated during ANI analysis)
df_sp <- read.csv('reorderd_speciesassignments.csv')
df_fabv <- read.csv('fabv_pa.csv')
## merge datasets
df_fabvx <- merge(df_fabv,df_sp,by='mod_species')
length(unique(df_sp$mod_species)) # check number of unique species groups in the dataset
#
## Modification being made to nitroreducens counts. it was found that one of the nitroreducens genomes was a clonal copy of an already existing strain, HBP1. This extra strain was removed
#
df_fabvy <- df_fabvx 
df_fabvy <- df_fabvy%>%filter(mod_species=='nitroreducens')%>%mutate(species_count.x=c(7,7))%>%mutate(category_count=c(2,5))%>%
  mutate(category_perc=c(20,80))%>%mutate(species_count.y=c(7,7))
df_fabvx <- df_fabvx%>%filter(mod_species!='nitroreducens')
df_fabvx <- rbind(df_fabvx,df_fabvy)
#
## plot the fabv counts within bars of species group log counts
df_fabvy <- df_fabvx%>%filter(mod_species!='antarctica')%>%mutate(category_log=log(category_count,base=10))%>%
  group_by(mod_species)%>%
  mutate(count_log=log(species_count.y,base=10))%>%mutate(perc_log=(category_perc/100)*count_log)
#
df_plotfabv <- df_fabvy%>%select(mod_species,clade.y,species_count.y,subclade,cladeorder,clade2,plot_order.y,genome_fabvcount,perc_log,plot_order.x,category_count,count_log)%>%unique()
#
df_plotfabv$genome_fabvcount <- as.factor(df_plotfabv$genome_fabvcount)
## df_subset is for the genome_count labels at the end of the graphs
df_subset <- df_plotfabv%>%select(mod_species,species_count.y,count_log)%>%unique()%>%mutate(genome_fabvcount=NA)
#
## plot the species counts as bars scaled with the log function. percentage of different fabv counts within the bars
#
p1 <- ggplot(df_plotfabv,
             aes(x=reorder(mod_species,plot_order.y),y=perc_log,fill=factor(genome_fabvcount,levels=c(2,1,0))))+
  coord_flip()+
  geom_bar(stat='identity',color='black')+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 0),
        axis.text.x = element_text(size = 30),axis.text.y = element_text(size = 15,hjust = 1),legend.position = 'none')+
  theme_transparent()+
  #theme_transparent()+
  scale_y_continuous(breaks=c(0,1,2,3,4),labels=c('0','1','2','3','4'),limits = c(0,4))+
  ylab(label=expression(Species~count~" "~"("~Log[10]~")"))+
  geom_text(data=df_subset,aes(x=mod_species,y=count_log,label=species_count.y),position=position_dodge(width=.1),hjust=-.1,size=5.5)+
  scale_fill_manual(values=c('#2c7bb6','#fdae61','#d7191c')) #0,1,2
#
save_plot('sec3_fabv_speciesgroup_counts.png',p1,base_height=20,base_aspect_ratio = .5,bg='white')


#### Section 4 ####
#### Plot the fabv gene tree for the expanded seq of refseq genomes
#
section_4_path <- paste(main_path,'section_4',sep='/')
setwd(section_4_path)
#
## Create a dataframe with the refseq genomes that demonstrate likely hgt fabv. reason column indicates either fabv duplication or unusual fabv presence given other species group members are fabv negative
df_interest <- data.frame('cleaned_filename'=c('91_AZRU01000001','91_JUEH01000001','1_LOHH01000001','9_CP025229','91_CP049140',
                                               '103_LFQO01000001','103_MOBN01000001','130_LHVG01000001'),
                          'reason'=c('dup','dup','dup','dup','dup',
                                     'unu','unu','unu'))
## read in the processed fabv homology search results. this list has been filtered for qcovs > 80 and pident > 30
df_fabv <- read.csv('fabv_homologs_processed.csv')
## merge datasets
df_fabv <- merge(df_fabv,df_interest,by='cleaned_filename',all.x=TRUE)
#
## Read in the fabV tree for plotting
tree_gene <- read.tree('RAxML_bipartitions.fabv_homologs_processed_trimmed')
## midpoint root because we are not using an outgroup, no cjaponicus here.
#tree_gene <- phytools::midpoint.root(tree_gene) #  cstack error for me on my machine. fortunately phangorn produces the same output 
tree_gene <- phangorn::midpoint(tree_gene)
#
## Create dataframe of newick tree and add fabv and species group metadata
df_treegene <- data.frame('sseqid'=tree_gene$tip.label,'order'=1:length(tree_gene$tip.label))
df_treegene <- merge(df_treegene,df_fabv,by='sseqid',all.x=TRUE)  
df_treegene <- df_treegene%>%mutate(xtiplab=case_when(cleaned_filename%in%df_interest$cleaned_filename~as.character(sseqid),TRUE~''))
#
#
## prune aeruginosa  genomes. there are ~3900 aeruginosa genomes and will skew the phylogeny while not providing additional info.
## keep only 10 aeruginosa total + 2 of interest. Also, remove nitroreducens 91_AZRU01000001 since it is a clone
#
df_aerug <- df_treegene%>%filter(mod_species=='aeruginosa')%>%arrange(desc(pident))
df_aerug <- df_aerug[1:10,]
df_treegene <- df_treegene%>%mutate(dontprune = case_when(mod_species!='aeruginosa'~1,
                                                          cleaned_filename%in%df_interest$cleaned_filename~1,
                                                          cleaned_filename%in%df_aerug$cleaned_filename~1,
                                                          TRUE~0)) 
#
## Drop list are all the genomes that will be removed from the final tree. pretty much is aeruginosa and one nitroreducens HBP1 clone. 
drop_list <- df_treegene%>%filter(dontprune==0|cleaned_filename=='91_AZRU01000001')%>%select(sseqid)
tree_gene1 <- drop.tip(tree_gene,as.character(drop_list$sseqid))
#
## Plot fabv gene tree, no alterations
#p1 <- ggtree(tree_gene1,branch.length='rate',layout='circular')
p1 <- ggtree(tree_gene1,branch.length='rate')+geom_treescale(x=0.1,y=250)
p2 <- p1%<+%df_treegene+
  #geom_tiplab(size=.5)+
  geom_tiplab(aes(label=mod_species),size=.5)+
  geom_tippoint(aes(color=as.factor(clade)),size=.3,alpha=.5)  
save_plot('sec4_large_fabv_genetree.png',p2,base_height=16,base_aspect_ratio = .8,dpi=1000)
#
#
## Plot the fabv gene tree but with non HGT fabV gene branches collapsed
## code used to get edge positions to help with plotting
## get edge positions
df_edge <- data.frame('edge'=tree_gene1$edge,'order'=1:length(tree_gene1$edge))
colnames(df_edge) <- c('parent','node','edge_num')
p <- ggtree(tree_gene1,ladderize=T)
p1 <- p%<+%df_edge+geom_label(aes(x=branch,label=node))
save_plot('sec4_node_positions.png',p1,base_height=16,base_aspect_ratio = .8,dpi=1000)
#
## Plotting the gene tree now. collapsed nodes are indicated to the left of the code
p1 <- ggtree(tree_gene1,branch.length='rate')+geom_treescale(x=0.1,y=250)
p2 <- p1 %>% collapse(node=683,'max')+geom_point2(aes(subset=(node==683)),shape=NA,fill='green') #flo marginalis, rhod, azoto, synxantha,lurida, paleroina, yama,paacis,tolasii, fluorescens,salo,poae
p2 <- p2 %>% collapse(node=894,'max')+geom_point2(aes(subset=(node==894)),shape=NA,fill='green') #flo aspleniii, agaricii
p2 <- p2 %>% collapse(node=958,'max')+geom_point2(aes(subset=(node==958)),shape=NA,fill='green') #flo chlororaphis
p2 <- p2 %>% collapse(node=902,'max')+geom_point2(aes(subset=(node==902)),shape=NA,fill='green') #flo protegens, sapo 
p2 <- p2 %>% collapse(node=1033,'max')+geom_point2(aes(subset=(node==1033)),shape=NA,fill='green') #otitdis, citronellolis_group, nitroreducens
p2 <- p2 %>% collapse(node=670,'max')+geom_point2(aes(subset=(node==670)),shape=NA,fill='blue') #aeruginosa
p2 <- p2 %>% collapse(node=543,'max')+geom_point2(aes(subset=(node==543)),shape=NA,fill='blue') #putida
p2 <- p2 %>% collapse(node=665,'max')+geom_point2(aes(subset=(node==665)),shape=NA,fill='blue') #some straminea
#
save_plot('sec4_collapsed_fabv_tree.png',p2,base_height=6,base_aspect_ratio = .6,dpi=300,bg='white')




















