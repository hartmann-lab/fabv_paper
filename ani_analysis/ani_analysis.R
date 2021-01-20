###
'ani_analysis.R contains the code used to process and manipulate ANI data '
###

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



main_path <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/manuscript/ismej_submission/github/fabv_paper/ani_analysis'


#### Section 0 ####
#### ANI analysis of 160 type strains, 30 study isolates, and C japonicus outgorup
#### Outputs are 4 data tables with the folliwing: all pair-wise ANI, highest pair-wise ANI for study isolate genomes vs type strain genomes, type strain genomes vs type strain genomes, and study isolate genomes vs study isolate genomes 
section_0_path <- paste(main_path,'section_0',sep='/')
setwd(section_0_path)
#
## Read in genome metadata and the mlsa_ani1.txt output from fastANI
df_meta <- read.csv('processed_genomes_metadata.csv',colClasses=c("character"))
df_sp <- read.csv('mlsa_ani1.txt',sep='\t',header=FALSE,colClasses=c("character"))
## Clean up column information
df_sp$V1 <- str_remove(as.character(df_sp$V1),'/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/assemblies/')
df_sp$V2 <- str_remove(as.character(df_sp$V2),'/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/assemblies/')
## need to match query and subject metadata 
df_sp <- merge(df_meta,df_sp,by.x='filename',by.y='V1')
df_sp <- df_sp%>%rename(speciesIDsubj = speciesID)%>%rename(filenameSubj=filename)
df_sp <- merge(df_meta,df_sp,by.x='filename',by.y='V2')
df_sp <- df_sp%>%rename(speciesIDquery = speciesID)%>%rename(filenameQuery=filename)
#
df_sp <- df_sp%>%rename(ani=V3)%>%select(speciesIDquery,speciesIDsubj,ani)
df_sp$ani <- as.numeric(df_sp$ani) 
#
## Each type strain filename is tied to a unique species ID
## Will now replace the filename with species ID
## The final df_sp will display all the hits '(SpeciesIDquery') against that match against a specific query ('SpeciesIDsubject'). 
## Basically, the pair-wise ANI of all combinations of subject and query
#
df_meta1 <- df_meta%>%select(speciesID,source)
df_sp <- merge(df_sp,df_meta1,by.x='speciesIDsubj',by.y='speciesID')
df_sp <- df_sp%>%rename(source_subj=source)
df_meta1 <- df_meta%>%select(speciesID,source)
df_sp <- merge(df_sp,df_meta1,by.x='speciesIDquery',by.y='speciesID')
df_sp <- df_sp%>%rename(source_query=source)
#
df_sp <- df_sp%>%filter(speciesIDsubj!=speciesIDquery) ## output 1
#
## extracting the highest matching ani for study isolate genomes against refseq genomes
df_sp1 <- df_sp%>%filter(source_subj=='USER')%>%filter(source_query=='REFSEQ')%>%group_by(speciesIDsubj)%>%filter(ani==max(ani)) ## output2 
write.csv(df_sp1,'study_type_ani.csv',row.names=FALSE)
## extracting refseq genomes with >=95 ANI with other refseq genomes
df_sp2 <- df_sp%>%filter(source_subj=='REFSEQ')%>%filter(source_query=='REFSEQ')%>%filter(ani>=95) ## output 3
## extracting the highest matching ani for study isolate genomes against other study isolate genomes
df_sp3 <- df_sp%>%filter(source_subj=='USER')%>%filter(source_query=='USER')%>%group_by(speciesIDsubj)%>%filter(ani==max(ani)) ## output 4


#### Section 1 ####
#### ANI analysis of 7,164 genomes used for expanded fabv analysis
#### Output is a grouping of non type strain genomes with type strain genomes according to ANI. Essentially assigning each genome a species group based on ANI threshold to best matching type strain
#### In instances where type strains have a high enough ANI to be the same species group, a species group'_group' suffix is added.
section_1_path <- paste(main_path,'section_1',sep='/')
setwd(section_1_path)
#
#
## read in all 160 type strain genomes. these will match the query column from the compiled fastANI results
#
df_rep <- read.csv('processed_genomes_metadata.csv')%>%filter(source=='REFSEQ')%>%
  select(filename,species,speciesID)
df_ani <-  read.csv('compiled_fastani_results.csv')
df_ani <- merge(df_ani,df_rep,by.x='query',by.y='filename')
#
## Remove subject genomes from fastANI results that did not pass checkM or QUAST cutoffs. This uses df_compiled that has a modified processed_genomes_metadata list of genomes
df_compiled <- read.csv('filtered_processed_genomes_metadata.csv')
df_ani <- df_ani%>%filter(subject%in%df_compiled$filename)
nrow(df_ani)
#
## make a new dataframe containing query to query ani results with identical columns to df_ani (pair-wise self ANI)
df_queryani <- df_ani%>%select(query,species,speciesID)%>%unique()%>%mutate(subject=query)%>%mutate(ani=100)%>%select(query,subject,ani,species,speciesID)
## Merge data sets
df_ani <- df_ani%>%select(query,subject,ani,species,speciesID)
df_ani <- rbind(df_ani,df_queryani)
#
## Phylogenetc context added: add clade designations to df_ani using those from pseudomonas_groups1.csv This uses the same metdata used for clade designations for the mlsa analysis
df_taxassignments <- read.csv('pseudomonas_groups1.csv',colClasses=c("character"))%>%select(species,clade,subclade)
df_taxassignments[is.na(df_taxassignments)] <- 'none'
df_ani <- merge(df_ani,df_taxassignments,by.x='speciesID',by.y='species')
## All type strain that have >=95% ANI to another type strain are merged into a species group_group
df_ani <- df_ani%>%mutate(mod_species=case_when(species%in%c('alcaliphila','sihuiensis','chengduensis','toyotomiensis')~'alcaliphila_group',
                                                species%in%c('amygdali','savastanoi','meliae','ficuserectae')~'amygdali_group',
                                                species%in%c('asplenii','fuscovaginae')~'asplenii_group', 
                                                species%in%c('brassicacearum','kilonensis')~'brassicacearum_group',
                                                species%in%c('chloritidismutans','kunmingensis')~'chloritidismutans_group',
                                                species%in%c('citronellolis','delhiensis')~'citronellolis_group',
                                                species%in%c('synxantha','libanensis')~'synxantha_group',
                                                species%in%c('luteola','zeshuii')~'luteola_group',
                                                species%in%c('oleovorans','pseudoalcaligenes')~'oleovorans_group',
                                                species%in%c('oryzihabitans','psychrotolerans')~'oryzihabitans_group',TRUE~(as.character(species))))
#
## Now matching all non type strain genomes to the best matching type strain ANI
df_ani <- df_ani%>%select(mod_species,ani,subject,clade,subclade)%>%group_by(mod_species,subject)%>%filter(ani==max(ani))%>%ungroup()
nrow(df_ani) # should result in a decrease of rows
# writing to file
write.csv(df_ani,'species_assignments.csv',row.names=FALSE)



#### Section 2 ####
#### Produce dataframe 'reorderd_speciesassignments.csv' which contains the plotting order used in phyloenetic_analysis.R, section 4
section_2_path <- paste(main_path,'section_2',sep='/')
setwd(section_2_path)

## Producing 'reorderd_speciesassignments.csv', which contains the plotting order used in phyloenetic_analysis.R, section 4
df_ani <- read.csv('species_assignments.csv')
df_pani <- df_ani%>%group_by(mod_species)%>%mutate(species_count=length(subject))%>%ungroup()%>%
  select(mod_species,species_count,clade,subclade)%>%unique()
# assigning clade plot ordering based on mlsa tree order:
df_cladeorder <- data.frame('clade'=c('fluorescens','syringae','lutea',
                                      'putida','anguilliseptica','straminea',
                                      'oleovorans','stutzeri','aeruginosa',
                                      'oryzihabitans','resinovorans',
                                      'linyingensis','pertucinogena','none'))
df_cladeorder$cladeorder <- seq(1:nrow(df_cladeorder))
df_pani <- merge(df_pani,df_cladeorder,by='clade')
df_pani <- df_pani%>%arrange(desc(cladeorder,species_count))%>%
  mutate(plot_order=seq(1:length(mod_species)))
df_pani <- df_pani%>%mutate(clade2=paste(clade,subclade,sep=''))%>%select(-plot_order)
df_pani[is.na(df_pani)] <- 0
df_pani <- df_pani%>%arrange(desc(cladeorder),desc(subclade),species_count)%>%mutate(plot_order=seq(1:length(mod_species)))
## Write reordered taxonomy to file
write.csv(df_pani,'reorderd_speciesassignments.csv',row.names=FALSE)




