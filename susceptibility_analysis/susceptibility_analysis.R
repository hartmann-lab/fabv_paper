###
'susceptiblity_analysis.R contains the code used to manipulate and analyze OD values from TCS susceptibility testing'
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
#library(egg)
library(ggpubr)
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

#abs(100*(df[column] - od_gc)/(od_gc-od_sc))

main_path <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/manuscript/genome_biology_submission/github/susceptibility_analysis'




#### Section 0 ####
#### Prepare metadata for study isolates ####
section_0_path <- paste(main_path,'section_0',sep='/')
setwd(section_0_path)
df_taxass<- read.csv('taxonomic_assignments.csv')%>%select(cleaned_filename,species)
#removing 205_ from cleaned_filename so it matches the sample column from df
df_taxass <- df_taxass%>%filter(str_detect(as.character(cleaned_filename),'205_'))%>%mutate(sample=str_replace(cleaned_filename,'205_',''))
#add fabv discovery data
df_fabv <- read.csv('pao1_fabv_homologs.csv')%>%select(cleaned_filename,pident)
#merge df_taxass and df_fabv
df_m <- merge(df_taxass,df_fabv,by='cleaned_filename',all.x=TRUE)
#create fabv presence absence by replacing na with fabi and numeric with fabv
df_m <- df_m%>%mutate(ENR=case_when(is.na(pident)==TRUE~'No FabV',is.numeric(pident)==TRUE~'FabV',TRUE~'NA'))%>%select(sample,ENR,species)
df_reference <- data.frame('sample'=c('KT2440','PAO1'),'ENR'=c('FabV'),'species'=c('putida','aeruginosa'))
df_m <- rbind(df_m,df_reference)
df_m <- df_m%>%rename(isolate=sample)%>%mutate(origin=case_when(isolate%in%c('HS_1','HS_2','HS_3','HS_4','HS_5')~'hospital',isolate%in%c('PAO1','KT2440')~'reference_strain',TRUE~'study_isolate'))
## Add metadata from phylogenetic analysis with clade info
df_clades <- data.frame('isolate'=c(
  '96A1','HS_2','HS_1','114A4',
  '45C2',
  '31A8','56A10','82B1','HS_5','88B1','34A1','4A7','6C6','20A1','119A3','66C3','KT2440',
  '109A1','69C1','57B2','99A1','95A6','62A4','8A1','115A1',
  'HS_4','HS_3','PAO1',
  '97C1','39A1','10A6','89C1'),
  'clade_assignment'=c(
    'fluorescens','fluorescens','fluorescens','fluorescens',
    'viridiflava',
    'putida','putida','putida','putida','putida','putida','putida','putida','putida','putida','putida','putida',
    'stutzeri','stutzeri','stutzeri','stutzeri','stutzeri','stutzeri','stutzeri','stutzeri',
    'aeruginosa','aeruginosa','aeruginosa',
    'oryzihabitans','oryzihabitans','oryzihabitans','oryzihabitans'))
df_m <- merge(df_m,df_clades,by='isolate')
write.csv(df_m,'isolate_metadata.csv',row.names=FALSE)


#### Section 1 ####
#### Visualize od curves ####
section_1_path <- paste(main_path,'section_1',sep='/')
setwd(section_1_path)
## read in OD data
df_sus <-  read.csv('all_isolate_od.csv')
## standardize decimals from the variable column. for example 0.25 and .25 are the same but are interpreted differently
df_sus <- df_sus%>%filter(variable!='SC')%>%mutate(mod_variable=case_when(variable=='GC'~'0',TRUE~variable))%>%
  mutate(mod_variable=as.numeric(mod_variable))

## calculate the mean, min, max, std error of 
df_sus <- df_sus%>%ungroup()%>%
  group_by(isolate,treatment,mod_variable)%>%
  mutate(mean_value=mean(value),min_value=min(value),max_value=max(value),stderr_value = sd(value)/sqrt(length(value)))%>%
  select(isolate,treatment,mod_variable,mean_value,min_value,max_value,stderr_value)%>%unique()
#@ add metadata
df_meta <- read.csv('isolate_metadata.csv')
df_sus <- merge(df_sus,df_meta,by='isolate')
#
plot_od_susceptibility <- function(input_df, input_title){
  plotout <- ggplot(input_df,aes(as.factor(mod_variable),y=mean_value))+
    #geom_text(aes(label=isolate))+
    geom_line(aes(group = isolate, color = origin, linetype = ENR), size = 1.5, alpha = .8)+
    geom_errorbar(aes(ymin = mean_value - stderr_value, ymax = mean_value + stderr_value), alpha = .1, width = .1)+
    geom_point(size = 1,alpha=.3)+
    theme(axis.text.y = element_text(size = 24, color = "black"),
          axis.text.x = element_text(size = 16, angle = 0, color = "black"),
          axis.title.y = element_text(size = 27),
          axis.title.x = element_text(size = 27),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(size = 25),
          #legend.title = element_text(size = 25, face = "bold"),
          legend.position = "none",
          #legend.position=c(.9,.9), #x,y
          legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25,vjust=-1.5))+
    scale_x_discrete(limits=c('0','0.0625','0.125','0.25','0.5','1','2','4','8','16','32','64','128'),
                     labels=c('0','0.0625    ','0.125','0.25','0.5','1','2','4','8','16','32','64','128'))+
    scale_y_continuous(limits=c(-0.05,1.6),breaks=c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "red", size = 1)+
    scale_color_manual(values=c('darkgray','#E69F00','#2b8cbe'))+
    ylab(label = "OD")+
    xlab(label = "Triclosan (mg/L)")+
    ggtitle(input_title)
  return(plotout)}
#
p1 <- plot_od_susceptibility(input_df=df_sus%>%filter(treatment=='con'),input_title='Growth inhibition')
p2 <- plot_od_susceptibility(input_df=df_sus%>%filter(treatment=='inh'), input_title='Growth inhibition with 40 µg/ml PAßN')
##
t <- ggarrange(p1,p2,nrow=1,ncol=2)
save_plot("sec0_od_curves.png", t, base_height = 10, base_aspect_ratio = 2)
 

#### T test (wilcoxon signed-rank test) of experimental vs growth control od for instances where expermental od>growth control od ####
## read in susceptibility info
df_mic <-  read.csv('all_isolate_od.csv')
## because variable column cotaining the treatment group has both numeric and characters, not all numbers are formattd the same. like 0.25 and .25. so I need to standardize this.
df_mic$variable <- as.character(df_mic$variable)
df_mic <- df_mic%>%mutate(numeric_variable=case_when(variable%in%c('SC','GC')~'0',TRUE~variable))
df_mic$numeric_variable <- as.numeric(df_mic$numeric_variable)
df_mic$numeric_variable <- as.character(df_mic$numeric_variable)
df_mic <- df_mic%>%mutate(variable_copy=variable)%>%mutate(variable=case_when(numeric_variable=='0'~variable,TRUE~numeric_variable))
## specifiy growth control for each sample
df_mic  <- df_mic%>%mutate(identifier=paste(isolate,treatment,replicate,sep='_'))
df_gc <- df_mic%>%filter(variable=='GC')%>%select(identifier,value)%>%rename(growth_control_od=value)
df_mic <- merge(df_mic,df_gc,by='identifier')%>%select(-identifier)
## subset all instances where an experimental value is greater than a growth control od
df_outofrange <- df_mic%>%filter(value>growth_control_od)%>%mutate(od_range=value-growth_control_od)
## calculate the wilcox signed rank test on the experimental od vs growth control for each sample
df_outofrange_summary <- df_outofrange%>%group_by(isolate,treatment,variable)%>%mutate(count=length(value))%>%
  mutate(wilcox_p_val=wilcox.test(unlist(value),unlist(growth_control_od),paired=TRUE)$p.value)%>%
  mutate(wilcox_statistic=wilcox.test(unlist(value),unlist(growth_control_od),paired=TRUE)$statistic)%>%
  ungroup()%>%
  select(isolate,treatment,variable,count,wilcox_p_val,wilcox_statistic)%>%unique()%>%
  mutate(identifier=paste(isolate,treatment,variable,sep='_'))%>%
  select(identifier,wilcox_p_val,wilcox_statistic)
#
## merge statistic info with df_mic
df_mic <- df_mic%>%mutate(identifier=paste(isolate,treatment,variable,sep='_'))
df_mic <- merge(df_mic,df_outofrange_summary,by='identifier',all.x=TRUE)
#
## for all instances where experimental od > growth control od, check whether p value is >0.05
## if this is true, that means the experimental od and growth control are not distinguishable. Set the experimental od equal to the growth control od
df_mic <- df_mic%>%mutate(value=abs(value))%>%mutate(od_value=case_when((value>growth_control_od & wilcox_p_val>0.05)~growth_control_od,
                                             TRUE~value))


#### calculate group percentage growth inhibition ####
df_mic <- df_mic%>%group_by(identifier)%>%mutate(mean_growth_control_od=mean(growth_control_od),mean_od_value=mean(od_value))%>%
  select(isolate,identifier,treatment,variable,mean_growth_control_od,mean_od_value)%>%unique()%>%
  mutate(group_growth_inhibition=(100*(mean_growth_control_od-mean_od_value)/mean_growth_control_od))

#
df_mic_summarize <- df_mic%>%group_by(identifier)%>%
  mutate(mean_perc_inh = mean(group_growth_inhibition))%>%
  ungroup()%>%
  select(identifier,isolate,treatment,variable,mean_perc_inh)%>%unique()
## modify variable column for plotting
df_mic_summarize <- df_mic_summarize%>%filter(!variable%in%c('GC','SC'))
df_mic_summarize$variable <- as.numeric(df_mic_summarize$variable)
## add metadata
df_meta <- read.csv('isolate_metadata.csv')
df_mic_summarize <- merge(df_mic_summarize,df_meta,by='isolate')

## plot growth inhibition curve
plot_inhibition_curve <- function(input_df, input_title){
  plotout <- ggplot(input_df,aes(as.factor(variable),y=mean_perc_inh))+
    #geom_text(aes(label=isolate))+
    geom_line(aes(group = isolate, color = origin, linetype = ENR), size = 1.5, alpha = .8)+
    geom_point(size = 1,alpha=.3)+
    theme(axis.text.y = element_text(size = 24, color = "black"),
          axis.text.x = element_text(size = 16, angle = 0, color = "black"),
          axis.title.y = element_text(size = 27),
          axis.title.x = element_text(size = 27),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(size = 25),
          #legend.title = element_text(size = 25, face = "bold"),
          legend.position = "none",
          #legend.position=c(.9,.9), #x,y
          legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25,vjust=-1.5))+
    scale_x_discrete(limits=c('0.0625','0.125','0.25','0.5','1','2','4','8','16','32','64','128'),
                     labels=c('0.0625    ','0.125','0.25','0.5','1','2','4','8','16','32','64','128'))+
    scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
    #geom_hline(yintercept = 0, linetype = "dotted", color = "red", size = 1)+
    geom_hline(yintercept = 95, linetype = "dotted", color = "red", size = 1)+
    scale_color_manual(values=c('darkgray','#E69F00','#2b8cbe'))+
    ylab(label = "Relative growth inhibition (%)")+
    xlab(label = "Triclosan (mg/L)")+
    ggtitle(input_title)
  return(plotout)}

## isolates '99A1','119A3','66C3','62A4' are removed from inhibitor curve because their growth control was strongly affected by the presence of the inhibitor itself
p1 <- plot_inhibition_curve(input_df=df_mic_summarize%>%filter(treatment=='con'),input_title='Growth inhibition')
p2 <- plot_inhibition_curve(input_df=df_mic_summarize%>%filter(treatment=='inh')%>%filter(!isolate%in%c('99A1','119A3','66C3','62A4')),input_title='Growth inhibition with 40 µg/mL PAβN')
t <- ggarrange(p1,p2,nrow=1,ncol=2)
save_plot("sec0_growth_inhibitions.png", t, base_height = 10, base_aspect_ratio = 2)


#### Extract MIC ####
## isolates '99A1','119A3','66C3','62A4' are removed from inhibitor MIC calculation because their growth control was strongly affected by the presence of the inhibitor itself
df_mic1 <- df_mic%>%filter(!variable%in%c('GC','SC'))%>%mutate(variable=as.numeric(variable))%>%filter(!isolate%in%c('99A1','119A3','66C3','62A4')|treatment!='inh')
## make two df, one with no MIC detected and one with MIC detected
df_noinh <- df_mic1%>%mutate(inhibited=case_when(group_growth_inhibition>=95~1,TRUE~0))%>%
  group_by(isolate,treatment)%>%mutate(no_inhibition=case_when(length(unique(inhibited))>1~0,TRUE~1))%>%filter(inhibited==0&no_inhibition==1)%>%ungroup()%>%
  group_by(isolate,treatment)%>%mutate(highest_mic=case_when(variable==max(variable)~1,TRUE~0))%>%filter(highest_mic==1)%>%select(-no_inhibition,-inhibited,-highest_mic)
## df_inh finds instances where perc_inh was greater than or equal to 95 and then takes the lowest value where this occurred as the MIC
df_inh <- df_mic1%>%mutate(inhibited=case_when(group_growth_inhibition>=95~1,TRUE~0))%>%filter(inhibited==1)%>%
  group_by(isolate,treatment)%>%mutate(lowest_inhibited=case_when(variable==min(variable)~1,TRUE~0))%>%
  filter(lowest_inhibited==1)%>%select(-inhibited,-lowest_inhibited)
df_mic2 <- rbind(df_noinh,df_inh)
## defining whether MICs are within the tested ranges
df_mic2 <- df_mic2%>%mutate(range=case_when(variable=='128'&group_growth_inhibition<95~'greater_than',variable=='0.0625'&group_growth_inhibition>=95~'less_than',TRUE~'within_range'))
#
## Add metadata for fabv p/a
df_m <- read.csv('isolate_metadata.csv')
df_mic2 <- merge(df_mic2,df_m,by='isolate')

plot_mic_values <- function(input_df,input_title){
  #
  plot_out <- ggplot(input_df,aes(x=ENR,y=variable))+
    geom_boxplot(outlier.fill=NULL,outlier.size=0,outlier.alpha=0,fill=NA)+
    geom_point(aes(fill=as.factor(clade_assignment)),shape=21,position=position_jitterdodge(jitter.width=1,dodge.width=.8),size=6,alpha=.8)+
    theme_classic()+
    theme(axis.text.y = element_text(size = 24, color = "black"),
          axis.text.x = element_text(size = 40, angle = 0, color = "black"),
          axis.title.y = element_text(size = 27),
          axis.title.x = element_text(size = 27),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25, face = "bold"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 25),
          panel.grid.major.y =element_line(size=1,linetype='dotted',colour='lightblue'))+
    scale_fill_brewer(palette='Dark2')+
    scale_y_continuous(breaks=c(0.0625,8,16,32,64,128),labels=c(0.0625,8,16,32,64,128))+
    scale_x_discrete(labels=c('+','-'))+
    ylab(label = "Triclosan (mg/L)")+
    xlab(label = expression(italic('fabV')~'detected in genome'))+
    ggtitle(input_title)
  }
#riclosan MIC with PAβN
p1 <- plot_mic_values(input_df=df_mic2%>%filter(treatment=='con'),input_title='Triclosan MIC alone')
p2 <- plot_mic_values(input_df=df_mic2%>%filter(treatment=='inh'),input_title='Triclosan MIC with 40 µg/ml PAβN')

#
t3 <- ggarrange(p1,p2,nrow=1,ncol=2)#,align='h',widths=c(1,1))
save_plot("sec0_MIC.png", t3, base_height = 10, base_aspect_ratio = 2)
write.csv(df_mic2,'compiled_mic_data.csv',row.names=FALSE)

#### Calculate median triclosan MICs ####
fabv_no_inhibitor <- df_mic2%>%select(isolate,treatment,variable,ENR,mean_od_value)%>%filter(treatment=='con'&ENR=='FabV')
fabv_with_inhibitor <- df_mic2%>%select(treatment,variable,ENR,mean_od_value)%>%filter(treatment=='inh'&ENR=='FabV')
nofabv_no_inhibitor <- df_mic2%>%select(treatment,variable,ENR,mean_od_value)%>%filter(treatment=='con'&ENR=='No FabV')
nofabv_with_inhibitor <- df_mic2%>%select(treatment,variable,ENR,mean_od_value)%>%filter(treatment=='inh'&ENR=='No FabV')

median(fabv_no_inhibitor$variable)
median(fabv_with_inhibitor$variable)
median(nofabv_no_inhibitor$variable)
median(nofabv_with_inhibitor$variable)

#### Calculate T-test for fabv genotype and triclosan MIC ####
#### Using wilcoxon rank-sum/mann-whitney U 
## fabv vs no fabv, no inhibitor present
wilcox.test(fabv_no_inhibitor$variable,nofabv_no_inhibitor$variable,paired=FALSE)
## fabv vs no fabv, no inhibitor present
wilcox.test(fabv_with_inhibitor$variable,nofabv_with_inhibitor$variable,paired=FALSE)

#### Calculate T-test for fabv genotype and triclosan growth inhibition %
#### Using wilcoxon rank-sum/mann-whitney U
wilcox.test(fabv_no_inhibitor$variable,nofabv_no_inhibitor$mean_od_value,paired=FALSE)
## fabv vs no fabv, no inhibitor present
wilcox.test(fabv_with_inhibitor$variable,nofabv_with_inhibitor$mean,paired=FALSE)




#### create data table with isolate info ####
df_com <- read.csv('compiled_mic_data.csv')
df_com <- df_com%>%mutate(range2=case_when(range=='greater_than'~'>',range=='less_than'~'≤',TRUE~''))%>%
  mutate(mic_value=paste(range2,as.character(variable),sep=''))%>%
  select(isolate,treatment,mic_value,clade_assignment,ENR)
## convert to long format
df_com <- dcast(df_com,isolate+clade_assignment+ENR~treatment,value.var='mic_value')
df_com <- df_com%>%rename(triclosan_MIC=con,triclosan_MIC_w_inhibitor=inh)%>%select(isolate,ENR,triclosan_MIC,triclosan_MIC_w_inhibitor,clade_assignment)
## order by isolate. dust isolates, hospital isolates, reference strains

isolate_order <- data.frame('isolate'=c('4A7','6C6','8A1','10A6','20A1','31A8','34A1','39A1','45C2','56A10','57B2','62A4','66C3','69C1','82B1','88B1','89C1','95A6','96A1','97C1','99A1','109A1','114A4','115A1','119A3','HS_1',
                                        'HS_2','HS_3','HS_4','HS_5','PAO1','KT2440'),'order'=seq(1,32))
#
df_com <- merge(df_com,isolate_order,by='isolate')
df_com <- df_com%>%arrange(order)%>%select(-order)%>%mutate(isolate=str_replace(isolate,'_',''),clade_assignment=paste('P. ',clade_assignment,sep=''))
write.csv(df_com,'isolate_mic_final.csv',row.names=FALSE)



