
'''
Phylogeny scripts and functions for MLSA sequence analysis and comparision with fabv sequence
'''

import os
from os.path import join as pjoin

import re

from multiprocessing import Pool
import subprocess

from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

import pandas as pd 
import numpy as np

import textwrap
from timeit import default_timer as timer


def print_HitStats(input_df):
	print('--------------------------HIT STATS--------------------------------------------')
	print('gene: ', input_df['qseqid'].unique().item(),'\n')
	print('genomes with hits: ', len(input_df['cleaned_filename'].unique()))
	if len(input_df['cleaned_filename'].unique()) < 191:
		print('Not all genomes have hits!\n\n')
	print('hits per genome ', )
	print('max: ', input_df.groupby('cleaned_filename').count()['speciesID'].max())
	print('min: ', input_df.groupby('cleaned_filename').count()['speciesID'].mean())
	print('mean: ', input_df.groupby('cleaned_filename').count()['speciesID'].min())
	print('quantile: \n', input_df.groupby('cleaned_filename').count()['speciesID'].quantile([.25,.5,.75]),'\n')
	print('qcovs: ')
	print('max: ', input_df.groupby(['cleaned_filename'])['qcovs'].transform(max).max())
	print('min: ', input_df.groupby(['cleaned_filename'])['qcovs'].transform(min).min())
	print('mean: ', input_df.groupby(['cleaned_filename'])['qcovs'].transform(max).mean())
	print('quantile: \n', input_df.groupby(['cleaned_filename'])['qcovs'].transform(max).quantile([.25,.5,.75]),'\n')
	print('pident')
	print('max: ', input_df.groupby(['cleaned_filename'])['pident'].transform(max).max())
	print('min: ', input_df.groupby(['cleaned_filename'])['pident'].transform(max).min())
	print('mean: ', input_df.groupby(['cleaned_filename'])['pident'].transform(max).mean())
	print('quantile: \n', input_df.groupby(['cleaned_filename'])['pident'].transform(max).quantile([.25,.5,.75]),'\n')
	print('bitscore')
	print('max: ', input_df.groupby(['cleaned_filename'])['bitscore'].transform(max).max())
	print('min: ', input_df.groupby(['cleaned_filename'])['bitscore'].transform(max).min())
	print('mean: ', input_df.groupby(['cleaned_filename'])['bitscore'].transform(max).mean())
	print('quantile: \n', input_df.groupby(['cleaned_filename'])['bitscore'].transform(max).quantile([.25,.5,.75]),'\n')
	print('-------------------------------------------------------------------------------')
	print('\n')


#### Section 0 ####
#### Extract reference sequences from the pao1 type strain faa list.
path_reference = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/mlsa_genes'
df_tracker = pd.read_csv(pjoin(path_reference,'mlsa_tracker.csv'))
#
locus_tracker = 0
for pao1_locus in df_tracker['pao1_locus_tag'].tolist():
	with open(pjoin(path_reference,'pao1_'+pao1_locus+'.txt'),'w') as outfile:
		for record in SeqIO.parse(pjoin(path_reference,'pao1_mlsareference.faa'),'fasta'):
			if record.id == pao1_locus:
				newseq = textwrap.wrap(str(record.seq),60,break_on_hyphens=False)
				outfile.write('>pao1_%s\n'%pao1_locus)
				[outfile.write('%s\n'%i) for i in newseq]
				print('MESSAGE: %s sequence found. Writing to file...'%pao1_locus)
				locus_tracker +=1
				break
print('Expected: ',locus_tracker)
print('Actual: ', len(df_tracker['pao1_locus_tag'].tolist()))

#### Section 1 ####
#### edit processed_genomes_metadata.csv so that there are individual columns for sseqid
#### and the cleaned_filename without the taxonomic identifier
#### Afterwards, merge all the hits tables into 1. Search for the best hit in each one and concatenate according
#### to the specific order given. 
i_path = '/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis'
# adding new columns to df_meta
df_meta = pd.read_csv(pjoin(i_path,'processed_genomes_metadata.csv'))
# replacing species names gathered from refseq headers and replacing them with the type_strains.csv identifiers
df_type = pd.read_csv(pjoin(i_path,'type_strains.csv'))
df_meta = df_meta.rename(columns={'species':'refseqSpecies'})
df_type = df_type[['species','id']]
df_meta = df_meta.merge(df_type,on='id',how='left')
# creating the genID and speciesID columns used for fasta delimiters
df_meta['species'] = df_meta['species'].replace(np.nan,'unknown',regex=True)
df_meta['genID'] = [s[s.find('_')+1:] for s in df_meta['cleaned_filename'].tolist()]
df_meta['speciesID'] = np.where(df_meta['species'] == 'unknown',df_meta['genID'],df_meta['species'])
df_meta.to_csv(pjoin(i_path,'processed_genomes_metadata.csv'),index=False)


#### Produce filteredhomolog outputs using GetHomologs.py


#### Section 2 ####
#### Extract sequences from hits table for a GetHomologs.py output
#### 
type_analysis_path = '/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis'
i_path = pjoin(type_analysis_path,'homologsearch_output')
o_path = pjoin(type_analysis_path,'singlegene_aln')
input_df = pjoin(type_analysis_path,'processed_genomes_metadata.csv')
# genename = 'pao1_fabi'
# keep_only_best = True
# species_delim = 'speciesID' #options: genID,speciesID
# gene_delim = 'homologID' #options: homologID
# delim_separator = '_' 
# gene_name_order = 1 #options: 1 = species,gene 2 = gene,species   NOTE: , indicates delimiter separator

def prep_ForAlignment(i_path,o_path,input_df,genename,seq_extension,keep_only_best,species_delim,gene_delim,delim_separator,gene_name_order,sequence_filter):
	## metadata containing species info, selecting the following columns
	df_meta = pd.read_csv(input_df)
	df_meta = df_meta[['cleaned_filename','source','id','species','genID', 'speciesID']]
	## Modifying GetHomologs.py blast output
	df_gene = pd.read_csv(pjoin(i_path,genename+'_homologs.csv'))
	df_gene = df_gene.merge(df_meta,on='cleaned_filename')
	df_gene['homologID'] = [i[i.find('_')+1:] for i in df_gene['sseqid']]
	df_gene['homologID'] = [i[i.find('_')+1:] for i in df_gene['homologID']]
	## Applying filters and reporting general statistics
	if sequence_filter[0] == True:
		print('MESSAGE: Removing hits with less than %s qcovs and %s pident'%(str(sequence_filter[1]),str(sequence_filter[2])))
		df_gene = df_gene[df_gene['qcovs']>=sequence_filter[1]]
		df_gene = df_gene[df_gene['pident']>=sequence_filter[2]]
	if keep_only_best == False:
		print('MESSAGE: Keeping all hits\n')
		print_HitStats(input_df=df_gene)
	if keep_only_best == True:
		df_gene = df_gene[df_gene.groupby(['cleaned_filename','qseqid'])['bitscore'].transform(max) == df_gene['bitscore']]
		if len(df_gene.groupby('cleaned_filename').filter(lambda x: len(x) > 1)) > 1:
			print('WARNING: bitscore filtering not sufficient, filtering by qcovs')
			df_gene = df_gene[df_gene.groupby('cleaned_filename')['qcovs'].transform(max) == df_gene['qcovs']]
			if len(df_gene.groupby('cleaned_filename').filter(lambda x: len(x) > 1)) > 1:
				print('WARNING: qcovs filtering not sufficient')
		print('MESSAGE: Keeping best hits\n')
		print_HitStats(input_df=df_gene)
	## Modifying gene dataframe to include specified delimters and separators
	# selecting the following columns	
	df_gene = df_gene[['sseqid','pident','qcovs','bitscore','qseqid','cleaned_filename','homologID','speciesID','genID']]
	# creating column with inputted delimiter IDs and separator
	if gene_name_order == 1:
		df_gene['inputtedDelimID'] = df_gene[species_delim]+delim_separator+df_gene[gene_delim]
	if gene_name_order == 2:
		df_gene['inputtedDelimID'] = df_gene[gene_delim]+delim_separator+df_gene[species_delim]
	## Exit function	
	# if len(df_gene['cleaned_filename'].unique()) < 191:
	# 	print('Not all genomes have hits, did not write any sequences to file!\n\n')
	# 	return(None)
	## writing selected sequences to file, making a csv copy with the same filename
	hit_search_list = df_gene['sseqid'].tolist()
	with open(pjoin(o_path,genename+'_filteredhomologs'+seq_extension),'w') as outfile:
		for record in SeqIO.parse(pjoin(i_path,genename+'_homologs'+seq_extension),'fasta'):
			if record.id in hit_search_list:
				newrecord = df_gene[df_gene['sseqid']==record.id]['inputtedDelimID'].item()
				newseq = textwrap.wrap(str(record.seq),60,break_on_hyphens=False)
				outfile.write('>%s\n'%newrecord)
				[outfile.write('%s\n'%i) for i in newseq]
	df_gene.to_csv(pjoin(o_path,genename+'_filteredhomologs.csv'),index=False)

# #for gene in ['pao1_fabi','pao1_fabv','']:
# gene_list = ['pao1_PA3297','pao1_PA3308','pao1_PA3257','pao1_PA1005','pao1_PA1294',
# 'pao1_PA2981','pao1_PA2961','pao1_PA0759','pao1_PA2964','pao1_PA3243']

# gene_list = ['pao1_PA3297',
# 'pao1_PA3308',
# 'pao1_PA3257',
# 'pao1_PA1005',
# 'pao1_PA1294',
# 'pao1_PA2981',
# 'pao1_PA2961',
# 'pao1_PA0759',
# 'pao1_PA2964',
# 'pao1_PA3243',
# 'pao1_gyrb',
# 'pao1_rpod',
# 'pao1_rpob']

## check output files and remove those that do not have 

gene_list = ['pao1_fabv']

##### using the mlsa_tracker.csv to generate gene lists
# path_reference = '/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis'
# df_tracker = pd.read_csv(pjoin(path_reference,'mlsa_tracker.csv'))
# gene_list =  ['pao1_%s'%i for i in df_tracker['pao1_locus_tag'].tolist()]

for gene in gene_list:
	prep_ForAlignment(i_path=i_path,
		o_path=o_path,
		input_df=input_df,
		genename=gene,
		seq_extension = '.ffn',#options :.ffn, .faa
		keep_only_best=True, #options: True, False
		species_delim='speciesID', #options: genID,speciesID
		gene_delim='homologID', #options: homologID
		delim_separator='_',
		gene_name_order=1, #options: 1 = species,gene ; 2 = gene,species
		sequence_filter = [False, 90, 50])




#### Align with mafft.
module load mafft 
cd /projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/singlegene_aln
for i in $(ls *.ffn)
do
echo  ${i%.*}
mafft --maxiterate 1000 --localpair --thread -1 ${i%.*}.ffn > ${i%.*}.aln
done

#### optional: trim each sequence with trimal 
module load anaconda3/2018.12
source activate trim 
cd /projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/singlegene_aln
for i in $(ls *_filteredhomologs.aln)
do
echo  ${i%.*}
trimal -in ${i%.*}.aln -out ${i%.*}_trimmed.aln -automated1
done


#### Batch individual trees
# cd /projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/singlegene_aln
# for i in $(ls *_trimmed.aln)
# do
# /home/agm9813/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -m GTRCAT -p 12345 -x 12345 -# 10 -s ${i%.*}.aln -n ${i%.*}
# done




#### Section 3 ####
#### Concatenated MLSA creation for a list of genes from section 1
type_analysis_path = '/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis'
homologsearch_path = pjoin(type_analysis_path,'homologsearch_output')
singlegenealn_path = pjoin(type_analysis_path,'singlegene_aln')
input_df = pjoin(type_analysis_path,'processed_genomes_metadata.csv')
seq_extension = '.aln'


gene_list = [i for i in os.listdir(singlegenealn_path) if i.endswith('trimmed.aln')]
gene_list = [i.replace('_filteredhomologs_trimmed.aln','') for i in gene_list]


# gene_list = ['pao1_PA3297',
# 'pao1_PA3308',
# 'pao1_PA3257',
# 'pao1_PA1005',
# 'pao1_PA1294',
# 'pao1_PA2981',
# 'pao1_PA2961',
# 'pao1_PA0759',
# 'pao1_PA2964',
# 'pao1_PA3243',
# 'pao1_gyrb',
# 'pao1_rpod',
# 'pao1_rpob']

# printing the max,min,mean of the max hit for each cleaned_filename
# for i in gene_list:
# 	df_gene = pd.read_csv(pjoin(o_path,'%s_filteredhomologs.csv'%i))
# 	print_HitStats(input_df=df_gene)


df_meta = pd.read_csv(input_df)
master_genetracker = 0 # check that all genes are being added
master_seqlentracker = 0 # check that all sequence lengths are the same
#
with open(pjoin(singlegenealn_path,'conctat34mlsa.aln'),'w') as outfile:
	for cf in df_meta['cleaned_filename'].tolist():
		# print(cf)
		species_id = df_meta[df_meta['cleaned_filename']==cf]['speciesID'].unique().item()
		genetracker = 0
		concatseq = ''
		#print(species_id)
		for x, gene in enumerate(gene_list):
			for record in SeqIO.parse(pjoin(singlegenealn_path,gene+'_filteredhomologs_trimmed'+seq_extension),'fasta'):
				if record.id.startswith(species_id):
					#print(gene)
					#print(str(record.seq))
					#concatseq  += 'newgene%s'%(x+1)
					#concatseq += gene
					concatseq += str(record.seq)
					genetracker += 1
					break
		if genetracker != len(gene_list):
			print('ERROR: expected %s genes, only found %s genes'%(len(str(gene_list))),str(genetracker))
			break
		print('MESSAGE: %s stats:\n number of genes added: %s\n total sequence length: %s \n' % (species_id,str(genetracker),str(len(concatseq))))
		master_genetracker += genetracker
		master_seqlentracker += len(concatseq)
		concatseq = textwrap.wrap(str(concatseq),60,break_on_hyphens=False)
		outfile.write('>%s\n'%species_id)
		[outfile.write('%s\n'%i) for i in concatseq]
print(master_genetracker/len(df_meta['cleaned_filename'].tolist()))
print(master_seqlentracker/len(df_meta['cleaned_filename'].tolist()))


#### Tree align concatenated sequence 
/home/agm9813/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s conctat34mlsa.aln -n conctat34mlsa












# #### Section 1 ####
# #### Concatenated MLSA creation for a list of genes from a GetHomologs.py output
# i_path = '/projects/b1042/HartmannLab/alex/fabv_data2/typeanalysis'
# prokka_path = pjoin(i_path,'prokka_output')

# gene_list = ['pao1_gyrb','pao1_rpob','pao1_rpod']
# input_delimiter = 'speciesID' #options: speciesID, genID


# # printing the max,min,mean of the max hit for each cleaned_filename
# for i in gene_list:
# 	df_gene = pd.read_csv(pjoin(i_path,'%s_homologs.csv'%i))
# 	print_HitStats(input_df=df_gene)


# ### Adding new ID tags for each genome input: genID, which is the cleaned_filename minus the species number code
# ### and speciesID, which is the species name. for isolates it uses the genID instead
# # format speciesID_geneID
# df_meta = pd.read_csv(pjoin(i_path,'processed_genomes_metadata.csv'))
# # selecting the following columns
# df_meta = df_meta[['cleaned_filename','source','id','species','genID', 'speciesID']]

# # extracting the best hit from each grouping
# df_phy = [pd.read_csv(pjoin(i_path,'%s_homologs.csv'%g)) for g in gene_list]
# df_phy = pd.concat(df_phy)
# df_phy = df_phy[df_phy.groupby(['cleaned_filename','qseqid'])['bitscore'].transform(max) == df_phy['bitscore']]
# # selecting the following columns
# df_phy = df_phy[['sseqid','pident','qcovs','bitscore','qseqid','cleaned_filename']]

# # merge data from df_meta with df_phy
# df_phy = df_phy.merge(df_meta,on='cleaned_filename')

# # check that there is a unique defline ID for each header
# df_phy.groupby(input_delimiter).filter(lambda x: len(x) > len(gene_list))


# # extracting best hits from the corresponding nucleotide sequences generated by prokka
# with open(pjoin(i_path,'mlsa_raw.txt'),'w') as outfile:
# 	for cf in df_phy['cleaned_filename'].unique().tolist():
# 		print(cf)
# 		df_t = df_phy[ df_phy['cleaned_filename']==cf ] 
# 		# add each sseqid to a list called temp_seq that is then turned into a dictionary
# 		temp_seq = []
# 		for record in SeqIO.parse(pjoin(prokka_path,cf+'.ffn'),'fasta'):
# 			if record.id in df_t['sseqid'].tolist():
# 				gene = df_t[df_t['sseqid']==record.id]['qseqid'].item()
# 				temp_seq.append([gene,str(record.seq)])
# 		temp_dict = {i[0]:i[1:] for i in temp_seq}
# 		# searches temp_dict to write the correct sequence in order
# 		concatenated_tempseq = ''
# 		for g in gene_list:
# 			concatenated_tempseq = concatenated_tempseq+temp_dict[g][0]
# 		concatseq = textwrap.wrap(concatenated_tempseq,60,break_on_hyphens=False)
# 		species_id = df_t[input_delimiter].unique().item()
# 		outfile.write('>%s\n'%species_id)
# 		[outfile.write('%s\n'%i) for i in concatseq]
# df_phy.to_csv(pjoin(i_path,'mlsa_raw.csv'),index=False)









































