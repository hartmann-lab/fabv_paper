'''
The list of genomes that were used for analysis in the fabv paper will be read and used to filter the fasta.gz files from NCBI. These will then be processed
by checkM at a later point. 
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

#### Section 0 ####
#### Downloaded all pseudomonas assemblies from NCBI. Keeping only the genomes that passed quality testing
#### in the first pass of the paper. Uploaded them to quest and then unzipped them.

path_i = '/Users/Owlex/Downloads/assemblies'
path_dfs = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/fabv_data2'
path_o = '/Users/Owlex/Downloads/assemblieskeep'

## Using the cleaned_filenames_pass.csv dataframe to filter the processed_genomes_metadata.csv
## This is becaused processed_genomes_metadata.csv contains the filenames

df_list = pd.read_csv(pjoin(path_df,'cleaned_filenames_pass.csv')) #passed quality testing
df_meta = pd.read_csv(pjoin(path_df,'processed_genomes_metadata.csv')) #matches cleaned_filename from study to filename from ncbi
# dropping the dust and hospital study genomes
df_meta = df_meta[~df_meta['source'].str.contains('USER')]
# subset for columns to keep and then mergge
df_meta = df_meta[['cleaned_filename','filename']]
df_meta = df_meta.merge(df_list,on='cleaned_filename')

# moving files
df_meta['gzname'] = [i+'.gz' for i in df_meta['filename']]

for i in df_meta['gzname'].tolist():
	os.system('cp %s %s'%(pjoin(path_i,i), pjoin(path_o,i)))

df_meta.to_csv(pjoin('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/fabv_data2','filtered_processed_genomes_metadata.csv'),index=False)	

# Uploaded the files from /Users/Owlex/Downloads/assemblieskeep to /projects/b1042/HartmannLab/alex/fabv_data2/assemblieskeep and unzipped

#### Section 1 ####
#### Preparing files from section 0. Need to make sure all genomes used in study are present.
#### Marking genomes from the original study that were not found in the newly downloaded genome set
#### Using the dataframe filtered_processed_genomes_metadata.csv generated in section 0.
#### keep repeating code in this section until i've found and confirmed the presence of all study genomes

path_main = '/projects/b1042/HartmannLab/alex/fabv_data2'
path_i = '/projects/b1042/HartmannLab/alex/fabv_data2/assemblies'

df_meta = pd.read_csv(pjoin(path_main,'filtered_processed_genomes_metadata.csv'))

# Marking genomes from the original study that were not found in the newly downloaded genome set
gfile_present = [1 if gfile in os.listdir(path_i) else 0 for gfile in df_meta['filename'].tolist()]
df_meta['gfile_present'] = gfile_present
# Writing text file to use for entrez dowload
df_subset = df_meta[df_meta['gfile_present']==0]
with open(pjoin(path_main,'entrez_submit.txt'),'w') as outfile:
	for gfile in df_subset['filename'].tolist():
		outfile.write(gfile[:gfile.find('.')+2]+'\n')

#### Section 2 ####
#### Preparing files from section 0 for chunked submission to quest to run checkm on them.
#### The structure of inputs for checkm is an input file directory, not files individually
#### Using the dataframe filtered_processed_genomes_metadata.csv generated in section 0.
#### Write a slurm submission script for each chunked submission
#### python2

#### This section was re-used to add another 1,000 genomes to checkM. These genomes passed contig and n50 filters

# path_main = '/projects/b1042/HartmannLab/alex/fabv_data2'
# path_i = '/projects/b1042/HartmannLab/alex/fabv_data2/assemblies'
# path_script = '/projects/b1042/HartmannLab/alex/fabv_data2/scripts'

# df_meta = pd.read_csv(pjoin(path_main,'filtered_processed_genomes_metadata.csv'))


path_main = '/projects/b1042/HartmannLab/alex/fabv_data3/large_analysis/checkm'
path_i = '/projects/b1042/HartmannLab/alex/fabv_data3/assemblies'
path_script = '/projects/b1042/HartmannLab/alex/fabv_data3/scripts'

df_meta = pd.read_csv(pjoin(path_main,'refseq_checkm_2.csv'))

# assigning every 500 genomes to a new processing folder
x = np.arange(0,(len(df_meta)/500)+1) # number of bin_id given number of gfiles per bin_id
x = np.repeat(x,500)
df_meta['bin_id'] = x[:len(df_meta)]
df_meta['bin_id'] = df_meta['bin_id'].astype(int).astype(str)

# making new folders for each processing folder group
[os.mkdir(pjoin(path_i,str(i)+'_group')) for i in df_meta['bin_id'].unique().tolist()]

# Copying files from assembly to the subfolder used for checkM
for b_id in df_meta['bin_id'].unique().tolist():
	df_subset = df_meta[df_meta['bin_id']==b_id]
	print ('bid' , b_id, '  ', len(df_subset))
	for gfile in df_subset['filename'].tolist():
		os.system('cp %s %s'%(pjoin(path_i,gfile), pjoin(path_i,str(b_id)+'_group',gfile)))


time_i = 3
for b_id in df_meta['bin_id'].unique().tolist():
	with open(pjoin(path_script,'%s_group.sh'%b_id),'w') as outfile:
		outfile.write('#!/bin/bash\n')
		outfile.write('#SBATCH -A b1042\n')
		outfile.write('#SBATCH -p genomics\n')
		outfile.write('#SBATCH -N 1\n')
		outfile.write('#SBATCH -n 24\n')
		outfile.write('#SBATCH -t %s:00:00\n'%time_i)
		outfile.write('#SBATCH --mem=0\n')
		outfile.write('#SBATCH --error=bin_%s.err\n'%b_id)
		outfile.write('#SBATCH --output=bin_%s.log\n'%b_id)
		outfile.write('#SBATCH	--job-name="bin_%s"\n'%b_id)
		outfile.write('#SBATCH --mail-type=BEGIN,END,FAIL\n')
		outfile.write('#SBATCH --mail-user=alexandermcfarland2022@u.northwestern.edu\n')
		outfile.write('module load python/anaconda checkm\n')
		outfile.write('cd /projects/b1042/HartmannLab/alex/fabv_data3/assemblies/%s_group\n'%b_id)
		outfile.write('checkm taxonomy_wf -t 24 -x fna genus Pseudomonas . checkm_output -f checkm_output_%s.txt --tab_table\n'%b_id)

## Submit bash scripts to cluster
#for f in *_group*; do sbatch $f; done

#### Section 3 ####
#### Compile checkM outputs from section 2 and merge with metadata assignments. might have to be done in groups

path_main = '/projects/b1042/HartmannLab/alex/fabv_data3'
checkm_o = [i for i in os.listdir(path_main) if i.find('checkm_output')>-1]
df_full = pd.DataFrame([])
for i in checkm_o:
	df_temp = pd.read_csv(pjoin(path_main,i),sep='\t')
	df_full= df_full.append(df_temp,ignore_index=True)
df_full.to_csv(pjoin(path_main,'refseq_checkm_complete.csv'),index=False)



















