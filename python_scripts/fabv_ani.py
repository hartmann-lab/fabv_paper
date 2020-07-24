'''
Runs through fastANI by creating slurm submission scripts. Compiles results into a single large dataframe.
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
#### Only doing the mlsa genomes
path_mash = '/home/agm9813/fastANI'
path_out = '/projects/b1042/HartmannLab/alex/fabv_data3/ani_analysis'
path_genome = '/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/assemblies'


with open(pjoin(path_out,'refgenomes.txt'),'w') as outfile:
	for i in os.listdir(path_genome):
		if i.endswith('.fna'):
			outfile.write('%s\n'%pjoin(path_genome,i))
		if i.endswith('.fasta'):
			outfile.write('%s\n'%pjoin(path_genome,i))

with open(pjoin(path_out,'inputgenomes.txt'),'w') as outfile:
	for i in os.listdir(path_genome):
		if i.endswith('.fna'):
			outfile.write('%s\n'%pjoin(path_genome,i))
		if i.endswith('.fasta'):
			outfile.write('%s\n'%pjoin(path_genome,i))


# ./fastANI --ql [QUERY_LIST] --rl [REFERENCE_LIST] -o [OUTPUT_FILE]

querylist = pjoin(path_out,'inputgenomes.txt')
referencelist = pjoin(path_out,'refgenomes.txt')

os.system('%s --ql %s --rl %s -o %s/mlsa_ani.txt -t 24'%(path_mash,querylist,referencelist,path_out))

#-----------------------------------------------------------------------------------------------------------


#### Section 1 ####
#### Setting up to do the mlsa genomes vs the entire search database. 
path_mash = '/home/agm9813/fastANI'
path_script = '/projects/b1042/HartmannLab/alex/fabv_data3/scripts/ani'
path_largeanalysis = '/projects/b1042/HartmannLab/alex/fabv_data3/large_analysis'
path_anioutput = '/projects/b1042/HartmannLab/alex/fabv_data3/large_analysis/ani_output'
path_referencegenome = '/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/assemblies'
path_refseqgenome = '/projects/b1042/HartmannLab/alex/fabv_data3/assemblies'
# 
## Use df_ref to make the reference list, use df_type for query
## large_ani_analysis.csv generated in fabv_LargeAnalysis.R, section 0 
#
df_ref = pd.read_csv(pjoin(path_largeanalysis,'large_ani_analysis.csv'))
df_type = df_ref[df_ref['representative']==1]
type_list = df_type['filename'].tolist()
#
#
def splits_RefFiles(chunk_size,df_ref,path_largeanalysis,path_refseqgenome):
	'''
	Makes a reference file into n number of equal sized chunks
	'''	
	full_ref_list = [i for i in df_ref[df_ref['representative']==0]['filename'].tolist()]
	chunked_ref_array = np.array_split(full_ref_list,chunk_size)
	for i in range(chunk_size):
		chunk_delim = str(i+1)
		print('chunk %s length: '%chunk_delim, len(chunked_ref_array[i]))
		with open(pjoin(path_largeanalysis,'refgenomes%s.txt'%chunk_delim),'w') as outfile:
			for genome in chunked_ref_array[i]:
				outfile.write('%s\n'%pjoin(path_refseqgenome,genome))

splits_RefFiles(chunk_size=4,df_ref=df_ref,path_largeanalysis=path_largeanalysis,path_refseqgenome=path_refseqgenome)



def mark_ANICompletion(chunk_size,df_type,path_largeanalysis,path_anioutput,tofile):
	'''
	Checks off genomes that have the specified number of chunk files output generated. writes to df
	'''
	merged_output_files_list = []
	for i in range(chunk_size):
		chunk_delim = str(i+1)
		output_file_temp = [f.replace('_%s.txt'%chunk_delim,'') for f in os.listdir(path_anioutput) if f.endswith('%s.txt'%chunk_delim)]
		merged_output_files_list.append(output_file_temp)
	output_files_list = [y for x in merged_output_files_list for y in x]
	df_outputfiles = pd.Series(output_files_list).value_counts().rename_axis('filename').to_frame('counts')
	df_outputfiles = df_outputfiles[df_outputfiles['counts']==2]
	df_outputfiles.reset_index(level=0,inplace=True)
	#
	df_type['ani_completed'] = np.where(df_type['filename2'].isin(df_outputfiles['filename'].tolist()),1,0)
	if tofile == True:
		print('writing ani tracker results to file')
		df_type.to_csv(pjoin(path_largeanalysis,'completed_ani_runs.csv'),index=False)


mark_ANICompletion(chunk_size=2,df_type=df_type,path_largeanalysis=path_largeanalysis,path_anioutput=path_anioutput,tofile=True)



def make_SubmissionScripts(chunk_size,time_i, threads_i, path_mash, path_referencegenome, path_largeanalysis, path_anioutput):
	'''
	Gives fastani command for inputted number of chunks for a single query genome
	'''
	df_todo = pd.read_csv(pjoin(path_largeanalysis,'completed_ani_runs.csv'))
	df_todo = df_todo[df_todo['ani_completed']==0]
	files_todo = df_todo['filename'].tolist()
	for sub_id, f in enumerate(files_todo):
		print('Submission file for:  %s' %f)
		with open(pjoin(path_script,'bin_%s_%s.sh'%(str(sub_id),str(0))),'w') as outfile:
			outfile.write('#!/bin/bash\n')
			outfile.write('#SBATCH -A b1042\n')
			outfile.write('#SBATCH -p genomics\n')
			# outfile.write('#SBATCH --constraint=quest9\n')
			outfile.write('#SBATCH -N 1\n')
			outfile.write('#SBATCH -n %s\n'%str(threads_i))
			outfile.write('#SBATCH -t %s:00:00\n'%time_i)
			outfile.write('#SBATCH --mem=0\n')
			outfile.write('#SBATCH --error=bin_%s_%s.err\n'%(str(sub_id),str(0)))
			outfile.write('#SBATCH --output=bin_%s_%s.log\n'%(str(sub_id),str(0)))
			outfile.write('#SBATCH	--job-name="bin_%s_%s"\n'%(str(sub_id),str(0)))
			outfile.write('#SBATCH --mail-type=BEGIN,END,FAIL\n')
			outfile.write('#SBATCH --mail-user=alexandermcfarland2022@u.northwestern.edu\n')
			outfile.write('\ncd /projects/b1042/HartmannLab/alex/fabv_data3/scripts/ani\n\n')
			# add the genome / refslist pairs to be written
			query_genome_name = pjoin(path_referencegenome,f)
			for i in range(chunk_size):
				chunk_delim = str(i+1)
				base_filename = f.replace('.fna','')
				outfile_name = pjoin(path_anioutput,'%s_%s.txt'%(base_filename,chunk_delim))
				ref_genomeinput = pjoin(path_largeanalysis,'refgenomes%s.txt'%(chunk_delim))
				outfile.write('%s -q %s --rl %s -o %s -t %s\n'%(path_mash,query_genome_name,ref_genomeinput,outfile_name,str(threads_i)))

make_SubmissionScripts(chunk_size=4,time_i=2,threads_i=24,
	path_mash=path_mash,path_referencegenome=path_referencegenome,path_largeanalysis=path_largeanalysis,path_anioutput=path_anioutput)

# ## Submit! 
# cd /projects/b1042/HartmannLab/alex/fabv_data3/scripts/ani
# for f in bin*.sh; do chmod u+x ./$f; done
# for f in *.sh; do sbatch $f; done

#--------------------------------------------------------------------------------------------------
#### Section 2 ####
#### Compile dataframes
path_script = '/projects/b1042/HartmannLab/alex/fabv_data3/scripts/ani'
path_largeanalysis = '/projects/b1042/HartmannLab/alex/fabv_data3/large_analysis'
path_anioutput = '/projects/b1042/HartmannLab/alex/fabv_data3/large_analysis/ani_output'
path_genome = '/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/assemblies'


df_ref = pd.read_csv(pjoin(path_largeanalysis,'large_ani_analysis.csv'))
df_type = df_ref[df_ref['representative']==1]


output_files_1 = [f.replace('_1.txt','') for f in os.listdir(path_anioutput) if f.endswith('1.txt')]
output_files_2 = [f.replace('_2.txt','') for f in os.listdir(path_anioutput) if f.endswith('2.txt')]
matching_pairs = [i for i in output_files_1 if i in output_files_2]
df_type['ani_completed'] = np.where(df_type['filename2'].isin(matching_pairs),1,0)
len(matching_pairs)
len(df_type)

os.system('cat /projects/b1042/HartmannLab/alex/fabv_data3/large_analysis/ani_output/*.txt > /projects/b1042/HartmannLab/alex/fabv_data3/large_analysis/raw_compiled_fastani_results.txt')

# Start here if want to just examine df in memory
df_complete = pd.read_csv(pjoin(path_largeanalysis,'raw_compiled_fastani_results.txt'),sep='\t',header=None)
colnames = ['query','subject','ani','matches','total']
df_complete.columns = colnames

# removing path names from query and subject columns
df_complete['query'] = df_complete['query'].str.replace('/projects/b1042/HartmannLab/alex/fabv_data3/typeanalysis/assemblies/','')
df_complete['subject'] = df_complete['subject'].str.replace('/projects/b1042/HartmannLab/alex/fabv_data3/assemblies/','')
df_complete.to_csv(pjoin(path_largeanalysis,'ani_values.csv'),index=False)
# making a filtered output file
df_complete = df_complete[df_complete['ani']>=95]
df_complete.to_csv(pjoin(path_largeanalysis,'compiled_fastani_results.csv'),index=False)




