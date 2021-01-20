#
'''
All processing scripts related to microbiome selection analysis for goose and soil samples
'''


import os
from os.path import join as pjoin 

import pandas as pd
import numpy as np
from scipy import stats

from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature

import multiprocessing as mp
from multiprocessing import Pool


#### General functions ------------------------------------------------------------------------------------------------------------------

def change_ToWorkingDirectory(directory_name):
	'''
	Change directory to OUTPUT_FILES_DIR
	'''
	os.chdir(directory_name)
	print('\nGOOD MESSAGE: changing working directory to \n{}\n'.format(directory_name))


def make_OutputDirectory(new_directory):
	if os.path.exists(new_directory) == False:
		print('GOOD MESSAGE: making {} directory'.format(new_directory))
		os.mkdir(new_directory)


def multiprocesssing_Submission(function_name,submission_list,processor_n,output_type):

	pool = Pool(processes = processor_n)
	if output_type == 'array':
		parallel_output = pool.starmap(function_name,submission_list)
		pool.close()
		return(parallel_output)
	if output_type == 'file':
		pool.starmap(function_name,submission_list)
		pool.close()


def identify_OrfPacBio(l,input_exp_identifier):
	'''
	identify all orfs in metagenomic reads from pacbio sequencing
	'''	
	os.system('conda run -n prokka_env prokka {}.fasta --outdir {} --locustag {} --prefix {} --fast --norrna --notrna --noanno --cpus 1 --force'.format(l,input_exp_identifier,input_exp_identifier,input_exp_identifier,input_exp_identifier))



def count_ORFs(input_f):
	'''
	'''
	df = pd.read_csv(input_f+'.tsv',sep='\t')
	orf_n = len(df)
	return(orf_n)

def calculate_GC(input_f):
	'''
	'''
	gc_count = 0
	actg_count = 0
	for record in SeqIO.parse(input_f+'.fna','fasta'):
		gc_count += str(record.seq).count('G')
		gc_count += str(record.seq).count('C')
		actg_count += len(record.seq)
	gc_content = gc_count/actg_count
	return(gc_content)




# mp.set_start_method("fork")

input_assemblies_path = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/hiseq_genome_assemblies/assemblies'

output_prokka_path = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/isolate_genomes'

change_ToWorkingDirectory(output_prokka_path)

genome_input_list = [[pjoin(input_assemblies_path,i.replace('.fasta','')),i.replace('_wgs.fasta','')] for i in os.listdir(input_assemblies_path) if i.endswith('_wgs.fasta')]

multiprocesssing_Submission(function_name=identify_OrfPacBio,
	submission_list=genome_input_list,
	processor_n=8,
	output_type='file')



genome_features_array = []
for i in genome_input_list:
	f = i[1]
	input_f = pjoin(output_prokka_path,f,f)
	#
	orf_n = count_ORFs(input_f=input_f)
	gc_content = calculate_GC(input_f=input_f)
	#
	genome_features_array.append([f,orf_n,gc_content])

genome_features_array = np.asarray(genome_features_array)
df = pd.DataFrame({'genome_id':genome_features_array[:,0],
	'orf_count':genome_features_array[:,1],
	'gc_content':genome_features_array[:,2]})
df.to_csv('isolate_genome_features.csv',index=None)





