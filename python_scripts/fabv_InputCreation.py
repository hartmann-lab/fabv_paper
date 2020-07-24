
'''
Assemblies from either REFSEQ or USER are supplied. USER assemblies must have a file with a 
df containing the assembly filename, id, and genus_species assignment
Outputs are directed to the user-defined DATA_PATH. 
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


def cat_Files(file_list, extension, DATA_PATH_, out_filename_path):
	'''
	General purpose function to cat thousands of files and delete the intermediate chunking files
	file_list must have the path joined to the name and also the file extension
	extension should be written without the <.> examples: <txt> <csv>
	DATA_PATH_ is the output directory for the chunks
	out_filename_path must have the path joined to the name but not the extension
	'''
	major_length = len(file_list) - len(file_list)%100  #get length to chunk that is divisble by 100
	i = 0
	while i < len(file_list):
		i2 = i
		i += 100
		if i == major_length:
			print('final', '   ', len(file_list[i2:])) #for tracking
			chunk = ' '.join(file_list[i2:])
			chunk_path = pjoin(DATA_PATH_, 'final_chunk')
			os.system('cat %s > %s.%s' % (chunk, chunk_path, extension))
			break
		chunk = ' '.join(file_list[i2:i])
		print(i, '   ', len(file_list[i2:i])) #for tracking
		chunk_path = pjoin(DATA_PATH_, '%s_chunk' % str(i))
		os.system('cat %s > %s.%s' % (chunk, chunk_path, extension))
	wild_chunk_path = pjoin(DATA_PATH_,'*_chunk.%s' % extension)
	os.system('cat %s > %s.%s' % (wild_chunk_path, out_filename_path, extension))
	os.system('rm %s' % wild_chunk_path)


def get_GenInputsPATHS(INPUT_DATA_PATH_, DATA_PATH_, INPUT_FILE_DF_):
	'''
	Makes all inputs global variables. Checks that the output DATA_PATH exists and if not creates it
	'''
# INPUT_DATA_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/pdua_project/tests_current/region_extraction/test_files'
# DATA_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/pdua_project/tests_current/region_extraction/test1'
# INPUT_FILE_DF = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/pdua_project/tests_current/region_extraction/test_files/isolate_info.csv'
	global INPUT_DATA_PATH, DATA_PATH, INPUT_FILE_DF
	INPUT_DATA_PATH = INPUT_DATA_PATH_
	DATA_PATH = DATA_PATH_
	INPUT_FILE_DF = INPUT_FILE_DF_
	if os.path.exists(DATA_PATH) == False: os.mkdir(DATA_PATH)



def parallel_BuildInputDataFrame(INPUT_DATA_PATH,raw_fasta):
	'''
	For one raw fasta input, extract filename, id, genus_species, and description from each raw refseq assembly
	'''
	#using SeqIO to parse files. Only need information from one delimiter that is not of a plasmid. Which is why 'if' and 'break' are used here
	#raw_fasta='GCF_000003215.1_ASM321v1_genomic.fna'
	for record in SeqIO.parse(pjoin(INPUT_DATA_PATH,raw_fasta),'fasta'):
		if record.description.find('lasmid') > -1: continue
		raw_description = record.description.replace('%s ' % record.id, '')
		#description_extract = re.findall(r' (.*?) ', raw_description)
		genus,species = raw_description.split()[0].lower(),raw_description.split()[1].lower()
		genus_species = genus+'_'+species
		genome_features = [raw_fasta, record.id, genus_species, raw_description]
		break
	return(genome_features)



def build_InputDataframe(INPUT_DATA_PATH,DATA_PATH,INPUT_FILE_DF,processes):
	'''
	Creates the main processed_genomes_metadata.csv file from a directory containing both user and refseq genomes
	Step1:Refseq genomes have a dataframe buit by accessing the defline of each raw fasta assembly
	Step2:User genomes have a dataframe supplied by the user
	Step3:Merges both dataframes after processing
	Step4:Standardizs species by assigning a code
	Step5:Makes a cleaned_filename for each assembly that is then used as the principal genome identifier
	'''
	###Step1:Refseq df creation
	##Using parallel function to extract filename, id, species, and description from each raw refseq assembly
	#getting list of refseq assemblies
	raw_fasta_inputs = [x for x in os.listdir(INPUT_DATA_PATH) if x.startswith('GCF')]
	#making parallel inputs
	parallel_input = [[INPUT_DATA_PATH,raw_fasta] for raw_fasta in raw_fasta_inputs]
	pool = Pool(processes = processes)
	#genome features stores each parallel output as a list
	genome_features = pool.starmap(parallel_BuildInputDataFrame, parallel_input[:])
	pool.close()
	#genome_features is turned into a pandas df: df_refseq
	genome_features = np.asarray(genome_features)
	df_refseq = pd.DataFrame({'filename' : genome_features[:,0],'id' : genome_features[:,1],'species' : genome_features[:,2],'description' : genome_features[:,3]})
	df_refseq['source'] = 'REFSEQ'
	# checking whether input data was supplied
	if os.path.exists(INPUT_FILE_DF) == True:
		###Step2:User df creation
		##read in user-inputted dataframe
		df_user = pd.read_csv(INPUT_FILE_DF)
		df_user['source'] = 'USER'
		####Step3:Merging both df and then making species code and cleaned_filename columns
		df_merge = df_refseq.append(df_user)
	else:
		df_merge = df_refseq
	##Step4:getting unique species and assigning it a number using dictionary
	unique_species_list = df_merge['species'].unique().tolist()
	species_key = dict(zip(unique_species_list, list(range(len(unique_species_list)))))
	#renaming species in species_code with number from species_key dictionary
	df_merge['species_code'] = df_merge['species']
	df_merge = df_merge.replace({'species_code':species_key})
	##Step5:making cleaned_filename which will be the designated genome name from now until forever
	#storing cleaned_filename in list that will then be added to the df_merge df
	cleaned_ids = []
	#cleaned filename is <species_ass>_<genome_id[3:]> for refseq genomes
	#cleanedfilename is <species_ass>_<genome_id> for user genomes
	for genome_id in df_merge['id'].tolist():
		df_genome = df_merge[df_merge['id']==genome_id]
		if df_genome['source'].item()=='REFSEQ':
			cleaned_id = genome_id[3:] #removing the NZ_ prefix in refseq assembly names
			cleaned_id = cleaned_id[:cleaned_id.find('.')] #removing the .1 suffix in refseq assembly names
		else:
			cleaned_id = genome_id
		#combining species_ass with genome_id to make cleaned_filename
		species_ass = df_genome['species_code'].item()
		cleaned_ids.append('%s_%s' % (species_ass, cleaned_id))
	#adding cleaned_ids to df_merge
	df_merge['cleaned_filename'] = cleaned_ids
	#writing to file
	df_merge.to_csv(pjoin(DATA_PATH,'processed_genomes_metadata.csv'),index=False)
	print('extracted %s assemblies'%(len(df_merge.cleaned_filename.tolist())))


def parallel_ProcessRawGenomes(df_asy,cf,INPUT_DATA_PATH,prokka_output,quast_output,blast_output):
	'''
	For the each raw assembly input:
	1.Create a clean fasta file 2.faa,ffn,gff files from prokka 3. Get QUAST metrics 4. make a BlastDB
	'''
	#each cf will have its own df_cf that contains the variables needed for processing
	df_cf = df_asy[df_asy['cleaned_filename']==cf]
	##Step1: making a new cleaned_assembly_fasta with cleaned def lines. this assembly is located in prokka_output and will be used for prokka processing
	cleaned_assembly_fasta = pjoin(prokka_output,cf+'.fasta')
	#making new assembly
	with open(cleaned_assembly_fasta,'w') as outfile:
		raw_file = df_cf['filename'].item()
		for contig_count, record in enumerate(SeqIO.parse(pjoin(INPUT_DATA_PATH,raw_file),'fasta')):
			contig_count+=1 
			outfile.write(">%s_%s\n%s\n" % (cf, str(contig_count),str(record.seq)))
	##Step2: faa, ffn, gff generation with prokka. all outputs go to prokka_output. all other outputs not listed are removed using rm_prokka_out
	os.system('prokka %s --notrna --norrna --noanno --cpus 1 --force --quiet --locustag %s --prefix %s --outdir %s' % (cleaned_assembly_fasta, cf, cf, prokka_output))
	rm_prokka_out = pjoin(prokka_output,cf)
	os.system('rm %s.err %s.fna %s.fsa %s.gbk %s.log %s.sqn %s.tbl %s.tsv %s.txt' %(rm_prokka_out,rm_prokka_out,rm_prokka_out,rm_prokka_out,rm_prokka_out,rm_prokka_out,rm_prokka_out,rm_prokka_out,rm_prokka_out)) 
	##Step3: getting quast metrics. only keeping transposed_report.tsv output. everything else is removed using rm_quast_out
	#quast outputs to a directory. so a new directory with cf is made for each cf quast
	cf_quast_output = pjoin(quast_output,cf)
	os.system('quast.py %s --fast --silent --threads 1 -o %s'%(cleaned_assembly_fasta,cf_quast_output))
	#the transposed_report.tsv file needs to read in to a pd.df because it is all messed up
	#after it is read in and processed it will be saved to the base quast_output directory. no headers are saved so that it can be processed efficiently by cat_Files
	#from the cf_quast_output folder process only transposed_report.tsv
	col_names = ['cleaned_filename','contigs_greater_0_bp', 'contigs_greater_1000_bp','contigs_greater_5000_bp', 'contigs_greater_10000_bp',  'contigs_greater_25000_bp', 'contigs_greater_50000_bp', 'total_length_greater_0_bp', 'total_length_greater_1000_bp', 'total_length_greater_5000_bp', 'total_length_greater_10000_bp', 'total_length_greater_25000_bp', 'total_length_greater_50000_bp',  'total_contigs', 'largest_contig', 'total_length', 'N50', 'N75', 'L50', 'L75',  'number_Ns_per_100_kbp']
	df_cf_transposed = pd.read_csv(pjoin(cf_quast_output,'transposed_report.tsv'),sep='\t',skiprows=1,names=col_names,engine='python')
	df_cf_transposed.to_csv(pjoin(quast_output,cf+'.csv'),index=False,header=False)
	#remove the cf_quast_output directory
	os.system('rm -R %s'%cf_quast_output)
	##Step4: make blastDB
	#faa path for cf
	cf_blast_input = pjoin(prokka_output,cf)
	#blast database output
	cf_blast_output = pjoin(blast_output,cf)
	os.system('makeblastdb -in %s.faa -out %s -dbtype prot -title "%s_db" -parse_seqids' % (cf_blast_input, cf_blast_output, cf))
	os.system('echo %s finished'%cf)



def process_RawGenomes(DATA_PATH,INPUT_DATA_PATH,processes):
	'''
	For each cleaned_filename that is in processed_genomes_metadata.csv do the following:
	1.Create a clean fasta file 2.faa,ffn,gff files from prokka 3. Get QUAST metrics
	'''
	#making an outdir for prokka files, quast files, and blast database files
	prokka_output = pjoin(DATA_PATH,'prokka_output')
	quast_output = pjoin(DATA_PATH,'quast_output')
	blast_output = pjoin(DATA_PATH,'blast_output')
	if os.path.exists(prokka_output) == False: os.mkdir(prokka_output)
	if os.path.exists(quast_output) == False: os.mkdir(quast_output)
	if os.path.exists(blast_output) == False: os.mkdir(blast_output)
	#reading in genomes metadata
	df_asy = pd.read_csv(pjoin(DATA_PATH,'processed_genomes_metadata.csv'))
	#making input list
	cleaned_filename_list = [x for x in df_asy['cleaned_filename'].tolist()]
	#creating inputs for parallel_ProcessRawGenomes
	parallel_input = [[df_asy,cf,INPUT_DATA_PATH,prokka_output,quast_output,blast_output] for cf in cleaned_filename_list]
	pool = Pool(processes = processes)
	pool.starmap(parallel_ProcessRawGenomes, parallel_input[:])
	pool.close()


def add_QUASTMetrics(DATA_PATH):
	'''
	Compiles all QUAST csv outputs into a master file which is then merged with processed_genomes_metadata.csv
	Step1: merge quast files into df_concat_quast
	Step2: add headers to df_concat_quast
	Step3: merge df_concat_quast with processed_genomes_metadata.csv and save as processed_genomes_metadata
	Step4: remove quast_output directory, the concat_quast.csv file (quast data is present in processed_genomes_metadata.csv now)
	'''
	#getting directory for quast_output
	quast_output = pjoin(DATA_PATH,'quast_output')
	##Step1:concatenating all csv outputs. skipping the first line as it is the header
	#input quast files
	quast_files = [pjoin(quast_output,x) for x in os.listdir(quast_output) if x.endswith('.csv')]
	#making the dataframe output
	df_concat_quast = pjoin(DATA_PATH,'concat_quast')
	cat_Files(quast_files,'csv',DATA_PATH,df_concat_quast)
	##Step2: add headers to df_concat_quast and save
	col_names = ['cleaned_filename','contigs_greater_0_bp', 'contigs_greater_1000_bp','contigs_greater_5000_bp', 'contigs_greater_10000_bp',  'contigs_greater_25000_bp', 'contigs_greater_50000_bp', 'total_length_greater_0_bp', 'total_length_greater_1000_bp', 'total_length_greater_5000_bp', 'total_length_greater_10000_bp', 'total_length_greater_25000_bp', 'total_length_greater_50000_bp',  'total_contigs', 'largest_contig', 'total_length', 'N50', 'N75', 'L50', 'L75',  'number_Ns_per_100_kbp']
	df_concat_quast = pd.read_csv(pjoin(DATA_PATH,'concat_quast.csv'),names = col_names)
	#df_concat_quast.to_csv(pjoin(DATA_PATH,'concat_quast.csv'))
	##Step3: merge with processed_genomes_metadata.csv
	df_metadata = pd.read_csv(pjoin(DATA_PATH, 'processed_genomes_metadata.csv'))
	df_merge = df_metadata.merge(df_concat_quast, on = 'cleaned_filename')
	df_merge.to_csv(pjoin(DATA_PATH, 'processed_genomes_metadata.csv'),index=False)
	##Step4: removing quast_output_directory and concat_quast.csv
	os.system('rm -R %s'%quast_output)
	os.system('rm %s'%pjoin(DATA_PATH,'concat_quast.csv'))


def run_InputCreation(INPUT_DATA_PATH_,DATA_PATH_,INPUT_FILE_DF_,processes_):
	s = timer()
	get_GenInputsPATHS(INPUT_DATA_PATH_,DATA_PATH_,INPUT_FILE_DF_)
	print('build_InputDataframe')
	build_InputDataframe(INPUT_DATA_PATH,DATA_PATH,INPUT_FILE_DF,processes_)
	print('process_RawGenomes')
	process_RawGenomes(DATA_PATH,INPUT_DATA_PATH,processes_)
	print('add_QUASTMetrics')
	add_QUASTMetrics(DATA_PATH)
	e = timer()
	print('total time:')
	print(e-s)

















