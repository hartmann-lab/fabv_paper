
'''
Each number indicates an individual run of the script to generate the output for the next section.
'''


import os
from os.path import join as pjoin
import pandas as pd
import numpy as np
import argparse
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
import textwrap

##--------1. Generate input_df and assign groups for multiple script submission-------------------------------
'''
Usage: Run in python to generate a dataframe of the samples to be assembled
'''

INPUT_PATH = '/projects/b1042/HartmannLab/alex/broad_genomes'
bam_list = [x.replace('.bam','') for x in os.listdir(INPUT_PATH) if x.endswith('.bam')==True if x.find('qsort')==-1]
#process assignment pa
pa = [1,2,3,4,5]
pa = np.repeat(pa,6).tolist()
#process df
df = pd.DataFrame({'sample':bam_list,'group':pa})
df.to_csv('sample_processing.csv',index=False)
#df becomes --input_df
#id from df.group becomes --group_id

#----------2. Create assembly-----------------------------------------------------------------------------------

'''
Usage: python assemble_genomes.py --input_df sample_processing.csv --group_id 1
requires: samtools, fastp(local), bbmap(local), spades(local), python/anaconda3
Can be run as a for loop for each item in input_df.sample or as a batch submission script using group_id to chunk files
'''
parser = argparse.ArgumentParser()
parser.add_argument('--input_df', required=True)
parser.add_argument('--group_id', required=True)
args = parser.parse_args()

INPUT_PATH = '/projects/b1042/HartmannLab/alex/broad_genomes'
os.chdir(INPUT_PATH)

def assemble_FromBam(input_df, group_id):
	'''
	Assemble BAM file with megahit by first converting to fq and trimming with fastp.
	'''
	FASTP_PATH = '/projects/p30629/fastp'
	BBMAP_PATH = '/home/agm9813/bbmap/bbmerge.sh'
	SPADES_PATH = '/home/agm9813/SPAdes-3.13.0-Linux/bin/spades.py'
	#
	input_df = pd.read_csv(input_df)
	input_df = input_df[input_df['group']==group_id]
	for sample in input_df['sample'].tolist():
		print(sample)
		#convert to fq
		os.system('samtools bam2fq -n %s.bam > %s.fq'%(sample,sample))
		#read quality
		os.system('%s --thread 16 -i %s.fq  -o %s_trim.fq --html %s.html --json %s.json '%(FASTP_PATH,sample,sample,sample,sample))
		#merge paired ends
		os.system('%s k=8 in=%s_trim.fq out=%s_merged.fq outu=%s_unmerged.fq ihist=%s_hist.txt'%(BBMAP_PATH,sample,sample,sample,sample))
		# assemble (long, 21 minuts for 88B1)
		os.system('%s -t 24 -s %s_merged.fq -o %s_spades'%(SPADES_PATH,sample,sample))

assemble_FromBam(str(args.input_df),int(args.group_id))


##Example submission script
# #!/bin/bash
# #SBATCH -A b1042
# #SBATCH -p genomics
# #SBATCH -N 1
# #SBATCH -n 24
# #SBATCH -t 5:00:00
# #SBATCH --mem=0
# #SBATCH --error=assembly5.err
# #SBATCH --output=assembly5.log
# #SBATCH --job-name="assembly5"
# #SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=alexandermcfarland2022@u.northwestern.edu
#
# module load python/anaconda3 bedtools/2.25.0 samtools/1.6 megahit/1.0.6.1
# cd /projects/p30629/submission_scripts
# python assemble_genomes.py --input_df sample_processing.csv --group_id 5


#----------3. Extract finished assemblies, read info, and run checkM on assembly-------------------------------------------------------
'''
Makes a new folder called assemblies containing the following in order of process:
1. .html and .json of reads after sequence chechking with fastp
2. histogram of merged reads from bbamp
3. all assemblies
Each spades run generates a folder with name <sample> in it is 
Usage: run in python using 24 cores for checkm
requires: checkm, python/anaconda
'''
INPUT_PATH = '/projects/b1042/HartmannLab/alex/broad_genomes'
os.chdir(INPUT_PATH)
input_df = pd.read_csv(pjoin(INPUT_PATH,'sample_processing.csv'))

#move files to assemblies folder
if os.path.exists(pjoin(INPUT_PATH,'assemblies')) == False: os.mkdir('assemblies')
for sample in input_df['sample'].tolist():
	#copy .html, .json, and _hist.txt files to assemblies
	os.system('cp %s.html %s.json %s_hist.txt ./assemblies'%(sample,sample,sample))
	#copy contigs.fasta and rename to sample
	os.system('cp ./%s_spades/contigs.fasta ./assemblies/%s.fasta'%(sample,sample))

#run checkm. the checkm file is save in the INPUT_PATH folder. assemblies should only contain sample-specific files
os.system('checkm taxonomy_wf -t 24 -x fasta genus Pseudomonas %s/assemblies %s/checkm_output -f %s/checkm_output.txt --tab_table'%(INPUT_PATH,INPUT_PATH,INPUT_PATH))
#convert into isolate_info.csv input file for fabv_InputCreation.py
df_checkm = pd.read_csv(pjoin(INPUT_PATH,'checkm_output.txt'),sep='\t')
#rename filenames
cleaned_filename = ['109A1',
'10A6',
'114A4',
'115A1',
'119A3',
'20A1',
'31A8',
'34A1',
'39A1',
'45C2',
'4A7',
'56A10',
'57B2',
'62A4',
'66C3',
'69C1',
'6C6',
'82B1',
'88B1',
'89C1',
'8A1',
'95A6',
'96A1',
'97C1',
'99A1',
'HS_1',
'HS_2',
'HS_3',
'HS_4',
'HS_5']

#getting real filename
filenames = []
for i in df_checkm['Bin Id'].tolist():
	i = i +'.fasta'
	filenames.append(i)
##create the isolate_input.csv file to use in fabv_InputCreation.py
df_isolate_input = pd.DataFrame({'filename':filenames,'id':cleaned_filename,'species':'unknown','description':'empty'})
df_isolate_input.to_csv(pjoin(INPUT_PATH,'isolate_info.csv'),index=False)


#----------4. Get read count and coverage for each genome -------------------------------------------------------
## Generates a data table containing the isolate sequence id, insert ength * read count, and name of corresponding assembly file
## the assembly file has been modified to exclude contigs with less than 200 bp per NCBI standards
INPUT_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/hiseq_genome_assemblies/assemblies'
os.chdir(INPUT_PATH)
# will output to fabv_miscellaneous 
OUT_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/fabv_data/miscellaneous'

# list of isolate ids that will be used instead of the long ass filename
cleaned_filename = ['109A1',
'10A6',
'114A4',
'115A1',
'119A3',
'20A1',
'31A8',
'34A1',
'39A1',
'45C2',
'4A7',
'56A10',
'57B2',
'62A4',
'66C3',
'69C1',
'6C6',
'82B1',
'88B1',
'89C1',
'8A1',
'95A6',
'96A1',
'97C1',
'99A1',
'HS_1',
'HS_2',
'HS_3',
'HS_4',
'HS_5']

## get unique filenames minus extensions
filenames = [file.replace('.html','') for file in os.listdir(INPUT_PATH) if file.endswith('.html')]
read_info = []

## list containing sequences that need to be removed from the assembly because NCBI has found them to be contaminants
contaminants = {'66C3':'NODE_45_length_558_cov_0.529002','97C1':
['NODE_28_length_764_cov_0.576138','NODE_35_length_548_cov_0.498812','NODE_40_length_444_cov_0.583596',
'NODE_30_length_699_cov_0.541958','NODE_31_length_673_cov_0.622711','NODE_36_length_540_cov_0.733656','NODE_38_length_458_cov_0.873112'],
 '109A1':['NODE_82_length_830_cov_0.019915']}

## open histogram txt file and calculate total read counts, then open the assembly and keep only sequences that are >= 200bp
for input_filename in filenames[:]:
	## read info gathering and adding to metadata table
	df_hist = pd.read_csv(input_filename+'_hist.txt',sep='\t')
	df_hist = df_hist.iloc[5:]
	df_hist.columns=['insert_length','read_count']
	df_hist['ln'] = df_hist['insert_length'].astype(int)*df_hist['read_count'].astype(int)
	#break
	read_count = df_hist['ln'].sum()
	# find the id within the filename
	file_id = [i for i in cleaned_filename if input_filename.find(i)>-1][0]
	# save id, read length, and new filename
	read_info.append([file_id,read_count,'%s_wgs.fasta'%file_id])
	## open assembly fasta and keep sequences that are >=200 bp
	with open('%s_wgs.fasta'%file_id,'w') as outfile:
		for record in SeqIO.parse(input_filename+'.fasta','fasta'):
			if file_id in list(contaminants.keys()):
				contaminating_seq = contaminants[file_id]
				if record.id in contaminating_seq:
					continue
			if len(record.seq)>200:
				outfile.write('>%s\n'%record.id)
				wrapped_seq = textwrap.wrap(str(record.seq),60,break_on_hyphens=False)
				[outfile.write('%s\n'%i) for i in wrapped_seq]
	print(input_filename)
	print(file_id)

# converting to array and then pd df
y = np.asarray(read_info).flatten().reshape(30,3)
df_readinfo = pd.DataFrame({'id':y[:,0],'read_count':y[:,1],'filename':y[:,2]})
# saving
df_readinfo.to_csv(pjoin(OUT_PATH,'assembly_readcounts.csv'),index=False)











