#
'''
Searching for triclosan tolerance determinants
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

def make_BlastDB(input_fasta, output_db, input_dbtype):
	'''
	'''
	os.system('makeblastdb -blastdb_version 5 -in {} -out {} -dbtype {} -title "{}" -parse_seqids'.format(input_fasta,output_db,input_dbtype,output_db))


def blast_TheDB(query,db_name,output_filename, input_dbtype):
	'''
	'''
	if input_dbtype == 'prot':
		os.system('blastp -query {} -db {} -max_target_seqs 10000000 -evalue 1e-6 -outfmt "10 qseqid sseqid mismatch positive gaps ppos pident qcovs evalue bitscore qframe sframe sstart send slen qstart qend qlen" -num_threads 8  -out {}'.format(query,db_name,output_filename))
	else:
		os.system('tblastn -query {} -db {} -max_target_seqs 10000000 -evalue 1e-6 -outfmt "10 qseqid sseqid mismatch positive gaps ppos pident qcovs evalue bitscore qframe sframe sstart send slen qstart qend qlen" -num_threads 8  -out {}'.format(query,db_name,output_filename))


#### Section 0 set up files and blast for ENR resistance determinants####
base_dir = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/isolate_genomes'
homologs_dir = pjoin(base_dir,'homologs')
ortho_dir = pjoin(homologs_dir,'ortho')

make_OutputDirectory(new_directory=homologs_dir)
change_ToWorkingDirectory(directory_name=homologs_dir)
combined_pseud_aa = 'combined_pseud_aa'

# make blast dataset
#os.system('cat {}/*/*.faa > {}.faa'.format(base_dir,combined_pseud_aa)) # this is dangerous to run by accident so it is hashed out when not used
# make blast db
make_BlastDB(input_fasta=combined_pseud_aa+'.faa', output_db=combined_pseud_aa, input_dbtype='prot')
# search blast database with tcs resistance determinant fasta
tcs_det = 'tcs_determinants'
blast_TheDB(query=tcs_det+'.txt',db_name=combined_pseud_aa,output_filename=combined_pseud_aa+'_hits.csv',input_dbtype='prot')


# for some reason qseqid is changed to pdb format? need to change that
blast_columns = ['qseqid','sseqid','mismatch', 'positive','gaps', 'ppos','pident','qcovs','evalue','bitscore','qframe','sframe','sstart', 'send', 'slen', 'qstart', 'qend', 'qlen']
df_b = pd.read_csv(pjoin(homologs_dir,combined_pseud_aa+'_hits.csv'),header=None,names=blast_columns)
df_b['sseqid'] = df_b['sseqid'].str.replace('pdb\|','').astype(str)
df_b['sseqid'] = df_b['sseqid'].str.replace('|','_').astype(str)
# filter hits for >=80 alignment and >=50 identity
df_b = df_b[(df_b['qcovs']>=80)&(df_b['pident']>=40)]
# if there are qseqids that match the same sseqid, keep the best occurrence
df_b = df_b.sort_values(['sseqid','pident','qcovs'],ascending=False).drop_duplicates(subset=['sseqid'])
# add isolate id
df_b['isolate_id'] =[i[:-6] for i in df_b['sseqid'].tolist()]

df_b.to_csv(pjoin(base_dir,'tcs_determinants_blastresults.csv'),index=None)

# assign clade identity
df_tax = pd.read_csv(pjoin(homologs_dir,'isolate_clades.csv'))
df_tax['clade'] = [c.replace('P. ','') for c in df_tax['clade'].tolist()]
df_b = df_b.merge(df_tax,on='isolate_id')


# for each isolate, extract hit sequences to isoalte specific faa
# isolate specific
multi_extract_seqs_input = [[sp,df_b,ortho_dir,combined_pseud_aa,False] for sp in df_b['isolate_id'].tolist()]
# clade wise
#multi_extract_seqs_input = [[sp,df_b,ortho_dir,combined_pseud_aa,True] for sp in df_b['clade'].unique().tolist()]

def extract_Seqs(sp,df_b,ortho_dir,combined_pseud_aa,clade_wise):
	'''
	'''
	if clade_wise == False:
		df_temp = df_b[df_b['isolate_id']==sp]
	if clade_wise == True:
		df_temp = df_b[df_b['clade']==sp]
	#
	with open(pjoin(ortho_dir,sp+'_tcs.faa'),'w') as outfile:
		blast_hit_list = df_temp['sseqid'].tolist()
		for record in SeqIO.parse(combined_pseud_aa+'.faa','fasta'):
			if record.id in blast_hit_list:
				outfile.write('>{}\n'.format(record.id))
				outfile.write('{}\n'.format(str(record.seq)))

# make a new folder for ortho output
make_OutputDirectory(new_directory=ortho_dir)

# mp.set_start_method("fork")
multiprocesssing_Submission(function_name=extract_Seqs,submission_list=multi_extract_seqs_input,processor_n=8,output_type='file')

# use orthoclustering tool to make these badbois into clusters
# !might have to run outside of ipython if not doing clade-wise approach!
os.system('ulimit -n 1396')
os.system('ulimit -Sn')
os.system('echo orthofinder -1 -t 8 -a 8 -X -y -S diamond -I 6.0 -f {}'.format(ortho_dir))  #orthofinder -t 8 -a 8 -X -S diamond -I 1.5 -f ortho

# read in orthogroup assignments
df_o = pd.read_csv(pjoin(homologs_dir,'ortho/OrthoFinder/Results_Jan14/Orthogroups','Orthogroups.tsv'),sep='\t')
# reformat commas to tabs for multiple genes within an orthogroup
df_o = df_o.replace(', ','\t',regex=True)
# mimicking the output of roary
gene_pa_colnames=["Gene","Non unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc"]

# read in counts of orthogroups
df_oc = pd.read_csv(pjoin(homologs_dir,'ortho/OrthoFinder/Results_Jan14/Orthogroups','Orthogroups.GeneCount.tsv'),sep='\t')
df_oc = df_oc[['Orthogroup','Total']]

df_o = df_o.merge(df_oc,on='Orthogroup')

### MAKING gene presence absence tabele input for SCOARY ###

## counting the number of unique isolates that have a representative in an orthogroup
ortho_members_list = []
with open(pjoin(homologs_dir,'ortho/OrthoFinder/Results_Jan14/Orthogroups','Orthogroups.txt')) as infile:
	for line in infile:
		orthogroup_id = line[:line.find(':')]
		ortho_members = line[line.find(':')+2:]
		ortho_members = ortho_members.split(' ')
		ortho_members = [m[:-6].replace('_','') for m in ortho_members]
		unique_isolates = len(set(ortho_members))
		ortho_members_list.append([orthogroup_id,unique_isolates])

ortho_members_list = np.asarray(ortho_members_list)
ortho_members_list = ortho_members_list[:-1] #removing OG0000006 because it is not in the other files

no_isolates_col = ortho_members_list[:,1].astype(float) # number of unique isolates that have a representative in an orthogroup
no_sequences_col = np.asarray(df_o['Total'].tolist()).astype(float) # total number of genes found in that orthogroup
avg_seq_per_isolate = no_sequences_col/no_isolates_col

# setting up the standard part of the gene p/a table
df_pa = pd.DataFrame({'Gene':df_o['Orthogroup'],'Non unique Gene name':'','Annotation':['effluxpump1','effluxpump2','effluxpump3','effluxpump4','effluxpump5','fabi','fabv','effluxpump6'],
	'No. isolates':no_isolates_col,'No. sequences':no_sequences_col,'Avg sequences per isolate':avg_seq_per_isolate,
	"Genome Fragment":'',"Order within Fragment":'',"Accessory Fragment":'',"Accessory Order with Fragment":'',"QC":'',"Min group size nuc":'',"Max group size nuc":'',"Avg group size nuc":''})

# getting just gene p/a qualitative info from df_o
df_paisolate = df_o.drop(['Orthogroup','Total'],axis=1)

#put it all together
df_pa = pd.concat([df_pa,df_paisolate],axis=1)

# write to file
df_pa.to_csv(pjoin(homologs_dir,'makeshift_gene_pa.csv'),index=None)

# reference roary output
# df_e = pd.read_csv('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/pseud_comparisons_project/results/roary_output/final_annotation_i95/gene_presence_absence.csv')

##### MAKING trait table ###
df_mic = pd.read_csv(pjoin(homologs_dir,'isolate_mic.csv'))
df_mic['isolate_id']=[n+'_tcs' for n in df_mic['isolate_id'].tolist()]

mictrait=[]
for x in df_mic['mic'].tolist():
	if x.find('128') >-1:
		mictrait.append(1)
	else:
		mictrait.append(0)
df_mic['highmic'] = mictrait
df_mic = df_mic[['isolate_id','highmic']]
df_mic.columns = ['Name','highmic']

df_mic.to_csv(pjoin(homologs_dir,'makeshif_traittable.csv'),index=None)

#### running scoary
change_ToWorkingDirectory(homologs_dir)
os.system('scoary -t makeshif_traittable.csv -g makeshift_gene_pa.csv -n phylogenetic_tree.newick')


#### Extracting members of each orthogroup and writing to long table format ####

## counting the number of unique isolates that have a representative in an orthogroup
ortho_members_list = []
with open(pjoin(homologs_dir,'ortho/OrthoFinder/Results_Jan14/Orthogroups','Orthogroups.txt')) as infile:
	for line in infile:
		orthogroup_id = line[:line.find(':')]
		ortho_members = line[line.find(':')+2:]
		ortho_members = ortho_members.split(' ')
		ortho_members = [[orthogroup_id,i.replace('\n','')] for i in ortho_members]
		#ortho_members = np.asarray(ortho_members).flatten().flatten()
		ortho_members_list.append(ortho_members)


flat_list = []
for sublist in ortho_members_list:
	for item in sublist:
		flat_list.append(item)

ortho_members_list = np.asarray(flat_list)
#ortho_members_list = ortho_members_list[:-1] #removing OG0000008 because it is not in the other files
df_orthomembers = pd.DataFrame({'orthogroup_id':ortho_members_list[:,0],'detected_ortholog':ortho_members_list[:,1]})
df_orthomembers.to_csv(pjoin(homologs_dir,'orthogrop_members_long.csv'),index=None)




















