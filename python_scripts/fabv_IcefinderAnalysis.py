import os
from os.path import join as pjoin
import pandas as pd
import re
import numpy as np
import requests, sys
from requests_toolbelt import MultipartEncoder
from Bio import SeqIO
import textwrap


#### Section 0 #####
#### Create dataframe with the icefinder job submission id
#### Merge with exisiting metadata
path_largeanalysis = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/large_analysis'
filename_islandviewerprediction = pjoin(path_largeanalysis,'sec2_4_toivsubmit.csv')
path_iceanalysis = pjoin(path_largeanalysis,'ice_analysis')
filename_icesubmission = pjoin(path_iceanalysis,'ice_analysis_submission.txt')
#
df_meta = pd.read_csv(filename_islandviewerprediction)
df_meta = df_meta[['filename','mod_species','cleaned_filename','total_contigs','description','clade']]
#
df_subid = pd.read_csv(filename_icesubmission,sep='\t',header=None)
df_subid.columns = ['refseq_id','mod_species','filename','ice_id','gross_hit']
df_subid = df_subid[['refseq_id','filename','ice_id','gross_hit']]
#
df_subid = df_subid.merge(df_meta,on='filename')
df_subid = df_subid.drop_duplicates()
df_subid.to_csv(pjoin(path_iceanalysis,'processed_ice_analysis_submission.csv'),index=False)

#### Section 0.1 ####
#### make new edited fasta files from the raw region fasta files downloaded from icefinder 
path_rawresults = pjoin(path_iceanalysis,'raw_results')
#
df_subid = pd.read_csv(pjoin(path_iceanalysis,'processed_ice_analysis_submission.csv'))
#
for f in os.listdir(path_rawresults):
	if f.endswith('.fas'):
		# input and output filenames
		infile_f = pjoin(path_rawresults,f)
		outfile_f = pjoin(path_iceanalysis,'processed_'+f)
		# defline renaming variables
		sub_id = f[:f.find('_')] # the genome submission id
		region_id = f[f.find('_R')+1:f.find('.f')] # the icefinder region id
		# rewriting deflines and writing to output file
		with open(outfile_f,'w') as outfile:
			for x,record in enumerate(SeqIO.parse(infile_f,'fasta')):
				orf_position = str(x+1)
				mod_record = sub_id+'_'+region_id+'_'+orf_position
				mod_seq = textwrap.wrap(str(record.seq),60,break_on_hyphens=False)
				outfile.write('>%s\n'% mod_record)
				[outfile.write('%s\n'%i) for i in mod_seq]
# concatenate the processed filenames into a single file
filename_concatfasta = 'processed_concatenated_icefinder.faa'
#
os.chdir('%s'%path_iceanalysis)
os.system('cat *processed*.fas > %s'%filename_concatfasta)


#### Section 0.2 ####
#### blast search for fabv 
filename_concatfasta = 'processed_concatenated_icefinder.faa'
filename_concatdb = 'processed_concatenated_icefinder'
filename_fabvquery = 'pao1_fabv.txt'
#
os.chdir('%s'%path_iceanalysis)
os.system('makeblastdb -in %s -out %s -dbtype prot -title "%s_db" -parse_seqids'%(filename_concatfasta,filename_concatdb,filename_concatdb))
os.system('blastp -query %s -db %s -max_target_seqs 100 -evalue 1e-6 -outfmt "10 sseqid stitle mismatch positive gaps ppos pident qcovs evalue bitscore" -num_threads 1  -out %s.csv' % (filename_fabvquery, filename_concatdb, filename_concatdb))

#### Section 0.3 ####
#### associate hits with existing metadata
filename_outputconcatdb = pjoin(path_iceanalysis,'processed_concatenated_icefinder.csv')
#
df_blastout = pd.read_csv(filename_outputconcatdb,header=None)
df_blastout.columns = ['sseqid','stitle','mismatch','positive','gaps','ppos','pident','qcovs','evalue','bitscore']
ice_id_list = [x[:x.find('_')] for x in df_blastout['sseqid'].tolist()]
df_blastout['ice_id'] = ice_id_list
#
df_subid = pd.read_csv(pjoin(path_iceanalysis,'processed_ice_analysis_submission.csv'))
df_subid = df_subid.merge(df_blastout,on='ice_id',how='outer')

































