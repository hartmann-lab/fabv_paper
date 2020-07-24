'''
Submit genbanks to islandviewer
'''

import os
from os.path import join as pjoin
import pandas as pd
import re
import numpy as np
import requests, sys
from requests_toolbelt import MultipartEncoder
from Bio import SeqIO


#### Section 0.0 ####
#### Make txt file to use for batch entrez to download _refseq_ .gbk/.gbff from the assembly NCBI database
#### Read in fabv_searchlist.csv generated in fabv_IslandViewerAnalysis.R, section 0
#### entrez record input must look like this: GCF_900104365

#input path/file
path_largeanalysis = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/large_analysis'
filename_islandviewerprediction = pjoin(path_largeanalysis,'sec2_4_toivsubmit.csv')
#output path/file
path_islandvieweranalysis = pjoin(path_largeanalysis,'islandviewer_analysis')
filename_entrezinput = pjoin(path_islandvieweranalysis,'entrez_submission.txt')


def get_GenomeSearchList(df_filepath):
	'''
	Takes 'filename' column for entrez submission, modifies so that it looks like this GCF_900104365
	'''
	df = pd.read_csv(df_filepath)
	df['refseq_id'] = [r[:r.find('.')] for r in df['filename'].tolist()]
	return(df)


def write_BatchEntrezSubmit(df_input,output_filename):
	'''
	Makes refseq list into a .txt file that is then submitted to batch entrez for downloading
	'''
	refseq_list = df_input['refseq_id'].unique().tolist()
	with open(output_filename,'w') as outfile: 
		[outfile.write('%s\n'%i) for i in refseq_list]


df_submit = get_GenomeSearchList(df_filepath=filename_islandviewerprediction)
write_BatchEntrezSubmit(df_input=df_submit,output_filename=filename_entrezinput)


#### Section 0.1 #### 
### gunzip and move islandviewer gbks to desired folder
### Write file containing islandviewer submission information

def untar_GbkTar(gbktar,path_downloads):
	'''
	untar the gbk
	'''
	os.chdir(path_downloads)
	os.system('tar xopf %s'%gbktar)


def process_GbkFolder(path_gbkfolder,path_islandviewer):
	'''
	gunzip, rename .gbff to .gbk, and then move contents from gbkfolder to the islandviewer folder. 
	'''
	os.chdir(gbk_folder)
	os.system('gunzip *.gz')
	for i in os.listdir(gbk_folder):
		if i.endswith('.gbff'):
			i = i.replace('.gbff','')
			os.system('mv %s.gbff %s.gbk'%(i,i))
			os.system('cp %s.gbk %s/%s.gbk'%(i,path_islandviewer,i))


def make_SubmissionDataframe(df_submit, path_islandviewer, output_filename):
	'''
	writes a csv containing the refseq id, cleaned_filename, islandviewer id, species Id, and contig count
	'''
	# gather matching ids
	submission_filenames = []
	for refseq_id in df_submit['refseq_id'].tolist():
		for file in os.listdir(path_islandviewer):
			if file.startswith(refseq_id):
				submission_filenames.append([refseq_id,file])
	# merge ids with existing dataframe
	df_temp = pd.DataFrame(submission_filenames,columns=['refseq_id','gbk_filename'])
	df_submit = df_submit.merge(df_temp,on='refseq_id')
	df_submit = df_submit[['cleaned_filename','reason','id','description','refseq_id','gbk_filename','total_contigs','mod_species','clade']].drop_duplicates()
	df_submit.to_csv(output_filename,index=False)

## extract gbks for islandviewer
path_downloads = '/Users/owlex/Downloads' # manual
entrez_gbks_tar = pjoin(path_downloads,'genome_assemblies_genome_gb_.tar') # manual
path_islandvieweranalysis = path_islandvieweranalysis
untar_GbkTar(gbktar=entrez_gbks_tar,path_downloads=path_downloads)
## put gbks in correct format and in proper directory
gbk_folder = pjoin(path_downloads,'ncbi-genomes-2020-05-18') # manual	
process_GbkFolder(path_gbkfolder=gbk_folder,path_islandviewer=path_islandvieweranalysis)
## make a db that will be used by islandviewer submission to get filenames 
df_submission = get_GenomeSearchList(df_filepath=filename_islandviewerprediction)
islandviewersubmit_filename = pjoin(path_islandvieweranalysis,'largeanalysis_fabvgenomes.csv')
make_SubmissionDataframe(df_submit=df_submission, path_islandviewer=path_islandvieweranalysis,output_filename=islandviewersubmit_filename)


#### Section 1.0 ####
#### Submit genbank files to islandviewer api
#
## submission info
path_islandviewerprediction = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/large_analysis/islandviewer_analysis'
submission_file = pjoin(path_islandviewerprediction,'largeanalysis_fabvgenomes.csv')
df_submission = pd.read_csv(submission_file)
#

#
## server info
#
xauth1='7ca32d2a-57a1-ab24-f92e-2c5beb578b65' # manual
server = "http://www.pathogenomics.sfu.ca/islandviewer"
ext = "/rest/submit/"
address='alexandermcfarland2022@u.northwestern.edu'
reference_gbk = 'NC_002516.2'
#
## Submit 
#
for r_id in df_submission['refseq_id']:
	df_temp = df_submission[df_submission['refseq_id']==r_id]
	mygenome = pjoin(path_islandviewerprediction,df_temp['gbk_filename'].item())
	name = r_id
	print(name)
	if int(df_temp['total_contigs'].item()) > 1:
		os.system('curl -X POST -H %s -Fgenome_file=@%s -Fgenome_name="%s" -Fref_accnum="%s" -Femail_addr="%s" -Fformat_type="GENBANK" http://www.pathogenomics.sfu.ca/islandviewer/rest/submit/'%(xauth1,mygenome,name,reference_gbk,address))
	else:
		os.system('curl -X POST -H %s -Fgenome_file=@%s -Fgenome_name="%s" -Femail_addr="%s" -Fformat_type="GENBANK" http://www.pathogenomics.sfu.ca/islandviewer/rest/submit/'%(xauth1,mygenome,name,address))

#GCF_900113905 error



#### Section 2.0 ####
### Parse results. Integrates 
path_islandvieweranalysis = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/large_analysis/islandviewer_analysis'

stdout_file = pjoin(path_islandvieweranalysis,'islandviewer_submission_stdout.txt')

submission_dataframe = pjoin(path_islandvieweranalysis,'largeanalysis_fabvgenomes.csv')



def parse_IVStdoutText(stdout_file):
	'''
	converts the stdout text from islandviewer submission into a df
	'''
	stdout_raw = []
	with open(stdout_file,'r') as infile:
		for l in infile:
			if l.find('GCF_')>-1: stdout_raw.append(l)
			if l.find('token')>-1: stdout_raw.append(l)
	stdout_list = ','.join(stdout_raw)
	stdout_list = stdout_list.replace('\n','')
	stdout_list = stdout_list.replace('token','')
	stdout_list = stdout_list.replace(' ','')
	stdout_list = stdout_list.replace(':','')
	stdout_list = stdout_list.replace('"','')
	stdout_list = stdout_list.replace('}','')
	stdout_list = stdout_list.replace('{','')
	stdout_list = stdout_list.split(',')
	#
	df_stdout = np.asarray(stdout_list)
	df_stdout = df_stdout.reshape(int(len(stdout_list)/2),2)
	df_stdout = pd.DataFrame({'refseq_id':df_stdout[:,0], 'url_id':df_stdout[:,1]})
	return(df_stdout)


def merge_StdoutSubmissionDf(df_sub, df_stdout, output_filename):
	'''
	merge the submission dataframe with df_stdout and write to file
	'''
	df_sub = pd.read_csv(submission_dataframe)
	df_sub = df_sub.merge(df_stdout, on = 'refseq_id')
	df_sub.to_csv(output_filename,index=False)


df_stdout = parse_IVStdoutText(stdout_file=stdout_file)
output_filename = pjoin(path_islandvieweranalysis,'largeanalysis_fabvgenomes_islandviewerid.csv')
merge_StdoutSubmissionDf(df_sub=submission_dataframe, df_stdout=df_stdout, output_filename=output_filename)


#### Section 2.1 #### 
#### Query islandviewer website 
import webbrowser
path_islandvieweranalysis = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/large_analysis/islandviewer_analysis'

iv_id_filename = pjoin(path_islandvieweranalysis,'largeanalysis_fabvgenomes_islandviewerid.csv')


def query_IVBrowser(iv_id_filename, refseq_query):
	'''
	enter a refseq_id to open up the island viewer website and checkout the results
	'''
	df = pd.read_csv(iv_id_filename)
	df = df[df['refseq_id']==refseq_query]
	ivweb_url = pjoin('http://www.pathogenomics.sfu.ca/islandviewer/results',df['url_id'].item())
	webbrowser.open(ivweb_url,new=2)
	#
	print(refseq_query, '  ' ,df['mod_species'].item(),'  ', df['reason'].item(), '  ', df['total_contigs'].item())



refseq_query = 'GCF_001269725'
query_IVBrowser(iv_id_filename=iv_id_filename,refseq_query=refseq_query)
refseq_query = 'GCF_001238395'
query_IVBrowser(iv_id_filename=iv_id_filename,refseq_query=refseq_query)
refseq_query = 'GCF_003732295'
query_IVBrowser(iv_id_filename=iv_id_filename,refseq_query=refseq_query)


#### Section 3.0 ####
#### Aggregate tsv files
path_islandviewerprediction = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/fabv_data3/fabv_transmission/islandviewer_prediction'
path_refseqtsv = pjoin(path_islandviewerprediction,'refseq_tsv')
# adding filename to each tsv group
append_df = []
for file in os.listdir(path_refseqtsv):
	df_tmp = pd.read_csv(pjoin(path_refseqtsv, file),sep='\t')
	df_tmp['tsvfile_id'] = file.replace('.tsv','')
	append_df.append(df_tmp)
df_predictions = pd.concat(append_df,sort=True)
#
























