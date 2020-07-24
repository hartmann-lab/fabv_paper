'''
Scripts for the large fabv analysis
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
#### Extract homologs of interest from the fabv_GetHomologs .faa fabv output
#### fabv_homologs_processed.csv generated in fabv_largeanalysis.R section 2.0
#### Align with mafft and then construct tree


def extract_HomologsOfInterest(df_input,seq_input,seq_output):
	'''
	Extract homologs of interest fom a fabv_GetHomologs output

	'''
	with open(seq_output,'w') as outfile:
		homolog_keep_list = df_input['sseqid'].tolist()
		for record in SeqIO.parse(seq_input,'fasta'):
			if record.id in homolog_keep_list:
				outfile.write('>%s\n'%record.id)
				newseq = textwrap.wrap(str(record.seq),60,break_on_hyphens=False)
				[outfile.write('%s\n'%i) for i in newseq]

path_largeanalysis = '/projects/b1042/HartmannLab/alex/fabv_data3/large_analysis'
path_homologsearchoutput = pjoin(path_largeanalysis,'homologsearch_output')
#
df_fabv = pd.read_csv(pjoin(path_largeanalysis,'fabv_homologs_processed.csv'))
seq_fabv_name = pjoin(path_homologsearchoutput,'pao1_fabv_homologs.faa')
seq_output_name = pjoin(path_homologsearchoutput,'fabv_homologs_processed.faa')
#
extract_HomologsOfInterest(
	df_input = df_fabv,
	seq_input = seq_fabv_name,
	seq_output = seq_output_name)
#
## Alignment, tree construction done on cluster
cd /projects/b1042/HartmannLab/alex/fabv_data3/large_analysis/homologsearch_output
mafft --maxiterate 1000 --localpair --thread -1 fabv_homologs_processed.faa > fabv_homologs_processed.aln
module load anaconda3/2018.12
source activate trim 
trimal -in fabv_homologs_processed.aln -out fabv_homologs_processed_trimmed.aln -automated1
/home/agm9813/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s fabv_homologs_processed_trimmed.aln -n fabv_homologs_processed_trimmed
#

