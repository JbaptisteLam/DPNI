"""
Mother homozigous for reference allele and foetus got alternate allele from the father 
VAFp = 1/2 FF

Mother homozygous for a alternate allele and fetus got reference allele
VAFm = 1 - (0,5 * FF)

FF = VAFpaverage+(1-VAFm average)
https://doi.org/10.3390/biotech10030017

"""


from functools import lru_cache
from os.path import join as osj
import argparse
import os
import pandas as pd
import pysam
import re
import subprocess
import sys

def print_progress_bar(i, max):
	'''
	i is your current iteration max should be equal to your maximum number of iterations - 1 (otherwise the bar stop at 99%)
	50 and 100 are arbutrary. 50 defines the size of the bar in the terminal. 100 makes this a percentage.
	The block using this function shoudl be followed by a print("\n") to skip to the next line.
	'''
	j = int(i*50/max)
	sys.stdout.write("\r")
	sys.stdout.write("[%-50s]%d%%"%('='*j, 100/max*i))
	sys.stdout.flush()

def systemcall(command):
	"""
	Passing command to the shell, return list containing stdout lines
	"""
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def average(iterable):
	return round(sum(iterable) / len(iterable), 4)

def preprocess(mum, dad, foetal, filetype):
	'''
	input: tsv of variants from trio family
	output: variant inside VAFp and VAFm statements in tsv format
	'''
	print(locals())
	for j, values in locals().items():
		filename = values+'.pickle'
		if not os.path.exists(filename) and j != 'filetype':
			if filetype == 'vcf':
				vcfTodataframe(filename)
			else:
				print(values)
				tsvTodataframe(values)
	
	#FOR TSV only
	mother = pd.read_pickle(mum+'.pickle')
	father = pd.read_pickle(dad+'.pickle')
	foetus = pd.read_pickle(foetal+'.pickle')
	
	filter_foetus, filter_foetus_homo = processtsv(mother, father, foetus)

	return filter_foetus, filter_foetus_homo

def processvcf(mother, father, foetus):
	dataframe_list = [ mother, father, foetus]
	for datafs in dataframe_list:
		datafs = datafs.loc[(datafs['REF'] == '.') | (datafs['ALT'] == '.') | (datafs['REF'].str.len() > 1) | (datafs['ALT'].str.len() > 1)]
	##foetus heterozygote avec alternate allele provenant du père
	filter_foetus = foetus.loc[(~foetus['POS'].isin(mother['POS'].to_list()))]
	print("#[INFO] Length alternate variant provide by father (mother is homozygote for reference allele)", len(filter_foetus))
	
	#foetus heterozygote avec alternate allele provenant de la mère sachant variant homozygote et père ayant donné un allèle de reference
	homotmp = mother.loc[mother['ASG2104747'].str.partition(':')[0] == "1/1"]
	homo = homotmp['POS'].to_list()
	print("#[INFO] Length allele from homozygote alternate variant provide by mother (father gave ref allele)", len(homotmp))

	filter_foetus_homo = foetus.loc[(foetus['POS'].isin(homo)) & ((foetus['FCL2104751'].str.partition(':')[0] == "1/0") | (foetus['FCL2104751'].str.partition(':')[0] == "0/1"))]
	return filter_foetus, filter_foetus_homo

def processtsv(mother, father, foetus):
	dataframe_list = [ mother, father, foetus]
	for datafs in dataframe_list:
		datafs = datafs.loc[datafs['varType'] == 'substitution']
	##foetus heterozygote avec alternate allele provenant du père
	filter_foetus = foetus.loc[(~foetus['start'].isin(mother['start'].to_list()))]
	print("#[INFO] Length alternate variant provide by father (mother is homozygote for reference allele)", len(filter_foetus))
	
	#foetus heterozygote avec alternate allele provenant de la mère sachant variant homozygote et père ayant donné un allèle de reference
	homotmp = mother.loc[mother['zygosity'] == "hom"]
	homo = homotmp['start'].to_list()
	print("#[INFO] Length allele from homozygote alternate variant provide by mother (father gave ref allele)", len(homotmp))

	filter_foetus_homo = foetus.loc[(foetus['start'].isin(homo)) & (foetus['zygosity'] == 'het')]
	return filter_foetus, filter_foetus_homo
	


def globalfilter(dataframe):
	"""
	Prenatal Testing. BioTech 2021, 10, 17.https://doi.org/10.3390/biotech10030017, parameters read depth and MAF SNP should be common enough to be detected
	Sims, D.; Sudbery, I.; Ilot, N.E.; Heger, A.; Ponting, C.P. Sequencing depth and coverage: Key considerations in genomic analyses.
Nat. Rev. Genet. 2014, 15, 121–132.
	MAF > 5% and got dbSNP ID (could also use gnomAD stats)
	"""
	
	df = dataframe['rsMAF'].str.replace(',', '.').astype('float')
	#MAF dbsnp 5% so common variant in population, and got dbsnp variant btw
	filter = df.loc[(df['rsMAF'] > 0.05) & (~df['rsId'].isnull())]
	return filter
	
#def estimateFF(filter, foetus):
#	#FOR VCF #DEPRECATED
#	AF_list = {} 
#	VAF_list = []
#	for j, var in filter.iterrows():
#		AF_list[var['POS']] = {}
#		for i, fields in enumerate(var['FORMAT'].split(':')):
#			AF_list[var['POS']][fields] = var[foetus].split(':')[i]
#		for  features in var['INFO'].split(';'):
#			if '=' in features:
#				keys = features.split('=')[0]
#				value = features.split('=')[1]
#				AF_list[var['POS']][keys] = value
#			else:
#				AF_list[var['POS']]['infos'] = features
#
#
#	for var in AF_list.values():
#		if 'VAF' in var:
#			VAF_list.append(float(var.get('VAF')))
#		else:
#			VAF_list.append(float(var.get('AF')))
#	VAF = average(VAF_list)
#
#	VAF_list = []
#	for var in AF_list.values():
#		if 'VAF' in var:
#			VAF_list.append(float(var['VAF']))
#	VAF = average(VAF_list)
#	print("#[INFO] VAF average: ", VAF)
#	return VAF

def estimateFF(filter, foetus):
	VAF = filter['alleleFrequency'].str.replace(',', '.').astype('float').mean()
	return VAF

def dataframetoregions(dataframe, save=False):
	bed = dataframe.get(['chr', 'start', 'end'])
	if len(bed.index) == 0:
		print("ERROR col are missing from  file exit")
		exit()
	#bed.set_axis(['#CHROM', 'START', 'END'], axis=1, inplace=True)
	if save:
		bed.to_csv(save, index=False, header=False, sep='\t')
	return bed 

def getUMI(bamfile, position):
	'''
	get number of different UMI carring alternate base for a given position #TODO
	
	'''
	zscore = {}
	#bam = pysam.AlignmentFile(bamfile, 'rb')
	for i, variants in position.iterrows(): #"chr1:876499-876499"
		print(variants)
		zscore[variants['START']] = []
		pos = str(variants['#CHROM'])+':'+str(variants['START'])+'-'+str(variants['END'])
		read = pysam.view(bamfile, pos)
		for reads in read.split('\n'):
			fields = reads.split('\t')
			#print(type(fields[0]))
			for items in fields:
				if items.startswith('RX') and items not in zscore[variants['START']]:
					print(items)
					zscore[variants['START']].append(items.split(':')[-1])
			
		print(zscore)
		exit()
	return 
	
@lru_cache
def vcfTodataframe(file, rheader=False):

	'''
	Take in input vcf file, or tsv and return a dataframe
	I"m gonna build my own vcf parser et puis c'est tout
	return 3 Dataframe, full, only sample, only info
	'''
	name, extension = os.path.splitext(file)
	
	header = []
	variants_tmp = []
	variants = []
	with open(file) as f:
		for lines in f:
			if lines.startswith('#'):
				header.append(lines.strip())
			else:
				variants_tmp.append(lines)
	
	col = header[-1].strip().split('\t')
	for v in variants_tmp:
		variants.append(v.strip().split('\t'))
	
	dfVar = pd.DataFrame(columns=col)

	print("#[INFO] Whole VCF to Dataframe")
	for i, var in enumerate(variants):
		print_progress_bar(i, len(variants)-1)
		rows = pd.Series(var, index=dfVar.columns)
		dfVar = dfVar.append(rows, ignore_index=True)
	
	print('\n')
	if rheader:
		return dfVar, header
	else:
		name = os.path.basename(file).split('.')[0]
		dfVar.to_pickle(osj(os.path.dirname(file), name+'.pickle'))
		return dfVar

def tsvTodataframe(file):
	df = pd.read_csv(file, sep='\t', skiprows=2, header=0)
	df.to_pickle(file+'.pickle')
	return df

def parseargs():
	parser = argparse.ArgumentParser(description="TODO")
	parser.add_argument("-d", "--dad", type=str, help="Absolute path of vcf variant from father")
	parser.add_argument("-f", "--foetus", type=str, help="Absolute path of vcf variant from cell free DNA, (maternal blood)")
	parser.add_argument("-m", "--mum", type=str, help="Absolute path of vcf variant from mother")
	parser.add_argument("-t", "--type", default='tsv', type=str, help="vcf or tsv, default tsv")
	args = parser.parse_args()
	return args

def main():
	args = parseargs()
	ffname = os.path.basename(args.foetus).split('.')[0]
	filter_father, filter_mother = preprocess(args.mum, args.dad, args.foetus, args.type)
	VAFp = estimateFF(filter_father, ffname)
	VAFm = estimateFF(filter_mother, ffname)
	print(dataframetoregions(filter_father))
	#getUMI("/home1/data/STARK/data/DPNI/trio/TWIST/FCL2104751.bwamem.bam", pos)
	#FF = VAFp + (1 - VAFm)
	#print("#[INFO] Estimation of Foetal fraction : ", FF)

if __name__ == "__main__":
	main()
	#getUMI()