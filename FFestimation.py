"""
Mother homozigous for reference allele and foetus got alternate allele from the father 
VAFp = 1/2 FF

Mother homozygous for a alternate allele and fetus got reference allele
VAFm = 1 - (0,5 * FF)

FF = VAFpaverage+(1-VAFm average)
https://doi.org/10.3390/biotech10030017

"""


from calendar import monthrange
import faulthandler
from functools import lru_cache
from os.path import join as osj
import argparse
import os
from unicodedata import numeric
import pandas as pd
import pysam
import re
import subprocess
import sys

#git combo FF, dossier TEST TODO

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

class Process:
	def __init__(self, mother, father, foetus):
		self.mother = self.filterquality(mother)
		self.father = self.filterquality(father)
		self.foetus = foetus

	def preprocess(self, filetype):
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
		mother = pd.read_pickle(self.mother+'.pickle')
		father = pd.read_pickle(self.father+'.pickle')
		foetus = pd.read_pickle(self.foetus+'.pickle')
		
		filter_foetus, filter_foetus_homo = self.processtsv(mother, father, foetus)
	
		return filter_foetus, filter_foetus_homo

	def filterquality(df):
		filter = df.loc[(df['varReadPercent'] > 30) & (df['varReadDepth'] > 30) & (df['QUALphred'] >300)]
		return df

	def processvcf(self, mother, father, foetus):
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
	
	def processtsv(self, mother, father, foetus):
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
		print(filter_foetus.head())
		print(filter_foetus_homo.head())
		return filter_foetus, filter_foetus_homo

class Paternalidentification(Process):
	def __init__(self, mother, father , foetus):
		super().__init_(mother, father, foetus)

	def subtract_maternal(foetus, mother):
		'''
		input: dataframe of maternal blood (unknown FF) and mother as Index Case --> 100 %
		output: dataframe keep only variants from Father, foetus and potentially denovo
		'''
		maternal_id_SNV = mother['variantID'].to_list()
		maternal_id_InDels = mother['cNomen'].to_list()
		df_subtract = foetus.loc[(~foetus['variantID'].isin(maternal_id_SNV)) | (~foetus['cNomen'].isin(maternal_id_InDels))]

		#foetus_from_m = foetus.loc[(foetus['variantID'].isin(maternal_id_SNV)) | (foetus['cNomen'].isin(maternal_id_InDels))]
		return df_subtract

	def identify_paternal(foetus_filter, father):
		'''
		input: fetal dataframe where commom variant with mother have been discarded
		output: try to estimate denovo variant and FF
		'''
		if not foetus_filter.dtypes['alleleFrequency'] == 'float64':
			foetus_filter['alleleFrequency'] = foetus_filter['alleleFrequency'].str.replace(',', '.').astype('float')
		paternal = foetus_filter.loc[(foetus_filter['alleleFrequency'] >= 3.5) & (foetus_filter['alleleFrequency'] <= 7.5)]
		#Probably denovo
		denovo = foetus_filter.loc[(foetus_filter['alleleFrequency'] > 7.5)]
		
		return paternal, denovo

	def getFF(foetus_from_p, foetus_filter):
		#df = pd.concat([foetus_from_m, foetus_from_p], axis=1, ignore_index=True)
		var = []
		for i, variants in foetus_from_p.iterrows():
			if variants['variantID'] in foetus_filter['variantID']:
				if variants['zygosity'] == 'hom' and variants['alleleFrequency'] >:
				df = foetus_from_p.loc[foetus_from_p['variantID']]
		#remove double line for homozygote variant 
		df_filter = df.drop_duplicates(subset=['variantID', 'cNomen'], keep='First')




#PURE FETAL FRACTION ESTIMATION based on publciation below
def globalfilter(df, rmasker, output, pattern):
	"""
	Prenatal Testing. BioTech 2021, 10, 17.https://doi.org/10.3390/biotech10030017, parameters read depth and MAF SNP should be common enough to be detected
	Sims, D.; Sudbery, I.; Ilot, N.E.; Heger, A.; Ponting, C.P. Sequencing depth and coverage: Key considerations in genomic analyses.
Nat. Rev. Genet. 2014, 15, 121–132.
	MAF > 5% and got dbSNP ID (could also use gnomAD stats)
	"""
	print("#[INFO] dtypes ", df.dtypes['rsMAF'])
	if not df.dtypes['rsMAF'] == 'float64':
		df['rsMAF'] = df['rsMAF'].str.replace(',', '.').astype('float')
	df['totalReadDepth'] = df['totalReadDepth'].astype('int')

	#MAF dbsnp 5% so common variant in population, and got dbsnp variant btw
	tmp = df.loc[(df['rsMAF'] > 0.05) & (~df['rsId'].isnull()) & (df['totalReadDepth'] > 30)]
	tmp.loc[:, 'chr'] = 'chr'+tmp.chr

	bedname = osj(output, pattern)
	foetusfilter = osj(output, pattern+'out')

	#repeatmasker
	dataframetoregions(tmp, bedname)
	print("#[INFO] BEDTOOLS Processing ... ")
	systemcall("bedtools intersect -v -a "+bedname+" -b "+rmasker+" -wa > "+foetusfilter)
	foetus_df = pd.read_csv(foetusfilter, sep='\t')
	foetus_df = foetus_df.set_axis(['chr','start', 'stop'], axis='columns') 
	print(tmp.head())
	print(foetus_df.head())
	filter = tmp.loc[tmp['start'].isin(foetus_df['start'].to_list())]

	return filter
	
def getsensibility(dfoetal, dfparents):
	return

def estimateFF(filter):
	VAF = filter['alleleFrequency'].str.replace(',', '.').astype('float').mean()
	return VAF

def dataframetoregions(dataframe, save):
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
	parser.add_argument("-r", "--rmasker", default='/home1/data/STARK/data/DPNI/trio/repeatmasker.bed', type=str, help="repeatmaskerfile")
	parser.add_argument("-o", "--output", default='/home1/data/STARK/data/DPNI/trio/TWIST', type=str, help='name of outputfolder')
	
	args = parser.parse_args()
	return args

def main():
	args = parseargs()
	ffname = os.path.basename(args.foetus).split('.')[0]
	
	
	filter_father, filter_mother = Process(args.mum, args.dad, args.foetus, args.type)
	

	VAFp = estimateFF(globalfilter(filter_father, args.rmasker, args.output, 'filter_father'))
	print("#[INFO] VAFp ", VAFp)

	VAFm = estimateFF(globalfilter(filter_mother, args.rmasker, args.output, 'filter_mother'))
	print("#[INFO] VAFm ", VAFm)

	#VAFp = estimateFF(filter_father, ffname)
	#VAFm = estimateFF(filter_mother, ffname)
	#print(dataframetoregions(filter_father))
	#getUMI("/home1/data/STARK/data/DPNI/trio/TWIST/FCL2104751.bwamem.bam", pos)
	FF = VAFp + (1 - VAFm)
	print("#[INFO] Estimation of Foetal fraction : ", FF)

if __name__ == "__main__":
	main()
	#getUMI()