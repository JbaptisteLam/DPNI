#!/usr/bin/python

"""
Mother homozigous for reference allele and foetus got alternate allele from the father 
VAFp = 1/2 FF

Mother homozygous for a alternate allele and fetus got reference allele
VAFm = 1 - (0,5 * FF)

FF = VAFpaverage+(1-VAFm average)
https://doi.org/10.3390/biotech10030017

"""


from calendar import monthrange
from functools import lru_cache
from os.path import join as osj
import argparse
import os
import pandas as pd
import pysam
import re
import subprocess
import sys

#git combo FF, dossier TEST TODO

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
	'''
	*passing command to the shell*
	*return list containing stdout lines*
	command - shell command (string)
	'''
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def average(iterable):
	return round(sum(iterable) / len(iterable), 4)

def getheader(file):
	with open(file, 'r') as f:
		header = []
		for lines in f:
			line = lines.strip()
			header.append(line)
	return header

class Process:
	def __init__(self, mother, father, foetus, filetype, output, filterqual):
		m, d, f = self.preprocess(mother, father, foetus, filetype, output, filterqual)
		self.mother = m
		self.father = d
		self.foetus = f
		self.filetype = filetype
		print("#[INFO] Length mother variants ", len(self.mother.index))
		print("#[INFO] Length father variants ", len(self.father.index))
		print("#[INFO] Length foetus variants ", len(self.foetus.index))
		self.header = {}
		fam = ['mother', 'father', 'foetus']
		for i, path in locals().items():
			if i in fam:
				name = os.path.basename(path).split('.')[0]
				self.header[name] = getheader(path)
				
	
	def preprocess(self, mother, father, foetus, filetype, output, filterqual):
		'''
		input: tsv of variants from trio family
		output: pickle object for whole family
		'''
		print("#[INFO] Locals init items ", locals())
		#make a pickle object in output folder if it not exists, even where the vcf are located
		for j, values in locals().items():
			if j != 'filetype' and j != 'self':
				path = os.path.basename(values)
				filename = osj(output, path+'.pickle')
				if not os.path.exists(filename):
					if filetype == 'vcf':
						vcfTodataframe(filename)
					else:
						tsvTodataframe(values)
		
		#FOR TSV only
		if filterqual:
			mother = self.filterquality(pd.read_pickle(mother+'.pickle'))
			father = self.filterquality(pd.read_pickle(father+'.pickle'))
		else:
			mother = pd.read_pickle(mother+'.pickle')
			father = pd.read_pickle(father+'.pickle')

		foetus = pd.read_pickle(foetus+'.pickle')
		
		#filter_foetus, filter_foetus_homo = self.processtsv(self.filterquality(mother), self.filterquality(father), foetus)
		return mother, father, foetus

	def filterquality(self, df):
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
	def __init__(self, mother, father , foetus, filetype, output, filterqual):
		super().__init__(mother, father , foetus, filetype, output, filterqual)

	def subtract_maternal(self):
		'''
		input: dataframe of maternal blood (unknown FF) and mother as Index Case --> 100 %
		output: dataframe keep only variants from Father, foetus and potentially denovo
		'''
		
		maternal_id_SNV = self.mother['variantID'].to_list()
		maternal_id_InDels = self.mother['cNomen'].to_list()
		df_subtract = self.foetus.loc[(~self.foetus['variantID'].isin(maternal_id_SNV)) | (~self.foetus['cNomen'].isin(maternal_id_InDels))]

		#foetus_from_m = foetus.loc[(foetus['variantID'].isin(maternal_id_SNV)) | (foetus['cNomen'].isin(maternal_id_InDels))]
		print("#[INFO] Length after removing mother variations ", len(df_subtract.index))
		return df_subtract

	def identify_paternal(self, foetus_filter):
		'''
		input: fetal dataframe where commom variant with mother have been discarded
		output: try to estimate denovo variant and FF
		'''
		if not foetus_filter.dtypes['alleleFrequency'] == 'float64':
			foetus_filter.loc[:, 'alleleFrequency'] = foetus_filter['alleleFrequency'].str.replace(',', '.').astype('float')
		paternal = foetus_filter.loc[(foetus_filter['alleleFrequency'] >= 0.035) & (foetus_filter['alleleFrequency'] <= 0.075)]
		#Probably denovo
		denovo = foetus_filter.loc[(foetus_filter['alleleFrequency'] > 0.075)]
		print("#[INFO] Variants comming from father between 4 and 15percent of FF ", len(paternal.index))
		print("#[INFO] After all filter, probably denovo variants (normally high for now cuz using 100percent ES foetus ", len(denovo.index))
		
		return paternal, denovo

	def getFF(self, paternal, foetus_filter):
		#equal to paternal if above func TODO
		#df = pd.concat([foetus_from_m, foetus_from_p], axis=1, ignore_index=True)
		#var = []
		#for i, variants in paternal.iterrows():
		#	if variants['variantID'] in foetus_filter['variantID'].to_list():
		#		if variants['alleleFrequency'] < 0.075:
		#			var.append(variants)
		#
		#print("#[INFO] LEN VAR ",len(var))
		##print(var[0:5])
		##df = foetus_from_p.loc[foetus_from_p['variantID']]
		#df = pd.DataFrame(var)
		#remove double line for homozygote variant 
		#df_filter = df.drop_duplicates(subset=['variantID', 'cNomen'], keep='Last')

		print("#[INFO] FF estimation: ", average(paternal['alleleFrequency']))
		paternal.to_csv('test_ff.tsv', sep='\t', columns=paternal.columns, index=False)
	
	def main_paternal(self):
		sub = self.subtract_maternal()
		paternal_var, denovo = self.identify_paternal(sub)
		self.getFF(paternal_var, sub)

class homozygotebased(Process):
	def __init__(self, mother, father , foetus, filetype, output, filterqual):
		super().__init__(mother, father , foetus, filetype, output, filterqual)


	#PURE FETAL FRACTION ESTIMATION based on publciation below
	def globalfilter(self, df, rmasker, output, pattern):
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
		self.dataframetoregions(tmp, bedname)
		print("#[INFO] BEDTOOLS Processing ... ")
		systemcall("bedtools intersect -v -a "+bedname+" -b "+rmasker+" -wa > "+foetusfilter)
		foetus_df = pd.read_csv(foetusfilter, sep='\t')
		foetus_df = foetus_df.set_axis(['chr','start', 'stop'], axis='columns') 
		print(tmp.head())
		print(foetus_df.head())
		filter = tmp.loc[tmp['start'].isin(foetus_df['start'].to_list())]

		return filter

	def getsensibility(self, dfoetal, dfparents):
		#TODO
		return

	def estimateFF(self, filter):
		VAF = filter['alleleFrequency'].str.replace(',', '.').astype('float').mean()
		return VAF

	def dataframetoregions(self, dataframe, save):
		bed = dataframe.get(['chr', 'start', 'end'])
		if len(bed.index) == 0:
			print("ERROR col are missing from  file exit")
			exit()
		#bed.set_axis(['#CHROM', 'START', 'END'], axis=1, inplace=True)
		if save:
			bed.to_csv(save, index=False, header=False, sep='\t')
		return bed 

	def getUMI(self, bamfile, position):
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

def parseargs(): #TODO continue subparser and add ML docker in script
	parser = argparse.ArgumentParser(description="TODO")
	subparsers = parser.add_subparsers()
	parser_a = subparsers.add_parser('standard',type=str, help='Standard use of FFestimation.py TRIO', dest='command')
	parser_b = subparsers.add_parser('paternal', type=str, help='Estimation of FF based on variant comming from father TRIO', dest='command')

	parser_c = subparsers.add_parser('seqff', type=str, help='ML model to estimate FF based on seqFF and regression model using read depth profile Maternal BLOOD only', dest='command')

	parser.add_argument("-q", "--quality", action='store_false', type=bool, help="Activate filtering on parents variants file, discarded variants with varReadDepth < 30, varReadPercent < 30 and qual Phred < 300, default True, set arg to remove filtering")
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
	#ffname = os.path.basename(args.foetus).split('.')[0]
	if args.command == 'standard':
		#1) From CDC of DPNI study
		pi = Paternalidentification(args.mum, args.dad, args.foetus, args.type, args.output, args.quality)
		pi.main_paternal()
	
	#2) From combo_FF seqFF modele
	elif args.command == 'seqff':
		print("In developpement exit !")
		exit()
		
	#3) From publication with UMI standard deviation
	elif args.command == 'paternal':
		p = homozygotebased(args.mum, args.dad, args.foetus, args.type, args.output, args.quality)
		VAFp = p.estimateFF(p.globalfilter(args.dad, args.rmasker, args.output, 'filter_father'))
		VAFm = p.estimateFF(p.globalfilter(args.mum, args.rmasker, args.output, 'filter_mother'))
		FF = VAFp + (1 - VAFm)
		print("#[INFO] Estimation of Foetal fraction : ", FF)

	#getUMI("/home1/data/STARK/data/DPNI/trio/TWIST/FCL2104751.bwamem.bam", pos)
	#FF = VAFp + (1 - VAFm)
	#print("#[INFO] Estimation of Foetal fraction : ", FF)

if __name__ == "__main__":
	main()