"""
Mother homozigous for reference allele and foetus got alternate allele from the father 
VAFp = 1/2 FF

Mother homozygous for a alternate allele and fetus got reference allele
VAFm = 1 - (0,5 * FF)

FF = VAFpaverage+(1-VAFm average)

"""


from functools import lru_cache
from os.path import join as osj
import argparse
import os
import pandas as pd
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

def preprocess(mum, dad, foetal):
	'''
	input: tsv of variants from trio family
	output: variant inside VAFp and VAFm statements in tsv format
	'''
	for j, values in locals().items():
		filename = osj(os.path.dirname(values), os.path.basename(values).split('.')[0]+'.pickle')
		if not os.path.exists(filename):
			vcfTodataframe(filename)

	mother = pd.read_pickle(osj(os.path.dirname(mum), os.path.basename(mum).split('.')[0]+'.pickle'))
	father = pd.read_pickle(osj(os.path.dirname(dad), os.path.basename(dad).split('.')[0]+'.pickle'))
	foetus = pd.read_pickle(osj(os.path.dirname(foetal), os.path.basename(foetal).split('.')[0]+'.pickle'))
	dataframe = [ mother, father, foetus]
	for datafs in dataframe:
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
	
def estimateFF(filter, foetus):
	AF_list = {} 
	VAF_list = []
	for j, var in filter.iterrows():
		AF_list[var['POS']] = {}
		for i, fields in enumerate(var['FORMAT'].split(':')):
			AF_list[var['POS']][fields] = var[foetus].split(':')[i]
		for  features in var['INFO'].split(';'):
			if '=' in features:
				keys = features.split('=')[0]
				value = features.split('=')[1]
				AF_list[var['POS']][keys] = value
			else:
				AF_list[var['POS']]['infos'] = features


	for var in AF_list.values():
		if 'VAF' in var:
			VAF_list.append(float(var.get('VAF')))
		else:
			VAF_list.append(float(var.get('AF')))
	VAF = average(VAF_list)

	VAF_list = []
	for var in AF_list.values():
		if 'VAF' in var:
			VAF_list.append(float(var['VAF']))
	VAF = average(VAF_list)
	print("#[INFO] VAF average: ", VAF)
	return VAF

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

def parseargs():
	parser = argparse.ArgumentParser(description="Filter tsv in DPNI and POOL context, basically tsv have to come from varank analysis ")
	parser.add_argument("-d", "--dad", type=str, help="Absolute path of vcf variant from father")
	parser.add_argument("-f", "--foetus", type=str, help="Absolute path of vcf variant from cell free DNA, (maternal blood)")
	parser.add_argument("-m", "--mum", type=str, help="Absolute path of vcf variant from mother")
	args = parser.parse_args()
	return args

def main():
	args = parseargs()
	ffname = os.path.basename(args.foetus).split('.')[0]
	filter_father, filter_mother = preprocess(args.mum, args.dad, args.foetus)
	VAFp = estimateFF(filter_father, ffname)
	VAFm = estimateFF(filter_mother, ffname)

	FF = VAFp + (1 - VAFm)
	print("#[INFO] Estimation of Foetal fraction : ", FF)

if __name__ == "__main__":
	main()