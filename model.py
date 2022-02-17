#!/usr/bin/python

"""
Python 3.8
SVM model to predict foetal fraction in DPNI analysis
"""

from scipy.stats import pearsonr
from sklearn import svm
from sklearn.externals import joblib
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
import pandas as pd
import argparse
import numpy as np
import os

class model():
	def __init__(self, kernel, verbose):
		self.model = svm.SVR(
        kernel=kernel,
        C=1.0,
        epsilon=0.01,
        tol=0.001,
        verbose=verbose
    )










###################### MAIN and ARGS

def main():
	args = parseargs()
	dataset_df = pd.read_table(args.dataset, index_col=args.index_col)
	"""
	return nothing
	"""
	verbose = False
	svm_model = model(kernel='linear', verbose=verbose)

	return

def parseargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--dataset", )

	args = parser.parse_args()
	return args



if __name__ == '__main__':
	main()


