#!/usr/bin/python

"""
Python 3.8
SVM model to predict foetal fraction in DPNI analysis
"""

from scipy.stats import pearsonr
from sklearn import preprocessing
from sklearn import svm
from sklearn.datasets import make_classification
from sklearn.externals import joblib
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
import pandas as pd
import argparse
import numpy as np
import os

INDEX_COL = 'id'
SEQFF_COL = 'seqff'
TARGET_COL = 'target'

class model():
	def __init__(self, kernel, dataset, verbose, kfold):
		X, y = make_classification(random_state=42)

		self.model = make_pipeline(StandardScaler(), SVR(
        kernel=kernel,
        C=1.0,
        tol=0.001,
        verbose=verbose,
		class_weight='balanced'
    ))
		self.dataset = pd.read_csv(dataset, sep='\t', header=0, index_col=None)
		self.kfold = kfold

	def preprocess(self):
		i = 0
		kf = KFold(self.kfold)
		for train_index, test_index in kf.split(self.dataset):
			train_df = self.dataset.iloc[train_index]
			test_df = self.dataset.iloc[test_index]
			print('#[INFO] kfold %d training set shape: %s' % (i, str(train_df.shape)))
			print('#[INFO] kfold %d testing set shape: %s' % (i, str(test_df.shape)))

			train_seqff = train_df.pop('seqff')
			#test_seqff = test_df.pop('seqff')

			train_target = train_df.pop('target')
			#test_target = test_df.pop('target')
			
			#Standardization
			scaler = preprocessing.StandardScaler().fit(train_df)
			X_scaled = scaler.transform(train_df)


			print(scaler.mean_)
			i += 1
			return scaler



###################### MAIN and ARGS

def main():
	args = parseargs()
	#dataset_df = pd.read_table(args.dataset, index_col=args.index_col)
	"""
	return nothing
	"""
	verbose = False
	svm_model = model('linear', args.dataset, verbose, args.kfold)
	svm_model.preprocess()
	

def parseargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--dataset", type=str, required=True, help="Absolute path of dataset")
	parser.add_argument("-k", "--kfold", type=int, default=5, help="kfold for cross validation in testing step")
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	main()


