import numpy as np
import sys

def score(vector1,vector2, normalize):
	"""Calculate Kullback-Leibler divergence for two vectores
	Arrays can be passed as input .It's converted to numpy arrays expplicitly
	before starting calculations"""

	p = np.asarray(vector1, dtype=np.float)
	q = np.asarray(vector2, dtype=np.float)

	p = p/np.sum(p) if np.sum(p)!=1.0 else p
	q = q/np.sum(q) if np.sum(q)!=1.0 else q

	KL_score = np.sum(np.where(p*q!=0 , p * np.log(p / q), 0))
	KL_score = KL_score/np.sqrt(len(p)) if normalize else KL_score
	
	return KL_score

def symmetric_score(vector1,vector2, normalize=True):
	"""Gives (KL(A,B)+KL(B,A))/2 . Where A and B are input vectors (arrays)"""
	return (score(vector1,vector2,normalize)+score(vector2,vector1,normalize))/2.0