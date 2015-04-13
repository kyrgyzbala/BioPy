#   stmo 
import KL
import numpy as np
import os
import sys
import operator

class DbFreqItem:	
	def __init__(self, name,score):
		self.name = name
		self.score = score

def get_genome_KL_scores(query_freq, sort_results=True, db_freqs_path='/home/sanjarbek/Data/signatures/tetra/genomes/'):
	
	db_items = []
	db_freq_files = os.listdir(db_freqs_path)
	for f in db_freq_files:
		item_name = f.split('.')[0]
		db_freq = [l.split(':')[1].strip() for l in open(os.path.join(db_freqs_path,f)).readlines()]
		db_items.append(DbFreqItem(item_name,KL.symmetric_score(query_freq,db_freq)))

	if sort_results:
		db_items.sort(key=operator.attrgetter('score'))
	return db_items

