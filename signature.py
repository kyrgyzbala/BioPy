from collections import OrderedDict
import sys
import numpy as np

def generate_dict(word_length):
	nt_set = ['A','T','C','G']

	if word_length==4:
		out_dict = {n1+n2+n3+n4:0 for n1 in nt_set for n2 in nt_set for n3 in nt_set for n4 in nt_set}
	elif word_length==3:
		out_dict = {n1+n2+n3:0 for n1 in nt_set for n2 in nt_set for n3 in nt_set}
	elif word_length==2:
		out_dict = {n1+n2:0 for n1 in nt_set for n2 in nt_set }
	elif word_length==1:
		out_dict = {n1:0 for n1 in nt_set}
	out_dict = OrderedDict(sorted(out_dict.items()))
	return out_dict

def get_frequencies(seq,word_length=4,log_file=[],chrom_name=[]):
	word_dict = generate_dict(word_length)
	for i in range(len(seq)-word_length-1):
		try:
			word_dict[str(seq[i:i+word_length].seq)]+=1
		except Exception, msg:
			if log_file:
				log_file.write('%s\n'%chrom_name)
				log_file.write(msg)
				log_file.write('\n')
			else:
				pass
	word_dict = OrderedDict(sorted(word_dict.items()))
	return word_dict

def read_frequency(f):
	out_freq = {}
	for l in open(f).readlines():
		parts = l.strip().split(':')
		out_freq[parts[0]]=parts[1]
	out_freq = OrderedDict(sorted(out_freq.items()))
	return out_freq

def read_frequencies(f):
	out_array = []
	for l in open(f).readlines():
		parts = l.strip().split()
		if parts==[]:
			continue
		out_array.append(parts[1:])
	return out_array

def read_window_frequencies(f):
	key_map = {}
	in_file = open(f)
	l = in_file.readline()
	cnt = 0
	freq_matrix = None
	while l:
		parts = l.strip().split()
		if l.strip()=='' or parts==[]:
			continue
		key_map[parts[0]]=cnt
		freqs = np.asarray(parts[1:], dtype=np.float)
		if freq_matrix==None:
			freq_matrix = np.matrix(freqs)
		else:
			freq_matrix = np.vstack([freq_matrix, freqs])
		cnt+=1
		l = in_file.readline()
	key_map = OrderedDict(sorted(key_map.items()))
	return key_map, freq_matrix

def merge_frequencies(freqs):
	out_freq = {}
	l = len(freqs[0])
	m = len(freqs)
	for i in range(l):
		k = freqs[0].keys()[i]
		out_freq[k]=0
		for j in range(m):
			out_freq[k] += int(freqs[j].values()[i])
	out_freq = OrderedDict(sorted(out_freq.items()))
	return out_freq

def write_freq_file(freqs, file_name):
	fout = open(file_name,'w')
	out_fmt = "%s:%f\n"
	[fout.write(out_fmt%(k,v)) for (k,v) in freqs.items()]
	fout.close()

def write_multiple_freqs_to_file(freqs, file_name, write_words_column=True):	
	fout = open(file_name,'w')
	freqs_size = len(freqs)
	vector_size = len(freqs[0].items())
	
	out_fmt = "%s\t" if write_words_column else ""
	out_fmt += '\t'.join(['%f' for i in range(freqs_size)])
	out_fmt += '\n'
	for i in range(vector_size):
		raw_values = [freqs[j].values()[i] for j in range(freqs_size)]
		if write_words_column:
			key=freqs[0].keys()[i]
			raw_values = [key]+raw_values
			fout.write(out_fmt%tuple(raw_values))
		else:
			fout.write(out_fmt%tuple(raw_values))
	fout.close()

def get_window_frequencies_old(gSource, word_size=4, window_size=5000, step_length=500):
	total_length = len(gSource.seq)
	freqs = []
	for i in range(0,total_length,step_length):
		wSeq = gSource[i:i+window_size]
		wFreqs = get_frequencies(wSeq, word_size)
		freqs.append(wFreqs)
	return freqs

def get_window_frequencies(gSource, word_size=4, window_size=5000, step_length=500):
	"""This is the faster version of get_window_frequencies_old.
	The idea came from Lorinc Pongor.
	The drawback is that you allocate dictionary instances for all windows at once,
	making it memory inefficient. But speed is faster in an order of magnitude"""
	total_length = len(gSource.seq)
	freqs = [generate_dict(word_size) for i in range(0,total_length,step_length)]
	# indices to show which frequencies are subject to increase
	start, end = 0, 1
	for i in range(0,total_length-word_size):
		if i!=0 and i%step_length==0:
			if i<window_size:
				end+=1
			else:
				start+=1
				if end!=len(freqs):
					end+=1
		nt_key = str(gSource.seq[i:i+word_size])
		for j in range(start,end):
			if freqs[j].has_key(nt_key):
				freqs[j][nt_key]+=1
	return freqs

def get_relative_abundance(freq, sub_freq, word_size=None):
	out_dict = {}
	
	freqs_sum = np.sum(np.asarray(freq.values(), dtype=np.float))
	sub_freqs_sum = np.sum(np.asarray(sub_freq.values(), dtype=np.float))
	if not word_size:
		word_size = len(freq.keys()[0])
	
	for k,v in freq.items():
		sub_f1 = float(sub_freq[k[:-1]])
		sub_f2 = float(sub_freq[k[1:]])
		abundance = (float(v)/freqs_sum) / (sub_f1*sub_f2/sub_freqs_sum**2)
		out_dict[k]=abundance
	out_dict = OrderedDict(sorted(out_dict.items()))
	return out_dict

def get_window_relative_abundances(freq_file, sub_freq_file, word_size):
	freqs_key_map, freqs_data_mat = read_window_frequencies(freq_file)
	sub_freqs_key_map, sub_freqs_data_mat = read_window_frequencies(sub_freq_file)
	
	freqs_data_mat = freqs_data_mat/np.sum(freqs_data_mat, axis=0)
	sub_freqs_data_mat = sub_freqs_data_mat/np.sum(sub_freqs_data_mat, axis=0)

	result = [generate_dict(word_size) for i in range(freqs_data_mat.shape[1])]
	
	freqs_keys = freqs_key_map.keys()
	freqs_values = freqs_key_map.values()
	for i in range(freqs_data_mat.shape[1]):
		for j in range(freqs_data_mat.shape[0]):
			f_key = freqs_keys[j]
			s1_key = f_key[:-1]
			s2_key = f_key[1:]
			freq_value = freqs_data_mat[j,i]
			first_sub_value = sub_freqs_data_mat[sub_freqs_key_map[s1_key],i]
			second_sub_value = sub_freqs_data_mat[sub_freqs_key_map[s2_key],i]
			
			if freq_value==0:
				continue
			relative_abundance = freq_value/(first_sub_value*second_sub_value)
			result[i][f_key] = relative_abundance
	
	return result

def get_window_relative_abundances_2(freq_file, sub_freq_file, word_size):
	# This version normalizes the frequency with genome wise frequencies,
	# as opposed to window wise frequencies, which is the case in get_window_relative_abundances

	freqs_key_map, freqs_data_mat = read_window_frequencies(freq_file)
	sub_freqs = read_frequency(sub_freq_file)
	
	freqs_data_mat = freqs_data_mat/np.sum(freqs_data_mat, axis=0)

	s_sum = np.sum(np.asarray(sub_freqs.values(), dtype=np.float))
	for k in sub_freqs.keys():
		sub_freqs[k] = float(sub_freqs[k])/s_sum

	result = [generate_dict(word_size) for i in range(freqs_data_mat.shape[1])]
	
	freqs_keys = freqs_key_map.keys()
	freqs_values = freqs_key_map.values()
	for i in range(freqs_data_mat.shape[1]):
		for j in range(freqs_data_mat.shape[0]):
			f_key = freqs_keys[j]
			s1_key = f_key[:-1]
			s2_key = f_key[1:]
			freq_value = freqs_data_mat[j,i]
			first_sub_value = sub_freqs[s1_key]
			second_sub_value = sub_freqs[s2_key]
			
			if freq_value==0:
				continue
			relative_abundance = freq_value/(first_sub_value*second_sub_value)
			result[i][f_key] = relative_abundance
	
	return result
