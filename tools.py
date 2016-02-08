completes_dir = '/home/sanjarbek/Data/NCBI/complete/'
drafts_dir = '/home/sanjarbek/Data/NCBI/draft/'

completes_map = '/home/sanjarbek/Data/info/map_complete_genomes.txt'
drafts_map = '/home/sanjarbek/Data/info/map_draft_genomes.txt'
all_map = '/home/sanjarbek/Data/info/map_all_genomes.txt'


def get_source_map(genome_type='all'):
	if genome_type=='all':
		return {l.strip().split('\t')[1]:l.split('\t')[0] for l in open(all_map).readlines()[1:]}
	elif genome_type=='complete':
		return {l.strip().split('\t')[1]:l.split('\t')[0] for l in open(completes_map).readlines()[1:]}
	elif genome_type=='draft':
		return {l.strip().split('\t')[1]:l.split('\t')[0] for l in open(drafts_map).readlines()[1:]}
