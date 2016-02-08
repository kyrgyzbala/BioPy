import MySQLdb as mdb
connection = mdb.connect(user="djuser",passwd="djuserpwd",db="QSportal")
import BioClasses as cl
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys

def setup_cursor():
    try:
        cursor = connection.cursor()
        return cursor
    except ConnectionDoesNotExist:
        print "database is not configured"
        return None

def get_topology_info(topid):
	cursor = setup_cursor()
	sqlcmd = """select td.Code, td.Pattern, o.Name
				from Topology t
				inner join TopologyDefinition td on t.TopologyDefinition_WID = td.WID
				inner join Topology_of_Organism too on too.Topology_WID = t.WID
				inner join Organism o on too.Organism_WID = o.WID
				where t.WID = %d"""

	cursor.execute(sqlcmd%topid)
	return cursor.fetchall()[0]

def get_gene_type(gid):
	cursor = setup_cursor()
	sqlcmd = """select gt.GeneShortType
				from Topology t				
				inner join GeneType gt on t.GeneType_WID=gt.WID
				where t.GeneID = %d"""

	cursor.execute(sqlcmd%gid)
	return cursor.fetchall()[0][0]

def get_topology_pattern(gid,accession=False):

	cursor = setup_cursor()
	if accession:
		sqlcmd = """SELECT td.pattern
					FROM QSportal.GI_to_Ref gtr
					inner join Topology t on t.GeneID=gtr.gid
					inner join TopologyDefinition td on td.WID = t.TopologyDefinition_WID
					where ref='%s';"""
	else:
		sqlcmd = """SELECT td.pattern
					FROM Topology t 
					inner join TopologyDefinition td on td.WID = t.TopologyDefinition_WID
					where t.GeneID='%s';"""
	sqlcmd = sqlcmd%gid
	cursor.execute(sqlcmd)
	res = cursor.fetchall()[0][0]
	return res

def gid2gene_list(gids,as_map=False):
	cursor = setup_cursor()
	sqlcmd = """SELECT source, pFrom, pTo, PID
				FROM GenePool 
				WHERE PID in (%s)"""
	sqlcmd = sqlcmd%",".join([str(g) for g in gids])
	cursor.execute(sqlcmd)
	all_res=cursor.fetchall()
	if not as_map:
		genes = [cl.Gene(source=g[0],pFrom=g[1],pTo=g[2],gid=g[3]) for g in all_res]
	else:
		genes = {g[3]:cl.Gene(source=g[0],pFrom=g[1],pTo=g[2],gid=g[3]) for g in all_res}
	return genes

def gid2acc(gids):
	cursor = setup_cursor()
	sqlcmd = """SELECT accession 
				FROM QSportal.GI_to_Ref
				WHERE gid in (%s)"""
	sqlcmd = sqlcmd%",".join([str(g) for g in gids])
	cursor.execute(sqlcmd)
	return cursor.fetchall()


def gid2seq_list(gids,as_map=False):
	cursor = setup_cursor()
	sqlcmd = """SELECT PID,header, sequence
				FROM Sequence
				where PID in (%s)"""
	sqlcmd = sqlcmd%",".join([str(g) for g in gids])	
	cursor.execute(sqlcmd)
	all_res=cursor.fetchall()

	if not as_map:
		genes = [SeqRecord(Seq(g[2]),id=g[1]) for g in all_res]
	else:
		genes = {g[0]:SeqRecord(Seq(g[2]),id=g[1]) for g in all_res}

	return genes

def gid2organism(gid):
	cursor = setup_cursor()
	sqlcmd = """select gp.PID, gp.source, o.Name
				from Topology t
				inner join Topology_of_Organism too on too.Topology_WID=t.WID
				inner join Organism o on o.WID = too.Organism_WID
				inner join GenePool gp on gp.PID=t.GeneID
				where gp.PID=%s"""
	sqlcmd = sqlcmd%gid
	cursor.execute(sqlcmd)
	retval = cursor.fetchone()
	return [retval[1], retval[2]]


def gid2loci_list(gids):
	cursor = setup_cursor()
	sqlcmd = """SELECT gp.source,t.GeneID,s.pFrom, s.pTo
				FROM Topology t
				INNER JOIN (SELECT ti.WID, min(gpi.pFrom) as pFrom, max(gpi.pTo) as pTo
				            FROM Topology ti
				            INNER JOIN GenePool gpi on ti.GeneID=gpi.PID
				            GROUP BY ti.WID
				        ) s on s.WID = t.WID
				INNER JOIN GenePool gp on t.GeneID = gp.PID
				WHERE t.GeneID in (%s)"""
	sqlcmd = sqlcmd%",".join([str(g) for g in gids])
	cursor.execute(sqlcmd)
	genes = [cl.Gene(source=g[0],pFrom=g[2],pTo=g[3],gid=g[1]) for g in cursor.fetchall()]
	return genes

def get_soloRs(genome_type=None,use_ref=True):
	cursor = setup_cursor()
	sqlcmd = """select o.status, o.Name, gp.source, gr.Ref, gp.pFrom, gp.pTo, gp.strand, gp.PID
				from Topology t
				inner join TopologyDefinition td on t.TopologyDefinition_WID = td.WID
				inner join GenePool gp on t.GeneID = gp.PID
				inner join GI_to_Ref gr on gr.gid = t.GeneID
				inner join Topology_of_Organism too on too.Topology_WID = t.WID
				inner join Organism o on too.Organism_WID = o.WID
				where td.WID=5 and o.status IN (%s)"""
	if not genome_type:
		sqlcmd = sqlcmd%'1,2,3'
	elif genome_type==1:
		sqlcmd = sqlcmd%'1'
	else:
		sqlcmd = sqlcmd%'2,3'
	cursor.execute(sqlcmd)
	#Use accession number vs GID for ID
	if use_ref:
		genes = [cl.Gene(source=g[2],gid=g[3],pFrom=g[4],pTo=g[5],strand=g[6]) for g in cursor.fetchall()]
	else:
		genes = [cl.Gene(source=g[2],gid=g[7],pFrom=g[4],pTo=g[5],strand=g[6]) for g in cursor.fetchall()]
	return genes

def get_nonsoloRs(genome_type=None,use_ref=True):
	cursor = setup_cursor()
	sqlcmd = """select o.status, o.Name, gp.source, gr.Ref, gp.pFrom, gp.pTo, gp.strand, gp.PID
				from Topology t
				inner join TopologyDefinition td on t.TopologyDefinition_WID = td.WID
				inner join GenePool gp on t.GeneID = gp.PID
				inner join GI_to_Ref gr on gr.gid = t.GeneID
				inner join Topology_of_Organism too on too.Topology_WID = t.WID
				inner join Organism o on too.Organism_WID = o.WID
				where t.GeneType_WID=2 and td.WID<>5 and o.status in (%s)"""

	if not genome_type:
		sqlcmd = sqlcmd%'1,2,3'
	elif genome_type==1:
		sqlcmd = sqlcmd%'1'
	else:
		sqlcmd = sqlcmd%'2,3'
	
	cursor.execute(sqlcmd)
	#Use accession number vs GID for ID
	if use_ref:
		genes = [cl.Gene(source=g[2],gid=g[3],pFrom=g[4],pTo=g[5],strand=g[6]) for g in cursor.fetchall()]
	else:
		genes = [cl.Gene(source=g[2],gid=g[7],pFrom=g[4],pTo=g[5],strand=g[6]) for g in cursor.fetchall()]
	
	return genes

def get_fasta(ids):
	cursor = setup_cursor()
	sqlcmd = """select header, sequence
				from Sequence
				where PID in (%s)"""
	
	sqlcmd = sqlcmd%",".join([str(id) for id in ids]) if isinstance(ids,list) else sqlcmd%ids
	cursor.execute(sqlcmd)
	seqs = [SeqRecord(Seq(g[1]),id=g[0],description="") for g in cursor.fetchall()]
	return seqs

def get_genes(ids):
	cursor = setup_cursor()
	sqlcmd = """SELECT source, PID, strand, pFrom, pTo, Product 
				FROM GenePool
				where PID in (%s)"""
	
	sqlcmd = sqlcmd%",".join([str(id) for id in ids]) if isinstance(ids,list) else sqlcmd%ids
	cursor.execute(sqlcmd)
	genes = [cl.Gene(source=g[0],gid=g[1],strand=g[2],pFrom=g[3],pTo=g[4]) for g in cursor.fetchall()]
	return genes

def get_genes_accession(ids):
	cursor = setup_cursor()
	sqlcmd = """SELECT source, PID, strand, pFrom, pTo, Product 
				FROM GenePool 
				where PID in (
	            SELECT gid
	            FROM QSportal.GI_to_Ref
	            where accession in (%s))"""
	
	sqlcmd = sqlcmd%",".join(['\''+str(id)+'\'' for id in ids]) if isinstance(ids,list) else sqlcmd%ids
	cursor.execute(sqlcmd)
	genes = [cl.Gene(source=g[0],gid=g[1],strand=g[2],pFrom=g[3],pTo=g[4]) for g in cursor.fetchall()]
	return genes

def get_RR_topologies():
	cursor = setup_cursor()
	sqlcmd = """select gp.source, group_concat(t.GeneID SEPARATOR '\t')
				from Topology t
				inner join (
				            select t.WID
				            from Topology t
				            inner join GenePool gp on t.GeneID = gp.PID
				            inner join GeneType gt on t.GeneType_WID=gt.WID
				            where t.TopologyDefinition_WID=20
				            group by t.WID
				            having group_concat(gt.GeneShortType ORDER BY gp.pFrom DESC SEPARATOR '') = 'RR' ) t2 on t.WID = t2.WID
				inner join GI_to_Ref gr on gr.gid = t.GeneID
				inner join GenePool gp on t.GeneID = gp.PID
				where GeneType_WID=2
				group by t.WID"""
	cursor.execute(sqlcmd)
	res = [[g[0]]+g[1].split('\t') for g in cursor.fetchall()]
	return res

def get_nonB_RR_topologies():
	cursor = setup_cursor()
	sqlcmd = """select gp.source, group_concat(t.GeneID SEPARATOR '\t')
				from Topology t
				inner join (
				            select t.WID
				            from Topology t
				            inner join GenePool gp on t.GeneID = gp.PID
				            inner join GeneType gt on t.GeneType_WID=gt.WID
				            where t.TopologyDefinition_WID=20
				            group by t.WID
				            having group_concat(gt.GeneShortType ORDER BY gp.pFrom DESC SEPARATOR '') = 'RR' ) t2 on t.WID = t2.WID
				inner join GI_to_Ref gr on gr.gid = t.GeneID
				inner join GenePool gp on t.GeneID = gp.PID
				inner join Topology_of_Organism too on too.Topology_WID = t.WID
				inner join Organism o on o.WID = too.Organism_WID
				where GeneType_WID=2 and not o.Name like '%Burkholderia%'
				group by t.WID"""
	cursor.execute(sqlcmd)
	res = [[g[0]]+g[1].split('\t') for g in cursor.fetchall()]
	return res

def get_RXR_topologies():
	cursor = setup_cursor()
	sqlcmd = """select gp.source, group_concat(t.GeneID SEPARATOR '\t')
				from Topology t
				inner join (
				            select t.WID
				            from Topology t
				            inner join GenePool gp on t.GeneID = gp.PID
				            inner join GeneType gt on t.GeneType_WID=gt.WID
				            where t.TopologyDefinition_WID=20
				            group by t.WID
				            having group_concat(gt.GeneShortType ORDER BY gp.pFrom DESC SEPARATOR '') = 'RXR' ) t2 on t.WID = t2.WID
				inner join GI_to_Ref gr on gr.gid = t.GeneID
				inner join GenePool gp on t.GeneID = gp.PID
				where GeneType_WID=2
				group by t.WID"""
	cursor.execute(sqlcmd)
	res = [[g[0]]+g[1].split('\t') for g in cursor.fetchall()]
	return res

def get_nonB_RXR_topologies():
	cursor = setup_cursor()
	sqlcmd = """select gp.source, group_concat(t.GeneID SEPARATOR '\t')
				from Topology t
				inner join (
				            select t.WID
				            from Topology t
				            inner join GenePool gp on t.GeneID = gp.PID
				            inner join GeneType gt on t.GeneType_WID=gt.WID
				            where t.TopologyDefinition_WID=20
				            group by t.WID
				            having group_concat(gt.GeneShortType ORDER BY gp.pFrom DESC SEPARATOR '') = 'RXR' ) t2 on t.WID = t2.WID
				inner join GI_to_Ref gr on gr.gid = t.GeneID
				inner join GenePool gp on t.GeneID = gp.PID
				inner join Topology_of_Organism too on too.Topology_WID = t.WID
				inner join Organism o on o.WID = too.Organism_WID
				where GeneType_WID=2 and not o.Name like '%Burkholderia%'				
				group by t.WID"""
	cursor.execute(sqlcmd)
	res = [[g[0]]+g[1].split('\t') for g in cursor.fetchall()]
	return res

def get_long_short_Rs_of_RR():
	shorts, longs = [], []
	rrs = get_RR_topologies()
	for rr in rrs:
		seqs = get_fasta(rr[1:])
		first, second = seqs
		if len(first.seq)<len(second.seq):
			shorts.append(first)
			longs.append(second)
		else:
			shorts.append(second)
			longs.append(first)
	return shorts, longs

def get_long_short_Rs_of_RXR():
	shorts, longs = [], []
	rrs = get_RXR_topologies()
	for rr in rrs:
		seqs = get_fasta(rr[1:])
		first, second = seqs
		if len(first.seq)<len(second.seq):
			shorts.append(first)
			longs.append(second)
		else:
			shorts.append(second)
			longs.append(first)
	return shorts, longs

def acc2src_org(acc_list, dict_out=True):
	cursor = setup_cursor()
	sqlcmd = """select gr.accession, gp.source, o.Name
				from GI_to_Ref gr
				inner join GenePool gp on gr.gid = gp.PID
				inner join Topology t on t.GeneID = gp.PID
				inner join Topology_of_Organism too on too.Topology_WID = t.WID
				inner join Organism o on o.WID = too.Organism_WID
				where gr.accession in ('%s')"""
	accs_param = "\',\'".join(acc_list)
	
	cursor.execute(sqlcmd%accs_param)
	if dict_out:
		res = {g[0]:g[2] for g in cursor.fetchall()}
	else:
		res = [[g[0],g[2]] for g in cursor.fetchall()]
	return res

def get_non_solo_topologies():
	cursor = setup_cursor()
	sqlcmd="""select t.WID, group_concat(t.GeneID)
			from Topology t
			inner join TopologyDefinition td on t.TopologyDefinition_WID = td.WID
			inner join GenePool gp on t.GeneID = gp.PID
			inner join GI_to_Ref gr on gr.gid = t.GeneID
			inner join Topology_of_Organism too on too.Topology_WID = t.WID
			inner join Organism o on too.Organism_WID = o.WID
			where td.QSSystem_WID=1 and not td.Code like 's%'
			group by t.WID"""
	cursor.execute(sqlcmd)
	topologies = []
	for g in cursor.fetchall():
		top_wid = g[0]
		top_ids = g[1].split(',')
		top_genes = sorted(get_genes(top_ids))
		topologies.append([top_wid, top_genes])
	return topologies

if __name__=='__main__':
	# rrs = get_nonB_RR_topologies()
	# fiveT, threeT = [], []
	# for rr in rrs:
	# 	seqs = get_genes(rr[1:])
	# 	first, second = seqs[0], seqs[1]
	# 	if first.strand=='-' and second.strand=='-':
	# 		fiveT.append(second)
	# 		threeT.append(first)
	# 	elif first.strand=='+' and second.strand=='+':
	# 		fiveT.append(first)
	# 		threeT.append(second)
	# 	else:
	# 		print 'shouldnt be here'
	
	# threeT = get_fasta([s.gid for s in threeT])
	# fiveT = get_fasta([s.gid for s in fiveT])

	# SeqIO.write(fiveT,open('/home/sanjarbek/Projects/AHL/RR_project/prep_files/NonB_5T.fasta','w'),'fasta')
	# SeqIO.write(threeT,open('/home/sanjarbek/Projects/AHL/RR_project/prep_files/NonB_3T.fasta','w'),'fasta')
	
	solors = get_soloRs(genome_type=2)+get_soloRs()
	nonsolors = get_nonsoloRs(genome_type=2)+get_nonsoloRs()
	allrs = solors + nonsolors
	
	sys.exit()

	shorts, longs = get_long_short_Rs_of_RR()
	print len(shorts), len(longs)
	shorts = [s for s in shorts if not 'Burkholderia' in s.id]
	longs = [s for s in longs if not 'Burkholderia' in s.id]
	print len(shorts), len(longs)

	for s in shorts:
		new_id = s.id.split('|')[1]+'_RR_SH_'+"_".join(s.id.split('[')[1][:-1].split())
		s.id = new_id
	for s in longs:
		new_id = s.id.split('|')[1]+'_RR_LO_'+"_".join(s.id.split('[')[1][:-1].split())
		s.id = new_id
	SeqIO.write(shorts+longs,open('/home/sanjarbek/Projects/AHL/RR_project/prep_files/NonB_RR.fasta','w'),'fasta')

	shorts, longs = get_long_short_Rs_of_RXR()
	print len(shorts), len(longs)
	shorts = [s for s in shorts if not 'Burkholderia' in s.id]
	longs = [s for s in longs if not 'Burkholderia' in s.id]
	print len(shorts), len(longs)
	
	for s in shorts:
		new_id = s.id.split('|')[1]+'_RXR_SH_'+"_".join(s.id.split('[')[1][:-1].split())
		s.id = new_id
	for s in longs:
		new_id = s.id.split('|')[1]+'_RXR_LO_'+"_".join(s.id.split('[')[1][:-1].split())
		s.id = new_id
	SeqIO.write(shorts+longs,open('/home/sanjarbek/Projects/AHL/RR_project/prep_files/NonB_RXR.fasta','w'),'fasta')
	
	short_ids =  [seq.id for seq in shorts if 'Burkholderia' in seq.id ]
	short_ids = [s.split('|')[1] for s in short_ids]
	
	long_ids =  [seq.id for seq in longs if 'Burkholderia' in seq.id ]
	long_ids = [s.split('|')[1] for s in long_ids]
	
	# SeqIO.write(all_seqs,open('/home/sanjarbek/Projects/AHL/RR_project/RRs_labelled.fasta','w'),'fasta')
	