import time, resource, platform, sys, argparse, gzip

def parse_seqs(path):
	handle = gzip.open(path) if path.split('.')[-1] == 'gz' else open(path)
	id = next(handle).split()[0][1:]
	seq = ''
	for line in handle:
		if line[0] == '>':
			yield id, seq
			id = line.split()[0][1:]
			seq = ''
		else:
			seq += line.rstrip()
	yield id, seq
	handle.close()

def log_time(start):
	current_time = time.time()
	program_time = round(current_time - start, 2)
	peak_ram = round(max_mem_usage(), 2)


def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	if platform.system() == 'Linux':
		return round((max_mem_self + max_mem_child)/float(1e6), 2)
	else:
		return round((max_mem_self + max_mem_child)/float(1e9), 2)

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="Centroid based sequence clustering"
	)
	parser.add_argument('--fna', type=str, required=True, metavar='PATH',
		help="""Path to nucleotide sequences""")
	parser.add_argument('--ani', type=str, required=True, metavar='PATH',
		help="""Path to tab-delimited file with fields: [qname, tname, num_alns, ani, qcov, tcov]""")
	parser.add_argument('--out', type=str, required=True, metavar='BASENAME',
		help="""Path to output file""")
	parser.add_argument('--exclude', type=str, metavar='PATH',
		help="""Path to list of sequence ids to exclude from clustering""")
	parser.add_argument('--keep', type=str, metavar='PATH',
		help="""Path to list of sequence ids to keep from clustering""")
	parser.add_argument('--min_ani', type=float, metavar='FLOAT', default=95,
		help="""Minimum average nucleotide identity (0...100, default=95)""")
	parser.add_argument('--min_qcov', type=float, metavar='FLOAT', default=10,
		help="""Minimum alignment coverage of longer sequence (0...100, default=10)""")
	parser.add_argument('--min_tcov', type=float, metavar='FLOAT', default=70,
		help="""Minimum alignment coverage of shorter sequence (0...100, default=70)""")
	parser.add_argument('--min_length', type=float, metavar='INT', default=1,
		help="""Minimum sequence length (default=1)""")
	return vars(parser.parse_args())


start = time.time()
args = parse_arguments()
seqs = {}
exclude = set([_.rstrip() for _ in open(args['exclude'])]) if args['exclude'] else None
keep = set([_.rstrip() for _ in open(args['keep'])]) if args['keep'] else None
for index, r in enumerate(parse_seqs(args['fna'])):
	id, seq  = r
	if len(seq) < args['min_length']:
		continue
	elif exclude and id in exclude:
		continue
	elif keep and id not in keep:
		continue
	else:
		seqs[id] = len(seq)
seqs = [x[0] for x in sorted(seqs.items(), key=lambda x: x[1], reverse=True)]
log_time(start)

num_edges = 0
edges = dict([(x,[]) for x in seqs])
handle = gzip.open(args['ani']) if args['ani'].split('.')[-1] == 'gz' else open(args['ani'])
for index, line in enumerate(handle):
	qname, tname, num_alns, ani, qcov, tcov = line.split()
	if qname == tname:
		continue
	elif qname not in edges or tname not in edges:
		continue
	elif float(qcov) < args['min_qcov'] or float(tcov) < args['min_tcov'] or float(ani) < args['min_ani']:
		continue
	edges[qname].append(tname)
	num_edges += 1
handle.close()
log_time(start)

clust_to_seqs = {}
seq_to_clust = {}

for seq_id in seqs:

	if seq_id in seq_to_clust:
		continue
	else:
		clust_to_seqs[seq_id] = [seq_id]
		seq_to_clust[seq_id] = seq_id
		for mem_id in edges[seq_id]:
			if mem_id not in seq_to_clust:
				clust_to_seqs[seq_id].append(mem_id)
				seq_to_clust[mem_id] = seq_id
log_time(start)


out = open(args['out'], 'w')
for seq_id, mem_ids in clust_to_seqs.items():
	out.write(seq_id + '\t' + ','.join(mem_ids)+'\n')
log_time(start)

