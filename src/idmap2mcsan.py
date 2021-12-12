#coding: utf-8
import sys,os
import re
from collections import OrderedDict
def fmt_chr(species, chr_id, sps_list=set([]), chr_list=set([]), gs0=''):
	g,s = species.split('_', 1)[:2]
	g,s = [re.compile(r'[a-z]', re.I).findall(v) for v in [g,s]]#list(g), list(s)
	for a in s:
		gs = g[0]+a
		if not gs in sps_list:
			break
	if gs in sps_list:
		for b in g[1:]:
			gs = g[0]+b
			if not gs in sps_list:
				break
	if gs in sps_list:
		for b in 'abcdefghijklmnopqrstuvwxyz':
			gs = g[0]+b
			if not gs in sps_list:
				break
			
	if gs0:
		gs = gs0
	elif gs in sps_list:
		raise ValueError('%s get no uniq id %s' % (species, gs))
	try: cid = re.compile(r'(\d+)').search(chr_id).groups()[0]
	except AttributeError:
		cid = '1'
#		print >>sys.stderr, chr_id, species
	cid = str(int(cid))
	rid = gs+cid
	if rid in chr_list:
		for i in range(1,10000):
			rid = gs+ str(int(cid)+i)
			if not rid in chr_list:
				break
	return gs, rid
def main(inSd=sys.argv[1], refdir='reference', outGff=sys.stdout):
	idmap_string = '%s/%s/id_mapping.txt'
	d_gene_num = {}
	d_chr_num = {}
	d_chrs = OrderedDict()
	d_sps = OrderedDict()
	sps_list=set([])
	chr_list=set([])
	for line in open(inSd):
		species = line.rstrip().split()[0]
		idmap = idmap_string % (refdir, species,)
		d_length = {}
		if not os.path.exists(idmap):
			with open('prepare.failed', 'a') as f:
				print >> f, idmap
			continue
		for line in open(idmap):
			if line.startswith('#'):
				continue
			temp = line.rstrip().split('\t')
			gene_id,transcript_id,protein_id,chr,start,end,strand = temp[1:8]
			gene_length = int(end)-int(start)+1
			gene_id = '%s|%s' % (species, gene_id)
			if gene_id in d_length: # and gene_length < d_length[gene_id]: 只取第一个，不论长短
				continue	# 一个基因多个转录本，则只取第一个的坐标做基因坐标
			else:
				d_length[gene_id] = gene_length
			if (species, chr) not in d_chrs:
				if species not in d_sps:
					gs, rid = fmt_chr(species, chr, sps_list, chr_list, )
					sps_list.add(gs)
					d_sps[species] = gs
				else:
					gs, rid = fmt_chr(species, chr, sps_list, chr_list, d_sps[species])
				chr_list.add(rid)
				d_chrs[(species, chr)] = rid

			try: d_gene_num[(species, chr)] += 1
			except KeyError: d_gene_num[(species, chr)] = 1

			chr = d_chrs[(species, chr)]
			print >> outGff, '\t'.join([chr, gene_id, start, end, strand])

			try: d_chr_num[species].add(chr)
			except KeyError: d_chr_num[species] = set([chr])

	for (species, chr),chr_id in d_chrs.items():
		print >> sys.stderr, '\t'.join([chr_id,chr,species, str(d_gene_num[(species, chr)])])
	for species, gs in d_sps.items():
		print >> sys.stderr, '\t'.join([species, 'X', gs, str(len(d_chr_num[species])) ])
if __name__ == '__main__':
	inSd=sys.argv[1]
	try: refdir =sys.argv[2]
	except IndexError: refdir='reference'
	main(inSd, refdir)
