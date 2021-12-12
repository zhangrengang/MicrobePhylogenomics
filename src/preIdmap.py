#!/bin/env python
#coding: utf-8
import sys,os
import re
import argparse
import numpy as np
from collections import OrderedDict
from Bio import SeqIO
from Gff import GffGenes, AugustusGtfGenes, GffLines
from RunCmdsMP import logger, run_cmd
from small_tools import mkdirs, is_gz
from small_tools import open_file as open

# jgi
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict
import tarfile
import gzip
import time

GFF_PARSER = {
	'refseq': GffGenes,
	'genbank': GffGenes,
	'gwh': GffGenes,
	'jgi': GffGenes,
	'cndb': GffGenes,
	'coge': GffGenes,
	'hgd': GffGenes,
	'self': GffGenes,
	'standard': GffGenes,
#	'standard': AugustusGtfGenes,
}

GFF_FIELD = {	# 从gff中提取的信息的类型和key
	# pse_id gene_id transcript_id   protein_id
#	'refseq': (('gene', 'Dbxref'), ('gene', 'Name'), ('mRNA', 'transcript_id'), ('CDS', 'protein_id')),
#	'refseq': (('gene', 'locus_tag'), ('gene', 'Name'), ('CDS', 'protein_id'), ('CDS', 'protein_id')),
#	'refseq': (('gene', 'ID'), ('gene', 'Name'), ('CDS', 'protein_id'), ('CDS', 'protein_id')),	# bacteria
	'refseq': (('gene', 'ID'), ('gene', 'locus_tag'), ('CDS', 'protein_id'), ('CDS', 'protein_id')),
#	'genbank': (('gene', 'locus_tag'), ('gene', 'Name'), ('mRNA', 'orig_transcript_id'), ('CDS', 'protein_id')),
#	'genbank': (('gene', 'locus_tag'), ('gene', 'Name'), ('mRNA', 'ID'), ('CDS', 'protein_id')),
#	'genbank': (('gene', 'locus_tag'), ('gene', 'locus_tag'), ('mRNA', 'ID'), ('CDS', 'protein_id')),
#	'genbank': (('gene', 'ID'), ('gene', 'Name'), ('mRNA', 'ID'), ('CDS', 'protein_id')),
#	'genbank': (('gene', 'ID'), ('gene', 'Name'), ('CDS', 'protein_id'), ('CDS', 'protein_id')), # bacteria
	'genbank': (('gene', 'ID'), ('gene', 'locus_tag'), ('CDS', 'protein_id'), ('CDS', 'protein_id')),
	'gwh': (('gene', 'ID'), ('gene', 'Accession'), ('mRNA', 'Accession'), ('CDS', 'Protein_Accession')),
	'jgi': (('gene', 'ID'), ('gene', 'Name'), ('mRNA', 'Name'), ('mRNA', 'Name')),
	'cndb': (('mRNA', 'ID'), ('mRNA', 'ID'), ('mRNA', 'ID'), ('mRNA', 'ID')),
#	'coge': (('gene', 'ID'), ('gene', 'ID'), ('mRNA', 'ID'), ('CDS', 'ID')),
#	'coge': (('gene', 'ID'), ('gene', 'ID'), ('gene', 'ID'), ('gene', 'ID')),
	'coge': (('gene', 'ID'), ('gene', 'ID'), ('CDS', 'ID'), ('CDS', 'ID')),
#	'hgd': (('gene', 'ID'), ('gene', 'Name'), ('mRNA', 'Name'), ('mRNA', 'Name')),
#	'hgd': (('gene', 'ID'), ('gene', 'ID'), ('mRNA', 'ID'), ('mRNA', 'ID')),
	'hgd': (('gene', 'ID'), ('gene', 'Name'), ('mRNA', 'ID'), ('mRNA', 'ID')),	
	'self': (('gene', 'ID'), ('gene', 'ID'), ('mRNA', 'ID'), ('mRNA', 'ID')),
#	'standard': (('gene', 'ID'), ('gene', 'ID'), ('mRNA', 'ID'), ('mRNA', 'ID')),
	'standard': (('mRNA', 'ID'), ('mRNA', 'ID'), ('mRNA', 'ID'), ('mRNA', 'ID')),
#	'standard': (('gene', 'ID'), ('gene', 'ID'), ('gene', 'ID'), ('gene', 'ID')),
#	'standard': (('gene', 'ID'), ('gene', 'ID'), ('mRNA', 'ID'), ('gene', 'ID')),
#	'standard': (('gene', 'Name'), ('gene', 'Name'), ('gene', 'Name'), ('gene', 'Name')),
#	'standard': (('gene', 'ID'), ('gene', 'Name'), ('mRNA', 'ID'), ('mRNA', 'ID')),
#	'standard': (('gene', 'ID'), ('gene', 'Name'), ('mRNA', 'Name'), ('mRNA', 'Name')),
#	'standard': (('gene', 'ID'), ('gene', 'Name'), ('gene', 'Name'), ('mRNA', 'Name')),
#	'standard': (('mRNA', 'Parent'), ('mRNA', 'Parent'), ('mRNA', 'ID'), ('mRNA', 'ID')),
#	'standard': (('mRNA', 'Name'), ('mRNA', 'Name'), ('mRNA', 'ID'), ('mRNA', 'ID')),
#	'standard': (('gene', 'ID'), ('gene', 'ID'), ('transcript', 'ID'), ('transcript', 'ID')),	# augustus
#	'standard': (('transcript', 'Parent'), ('transcript', 'Parent'), ('transcript', 'ID'), ('transcript', 'ID')),  # augustus
#	'standard': (('transcript', 'ID'), ('transcript', 'ID'), ('transcript', 'ID'), ('transcript', 'ID')), 
#	'standard': (('gene', 'ID'), ('gene', 'ID'), ('mRNA', 'ID'), ('polypeptide', 'ID')),
#	'standard': (('gene', 'ID'), ('gene', 'Name'), ('CDS', 'protein_id'), ('CDS', 'protein_id')),
#	'standard': (('gene', 'ID'), ('gene', 'ID'), ('CDS', 'ID'), ('CDS', 'ID')),
}

CDS_FIELD = {	# cds序列对应的id
	'refseq': 'protein_id',
	'genbank': 'protein_id',
	'gwh': 'transcript_id',
	'jgi': 'transcript_id',
	'cndb': 'transcript_id',
	'coge': 'protein_id',
	'hgd': 'transcript_id',
	'self': 'transcript_id',
	'standard': 'transcript_id',
}
PEP_FIELD = {   # pep序列对应的id
	'refseq': 'protein_id',
	'genbank': 'protein_id',
	'gwh': 'protein_id',
	'jgi': 'protein_id',
	'cndb': 'protein_id',
	'coge': 'transcript_id',
	'hgd': 'transcript_id',
	'self': 'protein_id',
	'standard': 'protein_id',
}


def parse_ncbi_cds_seqid(id):
	id = id.split()[0]
	return re.compile(r'\S+_cds_([A-Z]+_?\d+\.?\d*)_\d+').match(id).groups()[0]
#	except AttributeError: return None
def parse_normal_id(id):
	id = id.split()[0]
	return id
def parse_jgi_pep_id(id):
	try: return re.compile(r'transcript=(\S+)').search(id).groups()[0]
	except AttributeError as e: pass
	id = id.split()[0]
	if id.endswith('.p'):
		return id[:-2]
	return id
def parse_jgi_gff_id(id):
	try: return re.compile(r'(\S+)\.v\d+\.?\d*$').match(id).groups()[0]
	except AttributeError as e: return id
def parse_prunus_pep_id(id):
	try: return re.compile(r'^(\S+)_\d+').match(id).groups()[0]
	except AttributeError as e: return id

def parse_coge_cds_pep_id(id):
	return id.split('||')[4]
def parse_hgd_pep_id(id):
    return id.split()[-1]

def get_ncbi_prefix_from_url(url):
	return url.strip('/').split('/')[-1]
def get_gwh_prefix_from_url(url):
	return url.strip('/').split('/')[-1].split('_')[-1]
def get_jgi_prefix_from_url(url):
	return
	
CDS_ID_PARSER = {
	'refseq': parse_ncbi_cds_seqid,
	'genbank': parse_ncbi_cds_seqid,
	'gwh': parse_normal_id,
	'jgi': parse_normal_id,
	'cndb': parse_normal_id,
	'coge': parse_coge_cds_pep_id,
	'hgd': parse_normal_id,
	'self': parse_normal_id,
	'standard': parse_normal_id,
}
PEP_ID_PARSER = {
	'refseq': parse_normal_id,
	'genbank': parse_normal_id,
	'gwh': parse_normal_id,
	'jgi': parse_jgi_pep_id,
	'cndb': parse_normal_id,
	'coge': parse_coge_cds_pep_id,
	'hgd': parse_hgd_pep_id,
	'self': parse_normal_id,
	'standard': parse_normal_id,
#	'standard': parse_prunus_pep_id,
}
GFF_ID_PARSER = {
    'refseq': parse_normal_id,
    'genbank': parse_normal_id,
    'gwh': parse_normal_id,
    'jgi': parse_jgi_gff_id,
    'cndb': parse_normal_id,
    'coge': parse_normal_id,
	'hgd': parse_normal_id,
    'self': parse_normal_id,
    'standard': parse_normal_id,
}

NCBI_SUFFIX = ['_genomic.gff.gz', '_cds_from_genomic.fna.gz', '_protein.faa.gz', '_genomic.fna.gz']
DL_SUFFIX = {	# 下载及匹配文件后缀
	'refseq': NCBI_SUFFIX,
	'genbank': NCBI_SUFFIX,
	'gwh': ['.gff.gz', '.CDS.fasta.gz', '.Protein.faa.gz', '.genome.fasta.gz'],
	'jgi': ['.gene.gff3.gz', '\.cds.fa.gz', '\.protein.fa.gz', '.fa.gz'],
	'cndb': ['.gene.gff.gz', '.gene.cds.fasta.gz', '.gene.pep.fasta.gz', '.fasta.gz'],	# 不通行
	'coge': ['.gid\d+.gff(.gz)?', '\d+-CDS.fasta(.gz)?', '\d+-CDS-prot.fasta(.gz)?', '.faa(.gz)?'],
	'hgd': ['.gff3.gz', '_cds.fa.gz', '_protein.fa.gz', '_genomic.fa.gz'],
	'self': ['gene.gff3', 'gene.cds.fa', 'gene.pep.faa', 'genome.fasta'],
	'standard': [],
}
GFF_MATCH = {	# 匹配GFF文件格式
	'refseq': 'GCF_\d+\S+_genomic.gff.gz',
	'genbank': 'GCA_\d+\S+_genomic.gff.gz',
	'gwh': 'GWH[A-Z]+\d+.*.gff.gz',
	'jgi': '[A-Z][a-z]+\S*_\d+_v\d+\S*.gene.gff3.gz',
	'cndb': '\w+\d+.gene.gff.gz',
	'coge': '\w+\S+cds\S+gid\d+\.gff(.gz)?',
	'hgd': '[A-Z][a-z]+_[a-z]+_\S+.gff3.gz',
	'self': 'final.gene.gff3(.gz)?',
}
NCBI_URL='(ftp|https)://ftp.ncbi.nlm.nih.gov/genomes/all/'
URL_MATCH = {	# 匹配网络链接格式
	'refseq': NCBI_URL + 'GCF',
	'genbank': NCBI_URL + 'GCA',
	'gwh': 'ftp://download.big.ac.cn/gwh',
	'cndb': 'ftp://ftp.cngb.org/pub/CNSA',
	'jgi': '[A-Z][a-z]+',
}
NCBI_PREFIX = get_ncbi_prefix_from_url
DL_PREFIX = {	# 由目录URL获取文件前缀
	'refseq': NCBI_PREFIX,
	'genbank': NCBI_PREFIX,
	'gwh': get_gwh_prefix_from_url,
	'cndb': lambda: '',
	'jgi': get_jgi_prefix_from_url,
}

NCBI_CHROM_MATCH = '(chromosome|linkage group|chr) ([^\s,]+),'
CHROM_MATCH = {	# 匹配染色体
	'refseq':  NCBI_CHROM_MATCH,
	'genbank': NCBI_CHROM_MATCH,
	'gwh': '(Chromosome |chr)(\S+)',
	'jgi': '(chr|LG)(\S+)', 
	'cndb': 'chr(\S+)',
	'coge': '(chr|LG)(\S+)',
	'hgd': '(chr|LG)(\S+)',
	'self': 'chr(\S+)',
	'standard': '(chr|LG|Hic_asm)[^\d]*([^\s,]+)',	# Hic_asm for ALLHIC
}

STD_SOURCE = 'standard'
IDMAP, CDS, PEP, GENOME = 'id_mapping.txt', 'cds.fa', 'pep.faa', 'genome.fasta.gz'



class Download:
	def __init__(self, url, source=None):
		self.url = url
		if source is None:
			source = guess_source_from_url(url)
		self.source = source
		self.prefix = DL_PREFIX[source](url)
		self.suffix = DL_SUFFIX[source]
	def download(self):	# 按链接下载指定文件
		if self.source in {'refseq', 'genbank', 'gwh'}:
			urls = self.format_urls()
			self.wget_urls(urls)
		elif self.source in {'jgi'}:
			DownloadJGI(self.url).run()
		else:
			raise ValueError('unknown source: {}'.format(self.source))
	def format_urls(self):	# 组合目录URL+文件前缀（由URL获取）+文件后缀（固定）
		urls = []
		for suffix in self.suffix:
			url = '{}/{}{}'.format(self.url, self.prefix, suffix)
			urls += [url]
		with open('url', 'a') as f:
			for url in urls:
				print >>f, url
		return urls
	def wget_urls(self, urls, retry=10):
		for url in urls:
			for i in range(retry):
				filename = os.path.basename(url)
				ckpfile = filename + '.ok'
				cmd = 'wget -c {}'.format(url)
				os.system(cmd)
#				if os.path.exists(ckpfile):
#					os.remove(ckpfile)
#					break
				if gzcheck(filename):
					break
				try: os.remove(filename)
				except OSError: pass
			else:
				raise ValueError('Failed to download {}'.format(url))
		
def is_gz(input_file):
    suffix = os.path.splitext(input_file)[-1]
    if suffix  == '.gz':
        return 'gzip'
    elif suffix  == '.bz2':
        return 'bzip2'
    else:
        return False
def gzcheck(gzfile, skip=True):
    gz_ok_file = gzfile+'.ok'
    if skip:
        if os.path.exists(gz_ok_file):
            return True #'%s:    OK (I guess)\n' % (gzfile,)
    cmd = is_gz(gzfile)
    cmd += ' -tv "%s"' % gzfile
    output = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True).communicate()
    result = output[1].rstrip()
    result = '\n'.join([ v for v in result.split('\n') if not v.endswith("extra field of 6 bytes ignored") ])
    ifok = result[-2:]
    if ifok in set(['OK','ok']):
#        try: os.mknod(gz_ok_file)
#        except OSError: pass
        result = True #result + '\n'
    else:
        result = False #'\033[1;31;40m' + result + '\033[0m\n'
    return result
	
class DownloadJGI:
	def __init__(self, release, 	# example: Phaseolus vulgaris UI111 v1.1 -> PulgarisUI111 + v1.1
				load_failed='jgi-query.log', configure=False, xml=None, retry_n=5, filter_files=False, 
				source='jgi'):
		self.release = release
		self.organism_abbreviation, self.version = self.parse_release()
		
		self.suffix = DL_SUFFIX[source]
#		self.prefix = '{}_\d+_{}'.format(self.organism_abbreviation, self.major_version)
#		self.prefix += '.{0,3}'
		self.prefix = '{}_\d+_\S*?[Vv]?.*\d'.format(self.organism_abbreviation)
		self.load_failed = load_failed
		self.configure = configure
		self.xml = xml
		self.retry_n = retry_n
		self.filter_files = filter_files
	def parse_release(self):
		words = self.release.split()
		org_abbr = ''.join([words[0][0]] + words[1:-1])
		version = words[-1]
		#major, minor = version.split('.')
		return org_abbr, version
	def validate_release(self):
		if not re.compile(r'[A-Z]').match(self.organism_abbreviation):
			raise ValueError('`{}` not starts with uppercase'.format(self.release))
		if not re.compile(r'v\d').match(self.version):
			raise ValueError('`{}` not has verion'.format(self.release))
	def format_regex(self):
		regex = []
		for suffix in self.suffix:
			regex += [self.prefix + suffix]
		return regex
	def run(self):
		logger.info('JGI query: {} {}'.format(self.organism_abbreviation, self.version))
		jgi_query_main(self)






def guess_source_from_url(url):
	for source, pattern in URL_MATCH.items():
		if re.compile(pattern).match(url):
			logger.info('guessed source: {}'.format(source))
			return source
	return STD_SOURCE

def guess_files(indir='.'):
	files = [fl for fl in os.listdir(indir) if os.path.isfile(fl)]
	# guess from gff
	gff_source = 0
	for source, pattern in GFF_MATCH.items():
		pattern = re.compile(pattern)
		for fl in files:
			if pattern.match(fl):
				logger.info('guessed source: {} from {}'.format(source, fl))
				gff_source = source
	# guess from all files
	for source, suffixes in DL_SUFFIX.items():
		print source
		tfiles = []
		used_files = set([])
#		got = False
		for suffix in suffixes:
			pattern = re.compile(suffix+'$')
			for fl in files:
#				if fl.endswith(suffix):
#					print suffix, fl
				if fl in used_files:
					continue
#				if fl.endswith(suffix):
				if pattern.search(fl):
					print suffix, fl
#					got = True
					used_files.add(fl)
					break
			else:
				fl = None
			tfiles += [fl]
		got = [f for f in tfiles if f is not None]
		if len(got) >= 3:
			logger.info('guessed input files: {}'.format(tfiles))
			source = gff_source if gff_source else source	# to fix bug
			logger.info('guessed source: {}'.format(source))
			return tfiles, source

class PreIdMap:
	def __init__(self, gff=None, cds=None, pep=None, genome=None, source=None):
		if gff is None:
			(gff, cds, pep, genome), source = guess_files('./')
		if source is None:
			source = guess_source_from_gff(gff)
		self.gff, self.cds, self.pep, self.source = gff, cds, pep, source
		self.genome = genome
		self.gff_parser = GFF_PARSER[source]
		self.gff_field = GFF_FIELD[source]
		logger.info('GFF_FIELD: {}'.format(self.gff_field))
		self.cds_id_parser = CDS_ID_PARSER[source]
		self.pep_id_parser = PEP_ID_PARSER[source]
		self.gff_id_parser = GFF_ID_PARSER[source]
		self.cds_field = CDS_FIELD[source]
		self.pep_field = PEP_FIELD[source]
		self.chrom_match = CHROM_MATCH[source]
		self.pre_genome()

	def pipe(self):
		# gff -> idmap
		logger.info('preparing {}'.format(IDMAP))
		lines = self.gff_to_idmap()
		idmap = IdMap(lines=lines)
		d_chrom, d_length = self.get_chrom()	# 无基因组时为空字典；d_chrom只匹配染色体
		nchr = len(d_chrom)
		stats = GenomeStats(d_length, d_chrom)
		ngene, nrna = idmap.stats()
		#logger.info('retrieved {} chromosomes'.format(len(d_chrom)))
		print >>sys.stderr, 'size\tN50\t#gene\t#mRNA\t#chrom\t%chrom'
		print >>sys.stderr, '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(
				stats.genome_size, stats.n50, ngene, nrna, nchr, stats.chrom_rate)
#		with open(IDMAP, 'w') as fout:
#			idmap.write(fout, d_chrom)

		# pep
		logger.info('preparing {}'.format(PEP))
		gp_map = idmap.one2many_map('gene_id', self.pep_field)	# gene -> pep (s)
		p_seqs = self.seq2dict(self.pep, self.pep_id_parser)
		with open(PEP, 'w') as fout:
			gp_map1_1 = self.get_gene_seq(p_seqs, gp_map, fout)	# gene -> pep

		with open(IDMAP, 'w') as fout:
			idmap.write(fout, d_chrom, gp_map=gp_map1_1)

		# cds
		logger.info('preparing {}'.format(CDS))
		pt_map = idmap.one2one_map(self.pep_field, self.cds_field)	# pep -> cds
		gt_map = self.convert_map(gp_map1_1, pt_map)				# gene -> cds
		c_seqs = self.seq2dict(self.cds, self.cds_id_parser)
		with open(CDS, 'w') as fout:
			gt_map1_1 = self.get_gene_seq(c_seqs, gt_map, fout, log=False)

	def gff_to_idmap(self):
		self.lines = []
		i = 0
#		print >>sys.stdout, self.gff_parser, self.gff
#		for line in GffLines(self.gff):
#			print >>sys.stdout, len(line), line
#			break
		for gene in self.gff_parser(self.gff):
			i += 1
		#	print >>sys.stdout, gene.id #, vars(gene)
		#	print >> sys.stderr, gene.id, len(gene.lines), 'lines'
		#	print >> sys.stderr, gene.id, gene.edges()
			if gene.type in {'pseudogene'}:
				continue
#			print >>sys.stdout, gene.id
			gene_line = [line for line in gene.lines if line.type == 'gene']
#			if not gene_line:  # 'mRNA' only
#				rnas = [gene]
#			else:
#				rnas = gene.rnas
			rnas = gene.rnas
			for rna in rnas:
				try: rna_type = rna.type
				except KeyError: continue
				if not rna.type in {'mRNA', 'transcript', 'CDS', 'gene'}:
					continue
				lines = gene_line + rna.lines
		#		if i < 2:
		#			for line in  lines:
		#				print line.ID, line.type, line
				try: values = self.get_fields(lines)
				except UnboundLocalError: continue
				except Exception: continue
				chrom, start, end, strand = rna.chrom, rna.start, rna.end, rna.strand
				description = rna.gene.attributes.get('product')
				line = values + [chrom, start, end, strand, description]
				line = IdMapLine(line)
				self.lines += [line]
		print >> sys.stderr, '{} records in {}'.format(i, self.gff)
		return self.lines
	def seq2dict(self, seqfile, id_parser):
		d_seqs = {}
		for rc in SeqIO.parse(open(seqfile), 'fasta'):
			try: id = id_parser(rc.description)
			except AttributeError as e:
				print >>sys.stderr, '[WARN] cannot be parsed:', rc.description
			#	raise AttributeError(e)
				continue
			rc.id = id
			d_seqs[id] = rc
		return d_seqs
	def get_gene_seq(self, d_seqs, id_map, fout, log=True):
		# get longest pep for each gene
		longest_ids = OrderedDict()	  # g -> p
		map_ids = []
		nofound = 0
		for gene, rnas in id_map.items():
			if isinstance(rnas, str):
				rnas = [rnas]
			records = []
			for rna in rnas:
				try: records += [d_seqs[rna]]
				except KeyError as e:
					if nofound < 10:
						print >>sys.stderr, rna, 'not found in seqfile'
					elif nofound < 11:
						print >>sys.stderr, '...'
					nofound += 1
			if not records:
				continue
			map_ids += [rc.id for rc in records]
			longest = max(records, key=lambda x: len(x.seq))
			longest_ids[gene] = longest.id
			seq = self.clip_abnormal_aa(longest)
			print >> fout, '>{}\n{}'.format(gene, seq)
		all_sids = [rc.id for rc in d_seqs.values()]
		diff_ids = set(all_sids) - set(map_ids)
		if diff_ids and log:
			print >>sys.stderr, '{} id in seqfile but not in mapfile: {} ...'.format(
					len(diff_ids), sorted(diff_ids)[:10])
		if nofound:
			print >>sys.stderr, '{} sequences not found in seqfile'.format(nofound)
		return longest_ids
	def convert_map(self, map_t, map_q):
		map_n = OrderedDict()
		for k,v in map_t.items():
			map_n[k] = map_q[v]
		return map_n
	def clip_abnormal_aa(self, aaseq):	# for orthofinder/ diamond
		seq = str(aaseq.seq)
		pattern = r'[^ACDEFGHIKLMNPQRSTVWXY\*]'
		abnormal = re.compile(pattern, re.I).findall(seq)
		if abnormal:
			print >>sys.stderr, 'warning: {} contains abnormal string {}'.format(aaseq.id, abnormal)
			seq = re.compile(pattern, re.I).sub('*', seq)
		return seq

	def get_fields(self, lines):
		values = []
		for type, key in self.gff_field:
			#print type, key
			for line in lines:
				#print line
				if line.type == type:
					try: value = line.attributes[key]
					except KeyError:
						print >> sys.stderr, 'Cannot get key `{}` from {}'.format(key, str(line))
						value = None
#						raise KeyError('Cannot get key `{}` from {}'.format(key, str(line)))
					break
#			else:
#				print >> sys.stderr, type, key, 'is not got from', line
			try: values += [value]
			except UnboundLocalError as e:
				print >> sys.stderr, type, key, 'is not got from', line
				raise UnboundLocalError(e) 
			del value
#		print >> sys.stderr, values
		try: values = map(self.gff_id_parser, values)
		except Exception as e :
			print >> sys.stderr, values
			raise Exception(e)
		return values
	def get_chrom(self):
		d_chrom = {}
		d_length = {}
		if not self.chrom_match or not self.genome:
			return d_chrom, d_length
		for rc in SeqIO.parse(open(self.genome), 'fasta'):
			d_length[rc.id] = len(rc.seq)
			try: 
				chrom = re.compile(self.chrom_match, re.I).search(rc.description).groups()[-1]
				d_chrom[rc.id] = chrom
			except AttributeError: continue
		return d_chrom, d_length
	def pre_genome(self):
		if self.genome is None:
			return
		if is_gz(self.genome):
			cmd = 'ln {} {}'.format(self.genome, GENOME)
		else:
			cmd = 'pigz {raw} && ln {raw}.gz {new}'.format(raw=self.genome, new=GENOME)
		run_cmd(cmd)

class GenomeStats:
	def __init__(self, d_length, d_chrom=None):
		self.d_length = d_length
		self.d_chrom = d_chrom
		self.size = sum(self.d_length.values())
		self.lengths = self.d_length.values()
	@property
	def count(self):
		return len(self.d_length)
	@property
	def mean(self):
		return np.mean(self.lengths)
	@property
	def chrom_rate(self):
		try: percent = 1e2 * sum([self.d_length[chrom] for chrom in self.d_chrom]) / self.size
		except ZeroDivisionError: percent = 0
		return '{:.1f}%'.format(percent)
	@property
	def n50(self, cutoff=50):
		accum = 0
		for length in sorted(self.d_length.values(), reverse=1):
			accum += length
			if 1e2*accum / self.size >= cutoff:
				return self.format_length(length)
	@property
	def genome_size(self):
		if self.size > 1e9:
			return '{:.1f}G'.format(self.size/1e9)
		if self.size > 1e6:
			return '{:.0f}M'.format(self.size/1e6)
		return self.size
	def format_length(self, length):
		if length > 1e6:
			return '{:.1f}M'.format(length/1e6)
		if length > 1e3:
			return '{:.0f}K'.format(length/1e3)
		return length

class IdMapLine:
	def __init__(self, line=None):
		if isinstance(line, str):
			self.line = line.strip().split('\t')
		else:
			self.line = line
		self.title = ['pse_id', 'gene_id', 'transcript_id', 'protein_id', 'chrom', 'start', 'end', 'strand', 'description']
		self.type = [str, str, str, str, str, int, int, str, str]
		self._set_attr()
	def _set_attr(self):
		for key, value, type in zip(self.title, self.line, self.type):
			setattr(self, key, type(value))
	def write(self, fout):
		line = [getattr(self, key) for key in self.title]
		print >> fout, '\t'.join(map(str, line))
	def write_title(self, fout):
		print >> fout, '#' + '\t'.join(map(str, self.title))
class IdMap:
	def __init__(self, idmap=None, lines=None):
		if lines is not None:
			self.lines = lines
		else:
			self.lines = self._parse(idmap)
	def __iter__(self):
		return iter(self.lines)
	def _parse(self, idmap):
		for line in open(idmap):
			if line.startswith('#'):
				continue
			yield IdMapLine(line)
	@property
	def gp_map(self):	# gene_id -> protein_id
		return self.one2many_map('gene_id', 'protein_id')
	@property
	def pt_map(self):	# protein_id -> transcript_id
		return self.one2one_map('protein_id', 'transcript_id')
	def one2one_map(self, kid_field, vid_field):
		d_map = OrderedDict()
		for line in self.lines:
			kid = getattr(line, kid_field)
			vid = getattr(line, vid_field)
			d_map[kid] = vid
		return d_map
	def one2many_map(self, kid_field, vid_field):
		d_map = OrderedDict()
		for line in self.lines:
			kid = getattr(line, kid_field)
			vid = getattr(line, vid_field)
			try: d_map[kid] += [vid]
			except KeyError: d_map[kid] = [vid]
		return d_map
	def write(self, fout, d_chrom=None, gp_map=None):	# gp_map: gene_id -> longest protein id
		chroms = self.chroms
		d_chrn = {}
		for chrom in chroms:
			try: num = int(re.compile(r'(\d+)').search(d_chrom[chrom]).groups()[0])
			except: num = 1000
			d_chrn[chrom] = num
		self.lines = sorted(self.lines, key=lambda x: (
				d_chrn[x.chrom], chroms.index(x.chrom), x.start, x.protein_id))
		for line in self.lines:
			try:
				if gp_map is not None and line.protein_id == gp_map[line.gene_id]:
					line.description = 'primary'
			except KeyError: pass
			if d_chrom is not None and line.chrom in d_chrom:
				line.chrom = '{}|{}'.format(d_chrom[line.chrom], line.chrom)
			line.write(fout)
	@property
	def chroms(self):
		chrs = []
		for line in self.lines:
			if line.chrom in set(chrs):
				continue
			chrs += [line.chrom]
		return chrs
	def stats(self):
		genes, rnas = set([]), set([])
		for line in self.lines:
			genes.add(line.gene_id)
			rnas.add(line.protein_id)
		return len(genes), len(rnas)
def guess_source_from_gff(gff):
	for source, pattern in GFF_MATCH.items():
		if re.compile(pattern).match(gff):
			logger.info('guessed source: {}'.format(source))
			return source
	source = STD_SOURCE
	logger.info('guessed source: {}'.format(source))
	return source

def main():
	logger.info(' '.join(sys.argv))
	try:
		args = makeArgparse()
		url = args.download
	except: url = None
	if url is not None:
		Download(url=url).download()
		gff, cds, pep, genome = [None]*4
	else:
		logger.info('input files: {}'.format(sys.argv[1:]))
		try:
			if len(sys.argv[1:]) >=4:
				gff, pep, cds, genome = sys.argv[1:5]
			elif len(sys.argv[1:]) >=3:
				gff, pep, cds =  sys.argv[1:4]
				genome = None
			else:
				raise ValueError('no enough inputs')
		except ValueError as e:
			logger.warn('{}. To guess input files'.format(e))
			gff, cds, pep, genome = [None]*4
	PreIdMap(gff, cds, pep, genome).pipe()
#	with open('standard.sh', 'a') as f:
#			print >> f, ' '.join(sys.argv)
def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("-dl", action="store",type=str, dest='download',default=None,
					help="download URL")
	args = parser.parse_args()
	return args

# jgi-query
def jgi_query_main(self):
	'''revised from https://github.com/glarue/jgi-query'''
	global LOCAL_XML
	LOCAL_XML = False
	self.regex = self.format_regex()
	DIRECT_REGEX = self.regex
	GET_ALL = False
	if GET_ALL or DIRECT_REGEX:
		INTERACTIVE = False
	else:
		INTERACTIVE = True
	# Get script location info
	SCRIPT_HOME = '{}/.config/jgi'.format(os.environ['HOME'])
	mkdirs(SCRIPT_HOME)
	# Config should be in same directory as script
	CONFIG_FILENAME = "jgi-query.config"
	CONFIG_FILEPATH = SCRIPT_HOME + "/{}".format(CONFIG_FILENAME)
	# Categories to store in default config file
	DEFAULT_CATEGORIES = ['ESTs',
				  'EST Clusters',
				  'Assembled scaffolds (unmasked)',
				  'Assembled scaffolds (masked)',
				  'Transcripts',
				  'Genes',
				  'CDS',
				  'Proteins',
				  'Additional Files']
	# Does config file exist?	# 读取或记录config，含账户密码
	if os.path.isfile(CONFIG_FILEPATH):  # use config file
		config_info = read_config(CONFIG_FILEPATH)
	else:  # no config present or configure flag used; run config dialog
		config_info = get_user_info()
		config_info["categories"] = DEFAULT_CATEGORIES
		make_config(CONFIG_FILEPATH, config_info)
		
	# Get user information for sign-on
	USER = config_info["user"]
	PASSWORD = config_info["password"]
	# New syntax	# 生成cookie
	global LOGIN_STRING
	LOGIN_STRING = (
		# "curl 'https://signon-old.jgi.doe.gov/signon/create' "
		"curl 'https://signon.jgi.doe.gov/signon/create' "
		"--data-urlencode 'login={}' "
		"--data-urlencode 'password={}' "
		"-s "  # suppress status output
		"-c cookies > /dev/null"
		.format(USER, PASSWORD)
		)
	if self.load_failed and os.path.isfile(self.load_failed):
		logfile = self.load_failed
		print("Reading URLs from \'{}\'".format(logfile))
		retry_from_failed(LOGIN_STRING, logfile)
		clean_exit("All files in log attempted.")
	org_input = self.organism_abbreviation
	if not org_input:
		if self.configure:
			sys.exit("Configuration complete. Script may now be used to query JGI. "
					 "Exiting now.")
		elif self.xml and self.xml != 1:
			# Use org_input because is already checked further down
			# and avoids re-writing this whole block
			org_input = get_org_name(self.xml)
			if not org_input:
				sys.exit("No organism specified. Exiting now.")
		else:
			sys.exit("No organism specified. Exiting now.")
	org_regex = re.compile(r'\.jgi.+\.(?:gov|org).+\/(.+)\/(?!\/)')
	try:  # see if it's in address form
		organism = org_regex.search(org_input).group(1)
	except AttributeError:  # not in address form, assume string is organism name
		organism = org_input
	
	# URL where remote XML file should be, if it exists
	org_url = ("https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?"
			   "organism={}".format(organism))
			   
	# Get xml index of files, using existing local file or curl API
	if self.xml:
		LOCAL_XML = True  # global referenced by clean_exit()
		xml_arg = self.xml
		if xml_arg == 1:  # --xml flag used without argument
			xml_index_filename = "{}_jgi_index.xml".format(organism)
		else:
			xml_index_filename = xml_arg
		print(
			'Retrieving information from JGI for query '
			'\'{}\' using local file \'{}\'\n'.format(organism, xml_index_filename))
	else:  # fetch XML file from JGI
		xml_index_filename = "{}_jgi_index.xml".format(organism)

		# Old syntax
		# xml_address = ("curl {} -b cookies -c cookies > {}"
		#                .format(org_url, xml_index_filename))

		# New syntax	# 生成xml
		xml_address = ("curl '{}' -L -b cookies > {}"
					   .format(org_url, xml_index_filename))
		try:  # fails if unable to contact server
			subprocess.check_output(LOGIN_STRING, shell=True)
		except subprocess.CalledProcessError as error:
			clean_exit("Couldn't connect with server. Please check Internet "
					  "connection and retry.")
		print(
			'Retrieving information from JGI for query \'{}\' using command '
			'\'{}\'\n'.format(organism, xml_address))
		subprocess.call(xml_address, shell=True)	# call for python2, run for python3
		print()  # padding
	# Parse xml file for content to download
	xml_root = None
	if os.path.getsize(xml_index_filename) == 0:  # happens if user and/or pw wrong
		clean_exit("Invalid username/password combination (or other issue).\n"
				  "Restart script with flag '-c' to reconfigure credentials.")
	try:
		xml_in = ET.ElementTree(file=xml_index_filename)
		xml_root = xml_in.getroot()
	except ET.ParseError:  # organism not found/xml file contains errors
		clean_exit("Cannot parse XML file or no organism match found.\n"
				  "Ensure remote file exists and has content at the "
				  "following address:\n{}".format(org_url))


	# Get categories from config (including possible user additions)
	# Will only be used if --filter_files flag
	DESIRED_CATEGORIES = config_info["categories"]


	# Choose between different XML parsers
	# if args.filter_files, user wants only those files in <desired_categories>
	file_list = get_file_list(xml_index_filename, filter_categories=self.filter_files)
	file_list = {self.version: file_list[self.version]}
	print 'file_list', file_list
	# Check if file has any categories of interest
	if not any(v["results"] for v in list(file_list.values())):
		print(("ERROR: no results found for '{}' in any of the following "
			   "categories:\n---\n{}\n---"
			   .format(organism, "\n".join(DESIRED_CATEGORIES))))
		clean_exit()


	# Decision tree depending on if non-interactive options given
	regex_filter = None
	user_choice = None
	display_info = True
	#print DIRECT_REGEX
	if GET_ALL:
		user_choice = 'a'
		display_info = False
	elif DIRECT_REGEX:
		user_choice = 'r'
		regex_filters = [re.compile(regex) for regex in DIRECT_REGEX]
		display_info = False

	url_dict = print_data(file_list, organism, display=display_info)

	#print url_dict
	# if not user_choice:
		# # Ask user which files to download from xml
		# user_choice = get_user_choice()
		# if user_choice == 'r':
			# regex_filter = get_regex()

	urls_to_get = set()
	# special case for downloading all available files
	# or filtering with a regular expression
	if user_choice in ('a', 'r'):
		for k, v in sorted(url_dict.items()):	# k: folder, v: NO.
			for u in v.values():
				fn = re.search('.+/([^\/]+$)', u).group(1)
				print fn
				for regex_filter in regex_filters:
					match = regex_filter.search(fn)
					if not match:
						continue
					urls_to_get.add(u)
	else:
		# Retrieve user-selected file urls from dict
		ids_dict = parse_selection(user_choice)
		for k, v in sorted(ids_dict.items()):
			for i in v:
				urls_to_get.add(url_dict[k][i])

	logger.info('urls: {}'.format(urls_to_get))

	# Calculate and display total size of selected data
	urls_to_get = sorted(urls_to_get)
	filenames = [u.split('/')[-1] for u in urls_to_get]
	file_sizes = get_sizes(file_list, sizes_by_url={})
	total_size = sum([file_sizes[url] for url in urls_to_get])
	size_string = byte_convert(total_size)
	num_files = len(urls_to_get)
	print(("Total download size for {} files: {}".format(num_files, size_string)))
	if INTERACTIVE:
		download = input("Continue? (y/n/[p]review files): ").lower()
		if download == "p":
			while download == "p":
				print('\n'.join(filenames))
				download = input("Continue with download? (y/n/[p]review files): ").lower()
		if download != "y":
			clean_exit("ABORTING DOWNLOAD")

	downloaded_files, failed_urls = download_list(
		urls_to_get, retries=self.retry_n)

	print("Finished downloading {} files.".format(len(downloaded_files)))

	if failed_urls and INTERACTIVE:
		n_broken = len(failed_urls)
		retry_broken = input(
			"{} files failed to download; retry them? (y/n): ".format(n_broken))
		if retry_broken.lower() in ('yes', 'y'):
			downloaded_files, failed_urls = download_list(
				failed_urls, retries=1)

	if failed_urls:
		log_failed(organism, failed_urls)

	# Kindly offer to unpack files, if files remain after error check
	if downloaded_files and INTERACTIVE:
		decompress = input(("Decompress all downloaded files? "
							"(y/n/k=decompress and keep original): "))
		if decompress != "n":
			if decompress == "k":
				keep_original = True
			else:
				keep_original = False
			decompress_files(downloaded_files, keep_original)
			print('Finished decompressing all files.')

	# Clean up and exit
	# "cookies" file is always created
	if INTERACTIVE:
		keep_temp = input("Keep temporary files ('{}' and 'cookies')? (y/n): "
						.format(xml_index_filename))
		if keep_temp.lower() not in "y, yes":
			clean_exit()
		else:
			print("Leaving temporary files intact and exiting.")
	else:
		clean_exit()
		
# JGI-query
# FUNCTIONS

def deindent(string):
    """
    Print left-justified triple-quoted text blocks

    """
    print(string)


def check_config(d, config_name):
    """
    Check filesystem for existence of configuration
    file, and return the full path of config file
    if found.

    """
    files = os.listdir(d)
    if config_name in files:
        config_path = d + "/{}".format(config_name)
        return config_path
    else:
        return None


def get_user_info():
    """
    Dialog with user to gather user information for
    use with the curl query. Returns a dict.

    """
    blurb = """
    === USER SETUP ===

    JGI access configuration:

    Before continuing, you will need to provide your JGI login credentials.
    These are required by JGI's curl api, and will be stored in a config
    file for future use (unless you choose to delete them).

    If you need to sign up for a JGI account, use the registration link at
    https://signon.jgi-psf.org/signon

    === CREDENTIALS ===
    """
    deindent(blurb)
    user_query = "JGI account username/email (or 'q' to quit): "
    pw_query = "JGI account password (or 'q' to quit): "
    user = raw_input(user_query)	# raw_input for python2, input for python3
    if user == "q":
        sys.exit("Exiting now.")
    pw = raw_input(pw_query)
    if pw == "q":
        sys.exit("Exiting now.")
    input_blurb = ("Proceed with USER='{}', PASSWORD='{}' to configure "
                   "script?\n([y]es, [n]o, [r]estart): ".format(user, pw))
    user_info = {"user": user, "password": pw}
    # while True:  # catch invalid responses
        # choice = raw_input(input_blurb)
        # if choice.lower() == "y":
            # return user_info
        # elif choice.lower() == "n":
            # sys.exit("Exiting now.")
        # elif choice.lower() == "r":
            # user_info = get_user_info()


def make_config(config_path, config_info):
    """
    Creates a config file <config_path> using
    credentials from dict <config_info>.

    """
    u = config_info["user"]
    p = config_info["password"]
    c = config_info["categories"]
    c = ",".join(c)
    header = ("# jgi-query.py user configuration information {}\n"
              .format("#" * 34))
    info = "user={}\npassword={}\ncategories={}".format(u, p, c)
    with open(config_path, 'w') as config:
        config.write(header)
        config.write(info)


def read_config(config):
    """
    Reads "user", "password" and "categories" entries
    from config file.

    """
    user, pw, categories = None, None, None
    with open(config) as c:
        for line in c:
            line = line.strip()
            if line.startswith("user"):
                user = line.split("=")[1]
            if line.startswith("password"):
                pw = line.split("=")[1]
            if line.startswith("categories"):
                cats = line.strip().split("=")[1]
                categories = [e.strip() for e in cats.split(",")]
    if not (user and pw):
        sys.exit("ERROR: Config file present ({}), but user and/or "
                 "password not found.".format(config))
    config_info = {"user": user, "password": pw, "categories": categories}
    return config_info


# /CONFIG

def xml_hunt(xml_file):
    """
    Gets list of all XML entries with "filename" attribute,
    and returns a dictionary of the file attributes keyed
    by a ":"-joined string of parent names.

    """
    root = ET.iterparse(xml_file, events=("start", "end"))
    parents = []
    matches = {}
    for event, element in root:
        if element.tag not in ["folder", "file"]:  # skip topmost categories
            continue
        if element.tag == "folder":
            if event == "start":  # add to parents
                parents.append(element.attrib["name"])
            elif event == "end":  # strip from parents
                del parents[-1]
            continue
        if event == "start" and element.tag == "file":
            parent_string = ":".join(parents)
            try:
                matches[parent_string].append(element.attrib)
            except KeyError:
                matches[parent_string] = [element.attrib]
    return matches


def format_found(d, filter_found=False):
    """
    Reformats the output from xml_hunt()

    """
    output = {}
    for p, c in sorted(d.items()):
        layers = [e for e in p.split(":") if e]
        if filter_found:
            if not any(cat in layers for cat in DESIRED_CATEGORIES):
                continue
        if len(layers) == 1:
            top = parent = layers[0]
        else:
            top = layers[-2]  # either -2 or -1 works well, != parent
            parent = layers[-1]  # either -2 or -1 works well, != top
        if top not in output:
            output[top] = defaultdict(dict)
        if parent not in output[top]:
            output[top][parent] = c
        else:
            output[top][parent].extend(c)
    return output


def get_file_list(xml_file, filter_categories=False):
    """
    Moves through the xml document <xml_file> and returns information
    about matches to elements in <DESIRED_CATEGORIES> if
    <filter_categories> is True, or all files otherwise

    """
    descriptors = {}
    display_cats = ['filename', 'url', 'size',
                    'label', 'sizeInBytes', 'timestamp']
    found = xml_hunt(xml_file)
    found = format_found(found, filter_categories)
    if not list(found.values()):
        return None
    category_id = 0
    for category, sub_cat in sorted(found.items()):
        c = category
        if c not in descriptors:
            category_id += 1
            descriptors[c] = defaultdict(dict)
            descriptors[c]["catID"] = category_id
        uid = 1
        for parent, children in sorted(sub_cat.items()):
            descriptors[c]["results"][parent] = defaultdict(dict)
            results = descriptors[c]["results"][parent]
            unique_children = uniqueify(children)
            for child in sorted(unique_children, key=lambda x: x['filename']):
                try:
                    results[uid]
                except KeyError:
                    results[uid] = {}
                for dc in display_cats:
                    try:
                        results[uid][dc] = child[dc]
                    except KeyError:
                        continue
                uid += 1

    return descriptors


def uniqueify(children):
    """
    Takes a list of child XML elements (dicts of attribs) as 
    returns a filtered list of only unique filenames for a given 
    month/year timestamp (e.g. duplicates are allowed if month/year 
    is different).
    
    """
    unique = {}
    for child in children:
        try:
            fn = child['filename']
            date = fmt_timestamp(child['timestamp'])
            date_string = (date.tm_mon, date.tm_year)
            uid = (fn, date_string)
        except KeyError:
            continue
        if fn not in unique:
            unique[uid] = child
        else:
            existing = unique[uid].get('fileType', None)
            if existing == 'Unknown':
                existing = None
            current = child.get('fileType', None)
            if current == 'Unknown':
                current = None
            if current is not None and existing is None:
                unique[uid] = child
        
    return unique.values()


def get_sizes(d, sizes_by_url=None):
    """
    Builds a dictionary of url:sizes from
    output of get_file_list()

    """
    for k, v in d.items():
        if isinstance(v, dict):
            if 'url' in v:
                address = v['url']
                size = int(v['sizeInBytes'])
                sizes_by_url[address] = size
            else:
                get_sizes(v, sizes_by_url)
    return sizes_by_url


def clean_exit(exit_message=None, remove_temp=True):
    """
    Perform a sys.exit() while removing temporary files and
    informing the user.

    """
    to_remove = ["cookies"]
    # don't delete xml file if supplied by user
    if not LOCAL_XML and remove_temp is True:
        try:
            to_remove.append(xml_index_filename)
        except NameError:
            pass
    for f in to_remove:
        try:
            os.remove(f)
        except OSError:
            continue
    if exit_message:
        print_message = "{}\n".format(exit_message)
    else:
        print_message = ""
    print >> sys.stderr, "{}Removing temp files and exiting".format(print_message)
#    sys.exit("{}Removing temp files and exiting".format(print_message))


def extract_file(file_path, keep_compressed=False):
    """
    Native Python file decompression for tar.gz and .gz files.

    TODO: implement .zip decompression

    """
    tar_pattern = 'tar.gz$'  # matches tar.gz
    gz_pattern = '(?<!tar)\.gz$'  # excludes tar.gz
    endings_map = {"tar": (tarfile, "r:gz", ".tar.gz"),
                   "gz": (gzip, "rb", ".gz")
                   }
    relative_name = os.path.basename(file_path)
    if re.search(tar_pattern, file_path):
        opener, mode, ext = endings_map["tar"]
        with opener.open(file_path) as f:
            file_count = len(f.getmembers())
            if file_count > 1:  # make sub-directory to unpack into
                dir_name = relative_name.rstrip(ext)
                try:
                    os.mkdir(dir_name)
                except FileExistsError:
                    pass
                destination = dir_name
            else:  # single file, extract into working directory
                destination = "."
            f.extractall(destination)
    elif re.search(gz_pattern, file_path):
        opener, mode, ext = endings_map["gz"]
        # out_name = file_path.rstrip(ext)
        out_name = relative_name.rstrip(ext)
        with opener.open(file_path) as f, open(out_name, "wb") as out:
            for l in f:
                out.write(l)
    else:
        print("Skipped decompression for '{}'"
              .format(file_path))
        return
    if not keep_compressed:
        os.remove(file_path)


def decompress_files(local_file_list, keep_original=False):
    """
    Decompresses list of files, and deletes compressed
    copies unless <keep_original> is True.

    """
    for f in local_file_list:
        extract_file(f, keep_original)


def fmt_timestamp(time_string):
    """
    Parses the timestamp string from an XML document
    of the form "Thu Feb 27 16:38:54 PST 2014"
    and returns a string of the form "2014".

    """
    # Remove platform-dependent timezone substring
    # of the general form "xxT"
    tz_pattern = re.compile("\s[A-Z]{3}\s")
    time_string = tz_pattern.sub(" ", time_string)

    # Get the desired time info
    time_info = time.strptime(time_string, "%a %b %d %H:%M:%S %Y")
    # year = str(time_info.tm_year)
    return time_info


def print_data(data, org_name, display=True):
    """
    Prints info from dictionary data in a specific format.
    Also returns a dict with url information for every file
    in desired categories.

    """
    print("\nQUERY RESULTS FOR '{}'\n".format(org_name))
    dict_to_get = {}
    for query_cat, v in sorted(iter(data.items()),
                               key=lambda k_v: k_v[1]["catID"]):
        print_list = []
        if not v["results"]:
            continue
        catID = v["catID"]
        dict_to_get[catID] = {}
        print_list.append(" {}: {} ".format(catID, query_cat).center(80, "="))
        results = v["results"]
        for sub_cat, items in sorted(iter(results.items()),
                                     key=lambda sub_cat_items:
                                     (sub_cat_items[0], sub_cat_items[1])):
            print_list.append("{}:".format(sub_cat))
            for index, i in sorted(items.items()):
                dict_to_get[catID][index] = i["url"]
                print_index = " {}:[{}] ".format(str(catID), str(index))
                date = fmt_timestamp(i["timestamp"])
                date_string = '{:02d}/{}'.format(date.tm_mon, date.tm_year)
                size_date = "[{}|{}]".format(i["size"], date_string)
                filename = i["filename"]
                margin = 80 - (len(size_date) + len(print_index))
                file_info = filename.ljust(margin, "-")
                print_list.append("".join([print_index, file_info, size_date]))
        if display is True:
            print('\n'.join(print_list))
            print()  # padding
    return dict_to_get


def get_user_choice():
    """
    Get user file selection choice(s)

    """
    choice = input(
        "Enter file selection ('q' to quit, "
        "'usage' to review syntax, 'a' for all, "
        "'r' for regex-based filename matching):\n> ")
    if choice == "usage":
        print()
        print(select_blurb)
        print()
        return get_user_choice()
    elif choice.lower() in ("q", "quit", "exit"):
        remove_temp = input("Remove index file? (y/n): ")
        remove_temp = remove_temp.lower() in ('y', 'yes', '')
        clean_exit(remove_temp=remove_temp)
    else:
        return choice


def parse_selection(user_input):
    """
    Parses the user choice string and returns a dictionary
    of categories (keys) and choices within each category
    (values).

    """
    selections = {}
    parts = user_input.split(";")
    for p in parts:
        if len(p.split(":")) > 2:
            clean_exit("FATAL ERROR: can't parse desired input\n?-->'{}'"
                      .format(p))
        category, indices = p.split(":")
        category = int(category)
        selections[category] = []
        cat_list = selections[category]
        indices = indices.split(",")
        for i in indices:
            try:
                cat_list.append(int(i))  # if it's already an integer
            except ValueError:
                try:
                    start, stop = list(map(int, i.split("-")))
                except:
                    clean_exit("FATAL ERROR: can't parse desired "
                              "input\n?-->'{}'".format(i))
                add_range = list(range(start, stop + 1))
                for e in add_range:
                    cat_list.append(e)
    return selections


def url_format_checker(u):
    """
    Checks the URL string and corrects it to the JGI Genome
    Portal format in cases where it is differently formatted,
    e.g. links listed in Phytozome.

    Such malformed links are prepended with a string which breaks
    normal parsing, for example:
    "/ext-api/downloads/get_tape_file?blocking=true&url=" is
    prepended to the standard Genome Portal URL format for (all?)
    Phytozome links and needs to be removed for cURL to use it.

    """
    if "url=" in u:
        u = u.split("url=")[-1]  # take the bit after the prepended string
    return u


def get_org_name(xml_file):
    """
    Checks an XML file for organism name information,
    for cases where an XML file is used without organism
    information supplied by the user. Returns None if
    no organism name is found.

    XML entry format is: <organismDownloads name="org_name">

    """
    name_pattern = r"name=\"(.+)\""
    org_line = None
    with open(xml_file) as f:
        for l in f:
            if "organismDownloads" in l:  # standardized name indicator
                org_line = l.strip()
                break  # don't keep looking, already found
    try:
        org_name = re.search(name_pattern, org_line).group(1)
        return org_name
    except TypeError:  # org_line still None
        return None


def is_xml(filename):
    """
    Uses hex code at the beginning of a file to try to determine if it's an
    XML file or not. This seems to be occasionally necessary; if pulling
    files from JGI tape archives, the server may error out and provide an
    XML error document instead of the intended file. This function should
    return False on all downloaded files, although false positives have not
    been thoroughly investigated.

    Adapted from http://stackoverflow.com/a/13044946/3076552

    """
    xml_hex = "\x3c"  # hex code at beginning of XML files
    read_length = len(xml_hex)
    with open(filename) as f:
        try:
            file_start = f.read(read_length)
        except UnicodeDecodeError:  # compressed files
            return False
        if file_start.startswith(xml_hex):  # XML file
            return True
        else:  # hopefully all other file types
            return False


def hidden_xml_check(file_list):
    """
    Checks a file list for any files that are actually XML error files,
    but which were intended to be of another format. Returns a list of
    all files not failing the test.

    """
    for f in list(file_list):  # iterate over copy
        if is_xml(f):
            if not f.lower().endswith("xml"):  # not recognized properly
                print("ERROR: '{}' appears to be malformed and will be left "
                      "unmodified.".format(f))
                file_list.remove(f)  # don't try to process downstream
    return file_list


def byte_convert(byte_size):
    """
    Converts a number of bytes to a more human-readable
    format.

    """
    # Calculate and display total size of selected data
    adjusted = byte_size / (1024 * 1024)  # bytes to MB
    if adjusted < 1:
        adjusted = byte_size / 1024
        unit = "KB"
    elif adjusted < 1024:
        unit = "MB"
    else:
        adjusted /= 1024
        unit = "GB"
    size_string = "{:.2f} {}".format(adjusted, unit)
    return size_string


def is_broken(filename, min_size_bytes=20):
    """
    Rudimentary check to see if a file appears to be broken.
    
    """
    if (
        not os.path.isfile(filename) or
        os.path.getsize(filename) < min_size_bytes or 
        (is_xml(filename) and not filename.lower().endswith('xml'))
    ):
        return True
    else:
        return False


def download_from_url(url, timeout=30, retry=0, min_file_bytes=20):
    """
    Attempts to download a file from JGI servers using cURL.

    Returns a tuple of (filename, cURL command used, success boolean)
    
    """
    success = True
    url = url.replace('&amp;', '&')
    filename = re.search('.+/(.+$)', url).group(1)
    url_prefix = "https://genome.jgi.doe.gov"
    download_command = (
        "curl '{}{}' -b cookies "
        "> {}".format(url_prefix, url, filename)
    )
    if not is_broken(filename, min_file_bytes):
        success = True
        print("Skipping existing file {}".format(filename))
    else:
        print("Downloading '{}' using command:\n{}"
            .format(filename, download_command))
        # The next line doesn't appear to be needed to refresh the cookies.
        #    subprocess.call(login, shell=True)
        subprocess.call(download_command, shell=True)
        if is_broken(filename, min_file_bytes):
            success = False
            if retry > 0:
                # success = False
                # this may be needed if initial download fails
                alt_cmd = download_command.replace(
                    'blocking=true', 'blocking=false')
                current_retry = 1
                while current_retry <= retry:
                    if current_retry % 2 == 1:
                        retry_cmd = alt_cmd
                    else:
                        retry_cmd = download_command
                    print(
                        "Trying '{}' again due to download error ({}/{}):\n{}"
                        .format(filename, current_retry, retry, retry_cmd)
                    )
                    subprocess.call(retry_cmd, shell=True)
                    if not is_broken(filename, min_file_bytes):
                        success = True
                        break
                    current_retry += 1
                    time.sleep(1)

    return filename, download_command, success


def get_regex():
    """
    Get regex pattern from user, compile and return.
    
    """
    #TODO make this exit gracefully if user can't
    # manage to get a working regex
    compile_success = False
    while compile_success is False:
        pattern = input('Regex pattern: ')
        try:
            pattern = re.compile(pattern)
            compile_success = True
        except:
            print('[!] ERROR: Regex pattern failed to compile.')

    return re.compile(pattern)


def retry_from_failed(login_cmd, fail_log, timeout=60, retries=3):
    """
    Try to download from URLs in a previously-generated log file.
    
    """
    organism = os.path.basename(fail_log).split('.')[0]
    fail_log = open(fail_log, 'r')
    url_list = fail_log.read().splitlines()
    try:  # fails if unable to contact server
        subprocess.check_output(login_cmd, shell=True)
    except subprocess.CalledProcessError as error:
        clean_exit("Couldn't connect with server. Please check Internet "
                  "connection and retry.")
    downloaded, failed = download_list(url_list)

    print('Finished downloading {} files'.format(len(downloaded)))
    if failed:
        log_failed(organism, failed)
    
    return downloaded, failed


def log_failed(organism, failed_urls):
    """
    Write failed URLs to a local log file.
    
    """
    fail_log = '{}.failed.log'.format(organism)
    print(
        '{} failed downloads logged to {}'.format(len(failed_urls), fail_log))
    # write failed URLs to local file
    with open(fail_log, 'w') as f:
        f.write('\n'.join(failed_urls))


def download_list(url_list, timeout=120, retries=3):
    """
    Attempts download command on a list of partial file
    URLs (completed by download_from_url()).

    Returns a list of successfully-downloaded files and a
    list of unsuccessful URLs
    
    """
    # Run curl commands to retrieve selected files
    # Make sure the URL formats conforms to the Genome Portal format
    downloaded_files = []
    broken_urls = []
    subprocess.call(LOGIN_STRING, shell=True)
    start_time = time.time()
    for url in url_list:
        current_time = time.time()
        # refresh the session cookie every 5 minutes
        if current_time - start_time > 300:
            subprocess.call(LOGIN_STRING, shell=True)
            start_time = time.time()
        fn, cmd, success = download_from_url(
            url, timeout=timeout, retry=retries)
        if not success:
            broken_urls.append(url)
        else:
            downloaded_files.append(fn)
    
    return downloaded_files, broken_urls

# /FUNCTIONS


if __name__ == '__main__':
	main()
