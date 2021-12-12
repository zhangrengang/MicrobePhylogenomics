import sys
from Bio import SeqIO

def main(sp=sys.argv[1], seqfile=sys.argv[2], fout=sys.stdout, sep='|'):
	for rc in SeqIO.parse(seqfile, 'fasta'):
		rc.id = sp+ sep + rc.id
		SeqIO.write(rc, fout, 'fasta')

if __name__ == '__main__':
	main()
