
### Quick install and start ###
```
git clone https://github.com/zhangrengang/MicrobePhylogenomics
cd MicrobePhylogenomics

# install
conda env create -f MicrobePhylogenomics.yaml
conda activate MicrobePhylogenomics

# start
cd test
../microbePhylogenomics
```

### Inputs ###
1. `species.design`, a table including two columns: the first is species name or id and the second is the URL of Genbank or Refseq; e.g.:
```
Devosia_insulae ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/970/465/GCF_000970465.2_ASM97046v2
Devosia_limi    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/128/975/GCF_900128975.1_IMG-taxon_2582580747_annotated_assembly
Maritalea_myrionectae   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/433/515/GCF_003433515.1_ASM343351v1
Cucumibacter_marinus    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/429/865/GCF_000429865.1_ASM42986v1
Pelagibacterium_halotolerans    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/004/005/GCF_013004005.1_ASM1300400v1
Devosia_sediminis       ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/411/825/GCF_016411825.1_ASM1641182v1
Devosia_faecipullorum   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/158/295/GCF_015158295.1_ASM1515829v1
Devosiaceae_bacterium   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/192/585/GCA_019192585.1_ASM1919258v1
Devosia_aurantiaca      ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/058/215/GCF_011058215.1_ASM1105821v1
Pelagibacterium_limicola        ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/694/405/GCF_015694405.1_ASM1569440v1
```
Note that the blank or other special characters are not supported in the first column.

2. Outgroup for the final phylogeny (defualt: the first row in `species.design`). Specify an outgroup:
```
../microbePhylogenomics species.design Cucumibacter_marinus
```

### Outputs ###
Check the output phylogeny `Species_Tree/singlecopy_cds.fa.treefile.rooted` by:
```
nw_display Species_Tree/singlecopy_cds.fa.treefile.rooted
```
