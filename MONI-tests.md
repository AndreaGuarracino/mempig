
```shell
cd $dir_base

PATH_READS=/home/guarracino/git/moni/data/SARS-CoV2/reads.fastq.gz
PATH_PANG=/home/guarracino/git/moni/data/SARS-CoV2/SARS-CoV2.1k.fa.gz
NAME=sars-cov2

PATH_READS=/scratch/03_08_22_R941_HG002_2_Guppy_6.1.2_5mc_cg_prom_sup.fq.gz
PATH_PANG=/scratch/chm13v2+grch38+hg002v101.fa.gz
NAME=chm13v2+grch38+hg002v101

# Build the index
\time -v moni build -r $PATH_PANG -o /scratch/$NAME -f

# Compute the matching statistics of reads against the pangenome
$MONI ms -i sars-cov2 -p $PATH_READS -o reads
# Output:
# - reads.lengths: length for the matching statistics of reads against the pangenome
# - reads.pointers: position //

# Compute the MEMs of reads against the pangenome
$MONI mems -i sars-cov2 -p $PATH_READS -o reads
# Output:
# - reads.mems: MEMs reposted as pairs of position and lengths in a FASTA-like format

$MONI mems -i sars-cov2 -p $PATH_READS -o reads -s

# Compute the MEM extension of reads against the pangenome
$MONI extend -i sars-cov2 -p $PATH_READS -o reads
```

```shell
alias moni=/home/guarracino/git/moni-align/build/moni
moni build -f scerevisiae8.fa.gz -o xxx/scerevisiae8
moni mems -i xxx/scerevisiae8 -p 1read.fq

# MONI
conda activate /lizardfs/guarracino/condatools/moni/0.2.2
mkdir -p /scratch/small_test_output/moni
moni build -r /scratch/small_test_output/impg/extracted.fasta -o /scratch/small_test_output/moni/extracted -f
# Compute the MEMs of reads against the pangenome
moni mems -i /scratch/small_test_output/moni/extracted -p $PATH_READS -o reads -s

## rvarkiR: Don't use the mems command. Use align with --report_mems option
```
