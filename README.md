# mempig
MEM Pangenome Injection Genotyping

## Paths

```shell
DIR_BASE=/lizardfs/guarracino/mempig
```

## Tools

### cosigt

Install `go`:

```shell
cd /lizardfs/guarracino/tools
wget -c https://go.dev/dl/go1.23.1.linux-amd64.tar.gz
tar -C /lizardfs/guarracino/tools -xzf go1.23.1.linux-amd64.tar.gz && rm go1.23.1.linux-amd64.tar.gz
#Add 'export PATH="/lizardfs/guarracino/tools/go/bin:$PATH"' to ~/.zshrc
```

Build `cosigt`:

```shell
cd /lizardfs/guarracino/git
git clone https://github.com/davidebolo1993/cosigt
cd cosigt
git checkout v0.1.0 # Newer versions have mandatory clustering
go mod init cosigt && go mod tidy && go build cosigt
#Add 'export PATH="/lizardfs/guarracino/git/cosigt:$PATH"' to ~/.zshrc
```

Create a `conda` environment for `cosigt` with all its dependencies but R:

```shell
conda create --prefix /lizardfs/guarracino/condatools/cosigt -c conda-forge -c bioconda -c anaconda -c vikky34v snakemake=7.32.4 cookiecutter=2.6.0 bwa-mem2=2.2.1 megadepth samtools=1.21 bedtools=2.31.1 python=3.9 pyyaml=6.0.2 pandas -y #r-base r-rjson=0.2.23 r-reshape2=1.4.4 r-nbclust=3.0.1 r-data.table r-ggplot2=3.5.1 r-dendextend=1.18.1 r-gggenes=0.5.1 bioconductor-rtracklayer time -y
```

It assumes that all `pggb` and its tools (`wfmash`, `seqwish`, `smoothxg`, `odgi`, `gfaffix`), `samtools`, `bedtools` are in `$PATH` are installed and included in system's `$PATH` environment variable so they can be executed from any directory.

### ropebwt3

```shell
cd /lizardfs/guarracino/git
git clone https://github.com/lh3/ropebwt3
cd ropebwt3
git checkout f77e399634456f6c0d3194a9c146878d79ab34c3
make
#Add 'export PATH="/lizardfs/guarracino/git/ropebwt3:$PATH"' to ~/.zshrc
```

### MONI (not fully functional)

```shell
conda create --prefix /lizardfs/guarracino/condatools/moni/0.2.2 -c conda-forge -c bioconda moni=0.2.2 -y

# cd /lizardfs/guarracino/git
# git clone https://github.com/maxrossi91/moni

# https://github.com/maxrossi91/moni/issues/2#issuecomment-2362544812
# cd moni
# mkdir build
# cd build
# cmake ..
# make -j 16 VERBOSE=1

# cd /lizardfs/guarracino/tools
# wget https://github.com/maxrossi91/moni/releases/download/v0.2.2/moni-0.2.2-Linux.tar.gz
# tar -xzvf moni-0.2.2-Linux.tar.gz && rm moni-0.2.2-Linux.tar.gz
# mv moni-0.2.2-Linux moni-0.2.2
```

## Data

### HGSVC3

```shell
cd /scratch
ls /lizardfs/guarracino/pangenomes/HGSVC3/*.fa.gz | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.' | while read SAMPLE; do URL=$(cat $DIR_BASE/data/1000G_*.index | grep $SAMPLE -m 1 | cut -f 1); echo $URL; wget -c $URL; wget -c $URL.crai; done
mv *.cram* $DIR_BASE/cram
```

### 39 samples

```shell
cd $DIR_BASE/cram
wget https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index

grep -Ff <(cat $DIR_BASE/impg/extracted.fasta.fai | cut -f 1 -d '#' | sort | uniq) $DIR_BASE/cram/1000G_698_related_high_coverage.sequence.index | cut -f 1  > $DIR_BASE/cram/to_download.txt
wget -i $DIR_BASE/cram/to_download.txt
ls $DIR_BASE/cram/*.cram | while read CRAM; do samtools index $CRAM; done
```

## MEM Pangenome Injection Genotyping

## `mempig` with 39 samples on the C4 region from HPRCy1 pangenome

Build indexes:

```shell
mkdir -p $DIR_BASE/ropebwt3/indexes

# Construct a BWT for both strands of the input sequences
ropebwt3 build $DIR_BASE/impg/extracted.fasta -do $DIR_BASE/ropebwt3/indexes/extracted.fmd

# Sampled suffix array
ropebwt3 ssa -o $DIR_BASE/ropebwt3/indexes/extracted.fmd.ssa -s8 -t 1 $DIR_BASE/indexes/extracted.fmd

# Sequence lengths
seqtk comp $DIR_BASE/impg/extracted.fasta | cut -f1,2 | gzip > $DIR_BASE/ropebwt3/indexes/extracted.fmd.len.gz
```

Run:

```shell
# Prepare clusters (they are mandatory for cosigt > v0.1.0)
Rscript /lizardfs/guarracino/git/cosigt/cosigt_smk/workflow/scripts/cluster.r $DIR_BASE/odgi/chopped.similarity.tsv $DIR_BASE/odgi/chopped.similarity.clusters.json

ROI_BED=$DIR_BASE/roi/roi.bed
REFERENCE_CRAM_FASTA=$DIR_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

mkdir -p $DIR_BASE/reads
cd $DIR_BASE/ropebwt3/
ls $DIR_BASE/cram/*.final.cram | while read READS_CRAM; do
    NAME=$(basename $READS_CRAM .final.cram)
    echo $READS_CRAM $INPUT_READ

    INPUT_READ=$DIR_BASE/reads/$NAME.roi.fa.gz
    samtools view -T $REFERENCE_CRAM_FASTA -@ 16 -L $ROI_BED -M -b $READS_CRAM | samtools fasta -@ 24 - | bgzip -l 9 -@ 8 > $INPUT_READ
    for l in `seq 5 100`; do
        sbatch -c 8 -p tux --job-name "${NAME}-l${l}" $DIR_BASE/scripts/run_ropebwt3.sh "$DIR_BASE" "$NAME" "$INPUT_READ" "$l" "$DIR_BASE/odgi/graph.gfa" "$DIR_BASE/odgi/chopped.similarity.tsv $DIR_BASE/odgi/chopped.similarity.clusters.json" "$DIR_BASE/impg/extracted.fasta"
    done
done

# To check that everything went well
grep skipping slurm-** -c | cut -f 2 -d ':' | uniq -c
#   3744 0
```

Collect results:

```shell
(echo -n "l p " | tr ' ' '\t'; find $DIR_BASE/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" -exec head -n 1 {} \; -quit | sed 's/#//g') > $DIR_BASE/ropebwt3/C4.test.tsv
find $DIR_BASE/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" | while read TSV; do
    l=$(basename $TSV .tsv | cut -f 3 -d '.' | sed 's/l//g')
    p=$(basename $TSV .tsv | cut -f 4 -d '.' | sed 's/p//g')
   # echo $TSV $NAME $l $p

    grep '^#' -v $TSV | awk -v OFS='\t' -v l=$l -v p=$p '{print(l,p,$0)}'
done | sort -k 1,1 -k 2,2n -k 3,3n -T /scratch >> $DIR_BASE/ropebwt3/C4.test.tsv
```

## `mempig` with 39 samples on the amylase locus from HPRCv2 pangenome (387 haplotypes from 193 individuals + chm13v2.0)

Prepare the pangenome:

```shell
cd /lizardfs/guarracino/robertsonian_translocation/partitioning/HPRCv2
cat */*.aln.paf > /scratch/HPRCv2-vs-chm13.aln.paf

mkdir -p $DIR_BASE/amylase
impg -p /scratch/HPRCv2-vs-chm13.aln.paf -r chr1:103304997-103901127 > $DIR_BASE/amylase/HPRCv2-vs-chm13.chr1_103304997_103901127.bed

bedtools sort -i HPRCv2-vs-chm13.chr1_103304997_103901127.bed | bedtools merge -d 1000000 | cut -f 1 -d '#' | sort | uniq -c | sort -k 1,1nr | awk '$1 == 2' | tr -s ' ' | cut -f 3 -d ' ' > HPRCv2.samples-with-two-haplotypes.txt

bedtools sort -i HPRCv2-vs-chm13.chr1_103304997_103901127.bed | bedtools merge -d 1000000 | grep -Ff HPRCv2.samples-with-two-haplotypes.txt > HPRCv2-vs-chm13.chr1_103304997_103901127.merged-and-two-haplitypes.bed

(samtools faidx /lizardfs/guarracino/robertsonian_translocation/assemblies/chm13v2.0.fa.gz chr1:103304997-103901127; ls /lizardfs/guarracino/pangenomes/HPRCv2/*.fa.gz | while read FASTA; do
    samtools faidx \
		-r <(grep -Ff <(cut -f 1 $FASTA.fai) HPRCv2-vs-chm13.chr1_103304997_103901127.merged-and-two-haplitypes.bed | awk '{{print $1":"$2+1"-"$3}}') \
        $FASTA
done) | bgzip -l 9 -@ 24 > amylase387.fa.gz

sbatch -c 96 -p tux --job-name pggb-amy --wrap "hostname; cd /scratch; pggb -i $DIR_BASE/amylase/amylase387.fa.gz -o $DIR_BASE/amylase/pggb.amylase387 -c 2 -D /scratch -G 1300 -t 96 -n 387"

odgi similarity -i $DIR_BASE/amylase/pggb.amylase387/amylase387.fa.gz.8fea6cf.11fba48.a2c905a.smooth.final.og > $DIR_BASE/amylase/pggb.amylase387/amylase387.similarity.tsv -t 24 -P
```

Build indexes:

```shell
mkdir -p $DIR_BASE/ropebwt3/indexes

# Construct a BWT for both strands of the input sequences
ropebwt3 build $DIR_BASE/amylase/amylase387.fa.gz -do $DIR_BASE/ropebwt3/indexes/extracted.fmd

# Sampled suffix array
ropebwt3 ssa -o $DIR_BASE/ropebwt3/indexes/extracted.fmd.ssa -s8 -t 1 $DIR_BASE/ropebwt3/indexes/extracted.fmd

# Sequence lengths
seqtk comp $DIR_BASE/amylase/amylase387.fa.gz | cut -f1,2 | gzip > $DIR_BASE/ropebwt3/indexes/extracted.fmd.len.gz
```

```shell
# Prepare clusters (they are mandatory for cosigt > v0.1.0)
Rscript /lizardfs/guarracino/git/cosigt/cosigt_smk/workflow/scripts/cluster.r /lizardfs/guarracino/mempig/amylase/pggb.amylase387/amylase387.similarity.tsv /lizardfs/guarracino/mempig/amylase/pggb.amylase387/amylase387.similarity.clusters.json

ROI_BED=$DIR_BASE/roi/roi.bed
REFERENCE_CRAM_FASTA=$DIR_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

mkdir -p $DIR_BASE/reads
cd $DIR_BASE/ropebwt3/
ls $DIR_BASE/cram/*.final.cram | while read READS_CRAM; do
    NAME=$(basename $READS_CRAM .final.cram)
    echo $READS_CRAM $INPUT_READ

    INPUT_READ=$DIR_BASE/reads/$NAME.roi.fa.gz
    samtools view -T $REFERENCE_CRAM_FASTA -@ 16 -L $ROI_BED -M -b $READS_CRAM | samtools fasta -@ 24 - | bgzip -l 9 -@ 8 > $INPUT_READ
    for l in `seq 5 30`; do
        sbatch -c 8 -p tux --job-name "${NAME}-l${l}" $DIR_BASE/scripts/run_ropebwt3.sh "$DIR_BASE" "$NAME" "$INPUT_READ" "$l" "$DIR_BASE/amylase/pggb.amylase387/amylase387.fa.gz.8fea6cf.11fba48.a2c905a.smooth.final.gfa" "$DIR_BASE/amylase/pggb.amylase387/amylase387.similarity.clusters.json" "$DIR_BASE/amylase/amylase387.fa.gz"
    done
done


(echo -n "l p " | tr ' ' '\t'; find $DIR_BASE/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" -exec head -n 1 {} \; -quit | sed 's/#//g') > $DIR_BASE/ropebwt3/amylase387.test.tsv
find $DIR_BASE/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" | while read TSV; do
    l=$(basename $TSV .tsv | cut -f 3 -d '.' | sed 's/l//g')
    p=$(basename $TSV .tsv | cut -f 4 -d '.' | sed 's/p//g')
   # echo $TSV $NAME $l $p

    grep '^#' -v $TSV | awk -v OFS='\t' -v l=$l -v p=$p '{print(l,p,$0)}'
done | sort -k 1,1 -k 2,2n -k 3,3n -T /scratch >> $DIR_BASE/ropebwt3/amylase387.test.tsv



# From https://github.com/sudmantlab/amylase_diversity_project/blob/main/pangenome/pggb/20231102_graph/haplotype_all_structures_plotinfo.chm13_colors.bed
chr1_chm13_103304997_103901127_0	99	115472	bundle0#CCCCCC
chr1_chm13_103304997_103901127_0	492840	596058	bundle1#666666


samtools faidx /lizardfs/guarracino/robertsonian_translocation/assemblies/chm13v2.0.fa.gz chr1
chr1:103305096-103420469
chr1:103797837-103901055
```