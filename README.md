# mempig
MEM Pangenome Injection Genotyping

## Paths

```shell
DIR_BASE=/scratch/small_test_output
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

```shell
cd $DIR_BASE/cram
wget https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index

grep -Ff <(cat $DIR_BASE/impg/extracted.fasta.fai | cut -f 1 -d '#' | sort | uniq) $DIR_BASE/cram/1000G_698_related_high_coverage.sequence.index | cut -f 1  > $DIR_BASE/cram/to_download.txt
wget -i $DIR_BASE/cram/to_download.txt
ls $DIR_BASE/cram/*.cram | while read CRAM; do samtools index $CRAM; done
```

## MEM Pangenome Injection Genotyping

Build indexes:

```shell
mkdir -p $DIR_BASE/ropebwt3

# Construct a BWT for both strands of the input sequences
ropebwt3 build $DIR_BASE/impg/extracted.fasta -do $DIR_BASE/ropebwt3/extracted.fmd

# Sampled suffix array
ropebwt3 ssa -o $DIR_BASE/ropebwt3/extracted.fmd.ssa -s8 -t 1 $DIR_BASE/extracted.fmd

# Sequence lengths
seqtk comp $DIR_BASE/impg/extracted.fasta | cut -f1,2 | gzip > $DIR_BASE/ropebwt3/extracted.fmd.len.gz
```

Test with 39 samples on the C4 region from HPRCy1 pangenome:

```shell
cd $DIR_BASE/ropebwt3

ROI_BED=$DIR_BASE/roi/roi.bed
REFERENCE_CRAM_FASTA=$DIR_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd $DIR_BASE/ropebwt3/
ls $DIR_BASE/cram/*.final.cram | while read READS_CRAM; do
    NAME=$(basename $READS_CRAM)

    INPUT_READ=$DIR_BASE/$NAME.roi.fa.gz
    samtools view -T $REFERENCE_CRAM_FASTA -@ 16 -L $ROI_BED -M -b $READS_CRAM | samtools fasta -@ 24 - | gzip > $INPUT_READ
    # echo $READS_CRAM $INPUT_READ
    for l in `seq 9 60`; do
        for p in `seq 1 50`; do
            ropebwt3 mem -l $l $DIR_BASE/ropebwt3/extracted.fmd $INPUT_READ -p $p -t 24 > $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.tsv

            python3 $DIR_BASE/ropebwt3-to-paf.py $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.tsv <(cut -f 1,2 $DIR_BASE/impg/extracted.fasta.fai) $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.paf

            gfainject --gfa $DIR_BASE/odgi/graph.gfa --paf $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.paf > $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.gaf

            gafpack \
                -g $DIR_BASE/odgi/graph.gfa \
                -a $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.gaf \
                --len-scale | \
                gzip > $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.gafpack.gz

            /scratch/cosigt/cosigt -p $DIR_BASE/odgi/paths_matrix.tsv.gz -g $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.gafpack.gz -o $DIR_BASE/ropebwt3/ #-i $NAME

            cp $DIR_BASE/ropebwt3/cosigt_genotype.tsv $DIR_BASE/ropebwt3/$NAME.cosigt_genotype.l$l.p$p.tsv
            cp $DIR_BASE/ropebwt3/sorted_combos.tsv $DIR_BASE/ropebwt3/$NAME.sorted_combos.l$l.p$p.tsv

            grep $NAME $DIR_BASE/ropebwt3/$NAME.cosigt_genotype.l$l.p$p.tsv | awk -v OFS='\t' -v name=$NAME -v l=$l -v p=$p '{print(name,l,p,$0)}'

            rm $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.tsv $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.paf $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.gaf $DIR_BASE/ropebwt3/$NAME-vs-extracted.mem.l$l.p$p.gafpack.gz
        done
    done
done > $DIR_BASE/ropebwt3/C4.test.tsv
```
