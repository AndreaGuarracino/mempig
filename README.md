# MEMPIG (MEM Pangenome Injection Genotyping)

## Paths

```shell
dir_base=/lizardfs/guarracino/mempig

export PATH="/lizardfs/guarracino/tools/bedtools2/bin:$PATH"
export PATH="/lizardfs/guarracino/tools/samtools-1.21:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/gafpack/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/GFAffix/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/impg/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/gfainject/target/release:$PATH"

export PATH="/lizardfs/guarracino/tools_for_cosigt/wfmash/build/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/seqwish/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/smoothxg/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/odgi/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/pggb:$PATH"

export PATH="/lizardfs/guarracino/git/seqtk:$PATH"
```

## Tools

### cosigt

Install `go`:

```shell
cd /lizardfs/guarracino/tools
wget -c https://go.dev/dl/go1.24.2.linux-amd64.tar.gz
tar -C /lizardfs/guarracino/tools -xzf go1.24.2.linux-amd64.tar.gz && rm go1.24.2.linux-amd64.tar.gz
#Add 'export PATH="/lizardfs/guarracino/tools/go/bin:$PATH"' to ~/.zshrc
```

Build `cosigt`:

```shell
cd /lizardfs/guarracino/git
git clone https://github.com/davidebolo1993/cosigt
cd cosigt
git checkout ed9f117a7e1dad23e262e9d78dd777a97a0fde74
go mod init cosigt && go mod tidy && go build cosigt
#Add 'export PATH="/lizardfs/guarracino/git/cosigt:$PATH"' to ~/.zshrc
```

## bwa-mem2

```shell
conda create --prefix /lizardfs/guarracino/condatools/bwa-mem2/2.2.1 -c conda-forge -c bioconda -c anaconda bwa-mem2=2.2.1 -y
```

<!-- Create a `conda` environment for `cosigt` with all its dependencies but R:

```shell
conda create --prefix /lizardfs/guarracino/condatools/cosigt -c conda-forge -c bioconda -c anaconda -c vikky34v snakemake=7.32.4 cookiecutter=2.6.0 bwa-mem2=2.2.1 megadepth samtools=1.21 python=3.9 pyyaml=6.0.2 pandas -y #bedtools=2.31.1  r-base r-rjson=0.2.23 r-reshape2=1.4.4 r-nbclust=3.0.1 r-data.table r-ggplot2=3.5.1 r-dendextend=1.18.1 r-gggenes=0.5.1 bioconductor-rtracklayer time -y
``` -->

### ropebwt3

```shell
cd /lizardfs/guarracino/git
git clone https://github.com/lh3/ropebwt3
cd ropebwt3
git checkout 36a64118322cdd5c20690ec05c28bf0b22c71a98
make
#Add 'export PATH="/lizardfs/guarracino/git/ropebwt3:$PATH"' to ~/.zshrc
```

<!-- ### MONI (not fully functional)

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
``` -->

### PGGB and Co.

It assumes that all `pggb` and its tools (`wfmash`, `seqwish`, `smoothxg`, `odgi`, `gfaffix`), `samtools`, `bedtools` are in `$PATH` are installed and included in system's `$PATH` environment variable so they can be executed from any directory.

```shell
export PATH="/lizardfs/guarracino/tools/bedtools2/bin:$PATH"
export PATH="/lizardfs/guarracino/tools/samtools-1.21:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/gafpack/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/GFAffix/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/impg/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/gfainject/target/release:$PATH"

export PATH="/lizardfs/guarracino/tools_for_cosigt/wfmash/build/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/seqwish/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/smoothxg/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/odgi/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/pggb:$PATH"
```

## Data

### Reference

```shell
mkdir -p $dir_base/reference
cd $dir_base/reference

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

# Remove alt contigs to avoid fragmented ROI-specific pangenome sequences
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa $(grep chr GRCh38_full_analysis_set_plus_decoy_hla.fa.fai  | grep '_' -v | cut -f 1) | bgzip -l 9 -@ 24 > GRCh38.fa.gz
samtools faidx GRCh38.fa.gz
```

### HGSVC3 (63 samples)

```shell
cd /scratch
ls /lizardfs/guarracino/pangenomes/HGSVC3/*.fa.gz | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.' | while read sample; do url=$(cat $dir_base/data/1000G_*.index | grep $sample -m 1 | cut -f 1); echo $url; wget -c $url; wget -c $url.crai; done
mv *.cram* $dir_base/cram
```

<!-- ### 39 samples

```shell
cd $dir_base/cram
wget https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index

grep -Ff <(cat $dir_base/impg/extracted.fasta.fai | cut -f 1 -d '#' | sort | uniq) $dir_base/cram/1000G_698_related_high_coverage.sequence.index | cut -f 1  > $dir_base/cram/to_download.txt
wget -i $dir_base/cram/to_download.txt
ls $dir_base/cram/*.cram | while read CRAM; do samtools index $CRAM; done
``` -->

## MEM Pangenome Injection Genotyping

## `mempig` with 63 HGSVC3's short read samples on the C4 region from the same 63 HGSVC3's assemblies

Run:

```shell
# Prepare clusters (they are mandatory for cosigt > v0.1.0)
Rscript /lizardfs/guarracino/git/cosigt/cosigt_smk/workflow/scripts/cluster.r $dir_base/odgi/chopped.similarity.tsv $dir_base/odgi/chopped.similarity.clusters.json

ROI_BED=$dir_base/roi/roi.bed
REFERENCE_CRAM_FASTA=$dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

mkdir -p $dir_base/reads
cd $dir_base/ropebwt3/
ls $dir_base/cram/*.final.cram | while read READS_CRAM; do
    NAME=$(basename $READS_CRAM .final.cram)
    echo $READS_CRAM $INPUT_READ

    INPUT_READ=$dir_base/reads/$NAME.roi.fa.gz
    samtools view -T $REFERENCE_CRAM_FASTA -@ 16 -L $ROI_BED -M -b $READS_CRAM | samtools fasta -@ 24 - | bgzip -l 9 -@ 8 > $INPUT_READ
    for l in `seq 5 100`; do
        sbatch -c 8 -p tux --job-name "${NAME}-l${l}" $dir_base/scripts/run_ropebwt3.sh "$dir_base" "$NAME" "$INPUT_READ" "$l" "$dir_base/odgi/graph.gfa" "$dir_base/odgi/chopped.similarity.tsv $dir_base/odgi/chopped.similarity.clusters.json" "$dir_base/impg/extracted.fasta"
    done
done

# To check that everything went well
grep skipping slurm-** -c | cut -f 2 -d ':' | uniq -c
#   3744 0
```

```shell
mkdir -p $dir_base/wfmash
cd $dir_base/wfmash
ls /lizardfs/guarracino/pangenomes/HGSVC3/*.fa.gz | grep -Ff <(ls $dir_base/cram/*.cram | while read f; do echo $(basename $f .final.cram); done) -w | while read fasta; do
    sample=$(basename $fasta .fa.gz)

    sbatch -p allnodes -c 24 --job-name $sample-vs-grch38 --wrap "hostmame; cd /scratch; wfmash $dir_base/reference/GRCh38.fa.gz $fasta -s 10k -p 95 -t 24 > $dir_base/wfmash/$sample-vs-grch38.aln.paf"
done

mkdir -p $dir_base/regions_of_interest
echo -e "chr6\t31972057\t32055418" > $dir_base/regions_of_interest/C4.bed

mkdir -p $dir_base/impg
ls $dir_base/wfmash/*.paf | while read paf; do
    sample=$(basename $paf .paf)

    impg query \
        -p $paf \
        -b $dir_base/regions_of_interest/C4.bed
done > $dir_base/impg/C4.projected.bedpe
bedtools sort -i $dir_base/impg/C4.projected.bedpe | \
    bedtools merge -d 100000 > $dir_base/impg/C4.merged.bed

ls $dir_base/wfmash/*.paf | while read paf; do
    sample=$(basename $paf "-vs-grch38.aln.paf")
    fasta=/lizardfs/guarracino/pangenomes/HGSVC3/$sample.fa.gz

    bedtools getfasta -fi $fasta -bed <(grep $sample -w $dir_base/impg/C4.merged.bed)
done | bgzip -@ 8 > $dir_base/impg/C4.extracted.fa.gz
samtools faidx $dir_base/impg/C4.extracted.fa.gz

mkdir $dir_base/pggb
pggb -i $dir_base/impg/C4.extracted.fa.gz -o  $dir_base/pggb/C4 -t 48 -D /scratch
mv $dir_base/pggb/C4/*smooth.final.og $dir_base/pggb/C4/C4.final.og # Rename the final ODGI graph in a more human-friendly way

mkdir $dir_base/odgi
odgi chop \
    -i $dir_base/pggb/C4/C4.final.og \
    -c 32 \
    -o $dir_base/odgi/C4.chopped.og
odgi paths \
    -i $dir_base/odgi/C4.chopped.og \
    -H | \
    cut -f 1,4- | \
    gzip > $dir_base/odgi/C4.paths_matrix.tsv.gz

######################################################################################################

conda activate /lizardfs/guarracino/condatools/bwa-mem2/2.2.1
bwa-mem2 index $dir_base/impg/C4.extracted.fa.gz
mkdir -p $dir_base/alignments/C4
cat impg/C4.extracted.fa.gz.fai | cut -f 1 -d '#' > samples.txt
ls $dir_base/cram/*cram | grep -Ff samples.txt | while read cram; do
    echo $cram
    sample=$(basename $cram .cram)

    # Extract reads covering the C4 region and then align them against the pangenome
    samtools view \
        -T $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        -L $dir_base/regions_of_interest/C4.bed \
        -M \
        -b \
        $cram | \
        samtools sort -n | \
        samtools fasta | \
            bwa-mem2 mem -t 6 $dir_base/impg/C4.extracted.fa.gz - | \
            samtools view -b -F 4 -@ 2 - \
            > $dir_base/alignments/C4/$sample.reads_vs_extracted.bam
done
conda deactivate

#=====================================================================================================

mkdir -p $dir_base/ropebwt3/C4/indexes
# Construct a BWT for both strands of the input sequences
ropebwt3 build $dir_base/impg/C4.extracted.fa.gz -do $dir_base/ropebwt3/C4/indexes/C4.extracted.fmd
# Sampled suffix array
ropebwt3 ssa -o $dir_base/ropebwt3/C4/indexes/C4.extracted.fmd.ssa -s8 -t 24 $dir_base/ropebwt3/C4/indexes/C4.extracted.fmd
# Sequence lengths
seqtk comp $dir_base/impg/C4.extracted.fa.gz | cut -f1,2 | gzip > $dir_base/ropebwt3/C4/indexes/C4.extracted.fmd.len.gz

ls $dir_base/cram/*cram | grep -Ff samples.txt | while read cram; do
    echo $cram
    sample=$(basename $cram .cram)

    # Extract reads covering the C4 region and MEME them against the pangenome
    samtools view \
        -T $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        -L $dir_base/regions_of_interest/C4.bed \
        -M \
        -b \
        $cram | \
        samtools fasta | \
            ropebwt3 mem -l 17 $dir_base/ropebwt3/C4/indexes/C4.extracted.fmd - -p 1 -t 6 > $dir_base/ropebwt3/C4/$sample.reads_vs_extracted.mem.tsv

    fasta=$dir_base/impg/C4.extracted.fa.gz
    python3 $dir_base/scripts/ropebwt3-to-paf.py $dir_base/ropebwt3/C4/$sample.reads_vs_extracted.mem.tsv <(cut -f 1,2 $fasta.fai) $dir_base/ropebwt3/C4/$sample.reads_vs_extracted.mem.paf
done

######################################################################################################

odgi view \
    -i $dir_base/odgi/C4.chopped.og \
    -g > $dir_base/odgi/C4.chopped.gfa
ls $dir_base/cram/*cram | grep -Ff samples.txt | while read cram; do
    echo $cram
    sample=$(basename $cram .cram)

    gfainject \
        --gfa $dir_base/odgi/C4.chopped.gfa \
        --paf $dir_base/ropebwt3/C4/$sample.reads_vs_extracted.mem.paf \
        > $dir_base/alignments/C4/$sample.injected.gaf
done
        #--bam $dir_base/alignments/C4/$sample.reads_vs_extracted.bam \

ls $dir_base/cram/*cram | grep -Ff samples.txt | while read cram; do
    echo $cram
    sample=$(basename $cram .cram)

    gafpack \
        --gfa $dir_base/odgi/C4.chopped.gfa \
        --gaf $dir_base/alignments/C4/$sample.injected.gaf \
        --len-scale | \
        gzip > $dir_base/alignments/C4/$sample.coverage.gafpack.gz
done

mkdir -p $dir_base/clusters
# First generate the similarity matrix using odgi
odgi similarity \
    -i $dir_base/odgi/C4.chopped.og \
    > $dir_base/odgi/C4.similarity.tsv
# Download the clustering script from cosigt repository
wget https://raw.githubusercontent.com/davidebolo1993/cosigt/16b18815cf9fdfcbf2afbf588a02740c27941ee3/cosigt_smk/workflow/scripts/cluster.r -P clusters
# Run the clustering script to generate the JSON
Rscript clusters/cluster.r $dir_base/odgi/C4.similarity.tsv $dir_base/clusters/C4.clusters.json

mkdir -p $dir_base/cosigt/C4
ls $dir_base/cram/*cram | grep -Ff samples.txt | while read cram; do
    echo $cram
    sample=$(basename $cram .cram)

    cosigt \
        -i $sample \
        -p $dir_base/odgi/C4.paths_matrix.tsv.gz \
        -g $dir_base/alignments/C4/$sample.coverage.gafpack.gz \
        -c $dir_base/clusters/C4.clusters.json \
        -o $dir_base/cosigt/C4

    mv $dir_base/cosigt/C4/cosigt_genotype.tsv $dir_base/cosigt/C4/$sample.cosigt_genotype.tsv
    mv $dir_base/cosigt/C4/sorted_combos.tsv $dir_base/cosigt/C4/$sample.sorted_combos.tsv
done

grep 'final' -h $dir_base/cosigt/C4/*.cosigt_genotype.tsv | column -t
```

## `mempig` with 39 short read samples on the C4 region from the HPRCy1 pangenome

Build indexes:

```shell
mkdir -p $dir_base/ropebwt3/indexes

# Construct a BWT for both strands of the input sequences
ropebwt3 build $dir_base/impg/extracted.fasta -do $dir_base/ropebwt3/indexes/extracted.fmd

# Sampled suffix array
ropebwt3 ssa -o $dir_base/ropebwt3/indexes/extracted.fmd.ssa -s8 -t 1 $dir_base/indexes/extracted.fmd

# Sequence lengths
seqtk comp $dir_base/impg/extracted.fasta | cut -f1,2 | gzip > $dir_base/ropebwt3/indexes/extracted.fmd.len.gz
```

Run:

```shell
# Prepare clusters (they are mandatory for cosigt > v0.1.0)
Rscript /lizardfs/guarracino/git/cosigt/cosigt_smk/workflow/scripts/cluster.r $dir_base/odgi/chopped.similarity.tsv $dir_base/odgi/chopped.similarity.clusters.json

ROI_BED=$dir_base/roi/roi.bed
REFERENCE_CRAM_FASTA=$dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

mkdir -p $dir_base/reads
cd $dir_base/ropebwt3/
ls $dir_base/cram/*.final.cram | while read READS_CRAM; do
    NAME=$(basename $READS_CRAM .final.cram)
    echo $READS_CRAM $INPUT_READ

    INPUT_READ=$dir_base/reads/$NAME.roi.fa.gz
    samtools view -T $REFERENCE_CRAM_FASTA -@ 16 -L $ROI_BED -M -b $READS_CRAM | samtools fasta -@ 24 - | bgzip -l 9 -@ 8 > $INPUT_READ
    for l in `seq 5 100`; do
        sbatch -c 8 -p tux --job-name "${NAME}-l${l}" $dir_base/scripts/run_ropebwt3.sh "$dir_base" "$NAME" "$INPUT_READ" "$l" "$dir_base/odgi/graph.gfa" "$dir_base/odgi/chopped.similarity.tsv $dir_base/odgi/chopped.similarity.clusters.json" "$dir_base/impg/extracted.fasta"
    done
done

# To check that everything went well
grep skipping slurm-** -c | cut -f 2 -d ':' | uniq -c
#   3744 0
```

Collect results:

```shell
(echo -n "l p " | tr ' ' '\t'; find $dir_base/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" -exec head -n 1 {} \; -quit | sed 's/#//g') > $dir_base/ropebwt3/C4.test.tsv
find $dir_base/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" | while read TSV; do
    l=$(basename $TSV .tsv | cut -f 3 -d '.' | sed 's/l//g')
    p=$(basename $TSV .tsv | cut -f 4 -d '.' | sed 's/p//g')
   # echo $TSV $NAME $l $p

    grep '^#' -v $TSV | awk -v OFS='\t' -v l=$l -v p=$p '{print(l,p,$0)}'
done | sort -k 1,1 -k 2,2n -k 3,3n -T /scratch >> $dir_base/ropebwt3/C4.test.tsv
```

## `mempig` with 39 samples on the amylase locus from HPRCv2 pangenome (387 haplotypes from 193 individuals + chm13v2.0)

Prepare the pangenome:

```shell
cd /lizardfs/guarracino/robertsonian_translocation/partitioning/HPRCv2
cat */*.aln.paf > /scratch/HPRCv2-vs-chm13.aln.paf

mkdir -p $dir_base/amylase
impg -p /scratch/HPRCv2-vs-chm13.aln.paf -r chr1:103304997-103901127 > $dir_base/amylase/HPRCv2-vs-chm13.chr1_103304997_103901127.bed

bedtools sort -i HPRCv2-vs-chm13.chr1_103304997_103901127.bed | bedtools merge -d 1000000 | cut -f 1 -d '#' | sort | uniq -c | sort -k 1,1nr | awk '$1 == 2' | tr -s ' ' | cut -f 3 -d ' ' > HPRCv2.samples-with-two-haplotypes.txt

bedtools sort -i HPRCv2-vs-chm13.chr1_103304997_103901127.bed | bedtools merge -d 1000000 | grep -Ff HPRCv2.samples-with-two-haplotypes.txt > HPRCv2-vs-chm13.chr1_103304997_103901127.merged-and-two-haplitypes.bed

(samtools faidx /lizardfs/guarracino/robertsonian_translocation/assemblies/chm13v2.0.fa.gz chr1:103304997-103901127; ls /lizardfs/guarracino/pangenomes/HPRCv2/*.fa.gz | while read FASTA; do
    samtools faidx \
		-r <(grep -Ff <(cut -f 1 $FASTA.fai) HPRCv2-vs-chm13.chr1_103304997_103901127.merged-and-two-haplitypes.bed | awk '{{print $1":"$2+1"-"$3}}') \
        $FASTA
done) | bgzip -l 9 -@ 24 > amylase387.fa.gz

sbatch -c 96 -p tux --job-name pggb-amy --wrap "hostname; cd /scratch; pggb -i $dir_base/amylase/amylase387.fa.gz -o $dir_base/amylase/pggb.amylase387 -c 2 -D /scratch -G 1300 -t 96 -n 387"

odgi similarity -i $dir_base/amylase/pggb.amylase387/amylase387.fa.gz.8fea6cf.11fba48.a2c905a.smooth.final.og > $dir_base/amylase/pggb.amylase387/amylase387.similarity.tsv -t 24 -P
```

Build indexes:

```shell
mkdir -p $dir_base/ropebwt3/indexes

# Construct a BWT for both strands of the input sequences
ropebwt3 build $dir_base/amylase/amylase387.fa.gz -do $dir_base/ropebwt3/indexes/extracted.fmd

# Sampled suffix array
ropebwt3 ssa -o $dir_base/ropebwt3/indexes/extracted.fmd.ssa -s8 -t 1 $dir_base/ropebwt3/indexes/extracted.fmd

# Sequence lengths
seqtk comp $dir_base/amylase/amylase387.fa.gz | cut -f1,2 | gzip > $dir_base/ropebwt3/indexes/extracted.fmd.len.gz
```

```shell
# Prepare clusters (they are mandatory for cosigt > v0.1.0)
Rscript /lizardfs/guarracino/git/cosigt/cosigt_smk/workflow/scripts/cluster.r /lizardfs/guarracino/mempig/amylase/pggb.amylase387/amylase387.similarity.tsv /lizardfs/guarracino/mempig/amylase/pggb.amylase387/amylase387.similarity.clusters.json

ROI_BED=$dir_base/roi/roi.bed
REFERENCE_CRAM_FASTA=$dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

mkdir -p $dir_base/reads
cd $dir_base/ropebwt3/
ls $dir_base/cram/*.final.cram | while read READS_CRAM; do
    NAME=$(basename $READS_CRAM .final.cram)
    echo $READS_CRAM $INPUT_READ

    INPUT_READ=$dir_base/reads/$NAME.roi.fa.gz
    samtools view -T $REFERENCE_CRAM_FASTA -@ 16 -L $ROI_BED -M -b $READS_CRAM | samtools fasta -@ 24 - | bgzip -l 9 -@ 8 > $INPUT_READ
    for l in `seq 5 30`; do
        sbatch -c 8 -p tux --job-name "${NAME}-l${l}" $dir_base/scripts/run_ropebwt3.sh "$dir_base" "$NAME" "$INPUT_READ" "$l" "$dir_base/amylase/pggb.amylase387/amylase387.fa.gz.8fea6cf.11fba48.a2c905a.smooth.final.gfa" "$dir_base/amylase/pggb.amylase387/amylase387.similarity.clusters.json" "$dir_base/amylase/amylase387.fa.gz"
    done
done


(echo -n "l p " | tr ' ' '\t'; find $dir_base/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" -exec head -n 1 {} \; -quit | sed 's/#//g') > $dir_base/ropebwt3/amylase387.test.tsv
find $dir_base/ropebwt3/ -name "cosigt_genotype.*.l*.p*.tsv" | while read TSV; do
    l=$(basename $TSV .tsv | cut -f 3 -d '.' | sed 's/l//g')
    p=$(basename $TSV .tsv | cut -f 4 -d '.' | sed 's/p//g')
   # echo $TSV $NAME $l $p

    grep '^#' -v $TSV | awk -v OFS='\t' -v l=$l -v p=$p '{print(l,p,$0)}'
done | sort -k 1,1 -k 2,2n -k 3,3n -T /scratch >> $dir_base/ropebwt3/amylase387.test.tsv



# From https://github.com/sudmantlab/amylase_diversity_project/blob/main/pangenome/pggb/20231102_graph/haplotype_all_structures_plotinfo.chm13_colors.bed
chr1_chm13_103304997_103901127_0	99	115472	bundle0#CCCCCC
chr1_chm13_103304997_103901127_0	492840	596058	bundle1#666666


samtools faidx /lizardfs/guarracino/robertsonian_translocation/assemblies/chm13v2.0.fa.gz chr1
chr1:103305096-103420469
chr1:103797837-103901055
```