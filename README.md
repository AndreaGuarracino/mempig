# MEMPIG (MEM Pangenome Injection Genotyping)

## Paths

```shell
dir_base=/lizardfs/guarracino/mempig

export PATH="/lizardfs/guarracino/tools/bedtools2/bin:$PATH"
export PATH="/lizardfs/guarracino/tools/samtools-1.22.1:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/impg/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/wfmash/build/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/seqwish/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/smoothxg/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/GFAffix/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/odgi/bin:$PATH" # v0.9.3
export PATH="/lizardfs/guarracino/tools_for_genotyping/pggb:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/gfainject/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/gafpack/target/release:$PATH" # f08c04b4596bd9c40b0e62641f37c2500b956a1a
export PATH="/lizardfs/guarracino/tools_for_genotyping/panplexity/target/release:$PATH" # 46d065c1b6669657479736ae127eb16d5eea454f
export PATH="/lizardfs/guarracino/tools_for_genotyping/bwa-mem2-2.2.1_x64-linux:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/cosigt:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/agc-1.1_x64-linux:$PATH"

export PATH="/lizardfs/guarracino/git/seqtk:$PATH"

#export PATH="/lizardfs/guarracino/condatools/emboss/6.6.0/bin/stretcher:$PATH" 
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
git checkout 312f69b645cf99ea669dfb6388dd15898f1ae4bd
go mod init cosigt && go mod tidy && go build cosigt
#Add 'export PATH="/lizardfs/guarracino/git/cosigt:$PATH"' to ~/.zshrc
```

### minimap2

```shell
conda create --prefix /lizardfs/guarracino/condatools/minimap2/2.28 -c conda-forge -c bioconda minimap2=2.28 -y
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

# NOT DOWNLOADED: NA24385 (from project PRJEB35491)
#cd /scratch
#wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR368/ERR3684866/CT21796_NA24385_NIST_UDI22_HiseqX_19122017_Proband_S1.bam

# MISSING: NA21487 is missing
```

<!-- ### 39 samples

```shell
cd $dir_base/cram
wget https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index

grep -Ff <(cat $dir_base/impg/extracted.fasta.fai | cut -f 1 -d '#' | sort | uniq) $dir_base/cram/1000G_698_related_high_coverage.sequence.index | cut -f 1  > $dir_base/cram/to_download.txt
wget -i $dir_base/cram/to_download.txt
ls $dir_base/cram/*.cram | while read CRAM; do samtools index $CRAM; done
``` -->

## Tests

### 63 HGSVC3's short read samples on regions from the same 63 HGSVC3's assemblies

<!-- Run:

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
``` -->

#### Assemblies vs reference

```shell
conda activate /lizardfs/guarracino/condatools/minimap2/2.28

mkdir -p $dir_base/minimap2
cd $dir_base/minimap2
ls /lizardfs/guarracino/pangenomes/HGSVC3/*.fa.gz | while read fasta; do
    sample=$(basename $fasta .fa.gz)

    sbatch -p allnodes -c 24 --job-name $sample-vs-grch38 --wrap "hostname; cd /scratch; minimap2 -x asm20 --eqx -c -t 24 $dir_base/reference/GRCh38.fa.gz $fasta | bgzip -l 9 -@ 4 > $dir_base/minimap2/$sample-vs-grch38.paf.gz"
done

conda deactivate
```

#### Genotyping

```shell
mkdir -p $dir_base/regions_of_interest
cat $dir_base/data/loci.bed | while read -r chrom start end name; do echo -e "$chrom\t$start\t$end\t$name" > $dir_base/regions_of_interest/$name.bed; done

#cat /lizardfs/guarracino/genotypify/data/loci/tableS1_locityper_paper.bed | grep HLA | cut -f 1-4 | while read -r chrom start end name; do echo -e "$chrom\t$start\t$end\t$name" > $dir_base/regions_of_interest/$name.bed; done
# echo -e "chr6\t31972057\t32055418\tC4" > $dir_base/regions_of_interest/C4.bed
# echo -e "chr22\t42077656\t42253758\tCYP2D6" > $dir_base/regions_of_interest/CYP2D6.bed
# echo -e "chr16\t70801860\t71230722\tHYDIN" > $dir_base/regions_of_interest/HYDIN.bed
# echo -e "chr2\t88691673\t88750929\tRPIA" > $dir_base/regions_of_interest/RPIA.bed
# echo -e "chr7\t141963337\t142192359\tMGAM" > $dir_base/regions_of_interest/MGAM.bed
# echo -e "chr20\t14638177\t14962407\tMACROD2" > $dir_base/regions_of_interest/MACROD2.bed
# echo -e "chr16\t83532778\t83743290\tCDH13" > $dir_base/regions_of_interest/CDH13.bed
# echo -e "chr1\t155089880\t155330966\tMTX1" > $dir_base/regions_of_interest/MTX1.bed
# #echo -e "chr8\t43199245\t43428244\tPOTEA" > $dir_base/regions_of_interest/POTEA.bed # BROKEN LOCUS
# echo -e "chr10\t133331967\t133667386\tCYP2E1" > $dir_base/regions_of_interest/CYP2E1.bed
echo -e "chr1\t103304997\t103901127\tAMY_andrea" > $dir_base/regions_of_interest/AMY_andrea.bed
echo -e "chr1\t103511079\t103832645\tAMY1A_SV_list" > $dir_base/regions_of_interest/AMY1A_SV_list.bed
echo -e "chr1\t103456064\t103863972\tAMY_davide" > $dir_base/regions_of_interest/AMY_davide.bed
echo -e "chr15\t82245249\t82528509\tGOLGA6L10" > $dir_base/regions_of_interest/GOLGA6L10.bed
#echo -e "chr15\t84069272\t84348216\tGOLGA6L4" > $dir_base/regions_of_interest/GOLGA6L4.bed
echo -e "chr15\t84029272\t84388216\tGOLGA6L4_ref" > $dir_base/regions_of_interest/GOLGA6L4-ref.bed
echo -e "chr16\t76429\t86429\tMPG" > $dir_base/regions_of_interest/MPG.bed
echo -e "chr16\t30748293\t30761591\tPHKG2" > $dir_base/regions_of_interest/PHKG2.bed

# panplexity parameters
w=100
k=16
d=100
t=auto

# odgi chop parameter
c=0 

ls $dir_base/regions_of_interest | cut -f 1 | cut -f 1 -d '.' | while read region2; do
    chrom=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 1)
    start=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 2)
    end=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 3)
    region=${chrom}_${start}_${end}

    echo $region

    mkdir -p $dir_base/impg/$chrom/$region
    ls $dir_base/minimap2/*.paf.gz | grep -Ff <(ls $dir_base/cram/*.cram | while read f; do echo $(basename $f .final.cram); done) -w | while read paf; do
        sample=$(basename $paf .paf.gz)

        impg query \
            -p $paf \
            -b $dir_base/regions_of_interest/$region2.bed
    done > $dir_base/impg/$chrom/$region/$region.projected.bedpe
    bedtools sort -i $dir_base/impg/$chrom/$region/$region.projected.bedpe | \
        > $dir_base/impg/$chrom/$region/$region.impg.bedpe
    #    bedtools merge -d 200000 | grep '#U#' -v > $dir_base/impg/$chrom/$region/$region.merged.bed

    bash /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/bedpe_merge.sh \
        $dir_base/impg/$chrom/$region/$region.impg.bedpe \
        200000 \
        > $dir_base/impg/$chrom/$region/$region.merged.bed
    bash /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/bedpe_filter.sh \
        $dir_base/impg/$chrom/$region/$region.merged.bed \
        $region \
        1000 \
        > $dir_base/impg/$chrom/$region/$region.merged.filtered.bed

    (bedtools getfasta -fi $dir_base/reference/GRCh38.fa.gz -bed $dir_base/regions_of_interest/$region2.bed | sed 's/>chr/>GRCh38#0#chr/g';
    ls $dir_base/minimap2/*.paf.gz | grep -Ff <(ls $dir_base/cram/*.cram | while read f; do echo $(basename $f .final.cram); done) -w | while read paf; do
        sample=$(basename $paf "-vs-grch38.paf.gz")
        fasta=/lizardfs/guarracino/pangenomes/HGSVC3/$sample.fa.gz

        bedtools getfasta -fi $fasta -bed <(grep $sample -w $dir_base/impg/$chrom/$region/$region.merged.filtered.bed)
    done) | bgzip -@ 8 > $dir_base/impg/$chrom/$region/$region.extracted.fa.gz
    samtools faidx $dir_base/impg/$chrom/$region/$region.extracted.fa.gz

    ######################################################################################################
    #conda activate /lizardfs/guarracino/condatools/bwa-mem2/2.2.1
    bwa-mem2 index $dir_base/impg/$chrom/$region/$region.extracted.fa.gz
    ls $dir_base/cram/*cram | while read cram; do
        echo $cram
        sample=$(basename $cram .cram)
        mkdir -p $dir_base/bwamem2/$sample/$chrom/$region
        # Extract reads covering the region and then align them against the pangenome
        samtools view \
            -T $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
            -L $dir_base/regions_of_interest/$region2.bed \
            -M \
            -b \
            $cram | \
            samtools sort -n | \
            samtools fasta | \
                bwa-mem2 mem -t 42 $dir_base/impg/$chrom/$region/$region.extracted.fa.gz - | \
                samtools view -b -F 4 -@ 6 - \
                > $dir_base/bwamem2/$sample/$chrom/$region/$region.realigned.bam
    done
    #conda deactivate
    #=====================================================================================================
    # dir_ropebwt3=$dir_base/bedtools/getfasta/$chrom/$region
    # mkdir -p $dir_ropebwt3
    # # Construct a BWT for both strands of the input sequences
    # ropebwt3 build $dir_base/impg/$chrom/$region/$region.extracted.fa.gz -do $dir_ropebwt3/$region.fasta.gz.fmd
    # # Sampled suffix array
    # ropebwt3 ssa -o $dir_ropebwt3/$region.fasta.gz.fmd.ssa -s8 -t 24 $dir_ropebwt3/$region.fasta.gz.fmd
    # # Sequence lengths
    # seqtk comp $dir_base/impg/$chrom/$region/$region.extracted.fa.gz | cut -f1,2 | gzip > $dir_ropebwt3/$region.fasta.gz.fmd.len.gz
    # fasta=$dir_base/impg/$chrom/$region/$region.extracted.fa.gz
    # ls $dir_base/cram/*cram | while read cram; do
    #     echo $cram
    #     sample=$(basename $cram .cram)
    #     mkdir -p $dir_base/ropebwt3/$sample/$chrom/$region
    #     # Extract reads covering the $region region and MEME them against the pangenome
    #     samtools view \
    #         -T $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    #         -L $dir_base/regions_of_interest/$region2.bed \
    #         -M \
    #         -b \
    #         $cram | \
    #         samtools sort -n | \
    #         samtools fasta -F 0x0 | \
    #             ropebwt3 mem -l 29 $dir_ropebwt3/$region.fasta.gz.fmd - -p 1 -t 6 > $dir_base/ropebwt3/$sample/$chrom/$region/$sample.reads_vs_extracted.mem.tsv
    #     python3 $dir_base/scripts/ropebwt3-to-paf.py $dir_base/ropebwt3/$sample/$chrom/$region/$sample.reads_vs_extracted.mem.tsv <(cut -f 1,2 $fasta.fai) $dir_base/ropebwt3/$sample/$chrom/$region/$region.realigned.paf
    # done
    ######################################################################################################

    mkdir -p $dir_base/pggb/$chrom/$region
    pggb -i $dir_base/impg/$chrom/$region/$region.extracted.fa.gz -o  $dir_base/pggb/$chrom/$region -t 48 -D /scratch -c 2 -k 101
    mv $dir_base/pggb/$chrom/$region/*smooth.final.og $dir_base/pggb/$chrom/$region/$region.final.og # Rename the final ODGI graph in a more human-friendly way
    odgi sort \
        -i $dir_base/pggb/$chrom/$region/$region.final.og \
        -Y -P -t 48 \
        -H <(odgi paths -i $dir_base/pggb/$chrom/$region/$region.final.og -L | grep GRCh38) \
        -o $dir_base/pggb/$chrom/$region/$region.chop0.og  \
        --temp-dir /scratch # | \
        #odgi chop -i - -o $dir_base/pggb/$chrom/$region/$region.chop$c.og -c $c
done

ls $dir_base/regions_of_interest | cut -f 1 | cut -f 1 -d '.' | while read region2; do
    chrom=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 1)
    start=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 2)
    end=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 3)
    region=${chrom}_${start}_${end}

    echo $region

    mkdir -p $dir_base/odgi/paths/matrix/$chrom
    odgi paths -i $dir_base/pggb/$chrom/$region/$region.chop$c.og \
        -H | \
        cut -f 1,4- | \
        gzip > $dir_base/odgi/paths/matrix/$chrom/$region.paths_matrix.tsv.gz

    mkdir -p $dir_base/odgi/view/$chrom
    odgi view \
        -i $dir_base/pggb/$chrom/$region/$region.chop$c.og \
        -g > $dir_base/odgi/view/$chrom/$region.gfa

    ls $dir_base/cram/*cram | while read cram; do
        echo $cram
        sample=$(basename $cram .cram)

        mkdir -p $dir_base/gfainject/$sample/$chrom
        # gfainject \
        #     --gfa $dir_base/odgi/view/$chrom/$region.gfa \
        #     --paf $dir_base/ropebwt3/$sample/$chrom/$region/$region.realigned.paf | \
        #     gzip > $dir_base/gfainject/$sample/$chrom/$region.gaf.gz
        gfainject \
            --gfa $dir_base/odgi/view/$chrom/$region.gfa \
            --bam $dir_base/bwamem2/$sample/$chrom/$region/$region.realigned.bam | \
            gzip > $dir_base/gfainject/$sample/$chrom/$region.gaf.gz
    done

    ls $dir_base/cram/*cram | while read cram; do
        echo $cram
        sample=$(basename $cram .cram)

        mkdir -p $dir_base/gafpack/$sample/$chrom
        gafpack \
            --gfa $dir_base/odgi/view/$chrom/$region.gfa \
            --gaf $dir_base/gfainject/$sample/$chrom/$region.gaf.gz \
            --len-scale \
            --weight-queries | \
            gzip > $dir_base/gafpack/$sample/$chrom/$region.gafpack.gz
    done

    # First generate the similarity matrix using odgi
    mkdir -p $dir_base/odgi/dissimilarity/$chrom
    odgi similarity \
        -i $dir_base/pggb/$chrom/$region/$region.chop$c.og --all --distances \
        > $dir_base/odgi/dissimilarity/$chrom/$region.tsv
  
    mkdir -p $dir_base/clusters/$chrom
    grep '^S' $dir_base/odgi/view/$chrom/$region.gfa | awk '{{print("node."$2,length($3))}}' OFS="\t" > $dir_base/odgi/view/$chrom/$region.node.length.tsv
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/filter.r \
        $dir_base/odgi/paths/matrix/$chrom/$region.paths_matrix.tsv.gz \
        $dir_base/odgi/view/$chrom/$region.node.length.tsv \
        no_filter \
        $dir_base/odgi/paths/matrix/$chrom/$region.shared.tsv
    region_similarity=$(cut -f 3 $dir_base/odgi/paths/matrix/$chrom/$region.shared.tsv | tail -1)

    # Run the clustering script to generate the JSON
    Rscript /lizardfs/guarracino/git/cosigt/cosigt_smk/workflow/scripts/cluster.r \
        $dir_base/odgi/dissimilarity/$chrom/$region.tsv \
        $dir_base/clusters/$chrom/$region.clusters.json \
        automatic \
        $region_similarity
done

ls $dir_base/regions_of_interest | cut -f 1 | cut -f 1 -d '.' | while read region2; do
    chrom=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 1)
    start=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 2)
    end=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 3)
    region=${chrom}_${start}_${end}

    echo $region
    mkdir -p /scratch/cosigt

    # # Generate random mask file with 221 nodes
    # # Each node has 50% chance of being 0 or 1
    # awk 'BEGIN {
    #     srand();  # Initialize random seed
    #     for (i = 1; i <= 221; i++) {
    #         print int(rand() * 2);  # Generates 0 or 1 with equal probability
    #     }
    # }' > /scratch/cosigt/$region.mask
    panplexity -i $dir_base/odgi/view/$chrom/$region.gfa -w $w -k $k -t $t -d $d -m /scratch/cosigt/$region.mask -b /scratch/cosigt/$region.linguistic.bed -c /scratch/cosigt/$region.linguistic.csv
    
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/coverage_outliers.r \
                $dir_base/odgi/paths/matrix/$chrom/$region.paths_matrix.tsv.gz \
                /scratch/cosigt/$region.final \
                $dir_base/odgi/view/$chrom/$region.node.length.tsv \
                /scratch/cosigt/$region.mask

    ls $dir_base/cram/*cram | \
        parallel -j 48 '
            echo {}
            sample=$(basename {} .cram)
            
            mkdir -p /scratch/cosigt/$sample/'$chrom'
            cosigt \
                -i $sample \
                -p '$dir_base'/odgi/paths/matrix/'$chrom'/'$region'.paths_matrix.tsv.gz \
                -g '$dir_base'/gafpack/$sample/'$chrom'/'$region'.gafpack.gz \
                -c '$dir_base'/clusters/'$chrom'/'$region'.clusters.json \
                -o /scratch/cosigt/$sample/'$chrom'/'$region' \
                -m /scratch/cosigt/'$region'.final.mask.tsv
        '
done

ls $dir_base/regions_of_interest | cut -f 1 | cut -f 1 -d '.' | while read region2; do
    chrom=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 1)
    start=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 2)
    end=$(cat $dir_base/regions_of_interest/$region2.bed | cut -f 3)
    region=${chrom}_${start}_${end}

    echo $region

    mkdir -p /scratch/benchmark/$chrom/$region

    # Step 1: Make TPR table
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/calc_tpr.r \
        $dir_base/odgi/dissimilarity/$chrom/$region.tsv \
        $dir_base/clusters/$chrom/$region.clusters.json \
        $dir_base/clusters/$chrom/$region.clusters.hapdist.tsv \
        /scratch/benchmark/$chrom/$region/tpr.tsv \
        /scratch/cosigt/*/$chrom/$region/*sorted_combos.tsv.gz

    # Step 2: Flip PGGB graph
    odgi flip \
        -i  $dir_base/pggb/$chrom/$region/$region.chop$c.og \
        -o  /scratch/benchmark/$chrom/$region/$region.flip.og -P -t 28 \
        --ref-flips <(grep '^P' $dir_base/odgi/view/$chrom/$region.gfa | cut -f 2 | grep "GRCh38#0")

    # Step 3: Convert OG to FASTA
    odgi paths \
        -i /scratch/benchmark/$chrom/$region/$region.flip.og \
        -f  -t 28 | sed 's/_inv$//g' > /scratch/benchmark/$chrom/$region/$region.flip.fasta

    # Step 4: Prepare combinations for QV
    mkdir -p /scratch/benchmark/$chrom/$region/qv_prep
    bash /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/prepare_qv.sh \
        /scratch/benchmark/$chrom/$region/tpr.tsv \
        /scratch/benchmark/$chrom/$region/$region.flip.fasta \
        /scratch/benchmark/$chrom/$region/qv_prep

    # Step 5: Calculate QV for each sample (parallelized)
    find /scratch/benchmark/$chrom/$region/qv_prep -mindepth 1 -maxdepth 1 -type d | \
    parallel -j 48 '
        sample=$(basename {})
        echo "Calculating QV for sample $sample"
        bash /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/calculate_qv.sh \
            {} \
            {}/qv.tsv
    '

    # Step 6: Combine QV results
    cat /scratch/benchmark/$chrom/$region/qv_prep/*/qv.tsv > /scratch/benchmark/$chrom/$region/bestqv.tsv

    # Step 7: Combine TPR and QV
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/combine_tpr_qv.r \
        /scratch/benchmark/$chrom/$region/tpr.tsv \
        /scratch/benchmark/$chrom/$region/bestqv.tsv \
        $region \
        /scratch/benchmark/$chrom/$region/tpr_qv.tsv
done

# Final step: Plot TPR results
cat $dir_base/regions_of_interest/*.bed > tmp.bed
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/plot_tpr.r \
    /scratch/benchmark/tpr \
    tmp.bed \
    /scratch/benchmark/*/*/tpr_qv.tsv
rm tmp.bed

# To reduce storage space
#pigz -9 /scratch/cosigt/*/*/*/sorted_combos.tsv
#pigz -9 /scratch/benchmark/*/*/*.og
rm /scratch/benchmark/*/*/*.fai
#(ls /scratch/benchmark/*/*/*.fasta; ls /scratch/benchmark/*/*/qv_prep/*/*.fasta) | while read f; do bgzip -l 9 -@ 48 $f; done
(ls /scratch/benchmark/*/*/*.fasta; ls /scratch/benchmark/*/*/qv_prep/*/*.fasta) | while read f; do rm $f; done


# ls _* -d | while read d; do echo $d; cat $dir_base/$d/benchmark/*/*/tpr_qv.tsv > $d.tsv; done
# cat $dir_base/regions_of_interest/*.bed > tmp.bed
# Rscript ./z4.R ciao tmp.bed $(ls *tsv | sort -V)

# cd /scratch/random-grid-search
# rm _random*.tsv ciao*
# ls _random* -d | while read d; do echo $d; cat /scratch/random-grid-search/$d/benchmark/*/*/tpr_qv.tsv > $d.tsv; done
# cat /lizardfs/guarracino/mempig/regions_of_interest/*.bed > tmp.bed
# Rscript /lizardfs/guarracino/mempig/z4.R ciao tmp.bed $(ls *tsv | sort -V)
```

### 39 short read samples on the C4 region from the HPRCy1 pangenome

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
