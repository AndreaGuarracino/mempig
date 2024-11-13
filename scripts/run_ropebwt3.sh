#!/bin/bash

# Arguments passed from the main script
DIR_BASE="$1"
NAME="$2"
INPUT_READ="$3"
L="$4"

# Load required modules or set up environment if necessary
hostname

# Loop through p values as specified
for p in $(seq 1 50); do
    mkdir -p /scratch/${NAME}.l${L}.p${p}
    cd /scratch/${NAME}.l${L}.p${p}
    
    ropebwt3 mem -l "$L" "$DIR_BASE/ropebwt3/indexes/extracted.fmd" "$INPUT_READ" -p "$p" -t 8 > "${NAME}-vs-extracted.mem.l${L}.p${p}.tsv"
    
    python3 "$DIR_BASE/ropebwt3-to-paf.py" "${NAME}-vs-extracted.mem.l${L}.p${p}.tsv" <(cut -f 1,2 "$DIR_BASE/impg/extracted.fasta.fai") "${NAME}-vs-extracted.mem.l${L}.p${p}.paf"
    
    gfainject --gfa "$DIR_BASE/odgi/graph.gfa" --paf "${NAME}-vs-extracted.mem.l${L}.p${p}.paf" > "${NAME}-vs-extracted.mem.l${L}.p${p}.gaf"
    
    gafpack -g "$DIR_BASE/odgi/graph.gfa" -a "${NAME}-vs-extracted.mem.l${L}.p${p}.gaf" --len-scale | pigz -9 -p 8 > "${NAME}-vs-extracted.mem.l${L}.p${p}.gafpack.gz"
    
    cosigt -p "$DIR_BASE/odgi/paths_matrix.tsv.gz" -g "${NAME}-vs-extracted.mem.l${L}.p${p}.gafpack.gz" -o . -i "$NAME" -c "$DIR_BASE/odgi/chopped.similarity.clusters.json"
    
    mv cosigt_genotype.tsv "$DIR_BASE/ropebwt3/cosigt_genotype.${NAME}.l${L}.p${p}.tsv"
    
    pigz -9 -p 8 sorted_combos.tsv
    mv sorted_combos.tsv.gz "$DIR_BASE/ropebwt3/sorted_combos.${NAME}.l${L}.p${p}.tsv.gz"
    
    # Cleanup intermediate files
    cd /scratch
    rm /scratch/${NAME}.l${L}.p${p} -r
done
