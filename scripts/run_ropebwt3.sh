#!/bin/bash

# Arguments passed from the main script
dir_base="$1"
NAME="$2"
INPUT_READ="$3"
L="$4"
GRAPH_GFA="$5"
CLUSTER_JSON="$6"
FASTA="$7"

# Load required modules or set up environment if necessary
hostname

# Loop through p values as specified
for p in $(seq 1 50); do
    mkdir -p /scratch/${NAME}.l${L}.p${p}
    cd /scratch/${NAME}.l${L}.p${p}
    
    ropebwt3 mem -l "$L" "$dir_base/ropebwt3/indexes/extracted.fmd" "$INPUT_READ" -p "$p" -t 8 > "${NAME}-vs-extracted.mem.l${L}.p${p}.tsv"
    
    python3 "$dir_base/ropebwt3-to-paf.py" "${NAME}-vs-extracted.mem.l${L}.p${p}.tsv" <(cut -f 1,2 "$FASTA.fai") "${NAME}-vs-extracted.mem.l${L}.p${p}.paf"
    
    gfainject --gfa "$GRAPH_GFA" --paf "${NAME}-vs-extracted.mem.l${L}.p${p}.paf" > "${NAME}-vs-extracted.mem.l${L}.p${p}.gaf"
    
    gafpack -g "$GRAPH_GFA" -a "${NAME}-vs-extracted.mem.l${L}.p${p}.gaf" --len-scale --weight-queries | pigz -9 -p 8 > "${NAME}-vs-extracted.mem.l${L}.p${p}.gafpack.gz"
    
    cosigt -p "$dir_base/odgi/paths_matrix.tsv.gz" -g "${NAME}-vs-extracted.mem.l${L}.p${p}.gafpack.gz" -o . -i "$NAME" -c "$CLUSTER_JSON"
    
    mv cosigt_genotype.tsv "$dir_base/ropebwt3/cosigt_genotype.${NAME}.l${L}.p${p}.tsv"
    
    pigz -9 -p 8 sorted_combos.tsv
    mv sorted_combos.tsv.gz "$dir_base/ropebwt3/sorted_combos.${NAME}.l${L}.p${p}.tsv.gz"
    
    # Cleanup intermediate files
    cd /scratch
    rm /scratch/${NAME}.l${L}.p${p} -r
done
