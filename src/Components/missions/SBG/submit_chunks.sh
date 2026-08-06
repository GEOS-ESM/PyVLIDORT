#!/bin/bash
# submit_chunks.sh

ISO_T1="2006-01-16T17:35:00"
ISO_T2="2006-01-16T17:36:00"
YAML="sbg_vlidort.yaml"
NCHUNKS=4
NCH=213

CHUNK_SIZE=$(( (NCH + NCHUNKS - 1) / NCHUNKS ))  # ceiling division

for i in $(seq 0 $((NCHUNKS - 1))); do
    ICH_START=$((i * CHUNK_SIZE))
    ICH_END=$(( (i + 1) * CHUNK_SIZE ))
    if [ $ICH_END -gt $NCH ]; then
        ICH_END=$NCH
    fi

    sbatch --job-name="vlidort_ch${ICH_START}-${ICH_END}" \
           --time=12:00:00 \
           --output="vlidort_ch${ICH_START}-${ICH_END}-%j.log" \
           sbg_vlidort_chunks.j "$ISO_T1" "$ISO_T2" "$YAML" \
           --ich_start $ICH_START --ich_end $ICH_END
done
