#!/bin/bash
# submit_and_monitor.sh
# Submits channel chunks and monitors until all complete.
# Reports failures per channel.

ISO_T1="2006-01-16T17:35:00"
ISO_T2="2006-01-16T17:36:00"
YAML="sbg_vlidort.yaml"
NCH=213
CHUNK_SIZE=1   # number of channels per job

# Calculate number of jobs from chunk size
NCHUNKS=$(( (NCH + CHUNK_SIZE - 1) / CHUNK_SIZE ))

# Arrays to track jobs
declare -A JOB_CH_START
declare -A JOB_CH_END
JOB_IDS=()

echo "============================================="
echo "Submitting $NCHUNKS jobs ($CHUNK_SIZE channels each)..."
echo "============================================="

for i in $(seq 0 $((NCHUNKS - 1))); do
    ICH_START=$((i * CHUNK_SIZE))
    ICH_END=$(( (i + 1) * CHUNK_SIZE ))
    if [ $ICH_END -gt $NCH ]; then
        ICH_END=$NCH
    fi

    JOBID=$(sbatch --parsable \
           --job-name="vlidort_ch${ICH_START}-${ICH_END}" \
           --output="vlidort_ch${ICH_START}-${ICH_END}-%j.log" \
           sbg_vlidort_chunks.j "$ISO_T1" "$ISO_T2" "$YAML" \
           --ich_start $ICH_START --ich_end $ICH_END)

    JOB_IDS+=("$JOBID")
    JOB_CH_START[$JOBID]=$ICH_START
    JOB_CH_END[$JOBID]=$ICH_END

    echo "  Submitted job $JOBID: channels ${ICH_START}-${ICH_END}"
done

echo ""
echo "============================================="
echo "All $NCHUNKS jobs submitted. Monitoring..."
echo "============================================="
echo ""

# Monitor loop
POLL_INTERVAL=60  # seconds between checks
COMPLETED=0
FAILED=0
SUCCEEDED=0
TOTAL=${#JOB_IDS[@]}

declare -A JOB_DONE

while [ $COMPLETED -lt $TOTAL ]; do
    sleep $POLL_INTERVAL

    for JOBID in "${JOB_IDS[@]}"; do
        # Skip already processed jobs
        if [ "${JOB_DONE[$JOBID]}" == "1" ]; then
            continue
        fi

        # Check job state
        STATE=$(sacct -j "$JOBID" --format=State --noheader -P | head -1 | tr -d ' ')

        case "$STATE" in
            COMPLETED)
                SUCCEEDED=$((SUCCEEDED + 1))
                COMPLETED=$((COMPLETED + 1))
                JOB_DONE[$JOBID]=1
                echo "[$(date '+%H:%M:%S')] ✓ Job $JOBID COMPLETED: channels ${JOB_CH_START[$JOBID]}-${JOB_CH_END[$JOBID]}"
                ;;
            FAILED|CANCELLED|CANCELLED+|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
                FAILED=$((FAILED + 1))
                COMPLETED=$((COMPLETED + 1))
                JOB_DONE[$JOBID]=1
                echo "[$(date '+%H:%M:%S')] ✗ Job $JOBID FAILED ($STATE): channels ${JOB_CH_START[$JOBID]}-${JOB_CH_END[$JOBID]}"
                echo "    Log: vlidort_ch${JOB_CH_START[$JOBID]}-${JOB_CH_END[$JOBID]}-${JOBID}.log"
                ;;
            *)
                # Still running/pending — do nothing
                ;;
        esac
    done

    # Progress update
    RUNNING=$((TOTAL - COMPLETED))
    echo "[$(date '+%H:%M:%S')] Progress: $SUCCEEDED succeeded, $FAILED failed, $RUNNING running/pending"
done

# Final summary
echo ""
echo "============================================="
echo "ALL JOBS COMPLETE"
echo "============================================="
echo "  Total:     $TOTAL"
echo "  Succeeded: $SUCCEEDED"
echo "  Failed:    $FAILED"
echo "  Chunk size: $CHUNK_SIZE channels/job"
echo ""

if [ $FAILED -gt 0 ]; then
    echo "FAILED JOBS:"
    echo "---------------------------------------------"
    for JOBID in "${JOB_IDS[@]}"; do
        STATE=$(sacct -j "$JOBID" --format=State --noheader -P | head -1 | tr -d ' ')
        if [ "$STATE" != "COMPLETED" ]; then
            echo "  Job $JOBID: channels ${JOB_CH_START[$JOBID]}-${JOB_CH_END[$JOBID]} ($STATE)"
            echo "    Log: vlidort_ch${JOB_CH_START[$JOBID]}-${JOB_CH_END[$JOBID]}-${JOBID}.log"
        fi
    done
    echo ""
    echo "To resubmit failed channels, run:"
    for JOBID in "${JOB_IDS[@]}"; do
        STATE=$(sacct -j "$JOBID" --format=State --noheader -P | head -1 | tr -d ' ')
        if [ "$STATE" != "COMPLETED" ]; then
            echo "  sbatch --job-name=\"vlidort_ch${JOB_CH_START[$JOBID]}-${JOB_CH_END[$JOBID]}\" --time=1:00:00 --output=\"vlidort_ch${JOB_CH_START[$JOBID]}-${JOB_CH_END[$JOBID]}-%j.log\" sbg_vlidort_chunks.j \"$ISO_T1\" \"$ISO_T2\" \"$YAML\" --ich_start ${JOB_CH_START[$JOBID]} --ich_end ${JOB_CH_END[$JOBID]}"
        fi
    done
    exit 1
else
    echo "All jobs succeeded!"
    exit 0
fi
