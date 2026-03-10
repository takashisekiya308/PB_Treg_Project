#!/bin/bash
set -euo pipefail

# ---- tools ----
STAR_BIN="PATH_to_STAR_reference"
RSEM_BIN="rsem-calculate-expression"

# ---- references ----
STAR_GENOME_DIR="PATH_to_STAR_genome_dir"
RSEM_REF="PATH_to_RSEM_reference"

# ---- resources ----
THREADS=4

# ---- sanity checks ----
if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: STAR not executable: $STAR_BIN" >&2
  exit 1
fi
if [[ ! -d "$STAR_GENOME_DIR" ]]; then
  echo "ERROR: STAR genomeDir not found: $STAR_GENOME_DIR" >&2
  exit 1
fi

if [[ ! -f "${RSEM_REF}.grp" ]]; then
  echo "ERROR: RSEM reference not found (missing ${RSEM_REF}.grp): $RSEM_REF" >&2
  exit 1
fi

echo "Using STAR: $STAR_BIN"
"$STAR_BIN" --version || true

shopt -s nullglob

for r1 in *_R1.fastq.gz; do
  base="${r1%_R1.fastq.gz}"
  r2="${base}_R2.fastq.gz"

  if [[ ! -f "$r2" ]]; then
    echo "ERROR: missing pair file: $r2 (for $r1)" >&2
    exit 1
  fi

  echo "==== Processing: $base ===="

  prefix="${base}."

  tx_bam="${prefix}Aligned.toTranscriptome.out.bam"
  if [[ ! -f "$tx_bam" ]]; then
    echo "ERROR: STAR did not produce transcriptome BAM: $tx_bam" >&2
    echo "Check STAR log: ${prefix}Log.final.out / ${prefix}Log.out" >&2
    exit 1
  fi

  "$RSEM_BIN" \
    --num-threads "$THREADS" \
    --paired-end \
    --bam "$tx_bam" \
    "$RSEM_REF" \
    "$base"

  echo "==== Done: $base ===="
done

echo "All done."
