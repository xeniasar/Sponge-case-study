#!/usr/bin/env bash

IN_DIR="results/bam_dedup_rmdup" 
OUT_DIR="results/gatk/gvcf"
REF="reference/genome.fna"
THREADS=40

IN_SUFFIX="_dedup.bam"

mkdir -p "${OUT_DIR}" logs

# Prepare reference indexes (safe to re-run)
samtools faidx "${REF}"

# GATK expects a sequence dictionary (*.dict) next to the reference
DICT="${REF%.*}.dict"
if [[ ! -f "${DICT}" ]]; then
  gatk CreateSequenceDictionary -R "${REF}"
fi

for BAM in "${IN_DIR}"/*"${IN_SUFFIX}"; do
  base="$(basename "${BAM}" "${IN_SUFFIX}")"
  OUT="${OUT_DIR}/${base}.g.vcf.gz"

  echo "[$(date)] HaplotypeCaller: ${base}"

  gatk HaplotypeCaller \
    -R "${REF}" \
    -I "${BAM}" \
    -O "${OUT}" \
    -ERC GVCF \
    --min-pruning 1 \
    --min-dangling-branch-length 1 \
    --native-pair-hmm-threads "${THREADS}" \
    2> "logs/${base}.haplotypecaller.stderr.log"

done
