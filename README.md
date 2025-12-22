# Sponge-case-study
BGE's WP11 Case study on Spongia officinalis pipeline used for the processing of whole-genome resequencing data.

This repository contains scripts used for SNP discovery from whole-genome resequencing (WGS) data.

The pipeline was developed for *Spongia officinalis* and can be adapted to other species with a reference genome.

---

## Pipeline overview

1. Low-quality read trimming and adapter removal (fastp)
3. Read mapping to reference genome (bwa)
4. BAM processing (sorting, marking duplicates)
5. Variant calling (GATK HaplotypeCaller)
6. Joint genotyping
7. Variant filtering

Each step is implemented as a separate script.

## Step 1: Read trimming (fastp)

Trim low-quality bases and adapter sequences using fastp.

```bash
bash scripts/01_fastp.sh

## Step 2: Read mapping (BWA-MEM)

Trimmed reads are mapped to the reference genome using BWA-MEM, and alignments are sorted with samtools.

```bash
bash scripts/02_bwa_mem.sh
