# Sponge-case-study
BGE's WP11 Case study on Spongia officinalis pipeline used for the processing of whole-genome resequencing data.

This repository contains scripts used for SNP discovery from whole-genome resequencing (WGS) data.

The pipeline was developed for *Spongia officinalis* and can be adapted to other species with a reference genome.

---

## Pipeline overview

1. Low-quality read trimming and adapter removal (fastp)
2. Read mapping to reference genome (bwa) and BAM sorting (samtools)
3. Mark duplicates
4. Remove duplicates and index BAMs (optional)
5. Variant calling (GATK HaplotypeCaller)
6. Joint genotyping
7. Variant filtering

Each step is implemented as a separate script.

## Step 1: Read trimming (fastp)

Trim low-quality bases and adapter sequences using fastp.

```bash
bash scripts/01_fastp.sh
```
## Step 2: Read mapping (BWA-MEM)

Trimmed reads are mapped to the reference genome using BWA-MEM, and alignments are sorted with samtools.

```bash
bash scripts/02_bwa_mem.sh
```
## Step 3: Mark duplicates (sambamba)

PCR and optical duplicates are marked using sambamba.

```bash
bash scripts/03_mark_duplicates.sh
```
## Step 4: Remove duplicates and index BAMs (optional)

Optionally remove reads marked as duplicates and index the resulting BAM files.

```bash
bash scripts/04_remove_duplicates_and_index.sh
```
## Step 5: Variant calling (GATK HaplotypeCaller, GVCF mode)

Variants are called per sample in GVCF mode using GATK HaplotypeCaller.

To enable parallelisation across chromosomes, the genome was split into
~6 Mb intervals. HaplotypeCaller was run independently for each interval, and
the resulting per-interval GVCFs were merged per sample using
`gatk GatherVcfs`.

```bash
bash scripts/05a_haplotypecaller_gvcf.sh
bash scripts/05b_gather_gvcfs.sh
```
## Step 6: Joint genotyping preparation (GenomicsDBImport per chromosome)

Per-sample GVCFs are imported into GenomicsDB workspaces per chromosome.

```bash
bash scripts/06_genomicsdbimport_per_chrom.sh
```
## Step 7: Joint genotyping (GenotypeGVCFs)

GenomicsDB workspaces are genotyped per chromosome using GATK GenotypeGVCFs.

```bash
bash scripts/07_genotypegvcfs_per_chrom.sh
```
