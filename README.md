# ProDia-HDZ-SI

Scripts used in the analysis pipeline for a publication paper. This repository stores utilities for processing Illumina sequencing reads through quality trimming, alignment, and variant calling on an HPC/SLURM cluster. 

!!!Important!!!
As this script is written it is designed to be run on the SWAN HPRC at University of Nebraska-Lincoln, running on other HPRC may require editting.

## Contents

- `prop_raw2variant_v5.sh` — Identifies paired-end (R1/R2) FASTQ files in a directory, then generates and submits a SLURM job per sample that runs the following pipeline:
  1. Quality trimming with **Trim Galore**
  2. Alignment to a reference genome with **bwa-mem**
  3. SAM → BAM conversion with **samtools**
  4. BAM collation and fixmate (**samtools**)
  5. BAM sorting (**samtools**)
  6. Read group assignment (**Picard** `AddOrReplaceReadGroups`)
  7. Duplicate marking (**Picard** `MarkDuplicates`)
  8. BAM indexing (**samtools**)
  9. Variant calling with **GATK4** `HaplotypeCaller` (GVCF output)

## Dependencies

The following tools must be available as modules on your cluster:

- SLURM (installed and configured)
- bwa-mem
- samtools
- Picard
- Cutadapt
- FastQC
- Trim Galore
- GATK4

## Usage

```bash
./prop_raw2variant_v5.sh -r <reads_directory> -o <output_directory> -d <reference_genome> [-e <email>] [-t <tasks_per_node>]
```

### Options

| Flag | Description |
|------|-------------|
| `-r` | Directory containing the Illumina reads |
| `-o` | Directory to save output files and the generated SLURM script |
| `-d` | Reference genome file for bwa-mem alignment |
| `-t` | Number of tasks per node (default: 8) |
| `-e` | Email address for SLURM notifications (optional) |
| `-h` | Show help message |

### Example

```bash
./prop_raw2variant_v5.sh -r /path/to/reads -o /path/to/output -d /path/to/reference_genome.fa -e your_email@example.com
```

### Notes

- R1/R2 FASTQ files must follow the naming pattern `"$READS_DIR"/*_R1_cat.fastq.gz` / `"$READS_DIR"/*_R2_cat.fastq.gz`. Update the script if your files are named differently.
- **Full paths** must be provided for the reads directory, output directory, and reference genome file, otherwise the generated SLURM scripts will be incorrect.
- The script automatically submits the generated SLURM job for each sample unless an error occurs.
