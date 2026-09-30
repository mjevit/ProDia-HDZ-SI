#!/bin/bash

# Usage function to display help
show_help() {
    echo "Usage: $0 -r <reads_directory> -o <output_directory> -d <reference_genome> [-e <email>]"
    echo
    echo "This script id's forward and reverse read fastas in a directory, then writes and submits a slurm which performs quality trimming using Trim Galore, aligns the reads to a specified reference genome using bwa-mem, converts the SAM file to BAM, sorts the BAM file (samtools), marks duplicates (picard), adds read group information (picard), indexes the sorted BAM file, and processes the BAM files through variant calling with GATK4."
    echo
    echo "Dependencies:"
    echo "  - SLURM (installed and configured on your cluster)"
    echo "  - bwa-mem (available as a module on your cluster)"
    echo "  - samtools (available as a module on your cluster)"
    echo "  - Picard (available as a module on your cluster)"
    echo "  - Cutadapt (available as a module on your cluster)"
    echo "  - FASTQC (available as a module on your cluster)"
    echo "  - Trim Galore (available as a module on your cluster)"
    echo "  - GATK4 (available as a module on your cluster)"
    echo
    echo "Options:"
    echo "  -r    Directory containing the Illumina reads"
    echo "  -o    Directory to save the output files and SLURM script"
    echo "  -d    Reference genome file for bwa-mem alignment"
    echo "  -t    Number of tasks per node (default: 8)"
    echo "  -e    Email address for SLURM notifications (optional)"
    echo "  -h    Show this help message"
    echo
    echo "IMPORTANT: This script assumes that the R1 and R2 files are named with the pattern *_R1_cat.fastq.gz and *_R2_cat.fastq.gz respectively."
    echo "           If your files are named differently, please modify the script accordingly."
    echo "           This script will submit the generated SLURM script in the specified output directory if no error occurs."
    echo "IMPORTANT: FULL PATHS must be provided for reads, output directory, and reference genome file or SLURMs will be incorrect"
    echo
    echo "Example usage:"
    echo "  $0 -r /path/to/reads -o /path/to/output -d /path/to/reference_genome.fa -e your_email@example.com"
}

# Default value for SLURM tasks per node
SLURM_NTASKS_PER_NODE=8

# Parse command-line arguments
while getopts "r:o:d:e:h" opt; do
  case $opt in
    r) READS_DIR="$OPTARG"
    ;;
    o) OUTPUT_DIR="$OPTARG"
    ;;
    d) REFERENCE_GENOME="$OPTARG"
    ;;
    e) EMAIL="$OPTARG"
    ;;
    t) SLURM_NTASKS_PER_NODE="$OPTARG"
    ;;
    h) show_help
       exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        show_help
        exit 1
    ;;
  esac
done

# Check if required directories and files are provided
if [ -z "$READS_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$REFERENCE_GENOME" ]; then
    echo "Reads directory, output directory, and reference genome file must be provided."
    show_help
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all pairs of Illumina reads in the directory
for SAMPLE in $(ls "$READS_DIR"/*_R1_cat.fastq.gz | sed 's/_R1_cat.fastq.gz//' | xargs -n 1 basename); do
  

  # Create a SLURM script for this sample
  SLURM_SCRIPT="$OUTPUT_DIR/${SAMPLE}_pipeline.slurm"
  
  cat <<EOT > "$SLURM_SCRIPT"
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_pipeline  # Job name
#SBATCH --output=$OUTPUT_DIR/${SAMPLE}_pipeline_out.%j   # Output file name with Job ID
#SBATCH --error=$OUTPUT_DIR/${SAMPLE}_pipeline_err.%j    # Error file name with Job ID
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$SLURM_NTASKS_PER_NODE                                 
#SBATCH --mem=40G                      
#SBATCH --time=140:00:00
#SBATCH --partition=guest              
EOT

  if [ -n "$EMAIL" ]; then
    cat <<EOT >> "$SLURM_SCRIPT"
#SBATCH --mail-type=END,FAIL           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=$EMAIL             # Where to send mail
EOT
  fi

  cat <<EOT >> "$SLURM_SCRIPT"

# Step 1: Quality trimming using Trim Galore
module purge
module load cutadapt
module load fastqc
module load trim_galore
trim_galore --paired --quality 30 --cores $SLURM_NTASKS_PER_NODE -o $OUTPUT_DIR $READS_DIR/${SAMPLE}_R1_cat.fastq.gz $READS_DIR/${SAMPLE}_R2_cat.fastq.gz && \\

# Step 2: Align reads to reference genome using bwa-mem
module purge
module load bwa
bwa mem -t $SLURM_NTASKS_PER_NODE $REFERENCE_GENOME $OUTPUT_DIR/${SAMPLE}_R1_cat_val_1.fq.gz $OUTPUT_DIR/${SAMPLE}_R2_cat_val_2.fq.gz > $OUTPUT_DIR/${SAMPLE}.sam && \\

# Step 3: Convert SAM to BAM
module purge
module load samtools
samtools view -bS --threads $SLURM_NTASKS_PER_NODE $OUTPUT_DIR/${SAMPLE}.sam > $OUTPUT_DIR/${SAMPLE}.bam && \\

#Step 4 - Collate bam
samtools collate --threads $SLURM_NTASKS_PER_NODE $OUTPUT_DIR/${SAMPLE}.bam $OUTPUT_DIR/${SAMPLE}_collated && \\

#Step 5 - Fixmate bam
samtools fixmate --threads $SLURM_NTASKS_PER_NODE -m $OUTPUT_DIR/${SAMPLE}_collated.bam $OUTPUT_DIR/${SAMPLE}fixmate.bam && \\

# Step 6: Sort BAM file
samtools sort --threads $SLURM_NTASKS_PER_NODE $OUTPUT_DIR/${SAMPLE}fixmate.bam -o $OUTPUT_DIR/${SAMPLE}_sorted.bam && \\

# Step 7: Add read group header
module purge
module load picard
picard -Xms512m -Xmx10g AddOrReplaceReadGroups I=$OUTPUT_DIR/${SAMPLE}_sorted.bam O=$OUTPUT_DIR/${SAMPLE}_sorted_rg.bam RGID=$SAMPLE RGPL=ILLUMINA RGLB=PE RGPU=1 RGSM=$SAMPLE && \\

# Step 8: Mark duplicates
picard -Xms512m -Xmx10g MarkDuplicates I=$OUTPUT_DIR/${SAMPLE}_sorted_rg.bam O=$OUTPUT_DIR/${SAMPLE}_sorted_dedup.bam M=$OUTPUT_DIR/${SAMPLE}_metrics.txt && \\

# Step 9: Index sorted BAM file
module purge
module load samtools
samtools index --threads $SLURM_NTASKS_PER_NODE -b $OUTPUT_DIR/${SAMPLE}_sorted_dedup.bam && \\

# Step 10: Call Variants using GATK
module purge
module load gatk4
gatk HaplotypeCaller --native-pair-hmm-threads $SLURM_NTASKS_PER_NODE -R $REFERENCE_GENOME -I $OUTPUT_DIR/${SAMPLE}_sorted_dedup.bam -O $OUTPUT_DIR/${SAMPLE}.g.vcf.gz --emit-ref-confidence GVCF


echo "Variant calling process complete for $SAMPLE. Filtered variants are saved in $OUTPUT_DIR/${SAMPLE}_filtered_variants.vcf.gz"
EOT

  echo "Generated SLURM script: $SLURM_SCRIPT"

  # Submit the SLURM script
  sbatch "$SLURM_SCRIPT"
  if [ $? -eq 0 ]; then
    echo "SLURM script submitted successfully: $SLURM_SCRIPT"
  else
    echo "Failed to submit SLURM script: $SLURM_SCRIPT"
  fi
done