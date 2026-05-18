#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=lmac_annotation
#SBATCH --output=annot.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8
cd /QRISdata/Q9140/lmac/annot_lmac
export TMPDIR=/QRISdata/Q9140/lmac/annot_lmac/tmp
mkdir -p $TMPDIR

echo "=================================================="
echo "L. maculans Target-Specific Annotation Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "=================================================="

# Check if databases exist, setup if needed
if [ ! -f "databases/uniprot_sprot.dmnd" ]; then
    echo "Setting up UniProt database (first time only)..."
    snakemake -s annotation.smk --profile profiles/bunya/ setup_databases
fi

if [ ! -d "databases/eggnog" ]; then
    echo "Setting up eggNOG database (first time only)..."
    conda activate eggnog-mapper
    download_eggnog_data.py --data_dir databases/eggnog
    conda activate snakemake8
fi

echo "✓ All databases ready"

# Run annotation pipeline
echo "Running annotation pipeline..."
snakemake -s annotation.smk --profile profiles/bunya/

echo "=================================================="
echo "Pipeline completed!"
echo "End time: $(date)"
echo "=================================================="

# Print results
echo "=== RESULTS ==="
find results -name "*.xlsx" -o -name "*.csv" | head -10
