#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --job-name=cont_report
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --output=cont_report.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

# Set up temp/cache directories
export TMPDIR=/scratch/user/uqctung/tmp
export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

# Unload java/11 (snakemake8 loads it; can conflict)
module unload java/11.0.27 2>/dev/null || true

echo "Starting contamination analysis (report only)..."
echo "Results will be in: contamination_screening/"
echo ""

# Run Snakemake
snakemake -s contamination_check.smk \
    --cores 8 \
    -pr \
    2>&1 | tee contamination_report_run.log

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo ""
echo "Key output files:"
echo ""
echo "1. contamination_screening/01_contamination_report.txt"
echo "   → Detailed breakdown by contamination category"
echo "   → Recommended thresholds for manual review"
echo ""
echo "2. contamination_screening/02_scaffold_summary.tsv"
echo "   → Excel-friendly: all scaffolds with contamination %"
echo "   → Sort/filter to find which to remove"
echo ""
echo "3. contamination_screening/03_species_breakdown.txt"
echo "   → What species are contaminating your assembly"
echo ""
echo "NEXT STEPS:"
echo "1. Review the reports above"
echo "2. Decide which scaffolds to remove"
echo "3. Create: contamination_screening/removal_list.txt"
echo "   (one scaffold ID per line)"
echo "4. Run: contamination_filtering.smk (when ready)"
echo ""