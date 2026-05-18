#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --job-name=lmac_annotation
#SBATCH --output=sm.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

cd /QRISdata/Q9140/lmac/annot_lmac

# ── Temp directories on QRISdata (visible to all compute nodes) ───────────────
export TMPDIR=/QRISdata/Q9140/lmac/annot_lmac/tmp
mkdir -p $TMPDIR
export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

echo "=================================================="
echo "L. maculans Target-Specific Annotation Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "=================================================="

# ── Step 1: Pre-pull all SIF images ──────────────────────────────────────────
SIF_DIR=/QRISdata/Q9140/lmac/annot_lmac/sifs
mkdir -p $SIF_DIR

echo "Pulling container images..."
[ ! -f $SIF_DIR/biopython.sif ] && \
    apptainer pull --tmpdir $APPTAINER_TMPDIR \
        $SIF_DIR/biopython.sif \
        docker://quay.io/biocontainers/biopython:1.81

[ ! -f $SIF_DIR/eggnog.sif ] && \
    apptainer pull --tmpdir $APPTAINER_TMPDIR \
        $SIF_DIR/eggnog.sif \
        docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0

[ ! -f $SIF_DIR/diamond.sif ] && \
    apptainer pull --tmpdir $APPTAINER_TMPDIR \
        $SIF_DIR/diamond.sif \
        docker://quay.io/biocontainers/diamond:2.1.8--h43eeafb_0

echo "✓ All container images ready"

# ── Step 2: Download UniProt/DIAMOND database ─────────────────────────────────
mkdir -p databases

if [ ! -f databases/uniprot_sprot.dmnd ]; then
    echo "Downloading UniProt Swiss-Prot..."
    wget -q https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
        -O databases/uniprot_sprot.fasta.gz
    gunzip databases/uniprot_sprot.fasta.gz

    echo "Building DIAMOND database..."
    apptainer exec \
        --bind /QRISdata \
        $SIF_DIR/diamond.sif \
        diamond makedb \
            --in databases/uniprot_sprot.fasta \
            -d databases/uniprot_sprot \
            --threads 4

    rm databases/uniprot_sprot.fasta
    echo "✓ DIAMOND database ready"
else
    echo "✓ DIAMOND database already exists"
fi

# ── Step 3: Download eggNOG database ─────────────────────────────────────────
if [ ! -d databases/eggnog ] || [ ! -f databases/eggnog/eggnog.db ]; then
    echo "Downloading eggNOG database (this takes ~20-30 min)..."
    mkdir -p databases/eggnog

    apptainer exec \
        --bind /QRISdata \
        $SIF_DIR/eggnog.sif \
        download_eggnog_data.py \
            --data_dir databases/eggnog \
            -y

    echo "✓ eggNOG database ready"
else
    echo "✓ eggNOG database already exists"
fi

echo "✓ All databases ready"

# ── Step 4: Run annotation pipeline ──────────────────────────────────────────
chmod +x profiles/bunya/status-sacct-robust.sh

echo "Running annotation pipeline..."

snakemake -s annotation.smk --unlock --profile profiles/bunya/

snakemake -s annotation.smk \
    --profile profiles/bunya/ \