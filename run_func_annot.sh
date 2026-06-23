#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=96:00:00
#SBATCH --job-name=func_annot
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --output=funcannot.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

export TMPDIR=/scratch/user/uqctung/tmp
export XDG_CACHE_HOME=/scratch/user/uqctung/.cache
export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p $TMPDIR $XDG_CACHE_HOME $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

# Add this line explicitly before snakemake
cd /QRISdata/Q9140/lmac/annot_lmac

snakemake -s annotation_tiberius.smk --profile profiles/bunya/