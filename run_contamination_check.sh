#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --job-name=contam_check
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --output=contam_check.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

cd /QRISdata/Q9140/lmac/annot_lmac

export TMPDIR=/scratch/user/uqctung/tmp
export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

snakemake -s contamination_check.smk --cores 4 --profile profiles/bunya/