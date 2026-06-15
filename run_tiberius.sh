#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --job-name=tiberius
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --output=tb.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

export TMPDIR=/scratch/user/uqctung/tmp
mkdir -p $TMPDIR

snakemake -s tiberius.smk --profile profiles/bunya/