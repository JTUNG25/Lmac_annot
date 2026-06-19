#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=96:00:00
#SBATCH --job-name=tb_evidence
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --output=tbe.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

module unload java/11.0.27
module load java/21.0.8
module load nextflow/25.04.6

export PATH=~/bin:$PATH 
export APPTAINER_TMPDIR=/scratch/user/uqctung/tmp
export APPTAINER_CACHEDIR=/scratch/user/uqctung/cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

export NXF_WORK=/scratch/user/uqctung/nextflow_work
mkdir -p $NXF_WORK

cd /QRISdata/Q9140/lmac/annot_lmac

python Tiberius/tiberius.py \
    --params_yaml profiles/bunya/nf_params.yaml \
    --nf_config profiles/bunya/nf.config
    -resume