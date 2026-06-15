#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --job-name=ti=tb_evidence
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --output=tbe.log

module load nextflow/25.04.6

export APPTAINER_TMPDIR=/scratch/user/uqctung/tmp
export APPTAINER_CACHEDIR=/scratch/user/uqctung/cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

export NXF_WORK=/scratch/user/uqctung/nextflow_work
mkdir -p $NXF_WORK

cd /QRISdata/Q9140/lmac/annot_lmac

apptainer exec \
    --nv \
    -B $PWD,$NXF_WORK,$APPTAINER_TMPDIR \
    /home/uqctung/containers/tiberius_latest.sif \
    python /opt/Tiberius/tiberius.py \
    --params_yaml profiles/bunya/nf_params.yaml \
    --nf_config profiles/bunya/nf.config