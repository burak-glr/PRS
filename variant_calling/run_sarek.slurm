#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=nextflow_sarek
#SBATCH --output=logs/sarek_%j.log
#SBATCH --error=logs/sarek_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=120:00:00

eval "$(conda shell.bash hook)"
conda activate env_nf

echo "Job started on $(hostname) at $(date)"

mkdir -p logs

nextflow run nf-core/sarek \
    -r 3.5.1 \
    -profile conda \
    -c sarek_slurm.config \
    --input sarek_samplesheet.csv \
    --outdir ./sarek_results \
    --igenomes_ignore \
    --fasta /cta/users/burak.guler/UV_seq/ref_gen/hg38/GRCh38.p14.genome.fa \
    --tools haplotypecaller \
    --joint_germline \
    --trim_fastq \
    --skip_tools baserecalibrator \
    --wes \
    -resume

echo "Job finished at $(date)"
