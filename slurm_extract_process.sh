#!/bin/bash
#SBATCH -J pre_loop       ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics                ## Queue
#SBATCH -t 4:00:00             ## Walltime/duration of the job
#SBATCH --mem-per-cpu=32G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=2            ## Number of processors for your task


module load bcftools/1.10.1
module load R/4.0.0

#Can extract all samples with for sample in `bcftools query -l <VCF> `;

for sample in AB1 CB4856 JU2007 QX1211 JU258 DL238 ;
do
  bcftools view -s $sample -x -Oz /projects/b1059/analysis/WI-20210121/isotype_only/WI.20210121.hard-filter.isotype.bcsq.vcf.gz | bcftools query -f'%CHROM %POS %REF %ALT [%SAMPLE %TBCSQ{*}]\n' > ${sample}_BCSQ.vcf ;
  Rscript reformating_sampleBCSQ.R $sample
done
