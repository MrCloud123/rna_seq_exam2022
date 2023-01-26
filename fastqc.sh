#!/bin/bash
#SBATCH --job-name="QC"
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=20GB
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=ji.yun@students.unibe.ch
#allocate ressources fastqc
module add UHTS/Quality_control/fastqc/0.11.9;
# link to fsdtq files contained in "READS" folder
ln -s /data/courses/rnaseq/toxoplasma_de/reads/*.fastq.gz .
# Run fastqc 
fastqc -t 10 /data/courses/rnaseq_course/toxoplasma_de/reads/*.fastq.gz -o /home/jyun/RNA-seq/out
#execute the script (sbatch fastqc.sh)

#After the executing, run the following code in the shell
#srun  --mem=20G --cpus-per-task=2 --time=01:00:00 --pty bash
#multiqc .

#produce stats for the fastq files with seqkit stat module

#srun  --mem=25G --cpus-per-task=2 --time=01:00:00 --pty bash
#module add UHTS/Analysis/SeqKit/0.13.2
#cd /data/courses/rnaseq_course/toxoplasma_de/reads/
#seqkit stat *.gz > /home/jyun/RNA-seq/QC/fastq_stats.txt

