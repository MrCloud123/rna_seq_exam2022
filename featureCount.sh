#!/bin/sh
#SBATCH --job-name="maping"
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH --nodes=7
#SBATCH --mem=12GB
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=ji.yun@students.unibe.ch
#SBATCH --array=0-15

module add UHTS/Aligner/hisat/2.2.1
#download the featureCounts in https://sourceforge.net/projects/subread/files/subread-2.0.3/ into /data/users/jyun/
/data/users/jyun/subread-2.0.2-Linux-x86_64/bin/featureCounts -B -a /data/courses/rnaseq_course/toxoplasma_de/reference/Mus_musculus.GRCm39.108.gtf.gz -p -T 8 -o matrix.txt shibai/SRR7821921.bam shibai/SRR7821922.bam shibai/SRR7821918.bam shibai/SRR7821919.bam shibai/SRR7821920.bam shibai/SRR7821937.bam shibai/SRR7821938.bam shibai/SRR7821939.bam shibai/SRR7821949.bam shibai/SRR7821950.bam shibai/SRR7821951.bam shibai/SRR7821952.bam shibai/SRR7821953.bam shibai/SRR7821968.bam shibai/SRR7821969.bam shibai/SRR7821970.bam
