#!/bin/bash

#SBATCH --job-name="genome_indexed"
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mem=12GB
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=ji.yun@students.unibe.ch
REF_DIR=/data/courses/rnaseq_course/toxoplasma_de/reference
module add UHTS/Aligner/hisat/2.2.1;
hisat2-build $REF_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz mg_indexed --ss $REF_DIR/mmus.ss --exon $REF_DIR/mmus.exon

#run in shell : sbatch genome_index