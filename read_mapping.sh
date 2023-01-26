#!/bin/sh
#SBATCH --job-name="maping"
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100GB
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=ji.yun@students.unibe.ch

# Add all the necessary modules for mapping
module add UHTS/Aligner/hisat/2.2.1

# Calling the files directories and names of fasta files as variables
INDEX_DIR=/data/courses/rnaseq_course/toxoplasma_de/reference/hisat2
FASTQ_DIR=/data/courses/rnaseq_course/toxoplasma_de/reads
NAMES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")
REF_DIR=/data/courses/rnaseq_course/toxoplasma_de/reference
# Run hisat2 mapping jobs

hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821922_1.fastq.gz -2 $FASTQ_DIR/SRR7821922_2.fastq.gz -S SRR7821922.sam --rna-strandness RF 
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821937_1.fastq.gz -2 $FASTQ_DIR/SRR7821937_2.fastq.gz -S SRR7821937.sam --rna-strandness RF 
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821938_1.fastq.gz -2 $FASTQ_DIR/SRR7821938_2.fastq.gz -S SRR7821938.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821939_1.fastq.gz -2 $FASTQ_DIR/SRR7821939_2.fastq.gz -S SRR7821939.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821949_1.fastq.gz -2 $FASTQ_DIR/SRR7821949_2.fastq.gz -S SRR7821949.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821950_1.fastq.gz -2 $FASTQ_DIR/SRR7821950_2.fastq.gz -S SRR7821950.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821921_1.fastq.gz -2 $FASTQ_DIR/SRR7821921_2.fastq.gz -S SRR7821921.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821951_1.fastq.gz -2 $FASTQ_DIR/SRR7821951_2.fastq.gz -S SRR7821951.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821952_1.fastq.gz -2 $FASTQ_DIR/SRR7821952_2.fastq.gz -S SRR7821952.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821953_1.fastq.gz -2 $FASTQ_DIR/SRR7821953_2.fastq.gz -S SRR7821953.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821968_1.fastq.gz -2 $FASTQ_DIR/SRR7821968_2.fastq.gz -S SRR7821968.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821920_1.fastq.gz -2 $FASTQ_DIR/SRR7821920_2.fastq.gz -S SRR7821920.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821918_1.fastq.gz -2 $FASTQ_DIR/SRR7821918_2.fastq.gz -S SRR7821918.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821919_1.fastq.gz -2 $FASTQ_DIR/SRR7821919_2.fastq.gz -S SRR7821919.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821969_1.fastq.gz -2 $FASTQ_DIR/SRR7821969_2.fastq.gz -S SRR7821969.sam --rna-strandness RF
hisat2 -x $INDEX_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -1 $FASTQ_DIR/SRR7821970_1.fastq.gz -2 $FASTQ_DIR/SRR7821970_2.fastq.gz -S SRR7821970.sam --rna-strandness RF