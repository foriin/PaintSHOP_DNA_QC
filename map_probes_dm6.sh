#!/bin/bash

# Check if the 'idx' directory exists
if [ ! -d "idx" ]; then
  mkdir -p idx

  # Check if 'fasta' folder exists, if not create one
  [ ! -d "fasta" ] && mkdir -p fasta

  # Navigate to the 'fasta' folder
  cd fasta

  # Download the latest Drosophila melanogaster genome
  latest_fasta=$(curl -l ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/ | grep 'dmel-all-chromosome-r' | sort -V | tail -n 1)
  wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/$latest_fasta

  # Unzip the fasta file
  gunzip $latest_fasta

  # Get the uncompressed fasta filename
  fasta_file=${latest_fasta%.gz}

  # Navigate back to the root directory
  cd ..

  # Run bowtie2-build to create the index in 'idx' folder
  bowtie2-build fasta/$fasta_file idx/dm6

else
  echo "The 'idx' folder already exists."
fi

fasta=$1

base_name=$(basename "$fasta" .fasta)

# Run bowtie2
bowtie2 -x idx/dm6 -f "$fasta" -a -p 12 -N 0 -S "${base_name}_dm6.sam"

# Convert SAM to BAM and sort
samtools view -b "${base_name}_dm6.sam" | samtools sort - > "${base_name}_dm6.srt.bam"

# Convert BAM to BED
mamba run -n process bedtools bamtobed -i "${base_name}_dm6.srt.bam" > "${base_name}_dm6.bed"

