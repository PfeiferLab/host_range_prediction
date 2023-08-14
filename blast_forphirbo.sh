#!/bin/bash
#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of "tasks" (default: 1 core per task)
#SBATCH -t 7-00:00:00   # time in d-hh:mm:ss
#SBATCH --mem=16000mb
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --export=NONE   # Purge the job-submitting shell environment

module load blast_plus/2.12.0
blastn -db ref_prok_rep_genomes -query phage.fasta -outfmt '6 scomnames' > phage.blastn.txt
awk '!x[$0]++' phage.blastn.txt > phage.blastn.uniq.txt

for file in phage_hosts/*
do
base=$(basename "$file")
blastn -db ref_prok_rep_genomes -query phage_hosts/${base} -outfmt '6 scomnames' > phage_hostsblast/${base}.blastn.txt
awk '!x[$0]++' phage_hostsblast/${base}.blastn.txt > phage_hostsblast/${base}.blastn.uniq.txt
done

