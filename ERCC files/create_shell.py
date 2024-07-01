import sys
files =  sys.argv[1:]

out = open("command.sh", "w+")
out.write("#!/bin/bash\n")
out.write("#SBATCH --partition=preempt --time=04:00:00 --output=h_rundata2.log -c 30 --mem=120G\n\n")

adapter = "TGGAATTCTCGGGTGCCAAGGA"

for j in files: 
    i = j.split(".")[0]
    out.write(f"cutadapt -j 30 -m 16 -a {adapter} -o {i}_trimmed.fastq {i}.fastq\n")
    out.write(f"~/.conda/envs/mirge38/bin/bowtie /scratch/public/apatil8/custom_data_cp/index/sequence -n 0 -5 4 -3 4 --norc -S --threads 30 {i}_trimmed.fastq {i}_umiHead.sam\n")
    out.write(f"samtools view {i}_umiHead.sam > {i}_reads_template.txt\n")
    out.write(f"python add_umiRead_andCount.py {i}_trimmed.fastq {i}_reads_template.txt\n")
    out.write(f"python prepare_Matt_panels.py {i}_trimmed_umi_bwtMapCounts.csv\n")
    out.write("\n")

